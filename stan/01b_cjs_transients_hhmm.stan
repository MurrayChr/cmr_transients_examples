// CJS model p(t)phi(t) i.e. time-dependent detection and survival, also 
// accounting for time-dependent proportion of residents to transients in each cohort
// 
// This code implements the 'multievent' or hierarchical hidden Markov model likelihood
//

data {
  int<lower=2> T;                   // number of years
  int<lower=1> N;                   // number of unique capture histories
  array[N,T] int y;                 // unique capture histories
  array[N] int<lower=1,upper=T-1> fc;     // first capture 
  array[N] int<lower=1> mult;             // capture history multiplicities
}

parameters {
  vector<lower=0,upper=1>[T] p;         // detection probability for 'residents'
  vector<lower=0,upper=1>[T-1] phi;     // survival probability
  vector<lower=0,upper=1>[T-1] pi_r;    // proportion 'residents' in each cohort
} 

model {
  // transition and emission matrices
  array[T-1] matrix[3,3] Gamma;  // transition
  array[T] matrix[3,2] Omega;    // emission/observation
  for (t in 1:(T-1)) {
    Gamma[t][1] = [phi[t], 0, 1-phi[t]];
    Gamma[t][2] = [     0, 0,        1];
    Gamma[t][3] = [     0, 0,        1];
  }
  for (t in 1:T) {
    Omega[t][1] = [p[t], 1-p[t]];
    Omega[t][2] = [   0,      1];
    Omega[t][3] = [   0,      1];
  }
  
  // priors
  phi ~ beta(1,1);
  pi_r ~ beta(1,1);
  p ~ beta(1,1);
  
  // likelihood
  for (n in 1:N) {
    // initialise forward algorithm
    array[T] row_vector[3] fwd;
    fwd[fc[n]] = [pi_r[fc[n]], 1-pi_r[fc[n]], 0];
    // recursion
    for (t in fc[n]:(T-1)) {
      fwd[t+1] = fwd[t]*diag_post_multiply(Gamma[t], Omega[t+1][:,y[n,t+1]]);
    }
    // increment log-likelihood 
    target += mult[n]*log(sum(fwd[T]));
  }
}
