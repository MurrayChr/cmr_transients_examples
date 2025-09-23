//
// Multisite-plus-transients model with time- and site-dependent survival, site-dependent 
// movement and time- and site-dependent detection. Here we treat the simplest 
// case of two sites, thus there are two states:
// 1 - alive at first site
// 2 - alive at second site
//
// This code implements the extended product-multinomial likelihood
//
// Important note: here the m-array data must be constructed from *only* those
// individuals that were recaptured at least once after first capture
//

functions {
  matrix get_multinomial_probs(
    data int T,
    array[] vector p,
    array[] vector phi,
    array[] vector m
  ) {
    // matrix to store multinomial probabilities
    matrix[2*(T-1), 2*(T-1)+1] pr = rep_matrix(0.0, 2*(T-1), 2*(T-1)+1);
    // auxillary matrices to define multinomial probabilities
    array[T-1] matrix[2,2] gamma;   // state transition probabilities
    array[T] matrix[2,2] omega;     // state-dependent detection probabilities
    for (t in 1:(T-1)) {
      gamma[t][1] = [phi[1][t]*m[1][1], phi[1][t]*m[1][2]];
      gamma[t][2] = [phi[2][t]*m[2][1], phi[2][t]*m[2][2]];
    }
    for (t in 1:T) {
      omega[t][1] = [p[1][t], 1-p[1][t]];
      omega[t][2] = [p[2][t], 1-p[2][t]];
    }
    // define entries of pr using matrix multiplication
    for (rI in 1:(T-1)) { // block row index, reverse order
      int I = T - rI;     // I runs from T-1 to 1
      // row indices to for 2x2 block
      int i1 = 1 + 2*(I-1);
      int i2 = 2*I;
      // diagonal block
      pr[i1:i2, i1:i2] = diag_post_multiply(gamma[I], omega[I+1][:,1]); 
      if (I < T-1) {
        matrix[2,2] temp = diag_post_multiply(gamma[I], omega[I+1][:,2]);
        for (J in (I+1):(T-1)) {
          // column indices for 2x2 block
          int j1 = 1 + 2*(J-1);
          int j2 = 2*J;
          pr[i1:i2, j1:j2] = temp * pr[(i1 + 2):(i2 + 2), j1:j2];
        }
      }
    }
    // enforce row-sum-to-one constraint
    for (i in 1:(2*(T-1))){
      pr[i, 2*(T-1)+1] = 1 - sum(pr[i, 1:2*(T-1)]);
    }
      return pr;
    }
}

data {
  int<lower=2> T;                                 // number of years
  array[2*(T-1), 2*(T-1)+1] int<lower=0> marr;    // m-array (see above)
  array[2, T-1] int<lower=0> Nmult;       // no. indiv per cohort and state that were ever recaptured
  array[2, T-1] int<lower=0> Nsingle;     // no. indiv per cohort and state that were never recaptured
}

parameters {
  array[2] vector<lower=0,upper=1>[T] p;        // detection probabilities
  array[2] vector<lower=0,upper=1>[T-1] phi;    // survival probabilities
  array[2] vector<lower=0,upper=1>[T-1] pi_r;   // residency probabilities
  array[2] simplex[2] m;                        // movement probabilities 
}

model {
  // calculate multinomial probabilities
  matrix[2*(T-1), 2*(T-1)+1] pr;
  pr = get_multinomial_probs( T, p, phi, m );

  // m-array capture-recapture likelihood...
  // ... for indiv recaptured at least once
  for (i in 1:(2*(T-1))) {
    marr[i] ~ multinomial(pr[i]');
  }
  for (t in 1:(T-1)) {
    for (s in 1:2) {
      target += Nmult[s,t] * log(pi_r[s,t]);
    }
  }
  // ... for indiv never recaptured 
  for (t in 1:(T-1)) {
    for (s in 1:2) {
      int i = s+2*(t-1);   // indexes row s in block t
      target += Nsingle[s,t] * log( pr[i][2*(T-1)+1] * pi_r[s,t] + (1 - pi_r[s,t]) );
    }
  }
}
