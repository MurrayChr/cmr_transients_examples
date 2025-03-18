//
// Multisite-plus-transients model with time- and site-dependent survival, site-dependent 
// movement and time- and site-dependent detection. Here we treat the simplest 
// case of two sites.
//
// This code implements the hierarchical hmm likelihood
//
// The model has four 'hidden' states for transients and residents at each of 
// the two sites:
// 1 - resident at site 1
// 2 - transient, first detected at site 1
// 3 - resident at site 2
// 4 - transient, first detected at site 2
//

functions {
  // create state transition matrices
  array[] matrix get_Gamma_matrices(
    data int T,
    array[] vector phi,
    array[] vector m
  ) {
    array[T-1] matrix[5,5] Gamma;
    for (t in 1:(T-1)) {
      Gamma[t][1] = [phi[1][t]*m[1][1], 0, phi[1][t]*m[1][2], 0, 1 - phi[1][t]];
      Gamma[t][2] = [                0, 0,                 0, 0,             1];
      Gamma[t][3] = [phi[2][t]*m[2][1], 0, phi[2][t]*m[2][2], 0, 1 - phi[2][t]];
      Gamma[t][4] = [                0, 0,                 0, 0,             1];
      Gamma[t][5] = [                0, 0,                 0, 0,             1];
    }
    return Gamma;
  }
  // create observation matrices
  array[] matrix get_Omega_matrices(
    data int T,
    array[] vector p
  ) {
    array[T] matrix[5,3] Omega;
    for (t in 1:T) {
      Omega[t][1] = [p[1][t],       0, 1-p[1][t]];
      Omega[t][2] = [      0,       0,         1];
      Omega[t][3] = [      0, p[2][t], 1-p[2][t]];
      Omega[t][4] = [      0,       0,         1];
      Omega[t][5] = [      0,       0,         1];
    }
    return Omega;
  } 
}

data {
  int<lower=2> T;                            // number of years
  int<lower=1> N;                            // number of unique capture histories
  array[N,T] int y;                          // unique capture histories
  array[N] int<lower=1,upper=T-1> fc;        // first capture occasion
  array[N] int<lower=1,upper=2> fc_state;  // first capture occasion
  array[N] int<lower=1> mult;                // capture history multiplicities
}

parameters {
  array[2] vector<lower=0,upper=1>[T] p;        // detection probabilities
  array[2] vector<lower=0,upper=1>[T-1] phi;    // survival probabilities
  array[2] vector<lower=0,upper=1>[T-1] pi_r;   // residency probabilities
  array[2] simplex[2] m;                        // movement probabilities 
}

model {
  // transition and emission matrices
  array[T-1] matrix[5,5] Gamma = get_Gamma_matrices(T, phi, m);
  array[T] matrix[5,3] Omega = get_Omega_matrices(T, p);
  
  // likelihood
  for (n in 1:N) {
    // initialise forward algorithm
    array[T] row_vector[5] fwd;
    if (fc_state[n]==1) {
      fwd[fc[n]] = [pi_r[fc_state[n]][fc[n]], 1-pi_r[fc_state[n]][fc[n]], 0, 0, 0];
    }
    if (fc_state[n]==2) {
      fwd[fc[n]] = [0, 0, pi_r[fc_state[n]][fc[n]], 1-pi_r[fc_state[n]][fc[n]], 0];
    }
    // recursion
    for (t in fc[n]:(T-1)) {
      fwd[t+1] = fwd[t]*diag_post_multiply(Gamma[t], Omega[t+1][:,y[n,t+1]]);
    }
    // increment log-likelihood 
    target += mult[n]*log(sum(fwd[T]));
  }
}
