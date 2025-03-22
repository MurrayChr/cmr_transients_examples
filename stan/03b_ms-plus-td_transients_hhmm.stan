//
// Multisite-plus-trap-dependence model incorporating transients. Trap-dependence 
// is modelled using the approach of Pradel and Sanz-Aguilar. Here we treat the 
// case of two sites, thus there are four states for residents:
//  1 - alive at first site, trap-aware
//  2 - alive at first site, trap-unaware
//  3 - alive at second site, trap-aware
//  4 - alive at second site, trap-unaware
//
// This code implements the hierarchical hmm likelihood
//
// The model has six 'hidden' states (and a 'dead' state):
// 1 - resident at site 1, trap-aware
// 2 - resident at site 1, trap-unaware
// 3 - transient, first detected at site 1
// 4 - resident at site 2, trap-aware
// 5 - resident at site 2, trap-unaware
// 6 - transient, first detected at site 2
// 7 - 'dead'
//

functions {
  // create state transition matrices
  array[] matrix get_Gamma_matrices(
    data int T,
    array[] vector pA,
    array[] vector pU,
    array[] vector phi,
    array[] vector m
  ) {
    // convenience parameters
    array[2] vector[T] qA;
    array[2] vector[T] qU;
    for (i in 1:2) {
      qA[i] = 1 - pA[i];
      qU[i] = 1 - pU[i];
    }
    // transition matrices (DEFINE THESE IN BLOCKS...)
    array[T-1] matrix[7,7] Gamma;
    for (t in 1:(T-1)) {
      // define 3x3 site block
      for (I in 1:2) {   // I indexes source site
        int first_i = 3*(I-1) + 1;
        int last_i = 3*I;
        for (J in 1:2) { // J indexes destination site
          int first_j = 3*(J-1) + 1;
          int last_j = 3*J;
          Gamma[t][first_i:last_i, first_j:last_j] = [
            [phi[I][t]*m[I][J]*pA[J][t+1], phi[I][t]*m[I][J]*qA[J][t+1], 0],
            [phi[I][t]*m[I][J]*pU[J][t+1], phi[I][t]*m[I][J]*qU[J][t+1], 0],
            [                           0,                            0, 0]
          ];
        }
      }
      // last row
      Gamma[t][7] = [0, 0, 0, 0, 0, 0, 1];
      // last column
      Gamma[t][:,7] = [1- phi[1][t], 1- phi[1][t], 1, 1- phi[2][t], 1- phi[2][t], 1, 1]';
    }
    return Gamma;
  }
  
  // create observation matrices
  array[] matrix get_Omega_matrices(
    data int T
  ) {
    array[T] matrix[7,5] Omega;
    for (t in 1:T) {
      Omega[t] = [
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1],
        [0, 0, 0, 0, 1],
        [0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1],
        [0, 0, 0, 0, 1],
        [0, 0, 0, 0, 1]
      ];
    }
    return Omega;
  }
}

data {
  int<lower=2> T;                      // number of years
  int<lower=1> N;                      // number of unique capture histories
  array[N,T] int y;                    // unique capture histories
  array[N] int<lower=1,upper=T-1> fc;  // first capture occasion
  array[N] int fc_state;               // first capture state
  array[N] int<lower=1> mult;          // capture history multiplicities
}

parameters {
  // detection probabilities
  array[2] vector<lower=0,upper=1>[T] pA;       // for trap-aware indiv
  array[2] vector<lower=0,upper=1>[T] pU;       // for trap-unaware indiv
  array[2] vector<lower=0,upper=1>[T-1] phi;    // survival probabilities
  array[2] vector<lower=0,upper=1>[T-1] pi_r;   // residency probabilities
  array[2] simplex[2] m;                        // movement probabilities 
}

model {
  // transition and emission matrices
  array[T-1] matrix[7,7] Gamma = get_Gamma_matrices(T, pA, pU, phi, m);
  array[T] matrix[7,5] Omega = get_Omega_matrices(T);
  
  // likelihood
  for (n in 1:N) {
    // initialise forward algorithm
    array[T] row_vector[7] fwd;
    if (fc_state[n]==1) { // first capture at site 1
      fwd[fc[n]] = [pi_r[1][fc[n]], 0, 1-pi_r[1][fc[n]], 0, 0, 0, 0];
    }
    if (fc_state[n]==3) { // first capture at site 2
      fwd[fc[n]] = [0, 0, 0, pi_r[2][fc[n]], 0, 1-pi_r[2][fc[n]], 0];
    }
    // recursion
    for (t in fc[n]:(T-1)) {
      fwd[t+1] = fwd[t]*diag_post_multiply(Gamma[t], Omega[t+1][:,y[n,t+1]]);
    }
    // increment log-likelihood 
    target += mult[n]*log(sum(fwd[T]));
  }
}
