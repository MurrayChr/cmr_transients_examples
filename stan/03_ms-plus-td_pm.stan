//
// Multi-site model with time- and site-dependent survival, site-dependent 
// movement and time- and site-dependent detection, and trap-dependence modelled
// using the approach of Pradel and Sanz-Aguilar. Here we treat the simplest 
// case of two sites, thus there are four states:
//  1 - alive at first site, trap-aware
//  2 - alive at first site, trap-unaware
//  3 - alive at second site, trap-aware
//  4 - alive at second site, trap-unaware
//
// This code implements the product-multinomial likelihood for data in the form 
// of an m-array. 
//
// Note 1: Pradel and Sanz-Aguilar frame their approach as a multievent model,
// but it can be formulated as a multistate model with unobservable 'trap-unaware' 
// states and we do so here.
//

functions {
  matrix get_multinomial_probs(
    data int T,
    data real eps,
    array[] vector pA,
    array[] vector pU,
    array[] vector phi,
    array[] vector m
  ) {
    // matrix to store multinomial probabilities
    matrix[4*(T-1), 4*(T-1)+1] pr = rep_matrix(0.0, 4*(T-1), 4*(T-1)+1); 
    
    // convenience parameters
    array[2] vector[T] qA;
    array[2] vector[T] qU;
    for (i in 1:2) {
      qA[i] = 1 - pA[i];
      qU[i] = 1 - pU[i];
    }
    
    // auxillary matrices to define multinomial probabilities
    array[T-1] matrix[4,4] gamma;   // state transition probabilities
    array[T] matrix[4,2] omega;     // state-dependent detection probabilities
    // define gamma
    for (t in 1:(T-1)) {
      // define 2x2 site blocks
      for (I in 1:2) {      // I indexes source site
        int first_i = 2*(I-1) + 1;
        int last_i = 2*I;
        for (J in 1:2) {    // J indexes destination site
          int first_j = 2*(J-1) + 1;
          int last_j = 2*J;
          gamma[t][first_i:last_i, first_j:last_j] = [
            [phi[I][t]*m[I][J]*pA[J][t+1], phi[I][t]*m[I][J]*qA[J][t+1]],
            [phi[I][t]*m[I][J]*pU[J][t+1], phi[I][t]*m[I][J]*qU[J][t+1]]
          ];
        }
      }
    }
    // define omega
    // note: here omega's are the same for all t, so one could just define
    // omega as a single matrix; we retain this code for ease of generalisation
    for (t in 1:T) {
      omega[t] = [
        [1-eps,   eps],
        [  eps, 1-eps],
        [1-eps,   eps],
        [  eps, 1-eps]
      ];
    }
    // define entries of pr using matrix multiplication
    for (rI in 1:(T-1)) { // block row index, reverse order
      int I = T - rI;     // I runs from T-1 to 1
      // row indices to for 4x4 block
      int i1 = 1 + 4*(I-1);
      int i2 = 4*I;
      // diagonal block
      pr[i1:i2, i1:i2] = diag_post_multiply(gamma[I], omega[I+1][:,1]); 
      if (I < T-1) {
        matrix[4,4] temp = diag_post_multiply(gamma[I], omega[I+1][:,2]);
        for (J in (I+1):(T-1)) {
          // column indices for 4x4 block
          int j1 = 1 + 4*(J-1);
          int j2 = 4*J;
          pr[i1:i2, j1:j2] = temp * pr[(i1 + 4):(i2 + 4), j1:j2];
        }
      }
    }
    // enforce row-sum-to-one constraint
    for (i in 1:(4*(T-1))){
      pr[i, 4*(T-1)+1] = 1 - sum(pr[i, 1:4*(T-1)]);
    }
    return pr;
  }
}

data {
  int<lower=2> T;                           // number of years
  array[4*(T-1), 4*(T-1)+1] int<lower=0> marr;    // m-array
}

transformed data {
  real<lower=0> eps = pow(10, -15);
}

parameters {
  // detection probabilities
  array[2] vector<lower=0,upper=1>[T] pA;     // for trap-aware indiv
  array[2] vector<lower=0,upper=1>[T] pU;     // for trap-unaware indiv
  array[2] vector<lower=0,upper=1>[T-1] phi;  // survival probabilities
  array[2] simplex[2] m;                      // movement probabilities 
}

model {
  // calculate multinomial probabilities
  matrix[4*(T-1), 4*(T-1)+1] pr;
  pr = get_multinomial_probs( T, eps, pA, pU, phi, m );
  // m-array capture-recapture likelihood
  for (i in 1:(4*(T-1))) {
    marr[i] ~ multinomial(pr[i]');
  }
}
