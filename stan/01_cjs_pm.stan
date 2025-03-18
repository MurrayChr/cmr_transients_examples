// CJS model p(t)phi(t) i.e. time-dependent detection and survival, without 
// transients
// 
// This code implements the product-multinomial likelihood
//

functions {
  matrix get_multinomial_probs(
    data int T,
    vector p,
    vector phi
  ) {
    // matrix to store multinomial probabilities
    matrix[T-1, T] pr;
    // convenience parameter
    vector[T] q = 1 - p;
    // define multinomial probabilities
    for (t in 1:(T-1)) {
      for (j in 1:(t-1)) {
        pr[t][j] = 0;
      }
      for (j in t:(T - 1)) {
        pr[t][j] = prod(phi[t:j]) * prod(q[(t+1):j]) * p[j+1] ;
      }
      pr[t][T] = 1 - sum(pr[t][1:(T-1)]);
    }
    return pr;
  }
}

data {
  int<lower=2> T;                      // number of years
   array[T-1, T] int<lower=0> marr;    // m-array 
}

parameters {
  vector<lower=0,upper=1>[T] p;         // detection probability for 'residents'
  vector<lower=0,upper=1>[T-1] phi;     // survival probability
} 

model {
  // multinomial probabilities
  matrix[T-1, T] pr;      
  pr = get_multinomial_probs(T, p ,phi);
  
 // priors
  p ~ beta(1,1);
  phi ~ beta(1,1);
  
  // likelihood...
  for (t in 1:(T - 1)) {
    marr[t] ~ multinomial(pr[t]');
  }
}
