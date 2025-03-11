// CJS model p(t)phi(t) i.e. time-dependent detection and survival, also 
// accounting for time-dependent proportion of residents to transients in each cohort
// multinomial formulation

// Important note: here the m-array data must be constructed from *only* those
// individuals that were recaptured at least once after first capture

data {
  int<lower=2> T;                      // number of years
   array[T-1, T] int<lower=0> marr;    // m-array (see above)
   array[T-1] int<lower=0> N_1;        // no. indiv per cohort that were ever recaptured
   array[T-1] int<lower=0> N_0;        // no. indiv per cohort that were never recaptured
}

parameters {
  vector<lower=0,upper=1>[T] p;         // detection probability for 'residents'
  vector<lower=0,upper=1>[T-1] phi;     // survival probability
  vector<lower=0,upper=1>[T-1] pi_r;    // proportion 'residents' in each cohort
} 

model {
  // declare object to store multinomial probabilities
  array[T-1] vector[T] pr;      // 'standard' probabilities...

  // convenience parameter
  vector[T] q=1-p;

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
 
 // priors
  p ~ beta(1,1);
  phi ~ beta(1,1);
  pi_r ~ beta(1,1);
  
  // likelihood...
  // ... for indiv recaptured at least once
  for (t in 1:(T - 1)) {
    marr[t] ~ multinomial(pr[t]);
    target += N_1[t] * log(pi_r[t]); 
  }
  // ... for indiv never recaptured
  for (t in 1:(T-1)) {
    target += N_0[t] * log( pr[t][T] * pi_r[t] + (1 - pi_r[t]) );
  }
}
