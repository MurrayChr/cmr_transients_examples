**R files for analyses**

**01_single_state_transients_epm.R** simulates and fits the CJS-plus-transients
model using the extended product-multinomial ('epm') likelihood, and compares estimates
to data-generating values. 

**02_single-state_transients_hhmm.R** simulates and fits the CJS-plus-transients
model using the hierarchical hidden markov model ('hhmm') likelihood, and compares 
estimates to data-generating values. 

**03_single-state_transients_compare.R** simulates and fits the CJS-plus-transients
model using both likelihoods, and compares the marginal posterior distributions. 