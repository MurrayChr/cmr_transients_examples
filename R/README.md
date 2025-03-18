**R files for analyses**

**01_single_state.R** simulates and fits a single-state CJS model (without transients) 
using the product-multinomial likelihood, and compares estimates to data-generating values. 

**01a_single_state_transients_epm.R** simulates and fits the CJS-plus-transients
model using the extended product-multinomial ('epm') likelihood, and compares estimates
to data-generating values. 

**01b_single-state_transients_hhmm.R** simulates and fits the CJS-plus-transients
model using the hierarchical hidden markov model ('hhmm') likelihood, and compares 
estimates to data-generating values. 

**01c_single-state_transients_compare.R** simulates and fits the CJS-plus-transients
model using both likelihoods, and compares the marginal posterior distributions. 

**02_multi-state.R** simulates and fits a simple multistate model (without transients) 
using the product-multinomial likelihood, and compares estimates to data-generating 
values.

**02a_multistate_transients_epm.R** simulates and fits a multistate-with-transients 
model using the extended product-multinomial ('epm') likelihood, and compares estimates
to data-generating values. 

**02b_multistate_transients_hhmm.R** simulates and fits the multistate-with-transients
model using the hierarchical hidden markov model ('hhmm') likelihood, and compares 
estimates to data-generating values. 

**02c_multistate_transients_compare.R** simulates and fits the multistate-with-transients
model using both likelihoods, and compares the marginal posterior distributions. 