**Stan model code**

**"01_cjs_pm.stan"** fits the CJS model (without transients) using the 
product-multinomial likelihood.

**"01a_cjs_transients_epm.stan"** fits the CJS-with-transients model using the 
extended product-multinomial likelihood.

**"01b_cjs_transients_hhmm.stan"** fits the CJS-with-transients model using the 
'multievent' or hierarchical hidden markov model likelihood.

**"02_ms_pm.stan"** fits a multi-state model (with two states, without transients) using the 
product-multinomial likelihood.

**"02a_ms_transients_epm.stan"** fits a multi-state-with-transients model 
using the extended product-multinomial likelihood. 

**"02b_ms_transients_hhmm.stan"** fits a multi-state-with-transients model 
using the 'multievent' or hierarchical hidden markov model likelihood.
