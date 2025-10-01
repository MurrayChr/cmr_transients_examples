## R files for analyses 

- `00_function_get_marray.R` contains the function to convert capture histories into a single- or multi-state m-array

- `04_compare_efficiency_epm_vs_hhmm.R` compares the computational efficiency of transient mixture models fitted with extended product-multinomial (EPM) or hiearachical hidden Markov model (HHMM) likelihoods (using functions in `04_data_sim_and_prep_functions.R`).

- Files starting with the prefices 01, 02, 03 all relate to a common *base* model (without transients), which is extended to include transients:

 | File prefix | Base model description |
 | :---: | :----: |
 | 01 | Single-state model |
 | 02 | Two-state model for two sites |
 | 03 | Four-state model for two sites incorporating trap-dependence |

- For each of the base models there are four files, 
    - the first (e.g. 01_*) simulates data from the base model and fits the model to that data using a product-multinomial likelihood
    - the second (e.g. 01a_*) extends the base model to include transients, simulates data from the transient mixture model, and fits the model to that data using the EPM likelihood
    - the third (e.g. 01b_*) is the same as the previous, but fits the model using the HHMM likelihood
    - the fourth (e.g. 01c_*) fits the model with EPM and HHMM likelihood to the same data generated from the model and compares marginal posterior distributions 
