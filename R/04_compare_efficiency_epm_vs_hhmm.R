#' Compare efficiency between epm and hhmm likelihoods for each of the three
#' models, across a range of study sizes
library(tidyverse)
library(cmdstanr)
source("R/00_function_get_marray.R")
source("R/04_data_sim_and_prep_functions.R")

# Compile Stan models 
mod_01_epm <- cmdstan_model("stan/01a_cjs_transients_epm.stan")
mod_01_hhmm <- cmdstan_model("stan/01b_cjs_transients_hhmm.stan")
mod_02_epm <- cmdstan_model("stan/02a_ms_transients_epm.stan")
mod_02_hhmm <- cmdstan_model("stan/02b_ms_transients_hhmm.stan")
mod_03_epm <- cmdstan_model("stan/03a_ms-plus-td_transients_epm.stan")
mod_03_hhmm <- cmdstan_model("stan/03b_ms-plus-td_transients_hhmm.stan")

# Arrange models, sim_data_* and prep_data_* functions in lists for looping
mods <- list(
  list("epm" = mod_01_epm, "hhmm" = mod_01_hhmm),
  list("epm" = mod_02_epm, "hhmm" = mod_02_hhmm),
  list("epm" = mod_03_epm, "hhmm" = mod_03_hhmm)
)
sim_data_fns <- list( sim_data_01, sim_data_02, sim_data_03 )
prep_data_fns <- list(
  list("epm" = prep_data_01_epm, "hhmm" = prep_data_01_hhmm),
  list("epm" = prep_data_02_epm, "hhmm" = prep_data_02_hhmm),
  list("epm" = prep_data_03_epm, "hhmm" = prep_data_03_hhmm)
)

# Set study size parameters
T_vals <- c(5, 10, 20, 40)
n_vals <- c(50, 100, 200, 400, 800)

# Sim, fit and save (so we don't have to do them all in one go)
for (T in T_vals[1:2]) {
  for (n in n_vals[1:3]) {
    # file_suffix (for saving model fits)
    fs <- str_c(str_pad(T, 2, "left", pad = "0"),"_", 
                str_pad(n, 3, "left", pad = "0"))
    for (m in 1:3) {
      # simulate data
      if (m == 1) {
        y <- sim_data_fns[[m]](T, n)
      }
      if (m %in% c(2,3)) {
        y <- sim_data_fns[[m]](T, n_per_site = as.integer(n/2)) # two sites
      }
      
      # fit with epm and add runtime to dataframe
      stan_data_epm <- prep_data_fns[[m]]$epm(y)
      fit_epm <- mods[[m]]$epm$sample(stan_data_epm, parallel_chains = 4)
      out_file <- str_c("outputs/sim/", m, "_epm_", fs, ".RDS")
      fit_epm$save_object(out_file)

      # fit with hhmm and add runtime to dataframe
      stan_data_hhmm <- prep_data_fns[[m]]$hhmm(y)
      fit_hhmm <- mods[[m]]$hhmm$sample(stan_data_hhmm, parallel_chains = 4)
      out_file <- str_c("outputs/sim/", m, "_hhmm_", fs, ".RDS")
      fit_hhmm$save_object(out_file)
    }
  }
}

# Assemble dataframe of runtimes
rt <- data.frame(
  model = rep(1:3, each = length(T_vals)*length(n_vals)),
  T = rep(rep(T_vals, each = length(n_vals)),3), 
  n = rep(n_vals, 3*length(T_vals)), 
  epm = NA,
  hhmm = NA
)

for (T in T_vals) {
  for (n in n_vals) {
    # file_suffix (for reading in saved model fits)
    fs <- str_c(str_pad(T, 2, "left", pad = "0"),"_", 
                str_pad(n, 3, "left", pad = "0"))
    for (m in 1:3) {
      # read in epm and hhmm model fits
      file_epm <- str_c("outputs/sim/", m, "_epm_", fs, ".RDS")
      fit_epm <- readRDS(file_epm)
      file_hhmm <- str_c("outputs/sim/", m, "_hhmm_", fs, ".RDS")
      fit_hhmm <- readRDS(file_hhmm)

      # add runtimes to dataframe
      filt <- (rt$T == T) & (rt$n == n) & (rt$model == m)
      rt$epm[filt] <- fit_epm$time()$total
      rt$hhmm[filt] <- fit_hhmm$time()$total
    }
  }
}

rt$speedup <- with(rt, hhmm / epm)
