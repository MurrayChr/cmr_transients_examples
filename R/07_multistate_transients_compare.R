#' Compare the posterior distributions for multistate-with-transient model fit 
#' using extended-product-multinomial and hierarchical hmm likelihoods
library(tidyverse)
library(cmdstanr)
library(bayesplot)
source("R/00_function_get_marray.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               ---- Simulate a single dataset ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set parameter values 
T <- 10                   # number of years
n <- 100                    # number of newly marked birds per site per year
N <- 2*n*(T-1)              # total number of individuals in all years
phi <- matrix(NA, nrow = 2, ncol = T-1)
phi[1,] <- rbeta(T-1, 7, 3)     # survival at site 1
phi[2,] <- rbeta(T-1, 8, 2)     # survival at site 2
pi_r <- matrix(NA, nrow = 2, ncol = T-1)
pi_r[1,] <- rbeta(T-1, 1, 1)  # residency probability at site 1
pi_r[2,] <- rbeta(T-1, 8, 2)  # residency probability at site 2
p <- matrix(NA, nrow=2, ncol = T)
p[1,] <- rbeta(T, 6, 4)         # detection at site 1
p[2,] <- rbeta(T, 8, 2)         # detection at site 2
m <- matrix(0, 2, 2)
m[1,] <- c(0.7, 0.3)          # movement prob from site 1 (given survival) 
m[2,] <- c(0.5, 0.5)          # movement prob from site 1 (given survival) 
# check following is TRUE
all(rowSums(m) == c(1,1)) 

# construct state transition matrices
Gamma <- list()
for (t in 1:(T-1)) {
  Gamma[[t]] <- matrix(
    c(phi[1,t]*m[1,1], phi[1,t]*m[1,2], 1-phi[1,t],
      phi[2,t]*m[2,1], phi[2,t]*m[2,2], 1-phi[2,t],
      0,                0,          1),
    nrow = 3, ncol = 3, byrow = TRUE)
}
# construct detection matrices
Omega <- list()
for (t in 1:T) {
  Omega[[t]] <- matrix(
    c(p[1,t],      0, 1-p[1,t],
      0, p[2,t], 1-p[2,t],
      0,      0,       1),
    nrow = 3, ncol = 3, byrow = TRUE)
}

# simulate survival process 
z <- matrix(NA, nrow=N, ncol=T)
fc <- rep(1:(T-1), each=2*n) # length(fc) should be N
fc_state <- rep(rep(c(1,2), each = n), T-1 ) # length(fc_state) should be N
res <- rep(NA, N) # residency indicator
for (i in 1:N) {
  res[i] <- rbinom(1, 1, pi_r[fc_state[i], fc[i]]) # simulate residency status
}
for (i in 1:N) {
  z[i,fc[i]] <- fc_state[i]
  if (res[i]==1) {
    for (t in fc[i]:(T-1)) {
      z[i,t+1] <- sample(1:3, 1, prob = Gamma[[t]][z[i,t],])
    }
  }
  if (res[i]==0) {
    z[i,(fc[i]+1):T] <- 3 # 'dead' if transient
  }
}

# simulate detection process
y <- matrix(0, nrow=N, ncol=T)       
for(i in 1:N) {
  y[i,fc[i]] <- fc_state[i]
  for (t in (fc[i]+1):T) {
    y[i,t] <- sample(1:3, 1, prob = Omega[[t]][z[i,t],])
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       ---- Fit model with extended product-multinomial likelihood ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# classify single- and multi-encounter capture histories
fc <- apply(y, 1, function(x) min(which(x!=0)))  # first capture
lc <- apply(y, 1, function(x) max(which(x!=3)))  # last capture
recapt <- lc > fc  # recaptured at least once

# create m-array for those individuals that were recaptured at least once
y_recapt <- y[recapt, ]
fc_recapt <- fc[recapt]
marr_recapt <- get_marray(y_recapt, nStates = 2)

# count the number of single- and multi-encounter capture histories per cohort
# and site
Nsingle <- matrix(NA, nrow = 2, ncol = T-1)
Nmult <- matrix(NA, nrow = 2, ncol = T-1)
for (t in 1:(T-1)) {
  for (i in 1:2) {
    Nsingle[i, t] <- sum((!recapt)&(fc==t)&(fc_state==i))
    Nmult[i, t] <- sum(recapt&(fc==t)&(fc_state==i))
  }
}

# compile and fit Stan model
file_epm <- "stan/04_ms_transients_epm.stan"
mod_epm <- cmdstan_model(file_epm)
stan_data_epm <- list(T = T, marr = marr_recapt, Nmult = Nmult, Nsingle = Nsingle)
fit_epm <- mod_epm$sample(stan_data_epm, parallel_chains = 4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               ---- Fit model with hierarchical hmm likelihood ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create reduced representation of data with the number of unique capture 
# histories and their multiplicities (how many times each occurs)
colnames(y) <- str_c("yr", 1:ncol(y))
reduced_representation <- as_tibble(cbind(y, fc)) %>%
  group_by_all() %>%
  summarise( n = n() ) %>%
  ungroup()
y_rr <- reduced_representation %>%
  select( -c("fc", "n") ) %>%
  as.matrix()
fc_rr <- reduced_representation$fc
mult_rr <- reduced_representation$n

# state at first capture
N_rr <- nrow(y_rr)
fc_state_rr <- rep(NA, N_rr)
for (n in 1:N_rr) {
  fc_state_rr[n] <- y_rr[n, fc_rr[n]]
}

# compile and fit Stan model
file_hhmm <- "stan/05_ms_transients_hhmm.stan"
mod_hhmm <- cmdstan_model(file_hhmm)
stan_data_hhmm <- list(T=T, N=N_rr, y=y_rr, fc=fc_rr, fc_state=fc_state_rr, mult=mult_rr)
fit_hhmm <- mod_hhmm$sample(stan_data_hhmm, parallel_chains = 4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   ---- Compare posteriors ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get all parameters common to two models (all parameters in this case)
get_common_params <- function(fit_1, fit_2) {
  params_1 <- fit_1$metadata()$model_params
  params_2 <- fit_2$metadata()$model_params
  intersect(params_1, params_2)
}
common_params <- get_common_params(fit_epm, fit_hhmm)
common_params <- common_params[common_params != "lp__"]

# format posterior draws and add model name
prepare_draws <- function(fit, params, name) {
  fit$draws(params, format = "df") %>%
    add_column(lk = name) %>%
    select(-starts_with(".")) %>%
    pivot_longer(-"lk")
}
post_epm <- prepare_draws(fit_epm, common_params, "epm")
post_hhmm <- prepare_draws(fit_hhmm, common_params, "hhmm")
post <- rbind(post_epm, post_hhmm) 

# plot comparing marginal posterior densities
post %>%
  ggplot(aes(x=value, colour=lk)) +
  geom_density(bounds=c(0,1)) +
  coord_cartesian(xlim=c(0,1)) +
  theme_light() +
  facet_wrap(vars(name), scales = "free")

# and lastly compare runtimes
max(fit_hhmm$metadata()$time$total) / max(fit_epm$metadata()$time$total)
