#' Compare the posterior distributions for CJS-with-transient model fit using
#'  extended-product-multinomial and hierarchical hmm likelihoods
library(tidyverse)
library(cmdstanr)
library(bayesplot)
source("R/00_function_get_marray.R")
source("R/00_function_plot_estimates.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               ---- Simulate a single dataset ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set parameter values 
T <- 5                       # number of years
n <- 50                      # number of newly marked birds each year
N <- n*(T-1)                  # total number of individuals in all years
phi <- rbeta(T-1, 7, 3)     # survival
p <- rbeta(T, 6, 4)         # detection
pi_r <- rbeta(T-1, 8, 2)    # proportion 'resident'

# manually set first phi and pi_r to examine identifiability
# phi[1] <- 0.5
# pi_r[1] <- 0.5

# simulate survival process 
# (we enforce transience in the detection process below)
z <- matrix(NA, nrow=N, ncol=T)
fc <- rep(1:(T-1), each=n) # length(fc) should be N
for (i in 1:N) {
  z[i,fc[i]] <- 1
  for (t in fc[i]:(T-1)) {
    z[i,t+1] <- rbinom( 1, 1, z[i,t]*phi[t] )
  }
}

# simulate detection process
y <- matrix(0, nrow=N, ncol=T)       
n_res <- round( n*pi_r )           # number 'resident' in each cohort 
n_tra <- n - n_res                 # number 'transient' in each cohort
res <- c()
for (i in 1:(T-1)) {
  res <- c(res, rep(TRUE, n_res[i]), rep(FALSE, n_tra[i]))   # !beware! running twice ! 
}
for(i in 1:N) {
  y[i,fc[i]] <- 1
  # draw subsequent detections for residents only
  if (res[i]) { 
    for (t in (fc[i]+1):T) {
      y[i,t] <- rbinom( 1, 1, z[i,t]*p[t] )
    }
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#          ---- Fit with extended-product-multinomial likelihood ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# classify individuals that were recaptured vs never recaptured
fc <- apply(y, 1, function(x) min(which(x!=0)))  # first capture
lc <- apply(y, 1, function(x) max(which(x!=0)))  # last capture
recapt <- lc > fc 

# create m-array for those individuals that were recaptured at least once
y_recapt <- y[recapt, ]
fc_recapt <- fc[recapt]
y_recapt_recoded <- y_recapt 
for (i in 1:nrow(y_recapt)) {
  y_recapt_recoded[i, fc_recapt[i]:T] <- 
    ifelse( y_recapt[i, fc_recapt[i]:T]==0, 2, y_recapt[i, fc_recapt[i]:T] )
}
marr_recapt <- get_marray(y_recapt_recoded, nStates = 1)

# count the number of individuals that were ever/never recaptured per cohort
n_never_recapt <- rep(NA, T-1)
n_recapt <- rep(NA, T-1)
for (t in 1:(T-1)) {
  n_never_recapt[t] <- sum( (!recapt)&(fc==t) )
  n_recapt[t] <- sum( (recapt)&(fc==t) )
}

# fit the model
file_epm <- "stan/01a_cjs_transients_epm.stan"
mod_epm <- cmdstan_model(file_epm)
stan_data_epm <- list( T=T, marr=marr_recapt, N_1=n_recapt, N_0=n_never_recapt )
fit_epm <- mod_epm$sample( stan_data_epm, parallel_chains = 4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#             ---- Fit with hierarchical HMM likelihood ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# recode '0' after first capture as '2'
fc <- apply(y, 1, function(x) min(which(x != 0)))
N <- nrow(y)
T <- ncol(y)
for (n in 1:N) {
  y[n, fc[n]:T] <- ifelse(y[n, fc[n]:T]==0, 2, y[n, fc[n]:T])
}

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

# fit the model
file_hhmm <- "stan/01b_cjs_transients_hhmm.stan"
mod_hhmm <- cmdstan_model(file_hhmm) 
stan_data_hhmm <- list(T = T, N = nrow(y_rr), y = y_rr, fc = fc_rr, mult = mult_rr)
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
post_epm <- prepare_draws(fit_epm, common_params, "EPM")
post_hhmm <- prepare_draws(fit_hhmm, common_params, "HHMM")
post <- rbind(post_epm, post_hhmm) 

# plot comparing marginal posterior densities
post %>%
  mutate(
    var = str_extract(name, "[a-z_]+(?=\\[)"),
    time = str_extract(name, "(?<=\\[)[1-9]+(?=\\])")
  ) %>%
  ggplot(aes(x=value, colour=lk)) +
  stat_density(geom="line", position="identity", bounds = c(0,1)) +
  coord_cartesian(xlim=c(0,1)) +
  scale_x_continuous(breaks = seq(0, 1, length.out = 3)) +
  theme_classic() +
  facet_grid(var ~ time) +
  labs(
    colour="Likelihood",
    x = "parameter value",
    y = "posterior density"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    aspect.ratio = 1.5
  )
ggsave("figs/01c_compare_marginal_posteriors.png", scale = 1.5)
  

