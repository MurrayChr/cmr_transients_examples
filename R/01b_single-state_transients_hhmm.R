#' Time-varying CJS model with time-varying transients
#' Here we simulate data from a Cormack-Jolly-Seber model with time-varying survival & detection,
#' in which a small, time-varying proportion of the newly caught individuals in
#' each occasion are transients: their subsequent probability of being detected is zero.
library(tidyverse)
library(cmdstanr)
library(bayesplot)
source("R/00_function_get_marray.R")
source("R/00_function_plot_estimates.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               ---- Simulate a single dataset ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set parameter values 
T <- 10                       # number of years
n <- 1000                      # number of newly marked birds each year
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
#                   ---- Fit CJS-with-transience model ----
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
file_tr <- "stan/02_cjs_transients_hhmm.stan"
mod_tr <- cmdstan_model(file_tr) 
stan_data_tr <- list(T = T, N = nrow(y_rr), y = y_rr, fc = fc_rr, mult = mult_rr)
fit_tr <- mod_tr$sample(stan_data_tr, parallel_chains = 4)

# compare estimates to truth 
plot_estimates(list("transients"=fit_tr), "phi", truth=phi, years=1:T )
plot_estimates(list("transients"=fit_tr), "p", truth=p, years=1:T )
plot_estimates(list("transients"=fit_tr), "pi_r", truth=pi_r, years=1:T )

#' The credible intervals for survival and residency probability in the first
#' year are much wider than the others, because they are not separately identifiable:
mcmc_scatter(fit_tr$draws(c("phi[1]", "pi_r[1]")), alpha=0.2) +
  coord_cartesian( xlim=c(0,1), ylim=c(0,1) ) +
  geom_point(
    mapping = aes(x=x,y=y), 
    data = tibble(x=phi[1], y=pi_r[1]),
    colour="red", size=2.5
  )


