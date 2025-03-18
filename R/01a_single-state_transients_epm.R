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
file_tr <- "stan/01_cjs_transients_epm.stan"
mod_tr <- cmdstan_model(file_tr)
stan_data_tr <- list( T=T, marr=marr_recapt, N_1=n_recapt, N_0=n_never_recapt )
fit_tr <- mod_tr$sample( stan_data_tr, parallel_chains = 4 )

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


