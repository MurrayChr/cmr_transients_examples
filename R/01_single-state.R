#' Time-varying CJS model without transients
#' Here we simulate data from a Cormack-Jolly-Seber model with time-varying 
#' survival & detection
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
n <- 100                     # number of newly marked birds each year
N <- n*(T-1)                  # total number of individuals in all years
phi <- rbeta(T-1, 7, 3)       # survival
p <- rbeta(T, 6, 4)           # detection

# simulate survival process 
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
for(i in 1:N) {
  y[i,fc[i]] <- 1
  for (t in (fc[i]+1):T) {
    y[i,t] <- rbinom( 1, 1, z[i,t]*p[t] )
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   ---- Fit CJS model ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create m-array 
marr <- get_marray(y, nStates = 1)

# fit the model
file <- "stan/01_cjs_pm.stan"
mod <- cmdstan_model(file)
stan_data <- list( T=T, marr=marr )
fit <- mod$sample( stan_data, parallel_chains = 4 )

# compare estimates to truth 
plot_estimates(list("cjs"=fit), "phi", truth=phi, years=1:T )
plot_estimates(list("cjs"=fit), "p", truth=p, years=1:T )


#' The credible intervals for survival in year T-1 and detection in year T
#' are wider than the others, because they are not separately identifiable:
pars <- c(str_c("phi[",T-1,"]"), str_c("p[",T,"]"))
mcmc_scatter(fit$draws(pars), alpha=0.2) +
  coord_cartesian( xlim=c(0,1), ylim=c(0,1) ) +
  geom_point(
    mapping = aes(x=x,y=y), 
    data = tibble(x=phi[T-1], y=p[T]),
    colour="red", size=2.5
  )
