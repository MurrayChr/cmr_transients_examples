#' Multi-state model without transients
#' Here we simulate data from a multi-site model with time- and site-dependent 
#' survival, site-dependent movement, and time- and site-dependent detection. 
#' We treat the simplest case of two sites, thus there are two states:
#'  1 - alive at first site
#'  2 - alive at second sitetwo-state
#'  In the simulation we also explicitly include the 'dead state.
library(tidyverse)
library(cmdstanr)
library(bayesplot)
source("R/00_function_get_marray.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               ---- Simulate a single dataset ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set parameter values 
T <- 10                     # number of years
n <- 1000                     # number of newly marked birds per site per year
N <- 2*n*(T-1)              # total number of individuals in all years
phi <- matrix(NA, nrow = 2, ncol = T-1)
phi[1,] <- rbeta(T-1, 7, 3)     # survival at site 1
phi[2,] <- rbeta(T-1, 8, 2)     # survival at site 2
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
for (i in 1:N) {
  z[i,fc[i]] <- fc_state[i]
  for (t in fc[i]:(T-1)) {
    z[i,t+1] <- sample(1:3, 1, prob = Gamma[[t]][z[i,t],])
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
#               ---- Fit multistate model ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# construct m-array
marr <- get_marray(y, nStates = 2)

# compile and fit Stan model
file <- "stan/03_ms_pm.stan"
mod <- cmdstan_model(file)
stan_data <- list(T = T, marr = marr)
fit <- mod$sample(stan_data, parallel_chains = 4)

# compare estimates to truth
plot_estimate <- function(fit, par, truth) {
  max_t <- ncol(truth)
  plt <- fit$summary(par) %>%
    mutate(
      site = as.integer( str_extract(variable, "(?<=\\[)[1-9](?=,)")  ),
      t = as.integer( str_extract(variable,"(?<=,)[0-9]+(?=\\])") )
    ) %>%
    add_column(
      truth = c(truth)
    ) %>%
    ggplot(aes(x=t)) +
    geom_pointrange(aes(y = median, ymin = q5, ymax = q95)) +
    geom_point(aes(y=truth), colour = "red") +
    theme_classic() +
    theme(
      panel.grid.major = element_line()
    ) +
    coord_cartesian(ylim = c(0,1)) +
    scale_x_continuous(breaks = 1:T) +
    facet_wrap(vars(site)) +
    labs(title = str_c("Estimates vs. truth for ", par))
  if (par == "m") {
    plt <- plt +
      labs(x = "site to")
  }
  plt
}

plot_estimate(fit, "p", p)
plot_estimate(fit, "phi", phi)
plot_estimate(fit, "m", m)
  

