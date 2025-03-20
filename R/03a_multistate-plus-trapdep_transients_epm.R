#' Multi-state model with trap-dependence and transients
#' Here we simulate data from a multi-site model with time- and site-dependent 
#' survival, site-dependent movement, and time- and site-dependent detection. 
#' We treat the simplest case of two sites, following Pradel and Sanz-Aguilar's
#' approach to modelling trap-dependence. Thus there are four states for residents
#'  1 - alive at first site, trap-aware
#'  2 - alive at first site, trap-unaware
#'  3 - alive at second site, trap-aware
#'  4 - alive at second site, trap-unaware
#'  In the simulation we also explicitly include the 'dead state.
library(tidyverse)
library(cmdstanr)
library(bayesplot)
source("R/00_function_get_marray.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               ---- Simulate a single dataset ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set parameter values 
T <- 20                     # number of years
n <- 1000                    # number of newly marked birds per site per year
N <- 2*n*(T-1)              # total number of individuals in all years
phi <- matrix(NA, nrow = 2, ncol = T-1)
phi[1,] <- rbeta(T-1, 7, 3)     # survival at site 1
phi[2,] <- rbeta(T-1, 8, 2)     # survival at site 2
pi_r <- matrix(NA, nrow = 2, ncol = T-1)
pi_r[1,] <- rbeta(T-1, 6, 4)  # residency probability at site 1
pi_r[2,] <- rbeta(T-1, 8, 2)  # residency probability at site 2
pA <- matrix(NA, nrow=2, ncol = T)
pA[1,] <- rbeta(T, 6, 4)         # detection at site 1 for trap-aware indiv
pA[2,] <- rbeta(T, 8, 2)         # detection at site 2 for trap-aware indiv
pU <- matrix(NA, nrow=2, ncol = T)
pU[1,] <- rbeta(T, 5, 5)         # detection at site 1 for trap-unaware indiv
pU[2,] <- rbeta(T, 6, 4)         # detection at site 2 for trap-unaware indiv 
m <- matrix(0, 2, 2)
m[1,] <- c(0.7, 0.3)          # movement prob from site 1 (given survival) 
m[2,] <- c(0.5, 0.5)          # movement prob from site 1 (given survival) 

# construct state transition matrices
Gamma <- list()
for (t in 1:(T-1)) {
  Gamma[[t]] <- matrix(
    c(phi[1,t]*m[1,1]*pA[1,t+1], phi[1,t]*m[1,1]*(1-pA[1,t+1]), phi[1,t]*m[1,2]*pA[2,t+1], phi[1,t]*m[1,2]*(1-pA[2,t+1]), 1-phi[1,t],
      phi[1,t]*m[1,1]*pU[1,t+1], phi[1,t]*m[1,1]*(1-pU[1,t+1]), phi[1,t]*m[1,2]*pU[2,t+1], phi[1,t]*m[1,2]*(1-pU[2,t+1]), 1-phi[1,t],
      phi[2,t]*m[2,1]*pA[1,t+1], phi[2,t]*m[2,1]*(1-pA[1,t+1]), phi[2,t]*m[2,2]*pA[2,t+1], phi[2,t]*m[2,2]*(1-pA[2,t+1]), 1-phi[2,t],
      phi[2,t]*m[2,1]*pU[1,t+1], phi[2,t]*m[2,1]*(1-pU[1,t+1]), phi[2,t]*m[2,2]*pU[2,t+1], phi[2,t]*m[2,2]*(1-pU[2,t+1]), 1-phi[2,t],
      0,                             0,                        0,                              0,          1),
    nrow = 5, ncol = 5, byrow = TRUE)
}
# check that rows sum to one (should return TRUE)
# lapply(Gamma, function(G) all(near(rowSums(G), rep(1, nrow(G))))) %>% unlist() %>% all()

# construct detection matrices
Omega <- list()
for (t in 1:T) {
  Omega[[t]] <- matrix(
    c(1, 0, 0, 0, 0,
      0, 0, 0, 0, 1,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 1,
      0, 0, 0, 0, 1),
    nrow = 5, ncol = 5, byrow = TRUE)
}

# simulate survival process 
z <- matrix(NA, nrow=N, ncol=T)
fc <- rep(1:(T-1), each=2*n) # length(fc) should be N
fc_state <- rep(rep(c(1,3), each = n), T-1 ) # length(fc_state) should be N
res <- rep(NA, N) # residency indicator
for (i in 1:N) {
  row_i <- case_when(
    fc_state[i] == 1 ~ 1,
    fc_state[i] == 3 ~ 2
  )
  res[i] <- rbinom(1, 1, pi_r[row_i, fc[i]]) # simulate residency status
}
for (i in 1:N) {
  z[i,fc[i]] <- fc_state[i]
  if (res[i]==1){
    for (t in fc[i]:(T-1)) {
      z[i,t+1] <- sample(1:5, 1, prob = Gamma[[t]][z[i,t],])
    }
  }
  if (res[i]==0) {
    z[i,(fc[i]+1):T] <- 5 # 'dead' if transient
  }
}

# simulate detection process
y <- matrix(0, nrow=N, ncol=T)       
for(i in 1:N) {
  y[i,fc[i]] <- fc_state[i]
  for (t in (fc[i]+1):T) {
    y[i,t] <- sample(1:5, 1, prob = Omega[[t]][z[i,t],])
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        ---- Fit multistate-plus-trapdep with transients model ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# classify single- and multi-encounter capture histories
fc <- apply(y, 1, function(x) min(which(x!=0)))  # first capture
lc <- apply(y, 1, function(x) max(which(x!=5)))  # last capture
recapt <- lc > fc  # recaptured at least once

# create m-array for those individuals that were recaptured at least once
y_recapt <- y[recapt, ]
fc_recapt <- fc[recapt]
marr_recapt <- get_marray(y_recapt, nStates = 4)

# count the number of single- and multi-encounter capture histories per cohort
# and site
Nsingle <- matrix(NA, nrow = 2, ncol = T-1)
Nmult <- matrix(NA, nrow = 2, ncol = T-1)
for (t in 1:(T-1)) {
  for (i in 1:2) {
    state_i <- 2*(i-1) + 1
    Nsingle[i, t] <- sum((!recapt)&(fc==t)&(fc_state==state_i))
    Nmult[i, t] <- sum(recapt&(fc==t)&(fc_state==state_i))
  }
}

# compile and fit Stan model
file <- "stan/03a_ms-plus-td_transients_epm.stan"
mod <- cmdstan_model(file)
stan_data <- list(T = T, marr = marr_recapt, Nmult = Nmult, Nsingle = Nsingle)
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

plot_estimate(fit, "pA", pA)
plot_estimate(fit, "pU", pU)
plot_estimate(fit, "phi", phi)
plot_estimate(fit, "pi_r", pi_r)
plot_estimate(fit, "m", m)
