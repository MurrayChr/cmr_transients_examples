#' Here we bundle code in 01c, 02c and 03c into functions to simulate and 
#' prepare data for epm and hhmm likelihoods, for each of the three models.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#             ---- Model 1: Single-state CJS with transients ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim_data_01 <- function(
    T,  # number of years (or sampling occasions)
    n   # number of newly marked individuals per year
) {
  N <- n*(T-1)                  # total number of individuals in all years
  phi <- rbeta(T-1, 7, 3)     # survival
  p <- rbeta(T, 6, 4)         # detection
  pi_r <- rbeta(T-1, 8, 2)    # proportion 'resident'
  
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
  
  # return capture history matrix
  return(y)
}

prep_data_01_epm <- function(y) {
  T <- ncol(y)
  
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
  
  # return list of data for Stan
  list( T=T, marr=marr_recapt, N_1=n_recapt, N_0=n_never_recapt )
}

prep_data_01_hhmm <- function(y) {
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
  
  # return list of data for Stan
  list(T = T, N = nrow(y_rr), y = y_rr, fc = fc_rr, mult = mult_rr)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         ---- Model 2: multi-state (multi-site) with transients ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim_data_02 <- function(
    T,           # number of years (or sampling occasions)
    n_per_site   # number of newly marked individuals per site per year
) {
  n <- n_per_site             # number of newly marked birds per site per year
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
  m[2,] <- c(0.5, 0.5)          # movement prob from site 2 (given survival) 
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
  
  # return capture history matrix
  return(y)
}

prep_data_02_epm <- function(y) {
  T <- ncol(y)
  
  # classify single- and multi-encounter capture histories
  fc <- apply(y, 1, function(x) min(which(x!=0)))  # first capture
  lc <- apply(y, 1, function(x) max(which(x!=3)))  # last capture
  recapt <- lc > fc  # recaptured at least once
  
  # state at first capture
  fc_state <- rep(NA, nrow(y))
  for (i in 1:nrow(y)) {
    fc_state[i] <- y[i, fc[i]]
  }
  
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
  
  # return data list for Stan
  list(T = T, marr = marr_recapt, Nmult = Nmult, Nsingle = Nsingle)
}

prep_data_02_hhmm <- function(y) {
  T <- ncol(y)
  fc <- apply(y, 1, function(x) min(which(x!=0)))  # first capture
  
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
  
  # return data list for Stan
  list(T=T, N=N_rr, y=y_rr, fc=fc_rr, fc_state=fc_state_rr, mult=mult_rr)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  ---- Model 3: multi-state (multi-site plus trap-dependence) with transients ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim_data_03 <- function(
    T,           # number of years (or sampling occasions)
    n_per_site   # number of newly marked individuals per site per year
) {
  n <- n_per_site
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
  m[2,] <- c(0.5, 0.5)          # movement prob from site 2 (given survival) 
  
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
  
  # return capture history matrix
  return(y)
}

prep_data_03_epm <- function(y) {
  T <- ncol(y)
  
  # classify single- and multi-encounter capture histories
  fc <- apply(y, 1, function(x) min(which(x!=0)))  # first capture
  lc <- apply(y, 1, function(x) max(which(x!=5)))  # last capture
  recapt <- lc > fc  # recaptured at least once
  
  # state at first capture
  fc_state <- rep(NA, nrow(y))
  for (i in 1:nrow(y)) {
    fc_state[i] <- y[i, fc[i]]
  }
  
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
  
  # return data list for Stan
  list(T = T, marr = marr_recapt, Nmult = Nmult, Nsingle = Nsingle)
}

prep_data_03_hhmm <- function(y) {
  T <- ncol(y)
  fc <- apply(y, 1, function(x) min(which(x!=0)))  # first capture
  
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
  
  # return data list for Stan
  list(T=T, N=N_rr, y=y_rr, fc=fc_rr, fc_state=fc_state_rr, mult=mult_rr)
}



