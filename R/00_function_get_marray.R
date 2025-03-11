# function to obtain m-array from capture histories
# can be multistate m-array (nStates > 1) or regular cjs m-array (nStates = 1)
# assumes the capture histories are encoded as follows
# 0 = not yet captured (before first capture)
# 1, ..., nStates = code the nStates states that an indiv can be captured in
# nStates + 1 = not observed (post first capture)
library( tidyverse )

get_marray <- function(y, nStates) {
  nYears <- ncol(y)
  
  # check that the capture history matrix has all the values one would expect
  # for nStates states
  expected_values <- 0:(nStates + 1)
  data_values <- sort(unique(as.vector(y)))
  if (!all(expected_values %in% data_values)) {
    warning <- str_c( 
      "Warning! For nStates = ", nStates, " we expect to see ", 
      str_flatten( expected_values, ", "), 
      " in the capture history matrix, but only ", 
      str_flatten( data_values, ", "), " are present."
    )
    print(warning)
  }

  # initialise m-array
  marr <- matrix(0, nrow = nStates*(nYears-1), ncol = nStates*(nYears-1) + 1)
  rownames(marr) <- str_c("rel", rep(1:(nYears-1), each = nStates),".",rep(1:nStates, nYears - 1))
  colnames(marr) <- c( str_c("rec", rep(2:nYears, each = nStates),".",rep(1:nStates, nYears - 1)) , "never")
  
  # function to increment m-array with contributions from a single capture history
  increment_marray <- function(marr, h, nStates) {
    nYears <- length(h)
    captures <- which(h %in% 1:nStates)
    n_captures <- length(captures)
    segments <- list()
    for (k in  1:n_captures) {
      if (k < n_captures) {
        segments[[k]] = c(captures[k],captures[k+1])
      }
      if (k == n_captures) {
        if (captures[k] < length(h)) {
          segments[[k]] = c(captures[k],NA)
        }
      }
    }
    for (s in segments) {
      # marray row index
      block_i <- s[1]
      state_i <- h[s[1]]
      marr_i <- (block_i-1)*nStates + state_i
      # marray column index
      if ( is.na(s[2]) ) {
        marr_j <- nStates*(nYears - 1) + 1
      } else {
        block_j <- s[2] - 1
        state_j <- h[s[2]]
        marr_j <- (block_j-1)*nStates + state_j
      }
      marr[ marr_i, marr_j ] <- marr[ marr_i, marr_j ] + 1
    }
    marr
  }
  
  # add each capture history's contribution to m-array
  for (i in 1:nrow(y)) {
    marr <- increment_marray(marr,y[i,],nStates = nStates)
  }
  
  # return
  marr
} 
