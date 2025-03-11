# Function to plot estimates from cmdstanr model fit(s)
library(tidyverse)

#' A function to plot estimates of time-varying parameters from one or more
#' cmdstanr fit objects. Optionally, the function can also plot 'true' values
#' when they are available (i.e. simulations).
plot_estimates <- function(
    fit_list,     # named list of cmdstanr model fits
    par,          # character of parameter to fit
    truth=NULL,   # vector of same length as model fit(s) estimates
    years=2013:2023 # study years
    ) {
  #' Checks:
  #' Every fit contains the parameter.
  #' Truth of same length as model fit(s) if not NULL.
  #' fit_list has names.
  n_mods <- length(fit_list)
  mod_names <- names(fit_list)
  if ( any(mod_names=="") | is.null(mod_names) ) {
    stop("Please provide names for the models in the fit_list.")
  }
  n_par_values <- rep(NA, n_mods)
  for (m in 1:n_mods) {
    n_par_values[m] <- fit_list[[m]]$metadata()$variables %>%
      sapply( function(s) str_split_1(s, "\\[")[1] ) %>%
      str_equal(par) %>%
      sum()
  }
  if ( any(n_par_values==0) ) {
    mod_miss_par <- which(n_par_values==0)
    # stop(str_c("Model(s) ",mod_miss_par," do not contain the variable ", par,"."))
    print(str_c("Model(s) ", str_flatten_comma(mod_miss_par)," do not contain the variable ", par,"."))
  }
  if ( length(unique(n_par_values))!=1 ) {
    print("Warning: fitted models have different number of values for ",par,".",
          "Make VERY sure the years for each value are correct.")
  }
  if ( any(n_par_values>length(years)) ) {
    stop("Too few years specified, please add more.")
  }
  if (!is.null(truth)) {
    if ( any( length(truth)!=n_par_values ) ) {
      cat( 
        "Warning: number of 'true' values, and estimated values for",par,
        "in at least one model, are not equal. 
        Make VERY sure the times for each value are correct."
      )
    }
  }
  # combine estimates if necessary
  summary_list <- list()
  for (m in 1:n_mods) {
    temp <- fit_list[[m]]$summary(par) %>%
      add_column( 
        year=years[1:n_par_values[m]],
        model=names(fit_list)[m]
      )
    if (!is.null(truth)) {
      if ( length(truth)>=n_par_values[m] ) {
        temp <- temp %>%
          add_column(truth=truth[1:n_par_values[m]])
      }
      if (length(truth)<n_par_values[m]) {
        truth_extra <- c(truth, rep(NA, n_par_values[m]-length(truth)))
        temp <- temp %>%
          add_column(truth=truth_extra)
      }
    }
    summary_list[[m]] <- temp
  }
  
  # plot
  xshift_values <- case_when(
    n_mods==1 ~ c(0,NA,NA),
    n_mods==2 ~ 0.1*c(-1,1,NA),
    n_mods==3 ~ 0.1*c(-1,0,1)
  )[1:n_mods]
  
  plot_title <- ifelse( is.null(truth), str_c("Estimates for ", par, "."),
                        str_c("Estimates vs. truth for ", par,".") )
  
  plot <- bind_rows(summary_list) %>%
    mutate( 
      xshift=case_when(
        model==names(fit_list)[1] ~ xshift_values[1],
        model==names(fit_list)[2] ~ xshift_values[2],
        model==names(fit_list)[3] ~ xshift_values[3],)
    ) %>%
    ggplot() +
    geom_pointrange( 
      aes(x=year+xshift,y=median,ymin=q5,ymax=q95,colour=model) ) +
    coord_cartesian( ylim=c(0,1) ) +
    theme_classic() +
    theme(
      panel.grid.major = element_line()
    ) +
    scale_x_continuous(breaks=years) +
    labs( x="year", y="estimate", title=plot_title )
  
  if (!is.null(truth)) {
    plot <- plot + geom_point( aes(x=year,y=truth) )
  }
  
  return(plot)
}









