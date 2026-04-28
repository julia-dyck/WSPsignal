#' model fitting parameter specification
#' 
#' 
#' @noRd


sim.setup_fit_pars = function(tte.dist = c("w", "dw", "pgw"),
                              prior.belief = c("none", "beginning", "middle", "end"),
                              prior.dist = c("fg", "gg", "fl", "ll"),
                              fit_pars_list = NULL) {
  # arg checks -----------------------------------------------------------------
  if (is.null(fit_pars_list)) {
    message("No 'fit_pars_list' was provided.\nPlease generate a prior parameter template using:\n\n  fit_pars_list <- sim.priors_template(tte.dist = c(",
            paste(tte.dist, collapse = ", "), "))\n\nThen fill in the required prior means and standard deviations.\nOnce complete, re-run sim.setup_simpars() with argument fitpars.list = your_filled_out_list.")
    return(invisible(NULL))
  }
  
  if (is.null(fit_pars_list$prior.sds)) {
    stop("fit_pars_list must contain element 'prior.sds'.")
  }
  
  for (dist_name in tte.dist) {
    df <- fit_pars_list[[dist_name]]
    
    if (!is.data.frame(df)) {
      stop(sprintf("fit_pars_list$%s must be a data frame.", dist_name))
    }
    
    if (!"prior.belief" %in% names(df)) {
      stop(sprintf("fit_pars_list$%s must contain a 'prior.belief' column.", dist_name))
    }
    
    if (!all(df$prior.belief %in% c("none", "beginning", "middle", "end"))) {
      stop(sprintf("Invalid 'prior.belief' values found in fit_pars_list$%s.", dist_name))
    }
    
    for (col in setdiff(names(df), "prior.belief")) {
      if (!is.numeric(df[[col]])) {
        stop(sprintf("All values in column '%s' of fit_pars_list$%s must be numeric.", col, dist_name))
      }
    }
    
    if (anyNA(df)) {
      stop(sprintf("Missing (NA) values detected in fit_pars_list$%s. Please ensure all parameters are filled in before proceeding.", dist_name))
    }
  }
  
  # fct body -------------------------------------------------------------------
  
  prior.sds <- fit_pars_list$prior.sds
  
  dist_pc = expand.grid(tte.dist = tte.dist,
                        prior.dist = prior.dist,
                        prior.belief = prior.belief,
                        prior.sd = prior.sds,
                        KEEP.OUT.ATTRS = FALSE,
                        stringsAsFactors = FALSE)
  
  fit_pc = list()
  
  for (dist in c("w", "dw", "pgw")) {
    if (dist %in% tte.dist) {
      df <- fit_pars_list[[dist]]
      df <- merge(dist_pc[dist_pc$tte.dist == dist, ], df, by = "prior.belief")
      df$tte.dist <- as.character(df$tte.dist)
      df$prior.dist <- as.character(df$prior.dist)
      df$prior.belief <- as.character(df$prior.belief)
      fit_pc[[dist]] <- df
    } else {
      fit_pc[[dist]] <- fit_pars_list[[dist]]
    }
  }
  
  return(fit_pc)
}




## END OF DOC