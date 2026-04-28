#' Template for prior specification in simulation study
#'
#' @description
#' Generates a structured template for specifying prior means and prior standard
#' deviations (sds) for the Weibull, double Weibull, or Power generalized Weibull
#' model parameters to be inserted into \code{\link{sim.setup_sim_pars}}.
#'
#'
#' @param tte.dist character vector specifying one or multiple modelling approaches; options are
#' \code{"w", "dw", "pgw"} (see \code{\link{bwsp_model}})
#' @param prior.sds numeric value setting the same prior sd for all scale
#' and shape parameters across all included model types; default is 10.
#'
#' @details
#' The returned list contains one data frame per time-to-event (tte) distribution
#' (\code{w}, \code{dw}, \code{pgw}). 
#' For each chosen tte distribution, rows corresponding to
#' different levels of prior belief about the the hazard function are provided, namely
#' \code{"none"}, \code{"beginning"}, \code{"middle"}, and \code{"end"}
#' \insertCite{dyck2024bpgwsppreprint}{WSPsignal}.
#' Given the template, prior means must be filled by the user before simulation.
#'
#' @return
#' A named list containing three data frames 
#' \code{$w}, \code{$dw}, and \code{$pgw} specifying prior beliefs and placeholder entries for prior means,
#' and one vector specifying the prior standard deviations to be considered per parameter in each prior belief
#' setting.
#'
#'
#' @examples
#' ####  prior elicitation --------------------------------------------------------
#' # try a few prior parameter combinations and see whether the resulting hazard 
#' # roughly matches the prior belief about the hazard form 
#'
#' # Expected event time can also be taken into account for some guidance, but 
#' # should not be prioritized.
#' # The reason is that we do not expect the model to accurately fit the hazard of the 
#' # data, but only catch the rough form by distinguishing the cases
#' # constant vs decreasing vs unimodal vs increasing hazard.
#'
#' # set prior means for Power generalized Weibull parameters:
#' plot_pgw(scale = 1, shape = 1, powershape = 1)     # under prior belief "none"
#' plot_pgw(scale = 20, shape = 5.5, powershape = 14) # under prior belief "beginning"
#' plot_pgw(scale = 180, shape = 1, powershape = 1)   # under prior belief "middle"
#' plot_pgw(scale = 300, shape = 4, powershape = 1)   # under prior belief "end"
#'
#'
#' #### specify parameter combinations for simulation study -----------------------
#'
#' fp_list = sim.priors_template(tte.dist = c("pgw"),
#'                               prior.sds = 10) # setup prior template
#' # fill in prior template with values chosen in prior elicitation
#' fp_list$pgw$scale.mean_pgw = c(1, 1, 20, 300)     # scale prior means
#' fp_list$pgw$shape.mean_pgw <- c(1, 0.207, 5.5, 4) # shape prior means
#' fp_list$pgw$powershape.mean_pgw = c(1, 1, 14, 1)  # powershape prior means
#'
#' fp_list # filled fitpars.list ready for sim.setup_sim_pars()

#'
#' @seealso
#' \code{\link{sim.setup_sim_pars}}
#'
#' @export



# add commentaries about support of prior means and sds (have a look at 
# gamprior_repar, logprior_repar

sim.priors_template = function(tte.dist = c("w", "dw", "pgw"),
                               prior.sds = 10){
  # argument checks ------------------------------------------------------------
  
  ### checks for tte.dist
  if (any(duplicated(tte.dist))) {
    warning("Duplicate entries removed from tte.dist.\n")
    tte.dist <- unique(tte.dist)
  }
  allowed_dists <- c("w","dw","pgw")
  if (any(is.na(match(tte.dist, allowed_dists))))
    stop(paste0("Argument tte.dist must be out of: ",
                paste(allowed_dists, collapse = ", "),
                ".\n"))
  
  
  ### checks for prior.sds
  if (!is.numeric(prior.sds) || any(is.na(prior.sds)) || any(prior.sds <= 0)) {
    stop("Argument prior.sds must be a numeric vector with all values > 0.\n")
  }
  
  if (any(duplicated(prior.sds))) {
    warning("Duplicate entries removed from prior.sds.\n")
    prior.sds <- unique(prior.sds)
  }
  
  # fct body -------------------------------------------------------------------
  
  pb = c("none", "beginning", "middle", "end")
  placeholder = rep(NA, length(pb))
  
  prior_list = list()
  
  # 1. templates for prior.means
  prior_list$w = data.frame(prior.belief = pb, 
                            scale.mean_w = placeholder, 
                            shape.mean_w = placeholder)
  # if "w" is not in tte.dist, remove rows 1:4 from dataframe not using dplyr
  if(!("w" %in% tte.dist)){
    prior_list$w = prior_list$w[0,]
  }
  
  # dw
  prior_list$dw = data.frame(prior.belief = pb, 
                            scale.mean_dw = placeholder,
                            shape.mean_dw = placeholder,
                            scale_c.mean_dw = placeholder,
                            shape_c.mean_dw = placeholder)
  if(!("dw" %in% tte.dist)){
    prior_list$dw = prior_list$dw[0,]
  }
  
  # pgw
  prior_list$pgw = data.frame(prior.belief = pb, 
                             scale.mean_pgw = placeholder,
                             shape.mean_pgw = placeholder,
                             powershape.mean_pgw = placeholder)
  if(!("pgw" %in% tte.dist)){
    prior_list$pgw = prior_list$pgw[0,]
  }
  
  # 2. provide sd vector in fourth list element
  prior_list$prior.sds = prior.sds
  
  return(prior_list)
}

# sim.priors_template()
# sim.priors_template(prior.sds = 100)
