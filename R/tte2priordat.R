#' Prior and data setup for Bayesian survival model fitting
#'
#' @description
#' Prepare time-to-event data and prior specifications for model fitting
#' with \code{\link{bwsp_model}}. 
#'
#' @param dat matrix or data frame with time in the first column and event status 
#' in the second column
#' @param tte.dist character indicating the modelling approach; options are \code{"w"}, 
#' \code{"dw"}, \code{"pgw"}; default is \code{"dw"}
#' @param prior.dist character indicating the prior distribution for the
#' parameters of the tte distribution; options are 
#' \code{"fg", "fl", "gg", "ll"}; default is \code{"ll"}
#' @param scale.mean prior mean of the scale parameter; default is 1
#' @param scale.sd prior standard deviation (sd) of the scale parameter; default is 10
#' @param shape.mean prior mean of the shape parameter; default is 1
#' @param shape.sd prior sd of the shape parameter; default is 10
#' @param scale_c.mean prior mean of the scale parameter for censored-at-half data 
#' (only for \code{tte.dist="dw"}); default is 1
#' @param scale_c.sd prior sd of the scale parameter for censored-at-half data 
#' (only for \code{tte.dist="dw"}); default is 10
#' @param shape_c.mean prior mean of the shape parameter for censored-at-half data 
#' (only for \code{tte.dist="dw"}); default is 1
#' @param shape_c.sd prior sd of the shape parameter for censored-at-half data 
#' (only for \code{tte.dist="dw"}); default is 10
#' @param powershape.mean prior mean of the power shape parameter (only for 
#' \code{tte.dist="pgw"})
#' @param powershape.sd prior sd of the power shape parameter (only for 
#' \code{tte.dist="pgw"})
#' 
#'
#' @details
#' The function prepares data to fit a Bayesian model to time-to-event data.
#' The distribution can be assumed a Weibull \code{("w")},
#' a double Weibull \code{("dw"}, estimating two Weibull models - one to the data as is and 
#' one to the data censored at mid of observation period), 
#' or a power generalized Weibull \code{("pgw")} model.
#' 
#' Only the parameters relevant to the chosen \code{tte.dist} need to be provided, that is:
#' \itemize{
#'   \item for \code{"w"}: \code{scale.mean}, \code{scale.sd}, \code{shape.mean}, \code{shape.sd}
#'   \item for \code{"dw"}: \code{scale.mean}, \code{scale.sd}, \code{shape.mean}, \code{shape.sd},
#'   \code{scale_c.mean}, \code{scale_c.sd}, \code{shape_c.mean}, \code{shape_c.sd}
#'   \item for \code{"pgw"}: \code{scale.mean}, \code{scale.sd}, \code{shape.mean}, \code{shape.sd}, 
#'   \code{powershape.mean}, \code{powershape.sd}
#' }
#' 
#' 
#' Implemented prior distributions for the scale and shape parameters are products 
#' of the following univariate distributional choices:
#' \tabular{lll}{
#' for scale parameter \tab for shape parameter(s) \tab  abbreviation \cr
#' fixed to prior mean \tab gamma  \tab fg \cr
#' gamma \tab gamma \tab gg \cr
#' fixed to prior mean \tab lognormal \tab fl \cr
#' lognormal \tab lognormal \tab ll \cr
#' }
#' 
#' 
#' Prior means suitable to reflect the prior belief can be worked out by plotting the
#' hazard and estimating the expected event time under different parameter combinations
#' using \code{\link{plot_pgw}}(\code{powershape = 1} reduces the power 
#' generalized Weibull distribution to Weibull) or \url{https://janoleko.shinyapps.io/pgwd/}.
#' 
#' Prior standard deviations should reflect the uncertainty about the prior belief
#' (i.e. set smaller standard deviation in case of high certainty about prior belief vs. larger 
#' standard deviation in case of low certainty).
#' 
#'
#' @return A named list in the format expected by \code{\link{bwsp_model}}.
#'
#' @examples
#' 
#' tte2priordat(dat = tte, tte.dist = "w", prior.dist = "ll", 
#'              scale.mean = 10, scale.sd = 2, 
#'              shape.mean = 1.5, shape.sd = 15)
#'              
#' tte2priordat(dat = tte, tte.dist = "dw", prior.dist = "ll",
#'              scale.mean = 10, scale.sd = 2, 
#'              shape.mean = 1.5, shape.sd = 15,
#'              scale_c.mean = 5, scale_c.sd = 1, 
#'              shape_c.mean = 1, shape_c.sd = 10)
#'              
#' tte2priordat(dat = tte, tte.dist = "pgw", prior.dist = "ll",
#'              scale.mean = 10, scale.sd = 2, 
#'              shape.mean = 1.5, shape.sd = 15,
#'              powershape.mean = 3, powershape.sd = 20)
#' 
#'
#' @export
#' 

tte2priordat = function(dat,
                        tte.dist = "dw",
                        prior.dist = "ll",
                        scale.mean = 1, scale.sd = 10,
                        shape.mean = 1, shape.sd = 10,
                        scale_c.mean = 1, scale_c.sd = 10,
                        shape_c.mean = 1, shape_c.sd = 10,
                        powershape.mean = 1, powershape.sd = 10) {
  
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
  
  ### checks for prior.dist
  if (any(duplicated(prior.dist))) {
    warning("Duplicate entries removed from prior.dist.\n")
    prior.dist <- unique(prior.dist)
  }
  allowed_priors <- c("fg","fl","gg","ll")
  if (any(is.na(match(prior.dist, allowed_priors))))
    stop(paste0("Argument prior.dist must be out of: ", paste(allowed_priors, collapse = ", "),".\n"))
  
  args_list <- args_list <- mget(names(formals()), environment())  # grabs only arguments

  
  # Helper: check required arguments
  req <- function(arg_names, args_list, dist) {
    miss <- arg_names[sapply(arg_names, function(x) is.null(args_list[[x]]))]
    if (length(miss)) {
      stop(sprintf(
        "For tte.dist='%s' additional arguments must be provided: %s",
        dist, paste(miss, collapse = ", ")
      ))
    }
  }
  
  # # Warn about irrelevant arguments
  # warn_irrelevant <- function(args_list, forbidden_args) {
  #   provided <- forbidden_args[sapply(forbidden_args, function(x) !is.null(args_list[[x]]))]
  #   if (length(provided)) {
  #     warning(sprintf(
  #       "For tte.dist='%s' the following arguments are ignored: %s",
  #       tte.dist,
  #       paste(provided, collapse = ", ")
  #     ))
  #   }
  # }
  # reformat information into standat format
  if (tte.dist == "w") {
    # warn_irrelevant(args_list, c("scale_c.mean","scale_c.sd","shape_c.mean","shape_c.sd",
    #                              "powershape.mean","powershape.sd"))
    req(c("scale.mean","scale.sd","shape.mean","shape.sd"), args_list, tte.dist)
    
    standat = tte2priordat_w(dat, scale.mean, scale.sd, shape.mean, shape.sd)
  }
  
  if (tte.dist == "dw") {
    # warn_irrelevant(args_list, c("powershape.mean","powershape.sd"))
    req(c("scale.mean","scale.sd","shape.mean","shape.sd",
          "scale_c.mean","scale_c.sd","shape_c.mean","shape_c.sd"),
        args_list, tte.dist)
    
    standat = tte2priordat_dw(dat, scale.mean, scale.sd, shape.mean, shape.sd,
                              scale_c.mean, scale_c.sd, shape_c.mean, shape_c.sd)
  }
  
  if (tte.dist == "pgw") {
    # warn_irrelevant(args_list, c("scale_c.mean","scale_c.sd","shape_c.mean","shape_c.sd"))
    req(c("scale.mean","scale.sd","shape.mean","shape.sd",
          "powershape.mean","powershape.sd"),
        args_list, tte.dist)
    
    standat = tte2priordat_pgw(dat, scale.mean, scale.sd, shape.mean, shape.sd,
                               powershape.mean, powershape.sd)
  }
  
  
  out = list(standat_list = standat, args_list = args_list)
  return(out)
}
