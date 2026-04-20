#' Construct ROPE boundaries under the null hypothesis
#'
#' @description
#' Internal helper to compute the Region Of Practical Equivalence (ROPE)
#' for shape parameters under the null hypothesis of constant hazard
#' (shape parameters equal to 1) for BWSP testing.
#'
#' @param tte.dist Character string specifying the time-to-event model:
#'   \code{"w"}, \code{"dw"}, or \code{"pgw"}.
#' @param prior.dist Prior distribution for shape parameters:
#'   \code{"ll"} (lognormal) or \code{"gg"} (gamma).
#' @param cred.level Numeric in (0,1). Credibility mass defining the ROPE.
#' @param shape.sd Prior standard deviation of the main shape parameter.
#' @param shape_c.sd Prior standard deviation of the censored-at-half shape
#'   parameter (only relevant for \code{"dw"}).
#' @param powershape.sd Prior standard deviation of the power-shape
#'   parameter (only relevant for \code{"pgw"}).
#'
#' @details
#' The ROPE is defined as the equal-tailed interval of the prior distribution
#' centered at mean 1 (null hypothesis of constant hazard).
#' Depending on \code{tte.dist}, one or two parameter-specific ROPEs are returned.
#'
#' @return
#' Named numeric vector containing lower and upper ROPE bounds.
#' Order depends on \code{tte.dist}:
#' \itemize{
#'   \item \code{"w"}: shape
#'   \item \code{"dw"}: shape and shape_c
#'   \item \code{"pgw"}: shape and powershape
#' }
#'
#' @noRd


setup_rope = function(tte.dist = "dw", 
                      prior.dist = "ll",
                      cred.level = 0.8,
                      shape.sd = 10, 
                      shape_c.sd = 10, 
                      powershape.sd = 10){
  
  # shape means under null hypothesis of constant hazard
  shape.mean = 1
  shape_c.mean = 1
  powershape.mean = 1
  
  ci_boundaries = c((1 - cred.level)/2, 1 - ((1 - cred.level)/2)) # percentage for equal-tailed confidence interval
  
  # calc rope
  if(prior.dist == "ll"){
    distpars_shape = logprior_repar(shape.mean, shape.sd) # from mean, sd to lognormal parameters
    distpars_shape_c = logprior_repar(shape_c.mean, shape_c.sd) 
    distpars_powershape = logprior_repar(powershape.mean, powershape.sd) 
    
    rope_shape = qlnorm(p = ci_boundaries, meanlog = distpars_shape[1], sdlog = distpars_shape[2])
    rope_shape_c = qlnorm(p = ci_boundaries, meanlog = distpars_shape_c[1], sdlog = distpars_shape_c[2])
    rope_powershape = qlnorm(p = ci_boundaries, meanlog = distpars_powershape[1], sdlog = distpars_powershape[2]) 
  }
  if(prior.dist == "gg"){
    distpars_shape = gamprior_repar(shape.mean, shape.sd) # from mean, sd to lognormal parameters
    distpars_shape_c = gamprior_repar(shape_c.mean, shape_c.sd) 
    distpars_powershape = gamprior_repar(powershape.mean, powershape.sd) 
    
    rope_shape = qgamma(p = ci_boundaries, shape = distpars_shape[1], rate = distpars_shape[2])
    rope_shape_c = qgamma(p = ci_boundaries, shape = distpars_shape_c[1], rate = distpars_shape_c[2])
    rope_powershape = qgamma(p = ci_boundaries, shape = distpars_powershape[1], rate = distpars_powershape[2]) 
  }
  
  # gather relevant ropes as vector
  if (tte.dist == "w") {
    rope_vect <- stats::setNames(rope_shape, c("shape_lower", "shape_upper"))
  }
  if (tte.dist == "dw") {
    rope_vect <- stats::setNames(c(rope_shape, rope_shape_c),c("shape_lower", "shape_upper",
                                                        "shape_c_lower", "shape_c_upper")
    )
  }
  if (tte.dist == "pgw") {
    rope_vect <- stats::setNames(c(rope_shape, rope_powershape),c("shape_lower", "shape_upper",
                                                           "powershape_lower", "powershape_upper")
    )
  }
  
  return(rope_vect)
  
}
