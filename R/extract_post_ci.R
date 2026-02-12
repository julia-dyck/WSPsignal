#' Extract credibility interval from Bayesian model output
#' 
#' Inner function to extract credibility interval from output returned by \link{\code{bwsp_model}}.
#'
#' @param mod.output A list returned by \code{bwsp_model()}.
#'   Must contain elements:
#'   \itemize{
#'     \item \code{$fit}: stan fit object (structure depends on \code{tte.dist})
#'     \item \code{$args_list}: list containing prior and distributional information
#'   }
#' @param ci.type character specifying interval type
#' @param cred.level numeric in (0,1). Credibility mass of the interval
#'
#' @details
#' Extracted parameters depend on \code{tte.dist}:
#'
#' \itemize{
#'   \item \code{"w"}: shape parameter \code{nu}
#'   \item \code{"dw"}: shape parameters \code{nu} (uncensored) and
#'         \code{nu} (censored-at-half)
#'   \item \code{"pgw"}: shape parameter \code{nu} and power-shape
#'         parameter \code{gamma}
#' }
#'
#' For \code{"HDI"}, intervals are computed using
#' \code{HDInterval::hdi()}.
#' For \code{"ETI"}, intervals are computed via
#' \code{stats::quantile()} using equal-tailed probabilities.
#'
#' @return
#' Numeric vector of interval bounds ordered as:
#'
#' \itemize{
#'   \item \code{"w"}: (lower, upper)
#'   \item \code{"dw"}: (shape_lower, shape_upper,
#'                      shape_c_lower, shape_c_upper)
#'   \item \code{"pgw"}: (shape_lower, shape_upper,
#'                       powershape_lower, powershape_upper)
#' }
#'
#'



extract_post_ci = function(mod.output,
                           ci.type = c("ETI", "HDI"), 
                           cred.level){
  
  # percentage for equal-tailed confidence interval
  ci_boundaries = c((1 - cred.level)/2, 1 - ((1 - cred.level)/2)) 
  
  # extract tte.dist from mod.output
  tte.dist = mod.output$args_list$tte.dist
  
  # ----------------------------------------------------------------------------
  
  # calculate ci depending on tte.dist and ci.type
  if(tte.dist == "w"){
    stan.output = mod.output$fit
    post.sample = rstan::extract(stan.output, pars = c("nu"))
    if(ci.type == "HDI"){
      ci = HDInterval::hdi(object = post.sample$nu, credMass = cred.level) 
    }
    if(ci.type == "ETI"){
      ci = stats::quantile(post.sample$nu, probs = ci_boundaries)
    }
  }
  
  if(tte.dist == "dw"){
    stan.output = mod.output$fit$uncens
    stan.output.c = mod.output$fit$cens
    post.sample = rstan::extract(stan.output, pars = c("nu"))
    post.sample.c = rstan::extract(stan.output.c, pars = c("nu"))
    if(ci.type == "HDI"){
      ci = c(HDInterval::hdi(object = post.sample$nu, credMass = cred.level),
             HDInterval::hdi(object = post.sample.c$nu, credMass = cred.level))
    }
    if(ci.type == "ETI"){
      ci = c(stats::quantile(post.sample$nu, probs = ci_boundaries),
             stats::quantile(post.sample.c$nu, probs = ci_boundaries))
    }
  }
  
  if(tte.dist == "pgw"){
    stan.output = mod.output$fit
    post.sample = rstan::extract(stan.output, pars = c("nu", "gamma"))
    if(ci.type == "HDI"){
      ci = c(HDInterval::hdi(object = post.sample$nu, credMass = cred.level),
             HDInterval::hdi(object = post.sample$gamma, credMass = cred.level))
    }
    if(ci.type == "ETI"){
      ci = c(stats::quantile(post.sample$nu, probs = ci_boundaries),
             stats::quantile(post.sample$gamma, probs = ci_boundaries))
    }
  }
  # add entry names
  return(ci)
}