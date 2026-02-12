#' Frequentist Weibull Shape Parameter Test
#' 
#' Frequentist hypothesis test based on the shape parameter of the Weibull distribution.
#' 
#' 
#' @param mod.output estimation output resulting from `fwsp_model(..., tte.dist = "w")`
#' @param cred.level vector of credibility levels for the tests to be performed
#' 
#' @return A vector containing the test results for each credibility level.
#' 
#' 
#' @details The model estimate is a summary of a `survival::Survreg` outcome which provides
#' \eqn{ln(1 / \nu)} as transform of the shape parameter estimate.
#' The transform \eqn{ln(1 / \nu) = 0} under the null hypothesis \eqn{\nu = 1}. 
#' The shape parameter test is performed on the transform leading the test outcome
#' equivalent to performing the test based on the shape parameter itself.
#' 
#' 
#' @noRd

fwsp_test_w = function(mod.output, cred.level = 1 - c(1:10/1000, 2:10/100)){ 
  
  # check whether mod.output is a survreg object
  if(!inherits(mod.output, "summary.survreg")){
    stop("Argument mod.output must be a summary.survreg object returned by fwsp_model(..., tte.dist = 'w').")
  }
  
  alphas = 1 - cred.level
  
  test.w.inside = function(alpha){
    rej.w = as.numeric(try(isTRUE(mod.output[[9]][2,4] < alpha)))
    return(rej.w)
  }
  rej.H0 = sapply(alphas, test.w.inside)
  return(rej.H0)
}

