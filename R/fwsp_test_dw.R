#' Frequentist double Weibull Shape Parameter Test
#' 
#' Frequentist hypothesis test based on the shape parameter of the Weibull distribution
#' fitted to the data as is and to the data censored at the middle of the observation period.
#' 
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


fwsp_test_dw = function(mod.output, cred.level = 1 - c(1:10/1000, 2:10/100)){ 
  
  alphas = 1 - cred.level
  
  # check whether mod.output is a list with names cens and uncens
  if(!is.list(mod.output) || 
     !all(c("uncens", "cens") %in% names(mod.output)) ||
     !inherits(mod.output$uncens, "summary.survreg") ||
     !inherits(mod.output$cens, "summary.survreg")){
    stop("Argument mod.output must be a list with elements 'uncens' and 'cens' containing \n summary.survreg objects returned by fwsp_model(..., tte.dist = 'dw').")
  }
  
  test.w.inside = function(alpha, survreg.obj){
    rej.w = as.numeric(try(isTRUE(survreg.obj[[9]][2,4] < alpha)))
    return(rej.w)
  }
  
  # single parameter test -> intermediate results
  rej.uncens = sapply(alphas, test.w.inside, survreg.obj = mod.output$uncens)
  rej.cens = sapply(alphas, test.w.inside, survreg.obj = mod.output$cens)
  
  # combine
  # rule: reject H0, if at least one shape differs significantly from null value
  rej.H0 = as.numeric(rej.uncens == 1 | rej.cens == 1)
  
  return(rej.H0) # here: 1 if H0 rejected, 0 if H0 not rejected 
}


