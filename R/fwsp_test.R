#' Frequentist Shape Parameter Test
#' 
#' Frequentist hypothesis test based on the shape parameter(s) of a distribution of the 
#' Weibull family.
#' 
#
#' 
#' @param mod.output model output resulting from \code{\link{fwsp_model}}
#' @param cred.level numeric or vector of confidence levels (i.e. 1 - significance level) 
#' for the test(s) to be performed 
#' 
#' @return binary vector, 0 if \eqn{H_0} is accepted, 1 if \eqn{H_1} is rejected
#' 
#' @details This function tests the null hypothesis that the shape parameter(s) of the
#' Weibull family distribution are equal to one.
#' The distribution specific definition of the null and alternative hypotheses 
#' can be seen in \insertCite{sauzet2022;textual}{WSPsignal}.
#' 
#' For the \code{"w"} and \code{"dw"} case, the model output is a summary of a 
#' `survival::Survreg` outcome which provides
#' \eqn{ln(1 / \nu)} as transform of the shape parameter estimate \eqn{\nu}.
#' The transform \eqn{ln(1 / \nu) = 0} under the null hypothesis \eqn{\nu = 1}. 
#' The shape parameter test is performed on the transform
#' equivalent to performing the test based on the shape parameter itself.
#' 
#' For the \code{"pgw"} case, the shape parameter test is performed on the logarithmized
#' parameter estimates, i.e.
#' \code{fwsp_test} tests the null hypothesis that the logarithm of the shape parameters 
#' of the power generalized Weibull distribution are equal to zero based on the 
#' shape estimates and their estimated standard errors extracted from the estimated
#' Hessian matrix.
#' Issues with standard error calculation from the estimated Hessian matrix
#' may lead to a \code{NA} test results which are then transformed to no signal (0) 
#' following \insertCite{sauzet2022;textual}{WSPsignal}.
#'
#' @references 
#' \insertAllCited{}
#' 
#' @examples 
#' # fit a model
#' mod = fwsp_model(dat = tte)
#' mod
#' # perform the shape parameter test at credibility level 0.95 
#' # or significance level 0.05
#' fwsp_test(mod.output = mod, cred.level = 0.95)
#' 
#' @export
#' 


fwsp_test = function(mod.output, cred.level = 0.96){
  
  # argument check mod.output
  
  # argument check cred.level
  if (!is.numeric(cred.level)) stop("Argument cred.level must be numeric.\n")
  if (any(cred.level < 0 | cred.level > 1)) stop("Argument cred.level must be between 0 and 1.\n")
  if (any(duplicated(cred.level))) {
    warning("Duplicate entries removed from cred.level.\n")
    cred.level <- unique(cred.level)
  }
  
  tte.dist = mod.output$tte.dist
  
  fit = mod.output$fit
  
  if(tte.dist == "w"){
    testres = fwsp_test_w(fit, cred.level)
  }
  
  if(tte.dist == "dw"){
    testres = fwsp_test_dw(fit, cred.level)
  }
  
  if(tte.dist == "pgw"){
    testres = fwsp_test_pgw(fit, cred.level)
  }
  
  # name each vector entry according to the credibility level
  names(testres) = paste0("fwsp_", tte.dist, "_", cred.level)
  return(testres)
}
