#' Reparametrization of gamma mean and sd
#' 
#' 
#' Can be used after specifying a gamma prior in terms of
#' mean and standard deviation (sd) to calculate the `shape` 
#' and `rate` as defined in \code{\link[stats]{qgamma}}.
#' 
#' @param mean mean of the gamma distribution
#' @param sd standard deviation of the gamma distribution
#' 
#' @return a vector with the `shape` 
#' and `rate` of the gamma distribution
#' 
#' @details
#' Mean \eqn{E(X)} and sd \eqn{sd(X)} of a gamma distributed random variable \eqn{X} 
#' parametrized in terms of the shape parameter \eqn{\alpha}
#' and rate parameter \eqn{\beta} are given by
#' 
#' \deqn{E(X) = \frac{\alpha}{\beta}}
#' \deqn{sd(X) = \frac{\alpha}{\beta^2}}
#' 
#' For \eqn{E(X), \;sd(X)>0} we can rearrange above equations to derive the shape and rate parameters:
#' 
#' \deqn{\alpha = \frac{\mu^2}{\sigma^2}}
#' 
#' \deqn{\beta = \frac{\mu}{\sigma^2}.}
#' 
#' 
#' This function can be used to reparametrize the gamma prior
#' mean and sd to the parameters used in \code{\link[stats]{qgamma}}, for instance to 
#' calculate a ROPE \insertCite{kruschke2018}{WSPsignal} based on the prior belief representing the 
#' null hypothesis in BWSP testing \insertCite{dyck2024bpgwsppreprint}{WSPsignal} or visualize the
#' prior distribution.
#' 
#' @examples
#' # obtain shape and rate for gamma distribution with mean = 1 and sd = 10
#' m = 1; s = 10
#' gampars = gamprior_repar(mean = m, sd = s)
#' shape = gampars[1] # shape parameter
#' rate = gampars[2] # rate parameter
#' 
#' # test: sample from the gamma distribution with the derived parameters
#' gamma_sample = rgamma(10000000, shape = shape, rate = rate)
#' m_emp = mean(gamma_sample); m_emp # estimated mean
#' s_emp = sd(gamma_sample); s_emp # estimated sd
#' 
#' # suppose, upper mean and sd reflect prior belief about the Weibull shape parameter and 
#' # calculate an 80% ROPE based on the parameters
#' rope = qgamma(p = c(0.1,0.9), shape = gampars[1], rate = gampars[2])
#' rope
#'
#' @references 
#' \insertAllCited{}
#'
#' @noRd

gamprior_repar = function(mean, sd){
  if(mean <= 0 | sd <= 0) stop("mean and sd must be positive.")
  
  # same reparametrization as in stan
  # given random variable t~Gamma with expectation mean and variance sd^2,
  # shape (also called shape in pgamma fct) and rate (also called rate in pgamma fct)
  # are decucted as
  shape = mean^2 / sd^2
  rate = mean / sd^2
  return(c(shape, rate)) # returned parameterization
  # can be used to calculate ROPEs with qgamma
}
