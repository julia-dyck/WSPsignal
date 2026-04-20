#' Reparametrization of lognormal mean and standard deviation
#' 
#' 
#' Can be used after specifying a lognormal prior in terms of
#' mean and standard deviation (sd) to calculate the location `meanlog` 
#' and scale `sdlog` as defined in \code{\link[stats]{qlnorm}}.
#' 
#' 
#' @param mean mean of the lognormal distribution
#' @param sd sd of the lognormal distribution
#' 
#' @return a vector with the location `meanlog` and scale `sdlog` of the 
#' lognormal distribution
#' 
#' @details
#' 
#' Mean \eqn{E(X)} and sd \eqn{sd(X)} of a lognormal random variable \eqn{X} 
#' parametrized in terms of the location `meanlog` \eqn{\mu} and scale `sdlog`
#' \eqn{\sigma} are given by
#' 
#' \deqn{E(X) = \exp\left( \mu + \frac{1}{2\sigma^2}  \right)}
#' \deqn{sd(X) = \sqrt{\left(\exp\left( 2\mu + \sigma^2 \right) \left( \exp(\sigma^2) - 1 \right) \right)}.}
#' 
#' If \eqn{|sd(X)|>|E(X)|>0}, rearranging above equations to
#' 
#' \deqn{\mu = \log(E(X)) - \frac{1}{2} \log\left(\frac{sd(X)^2}{E(X)^2} - 1\right)}
#' 
#' \deqn{\sigma = \sqrt{\log\left(\frac{sd(X)^2}{E(X)^2} + 1\right)}.}
#' 
#' is possible such that the location parameter `meanlog` and the 
#' scale parameter `logmean` of the lognormal distribution are 
#' obtained.
#' 
#' The application purpose is to reparametrize lognormal prior
#' mean and sd to the parameters used in \code{\link[stats]{qlnorm}}, for instance to 
#' calculate a ROPE \insertCite{kruschke2018}{WSPsignal} based on the prior belief representing the 
#' null hypothesis in BWSP testing \insertCite{dyck2024bpgwsppreprint}{WSPsignal}.
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @examples
#' # obtain location and scale for lognormal distribution with mean = 1 and sd = 10
#' m = 1; s = 10
#' logpars = logprior_repar(mean = m, sd = s)
#' mu = logpars[1] # location parameter
#' sigma = logpars[2] # scale parameter
#' 
#' # test: sample from the lognormal distribution with the obtained parameters
#' lnorm_sample = rlnorm(10000000, meanlog = mu, sdlog = sigma)
#' m_emp = mean(lnorm_sample); m_emp # estimated mean
#' s_emp = sd(lnorm_sample); s_emp # estimated sd
#' 
#' # suppose, upper mean and sd reflect prior belief about the Weibull shape parameter and 
#' # calculate an 80% ROPE based on the parameters
#' rope = qlnorm(p = c(0.1,0.9), meanlog = logpars[1], sdlog = logpars[2])
#' rope
#' 
#'@noRd

logprior_repar = function(mean, sd){
  if(mean <= 0 | sd <= 0) stop("mean and sd must be positive.")
  if(sd <= mean) stop("sd must be larger than the mean.")
  
  # same reparametrization as in stan
  # given random variable t~logN with expectation mean and variance sd^2,
  # location mu (or meanlog in plnorm fct) and scale (or sdlog in plnorm fct)
  # are decucted as
  mu = log(mean) - log(sd^2/mean^2 -1)*0.5
  sigma = sqrt(log(sd^2/mean^2 +1))
  return(c(mu, sigma)) # returned parametrization
  # can be used to calculate ROPEs with qlnorm
}
