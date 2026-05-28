#' Fit Bayesian double Weibull models to time-to-event data
#'
#' 
#' The function applies the \code{\link[rstan]{sampling}} command to fit two 
#' Weibull models to time-to-event (tte) data with Gamma or Lognormal priors for the
#' parameters of the Weibull distributions - one to the tte data as is, one to the
#' tte data censored at mid of observation period.
#'
#' @param datstan A named list of data for the stanmodel; a data frame can be 
#' reformated with \code{\link{tte2standat}}.
#' @param prior.dist A character string indicating the prior distribution for the
#' parameters of the pgW distribution. Options are 
#' \code{"fg", "fl", "gg", "ll"} (see details).
#' @param chains The number of Markov chains to run
#' @param iter The total number of iterations per chain (including warmup)
#' @param warmup The number of warmup iterations per chain
#' 
#' 
#' @return A list of two stanfit objects; \code{$uncens} contains the stanfit object
#' obtained from the model to the tte data as is, \code{$cens} the stanfit object
#' obtained from the model to the tte data censored at mid of observation period.
#' 
#' 
#' @details
#' The function applies the \code{\link[rstan]{sampling}} command to fit two 
#' Weibull models to time-to-event data with Gamma or Lognormal priors for the
#' parameters of the Weibull distributions - one to the tte data as is, one to the
#' tte data censored at mid of observation period.
#' 
#' The posterior of each Bayesian model is proportional to the likelihood times the prior. 
#' The likelihood is
#' \deqn{\mathcal{L}(t| \Theta) = \prod_{i=1}^N S(t_i)^{1-d_i}\cdot f(t_i)^{d_i}} 
#' with \eqn{S(t)} the survival function of the Weibull distribution and \eqn{f(t)} the
#' density function of the Weibull distribution. The pair \eqn{(t_i, d_i)} are the observed
#' tte observations.
#' 
#' The priors are either independent univariate Gamma or Lognormal distribution
#' for the parameters of the pgW distribution. 
#' Implemented distributional choices for the joint prior are products of the following:
#' \tabular{llll}{
#' for scales  \tab for shapes  \tab abbreviation \cr
#' fixed to prior mean \tab Gamma \tab fg \cr
#' Gamma \tab Gamma \tab ggg \cr
#' fixed to prior mean \tab Lognormal \tab fl \cr
#' Lognormal \tab Lognormal \tab ll \cr
#' }
#' 
#' 
#' @noRd



bwsp_model_dw = function(datstan, 
                     prior.dist = c("fg","fl","gg","ll"),
                     chains = 4,
                     iter = 11000,
                     warmup = 1000
){
  if(prior.dist == "fg"){
    output_uncens = rstan::sampling(
      object = stanmodels$w_tte_gammaprior_scalefixed,  # Stan model
      data = datstan$uncens,     # named list of data
      chains = chains,    # number of Markov chains
      warmup = warmup,    # number of warmup iterations per chain
      iter = iter         # total number of iterations per chain (including warmup)
    )
    output_cens = rstan::sampling(
      object = stanmodels$w_tte_gammaprior_scalefixed,  # Stan model
      data = datstan$cens,     # named list of data
      chains = chains,    # number of Markov chains
      warmup = warmup,    # number of warmup iterations per chain
      iter = iter         # total number of iterations per chain (including warmup)
    )
  }
  if(prior.dist == "fl"){
    output_uncens = rstan::sampling(
      object = stanmodels$w_tte_lognormalprior_scalefixed,  # Stan model
      data = datstan$uncens,     # named list of data
      chains = chains,    # number of Markov chains
      warmup = warmup,    # number of warmup iterations per chain
      iter = iter         # total number of iterations per chain (including warmup)
    )
    output_cens = rstan::sampling(
      object = stanmodels$w_tte_lognormalprior_scalefixed,  # Stan model
      data = datstan$cens,     # named list of data
      chains = chains,    # number of Markov chains
      warmup = warmup,    # number of warmup iterations per chain
      iter = iter         # total number of iterations per chain (including warmup)
    )
  }
  if(prior.dist == "gg"){
    output_uncens = rstan::sampling(
      object = stanmodels$w_tte_gammaprior,  # Stan model
      data = datstan$uncens,     # named list of data
      chains = chains,    # number of Markov chains
      warmup = warmup,    # number of warmup iterations per chain
      iter = iter         # total number of iterations per chain (including warmup)
    )
    output_cens = rstan::sampling(
      object = stanmodels$w_tte_gammaprior,  # Stan model
      data = datstan$cens,     # named list of data
      chains = chains,    # number of Markov chains
      warmup = warmup,    # number of warmup iterations per chain
      iter = iter         # total number of iterations per chain (including warmup)
    )
  }
  if(prior.dist == "ll"){
    output_uncens = rstan::sampling(
      object = stanmodels$w_tte_lognormalprior,  # Stan model
      data = datstan$uncens,     # named list of data
      chains = chains,    # number of Markov chains
      warmup = warmup,    # number of warmup iterations per chain
      iter = iter         # total number of iterations per chain (including warmup)
    )
    output_cens = rstan::sampling(
      object = stanmodels$w_tte_lognormalprior,  # Stan model
      data = datstan$cens,     # named list of data
      chains = chains,    # number of Markov chains
      warmup = warmup,    # number of warmup iterations per chain
      iter = iter         # total number of iterations per chain (including warmup)
    )
  }
  
  
  output = list(uncens = output_uncens, cens = output_cens)
  return(output)
}


