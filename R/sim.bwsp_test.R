#' Bayesian Weibull Shape Parameter Test in simulation study
#' 
#' @description 
#' Inner fct that performs bwsp_test based on posterior statistics that were 
#' summarized with\link{\code{sim.stanfit.to.poststats}}  (not on mod.output from \link{\code{bwsp_model}})
#' 
#'
#' @param credregion vector of length 2 or 4 with the lower and upper boundaries of the 
#' credibility interval (CI) reflecting the posterior distribution of the shape
#' parameter(s); required order: 1. lower CI boundary of first shape parameter, 
#' 2. upper CI boundary of first shape parameter; and if existent: 3. lower CI 
#' boundary of second shape parameter, 4. upper CI boundary of second shape parameter
#' @param nullregion vector of length 2 or 4 with the lower and upper boundaries of the 
#' region of practical equivalence (ROPE) reflecting the null hypothesis region of the shape
#' parameter(s); required order: 1. lower ROPE boundary of first shape parameter, 
#' 2. upper ROPE boundary of first shape parameter; and if existent: 3. lower ROPE 
#' boundary of second shape parameter, 4. upper ROPE boundary of second shape parameter
#' @param tte.dist character specifying the modelling 
#' approach used to obtain the posterior samples; options are \code{"w", "dw", "pgw"};
#' @param option numeric value out of \code{1,2,3}; rule to be used to deduct a 
#' binary outcome (signal/no signal) from the HDI+ROPE test results of each shape 
#' parameter (see details)
#' 
#' 
#' @return binary, 0 if \eqn{H_0} is accepted, 1 if \eqn{H_1} is rejected; see details for definition
#' of \eqn{H_0} and \eqn{H_1}
#'
#'
#'
#' @section Test concept: 
#' 
#' The Bayesian Weibull shape parameter (WSP) test is a hypothesis test 
#' for signal detection of adverse drug reactions.
#' It is based on the principle of non-constant hazard \insertCite{cornelius2012}{WSPsignal}
#' which associates a constant hazard function with the absence of a drug-event 
#' association and a non-constant hazard function with the presence of a drug-event
#' association.
#' 
#' This can be formalized as the following hypotheses
#' depending on the underlying model:
#' 
#' \tabular{lcc}{
#'         \tab \eqn{H_0} \tab \eqn{H_1} \cr
#'  hypothesis \tab constant hazard function \tab non-constant hazard function \cr
#'  under Weibull model \tab \eqn{\nu = 1} \tab \eqn{\nu \neq 1} \cr
#'  under double Weibull model \tab \eqn{\nu_1 = 1 \text{ and } \nu_2 = 1} \tab \eqn{\nu_1 \neq 1 \text{ or } \nu_2 \neq 1} \cr
#'  under PGW model \tab \eqn{\nu = 1 \text{ and } \gamma = 1} \tab \eqn{\nu \neq 1 \text{ or } \gamma \neq 1} \cr
#' }
#' 
#' 
#' 
#' @section Bayesian test components: 
#' 
#' Information on the Bayesian 
#' variant of the Power Generalized Weibull (PGW) shape 
#' parameter test can be found in \insertCite{dyck2024bpgwsppreprint;textual}{WSPsignal}.
#' The same concept applies to the construction of the Bayesian Weibull and double
#' Weibull shape parameter test.
#' 
#' The region of practical equivalence (ROPE) specified in the 
#' \code{nullregion} argument represents the expected parameter value under \eqn{H_0}.
#' The credibility region(s) specified in the \code{credregion} argument represent
#' the posterior distribution of each shape parameter. 
#' For the ROPE we recommend to set up an equal-tailed interval (ETI) with 
#' \deqn{[q_{(1-\alpha)/2}, q_{(1+\alpha)/2}]}
#' based on the quantiles \eqn{q} of the shape parameters' prior distributions under 
#' \eqn{H_0} at a chosen credibility level \eqn{1 - \alpha}.
#' For the posterior CI, we recommend either an ETI at the same credibility level obtained from
#' the empirical quantiles of the posterior distribution per shape parameter 
#' or a highest density interval (HDI, \insertCite{kruschke2015;textual}{WSPsignal})
#' \deqn{HDI(\nu) = \{\nu \; |\; p_1(\nu) \geq w\} \text{ with } w\in [0,1] 
#' \text{ such that} \int_{\nu \; | \; p_1(\nu) \geq w} p_1(\nu|t)\;  d\nu = 1 - \alpha}
#' at the same credibility level.
#' See Examples for exemplary ROPE and post CI setup in R.
#' Best choices for ROPE and CI can be obtained from a simulation study to tune the 
#' BWSP test specifications.
#' 
#' The HDI+ROPE test checks the 
#' relationship between ROPE and credibility region(s) leading to either acceptance,
#' rejection or no decision regarding the null hypothesis for a single shape parameter.
#' Options to generate a binary outcome, i.e. a signal or not, from HDI+ROPE test results 
#' based on one (in case of \code{"w"}) or two (in case of \code{"dw", "pgw"}) shape parameters are:
#' \tabular{ccccc}{
#' HDI+ROPE \tab HDI+ROPE \tab combination \tab combination \tab combination \cr
#' outcome \tab outcome \tab rule \tab rule \tab rule \cr
#' for shape_1\tab for shape_2 \tab (\code{option = 1}) \tab (\code{option = 2}) \tab (\code{option = 3}) \cr
#' rejection \tab (none) \tab signal \tab signal \tab signal \cr
#' acceptance \tab (none) \tab - \tab - \tab - \cr
#' no decision \tab (none) \tab signal \tab - \tab - \cr
#' rejection \tab rejection \tab signal \tab signal \tab signal \cr
#' acceptance \tab rejection \tab signal \tab - \tab - \cr
#' rejection \tab acceptance \tab signal \tab - \tab - \cr
#' acceptance \tab acceptance \tab - \tab - \tab - \cr
#' no decision \tab rejection \tab signal \tab signal \tab - \cr
#' no decision \tab acceptance \tab - \tab - \tab - \cr
#' rejection \tab no decision \tab signal \tab signal \tab - \cr 
#' acceptance \tab no decision \tab - \tab - \tab - \cr
#' no decision \tab no decision \tab signal \tab - \tab - \cr
#' }
#' 
#' The hypotheses as stated above (see test concept) are implemented in \code{option = 1} whereas
#' \code{option = 2} and \code{option = 3} lead to a signal in fewer cases.
#' 
#' More details on the HDI+ROPE test, recommendations for interval specifications
#' and the combination rules
#' can be found in \insertCite{kruschke2018}{WSPsignal} and
#' \insertCite{dyck2024bpgwsppreprint}{WSPsignal}.
#' 
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @examples
#' #### Exemplary conduction of a test from data and prior to test result:
#' 
#' # under weibull model:
#' 
#' # 1. specify ROPE reflecting the null hypothesis:
#' # we choose an 80% confidence interval around the 
#' # null value (1 for both shape parameters)
#' logpars = logprior_repar(1, 10) # get parameters of a 
#'                                 # lognormal distribution with 
#'                                 # mean 1 and sd 10
#' 
#' rope = qlnorm(p = c(0.1,0.9), meanlog = logpars[1], sdlog = logpars[2])
#' 
#' # 2. Prior specification and model fitting:
#' # we formalize a prior belief (here "no association
#' # between drug and event", therefore prior mean = 1 for shape parameter)
#' # and reformat our tte data to fit the model in the following
#' standat = tte2priordat(dat = tte,   # reformat the data
#'                       tte.dist = "w",
#'                       scale.mean = 1, 
#'                       scale.sd = 10,
#'                       shape.mean = 1, 
#'                       shape.sd = 10)
#' 
#' fit = bwsp_model(datstan = standat,      # fit the model
#'                   tte.dist = "w", 
#'                   prior.dist = "ll",       
#'                   chains = 4,              
#'                   iter = 110,             # (posterior sample
#'                   warmup = 10)            # is small for demo purpose)
#' 
#' # 3. HDI specification and extraction:
#' # extract 80% HDIs representing posterior samples of the shape parameters
#' post.sample = rstan::extract(fit, pars = c("nu"))
#' nu.hdi = HDInterval::hdi(object = post.sample$nu, credMass = 0.8)
#' 
#' # 4. conduct the BWSPtest
#' bwsp_test(credregion = nu.hdi, 
#'           nullregion = rope, 
#'           tte.dist = "w", 
#'           option = 1)
#' # returns a signal
#' 
#' # under pgw model:
#' 
#' # 1. specify ROPE reflecting the null hypothesis:
#' # we choose an 80% confidence interval around the 
#' # null value (1 for both shape parameters)
#' 
#' logpars = logprior_repar(1, 10) # get parameters of a 
#'                                 # lognormal distribution with 
#'                                 # mean 1 and sd 10
#' 
#' rope = qlnorm(p = c(0.1,0.9), meanlog = logpars[1], sdlog = logpars[2])
#' 
#' # 2. Prior specification and model fitting:
#' # we formalize a prior belief (here "no association
#' # between drug and event", therefore prior mean = 1 for both shape parameters)
#' # and reformat our tte data to fit the model in the following
#' standat = tte2priordat(dat = tte,          # reformat the data
#'                       tte.dist = "pgw",
#'                       scale.mean = 1, 
#'                       scale.sd = 10,
#'                       shape.mean = 1, 
#'                       shape.sd = 10,
#'                       powershape.mean = 1, 
#'                       powershape.sd = 10)
#' 
#' fit = bwsp_model(datstan = standat,     # fit the model
#'                   tte.dist = "pgw",
#'                   prior.dist = "ll",       
#'                   chains = 4,              
#'                   iter = 110,            # (posterior sample
#'                   warmup = 10)           # is small for demo purpose)
#' 
#' # 3. HDI specification and extraction:
#' # extract 80% HDIs representing posterior samples of the shape parameters
#' post.samples = rstan::extract(fit, pars = c("nu", "gamma")) 
#' nu.hdi = HDInterval::hdi(object = post.samples$nu, credMass = 0.8)
#' ga.hdi = HDInterval::hdi(object = post.samples$gamma, credMass = 0.8)
#' 
#' # 4. conduct the BWSP test
#' bwsp_test(credregion = c(nu.hdi, ga.hdi), 
#'           nullregion = rope, 
#'           tte.dist = "pgw", 
#'           option = 1)
#' 
#' # returns a signal
#' 
#' @noRd


sim.bwsp_test = function(credregion, 
                     nullregion, 
                     tte.dist = c("w", "dw", "pgw"), 
                     option = c(1,2,3)){
  
  # argument check for mod
  if(tte.dist != "w" & tte.dist != "dw" & tte.dist != "pgw"){
    stop("tte.dist must be one of 'w', 'dw' or 'pgw'")
  }
  # argument check for option 
  if(option != 1 & option != 2 & option != 3){
    stop("option must be one of 1, 2 or 3")
  }
  
  # test under Weibull model
  if(tte.dist == "w"){
    # argument check for nullregion
    if(length(nullregion) != 2){
      stop("Argument nullregion must be a vector of length 2.")
    }
    # argument check for credregion
    if(length(credregion) != 2){
      stop("Argument credregion must be a vector of length 2.")
    }
    shape1credregion = credregion
    shape1ropehdi_res = hdi_plus_rope(nullregion = nullregion, credregion = shape1credregion)
    res = shape1ropehdi_res
    # options for "w" model
    if(option == 1){ # HDI+ROPE test rejects H0 or undecided -> signal
      out = ifelse(is.na(res), 1, res)
      return(out)
    }
    if(option == 2 | option == 3){ # HDI+ROPE test rejects H0 -> signal
      out = ifelse(is.na(res), 0, res)
      return(out)
    }
    else{
      stop("option must be 1, 2 or 3")
    }
  }
  # test under double or pgw model
  if(tte.dist == "dw" | tte.dist == "pgw"){
    
    # argument check for nullregion
    if(length(nullregion) != 2 && length(nullregion) != 4){
      stop("Argument nullregion must be a vector of length 2 or 4.")
    }
    if(length(nullregion) == 2){
      warning("Argument nullregion is a vector of length 2. 
              Assuming that the same null region applies to both shape parameters.")
      nullregion = c(nullregion, nullregion)
    }
    # argument check for credregion
    if(length(credregion) != 4){
      stop("Argument credregion must be a vector of length 4.")
    }
    shape1nullregion = nullregion[1:2]
    shape2nullregion = nullregion[3:4]
    
    shape1credregion = credregion[1:2]
    shape2credregion = credregion[3:4]
    
    shape1ropehdi_res = hdi_plus_rope(nullregion = shape1nullregion, credregion = shape1credregion)
    shape2ropehdi_res = hdi_plus_rope(nullregion = shape2nullregion, credregion = shape2credregion)
    res = c(shape1ropehdi_res, shape2ropehdi_res)
    
    if(option == 1){
      if(sum(is.na(res)) == 2){ # both undecided
        return(1)
      }
      
      if(sum(is.na(res)) == 1){ # one undecided, other leads the combined result
        out = res[!is.na(res)]
        return(out)
      }
      
      if(sum(is.na(res)) == 0){ # both decided
        out = ifelse(sum(res) == 0, 0, 1)
        return(out)
      }
    }
    if(option == 2){
      if(sum(is.na(res)) == 2){ # both undecided
        return(0)
      }
      
      if(sum(is.na(res)) == 1){ # one undecided, other leads the combined result
        out = res[!is.na(res)]
        return(out)
      }
      
      if(sum(is.na(res)) == 0){ # both decided
        out = ifelse(sum(res) == 2, 1, 0)
        return(out)
      }
    }
    if(option == 3){
      if(sum(is.na(res)) == 2){ # both undecided
        return(0)
      }
      
      if(sum(is.na(res)) == 1){ # one undecided
        return(0)
      }
      
      if(sum(is.na(res)) == 0){ # both decided
        out = ifelse(sum(res) == 2, 1, 0)
        return(out)
      }
    }
  }
  
}