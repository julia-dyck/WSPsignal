#' Bayesian Weibull Shape Parameter Test
#' 
#' @description 
#' Bayesian Weibull Shape Parameter (BWSP) test of the constant hazard (null-)hypothesis,
#' based on the shape parameter(s) of Weibull family of distributions.
#' 
#'
#' @param mod.output model output resulting from \code{\link{bwsp_model}}
#' @param cred.level numeric or vector of credibility levels; default is 0.8
#' @param ci.type  character indicating whether to extract an equal tailed
#' interval (\code{"ETI"}) or highest posterior density interval (\code{"HDI"}) as
#' posterior credibility interval (CI) for the BWSP test; default is \code{"HDI"}
#' @param sensitivity.option numeric value out of \code{1,2,3}; combination rule to deduct a 
#' binary outcome (signal/no signal) from one or two shape parameter tests; 
#' default is \code{sensitivity.option = 2} (see details)
#' 
#' 
#' @return binary vector, 0 if \eqn{H_0} is accepted (no signal), 1 if \eqn{H_1} is accepted (signal)
#'
#'
#' @section Test concept: 
#' 
#' The BWSP test is based on the principle of (non-)constant hazard 
#' \insertCite{cornelius2012}{WSPsignal}
#' which associates a constant hazard function with the absence of a drug-adverse event
#' association 
#' and a non-constant hazard with the presence of a drug-adverse event
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
#'  under Power generalized Weibull model \tab \eqn{\nu = 1 \text{ and } \gamma = 1} \tab \eqn{\nu \neq 1 \text{ or } \gamma \neq 1} \cr
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
#' The region of practical equivalence (ROPE) represents the expected parameter value under \eqn{H_0}.
#' The posterior credibility interval(s) (CI) represent the posterior distribution 
#' of each shape parameter. 
#' For the ROPE, the function sets up an equal-tailed interval (ETI) 
#' \deqn{[q_{(1-\alpha)/2}, q_{(1+\alpha)/2}]}
#' based on the quantiles \eqn{q} of the shape parameters' prior distributions under 
#' \eqn{H_0} at a chosen credibility level \eqn{1 - \alpha}.
#' 
#' For the posterior CI, the function calculates either an ETI at the same credibility level obtained from
#' the empirical quantiles of the posterior distribution per shape parameter 
#' or a highest density interval (HDI, \insertCite{kruschke2015;textual}{WSPsignal})
#' \deqn{HDI(\nu) = \{\nu \; |\; p_1(\nu) \geq w\} \text{ with } w\in [0,1] 
#' \text{ such that} \int_{\nu \; | \; p_1(\nu) \geq w} p_1(\nu|t)\;  d\nu = 1 - \alpha}
#' 
#' where \eqn{\nu} is one of the shape parameters, \eqn{p_1} it's posterior density and \eqn{t} the time variable.
#' 
#' 
#' The CI+ROPE test \insertCite{kruschke2018}{WSPsignal} checks the 
#' relationship between ROPE and posterior CI leading to either acceptance,
#' rejection or no decision regarding the null hypothesis for a single shape parameter.
#' Sensitivity options to generate a binary outcome, i.e. a signal or not, from CI+ROPE test results 
#' based on one (in case of \code{"w"}) or two (in case of \code{"dw", "pgw"}) shape parameters are:
#' \tabular{ccccc}{
#' HDI+ROPE \tab HDI+ROPE \tab combination \tab combination \tab combination \cr
#' outcome \tab outcome \tab rule \tab rule \tab rule \cr
#' for shape_1\tab for shape_2 \tab (\code{sensitivity.option = 1}) \tab (\code{sensitivity.option = 2}) \tab (\code{sensitivity.option = 3}) \cr
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
#' The hypotheses as stated above (see test concept) are implemented in \code{sensitivity.option = 1} whereas
#' \code{sensitivity.option = 2} and \code{sensitivity.option = 3} lead to a signal in fewer cases.
#' 
#' More details on the CI+ROPE test, recommendations for interval specifications
#' and the combination rules
#' can be found in \insertCite{kruschke2018;textual}{WSPsignal} and
#' \insertCite{dyck2024bpgwsppreprint;textual}{WSPsignal}.
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
#' # 1. prior specification
#' # we formalize a prior belief (here "no association
#' # between drug and event", therefore prior mean = 1 for shape parameter)
#' # and reformat our tte data to fit the model in the following
#' dat_list = tte2priordat(dat = tte,   # reformat the data
#'                       tte.dist = "w",
#'                       prior.dist = "ll", 
#'                       scale.mean = 1, 
#'                       scale.sd = 10,
#'                       shape.mean = 1, 
#'                       shape.sd = 10)
#' 
#' # 2. model fitting
#' fit = bwsp_model(datstan = dat_list,      # fit the model       
#'                  chains = 4,              
#'                  iter = 110,             # (posterior sample is
#'                  warmup = 10)            # small for demo purpose)
#' fit$fit
#' 
#' # 3. BWSP test
#' bwsp_test(mod.output = fit,
#'           cred.level = 0.8,
#'           ci.type = "HDI",
#'           sensitivity.option = 2)
#' 
#' # under pgw model:
#' 
#' # 1. prior specification
#' # we formalize a prior belief (here "no association
#' # between drug and event", therefore prior mean = 1 for both shape parameters)
#' # and reformat our tte data to fit the model in the following
#' dat_list = tte2priordat(dat = tte,          # reformat the data
#'                       tte.dist = "pgw",
#'                       prior.dist = "ll",
#'                       scale.mean = 1, 
#'                       scale.sd = 10,
#'                       shape.mean = 1, 
#'                       shape.sd = 10,
#'                       powershape.mean = 1, 
#'                       powershape.sd = 10)
#' 
#' # 2. model fitting
#' fit = bwsp_model(datstan = dat_list,     # fit the model      
#'                   chains = 4,              
#'                   iter = 110,            # (posterior sample
#'                   warmup = 10)           # is small for demo purpose)
#' 
#' # 3. BWSP test
#' bwsp_test(mod.output = fit,
#'           cred.level = 0.8,
#'           ci.type = "HDI", 
#'           sensitivity.option = 2)
#' 
#' 
#' 
#' @export
#'


bwsp_test = function(mod.output, 
                     cred.level = 0.8,
                     ci.type = "HDI",
                     sensitivity.option = 2){
  
  # argument check mod.output
  valid_mod.output <-
    is.list(mod.output) &&
    all(c("fit", "args_list") %in% names(mod.output)) &&
    is.list(mod.output$args_list) &&
    all(c("tte.dist", "prior.dist") %in% names(mod.output$args_list))
  
  if (!valid_mod.output) {
    stop("Argument 'mod.output' must be a list generated by bwsp_model().")
  } 
  
  # argument check cred.level
  if (!is.numeric(cred.level)) stop("Argument cred.level must be numeric.\n")
  if (any(cred.level < 0 | cred.level > 1)) stop("Argument cred.level must be between 0 and 1.\n")
  if (any(duplicated(cred.level))) {
    warning("Duplicate entries removed from cred.level.\n")
    cred.level <- unique(cred.level)
  }
  
  # argument check ci.type
  if (any(duplicated(ci.type))) {
    warning("Duplicate entries removed from ci.type.\n")
    ci.type <- unique(ci.type)
  }
  allowed_ci <- c("ETI","HDI")
  if (any(is.na(match(ci.type, allowed_ci))))
    stop(paste0("Argument ci.type must be out of: ",paste(allowed_ci, collapse = ", "),".\n"))
  
  # argument check sensitivity.option
  if (any(duplicated(sensitivity.option))) {
    warning("Duplicate entries removed from sensitivity.option.\n")
    sensitivity.option <- unique(sensitivity.option)
  }
  allowed_so <- 1:3
  if (any(is.na(match(sensitivity.option, allowed_so))))
    stop(paste0("Argument sensitivity.option must be out of: ",paste(allowed_so, collapse = ", "),".\n"))
  
  # ----------------------------------------------------------------------------
# inner fct for one combination of cred.level, ci.type and sensitivity.option
  test_one = function(cred.level_ci.type_sensitivity.option_vect){
    cred.level = as.numeric(cred.level_ci.type_sensitivity.option_vect[1])
    ci.type = cred.level_ci.type_sensitivity.option_vect[2]
    sensitivity.option = as.numeric(cred.level_ci.type_sensitivity.option_vect[3])
    
    tte.dist = mod.output$args_list$tte.dist
    prior.dist = mod.output$args_list$prior.dist
    shape.sd = mod.output$args_list$shape.sd
    shape_c.sd = mod.output$args_list$shape_c.sd
    powershape.sd = mod.output$args_list$powershape.sd
    
    rope_vect = setup_rope(tte.dist = tte.dist, prior.dist = prior.dist, 
                           cred.level = cred.level, shape.sd = shape.sd, 
                           shape_c.sd = shape_c.sd, powershape.sd = powershape.sd)
    
    ci_vect = extract_post_ci(mod.output = mod.output, ci.type = ci.type, cred.level = cred.level)
    
    # test under Weibull model
    if(tte.dist == "w"){
      # cat("shape ROPE:", rope_vect, "\n")
      # cat("shape posterior CI:", ci_vect, "\n\n")
      
      res = hdi_plus_rope(nullregion = rope_vect, credregion = ci_vect)
      
      # sensitivity.options for "w" model
      if(sensitivity.option == 1){ # CI+ROPE test rejects H0 or undecided -> signal
        out = ifelse(is.na(res), 1, res)
      }
      if(sensitivity.option == 2 | sensitivity.option == 3){ # CI+ROPE test rejects H0 -> signal
        out = ifelse(is.na(res), 0, res)
      }
    }
    # test under double or pgw model
    if(tte.dist == "dw" | tte.dist == "pgw"){
      # cat("shape 1 ROPE:", rope_vect[1:2], "\nshape 2 ROPE:", rope_vect[3:4], "\n")
      # cat("shape posterior CI:", ci_vect[1:2], "\nshape 2 posterior CI:", ci_vect[3:4], "\n\n")
      
      res.shape = hdi_plus_rope(nullregion = rope_vect[1:2], credregion = ci_vect[1:2]) # single shape parameter CI+ROPE test
      res.shape2 = hdi_plus_rope(nullregion = rope_vect[3:4], credregion = ci_vect[3:4]) # single shape_c parameter CI+ROPE test
      res = c(res.shape, res.shape2)
      
      # sensitivity.options for "dw" or "pgw" model
      if(sensitivity.option == 1){
        if(sum(is.na(res)) == 2){ # both undecided
          out = 1
        }
        
        if(sum(is.na(res)) == 1){ # one undecided, other leads the combined result
          out = res[!is.na(res)]
        }
        
        if(sum(is.na(res)) == 0){ # both decided
          out = ifelse(sum(res) == 0, 0, 1)
        }
      }
      if(sensitivity.option == 2){
        if(sum(is.na(res)) == 2){ # both undecided
          out = 0
        }
        
        if(sum(is.na(res)) == 1){ # one undecided, other leads the combined result
          out = res[!is.na(res)]
        }
        
        if(sum(is.na(res)) == 0){ # both decided, 2 rejections
          out = ifelse(sum(res) == 2, 1, 0)
        }
      }
      if(sensitivity.option == 3){
        if(sum(is.na(res)) == 2){ # both undecided
          out = 0
        }
        
        if(sum(is.na(res)) == 1){ # one undecided
          out = 0
        }
        
        if(sum(is.na(res)) == 0){ # both decided
          out = ifelse(sum(res) == 2, 1, 0)
        }
      }
    }
    
    # name each vector entry according to the credibility level
    names(out) = paste0("bwsp_", tte.dist, "_", prior.dist, "_", cred.level, "_", ci.type, "_", sensitivity.option)
    return(out)
    
  }
  
  # apply for multiple cred.levels
  test_specifications = expand.grid(cred.level, ci.type, sensitivity.option,
                                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  test.res.list <- lapply(
    seq_len(nrow(test_specifications)),
    function(i) test_one(unlist(test_specifications[i, ]))
  )
  
  test.res <- do.call(c, test.res.list)
  return(test.res)
  
  
}

## END OF DOC