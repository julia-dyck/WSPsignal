#' Set up simulation parameters
#' 
#' Sets up parameters for a simulation study to tune the Weibull shape parameter (WSP) test.
#' Simulation parameters encompass data generating process (DGP) parameters (\code{N,..., study.period}),
#' tuning parameters for the WSP test (\code{est.approach,..., sensitivity.option}), 
#' and additional parameters (\code{reps,..., stanmod.warmup}).
#' 
#' @return A list containing all simulation specifications in the format required for
#' \code{\link{sim.run}} or \code{\link{sim.run_parallel}}.
#'
#' @param N vector of sample sizes
#' @param br vector of background rates (observed in population on average)
#' @param adr.rate vector of adverse drug reaction rates as proportions of the background rates
#' @param adr.relsd vector of relative standard deviations from the adverse drug reaction times
#' @param study.period scalar specifying the length of the study period
#' @param est.approach character vector specifying one or two estimation approaches; options are
#'        Bayesian \code{"b"} and frequentist \code{"f"}
#' @param tte.dist character vector specifying one or multiple modelling approaches; options are
#' \code{"w", "dw", "pgw"} (see \code{\link{bwsp_model}}, \code{\link{fwsp_model}})
#' @param prior.dist character indicating the prior distribution for the parameters 
#' of the Bayesian time-to-event (tte) distribution; options are "fg", "fl", "gg", "ll" (see \code{\link{bwsp_model}})
#' @param fitpars.list list with one dataframe per tte distribution containing
#' the prior specifications for BWSP model fitting; setup with \code{\link{sim.priors_template}}
#' @param post.ci.type character indicating whether to extract equal tailed
#' intervals (\code{"ETI"}) or highest posterior density intervals (\code{"HDI"}) as
#' credibility interval (CI) for BWSP testing (see \code{\link{bwsp_test}})
#' @param cred.level vector of credibility levels used for construction
#' of region of practical equivalence (ROPE) and posterior CI for the BWSP tests (see \code{\link{bwsp_test}}) or
#' vector of confidence levels for FWSP tests (see \code{\link{fwsp_test}}).
#' @param sensitivity.option vector of sensitivity options for the BWSP test 
#' (see \code{\link{bwsp_test}})
#' @param reps number of repetitions for each simulation scenario, default is 100
#' @param batch.size number of simulation repetitions to be saved in a batch (see details); 
#' default is 10
#' @param resultpath directory where output of this function and intermediate results of the simulation
#' are saved
#' @param stanmod.chains number of Markov chains (see \code{\link[rstan]{sampling}}); 
#' default is 4
#' @param stanmod.iter total number of iterations per chain including warmup 
#' (see \code{\link[rstan]{sampling}}); default is 11000
#' @param stanmod.warmup number of warmup (aka burn-in) iterations per chain 
#' (see \code{\link[rstan]{sampling}}); default is 1000
#' 
#' 
#' @details
#' The purpose of the simulation study is to evaluate the performance of different WSP 
#' tests for data scenarios of interest following the tuning
#' scheme developed in \insertCite{dyck2024bpgwsppreprint;textual}{WSPsignal}.
#' 
#' DGP parameters (\code{N,..., study.period}) should 
#' reflect the data characteristics of interest.
#' Within simulation, data are generated with \code{\link{sim.datagen_tte}}.
#' 
#' Tuning parameters for the WSP test (\code{est.approach,..., sensitivity.option})
#' lead to a range of tuning combinations for the WSP tests 
#' evaluated during the simulation study to find the best test specification. 
#' Argument \code{fitpars.list} needs to be prepared with 
#' \code{\link{sim.priors_template}}. Note, that the \code{tte.dist} argument 
#' in \code{\link{sim.priors_template}} and in \code{\link{sim.setup_sim_pars}} must
#' match.
#' 
#' 
#'
#' Additional parameters (\code{reps,..., stanmod.warmup}) specify simulation 
#' settings and specifications for posterior sampling.
#' Simulation settings encompass the number of repetitions per simulation scenario, 
#' the directory in which to save results
#' and batch saving.
#' Batch saving is done to prevent losing simulation results in case of an
#' interruption of simulation e.g. due to termination of the R session.
#' Posterior sampling specifications encompass the number of chains, iterations and 
#' the length of the warmup phase. See \code{\link[rstan]{sampling}} for 
#' more details on the posterior estimation function.
#' 
#' @references 
#' \insertAllCited{}
#' 
#' 
#' @examples
#' #### specify parameter combinations for simulation study 
#' 
#' # setup prior template
#' fp_list = sim.priors_template(tte.dist = c("w", "pgw"), prior.sds = 10) 
#' fp_list # fitpars.list template
#' 
#' # fill in prior template with values chosen in prior elicitation
#' # for weibull models:
#' fp_list$w[,2] = c(1, 1, 180, 300) # scale prior means
#' fp_list$w[,3] = c(1, 0.207, 1, 4) # shape prior means
#' 
#' # for pgw models:
#' fp_list$pgw[,2] = c(1, 1, 20, 300)   # scale prior means
#' fp_list$pgw[,3] = c(1, 0.207, 5.5, 4)# shape prior means
#' fp_list$pgw[,4] = c(1, 1, 14, 1)     # powershape prior means
#' 
#' fp_list # fitpars.list filled with means
#' 
#' # setup parameter combination list for simulation
#' pc_list = sim.setup_sim_pars(N = c(500, 3000, 5000),             # DGP parameters
#'                              br = 0.1,                           # |
#'                              adr.rate = c(0, 0.5, 1),            # |
#'                              adr.relsd = 0.05,                   # |
#'                              study.period = 365,                 # v
#'                              est.approach = c("f", "b"),         # Tuning parameters
#'                              tte.dist = c("w", "pgw"),           # |
#'                              prior.dist = c("ll", "gg"),         # |
#'                              fitpars.list = fp_list,             # |
#'                              post.ci.type = c("ETI", "HDI"),     # |
#'                              cred.level = c(seq(0.5,0.9, by = 0.05)), 
#'                              sensitivity.option = 1:3,           # v
#'                              reps = 100,                         # Additional parameters
#'                              batch.size = 10,                    # |
#'                              resultpath = paste0(getwd(),"/simulation_results"),
#'                              stanmod.iter = 11000,               # |
#'                              stanmod.warmup = 1000               # v
#'                              )
#'                              
#' 
#' pc_list 
#'                              
#'                              
#'
#' @export

sim.setup_sim_pars = function(N,                 # dgp parameters
                              br,                # |
                              adr.rate,          # |
                              adr.relsd,         # v
                              study.period,      # -
                              est.approach,      # tuning parameters
                              tte.dist,          # |
                              prior.dist,        # |
                              fitpars.list,      # |
                              post.ci.type = c("ETI", "HDI"),     
                              cred.level,        # v
                              sensitivity.option = 1:3,
                              reps = 100,        # additional parameters
                              batch.size = 10,
                              resultpath = paste0(getwd(), "/results_raw"),
                              stanmod.chains = 4,
                              stanmod.iter = 11000,
                              stanmod.warmup = 1000
){   
  
  # Argument checks ------------------------------------------------------------
  
  ## for DGP parameters ---
  if (!is.numeric(N)) stop("Argument N must be numeric.\n")
  if (any(N <= 0)) stop("Argument N must contain positive values.\n")
  if (any(N %% 1 != 0)) {
    warning("Argument N contains non-integer values that will be rounded.\n")
    N <- round(N)
  }
  if (any(duplicated(N))) {
    warning("Duplicate entries removed from N.\n")
    N <- unique(N)
  }
  
  if (!is.numeric(br)) stop("Argument br must be numeric.\n")
  if (any(br < 0 | br > 1)) stop("Argument br must be between 0 and 1.\n")
  if (any(duplicated(br))) {
    warning("Duplicate entries removed from br.\n")
    br <- unique(br)
  }
  
  if (!is.numeric(adr.rate)) stop("Argument adr.rate must be numeric.\n")
  if (any(adr.rate < 0 | adr.rate > 1)) stop("Argument adr.rate must be between 0 and 1.\n")
  if (any(duplicated(adr.rate))) {
    warning("Duplicate entries removed from adr.rate.\n")
    adr.rate <- unique(adr.rate)
  }
  if (!any(adr.rate == 0)){
    warning("Control case (adr.rate = 0) necessary for simulation study. Value 0 is added to specified vector adr.rate.\n")
  }
  
  if (!is.numeric(adr.relsd)) stop("Argument adr.relsd must be numeric.\n")
  if (any(adr.relsd <= 0)) stop("Argument adr.relsd must be positive.\n")
  if (any(duplicated(adr.relsd))) {
    warning("Duplicate entries removed from adr.relsd.\n")
    adr.relsd <- unique(adr.relsd)
  }
  
  if (!is.numeric(study.period) || length(study.period) != 1)
    stop("Argument study.period must be a single numeric value.\n")
  
  ## for fit parameters ---
  
  ### checks for est.approach
  if (any(duplicated(est.approach))) {
    warning("Duplicate entries removed from est.approach.\n")
    est.approach <- unique(est.approach)
  }
  allowed_apprs <- c("b", "f")
  if (any(is.na(match(est.approach, allowed_apprs))))
    stop(paste0("Argument est.approach must be out of: ",
                paste(allowed_apprs, collapse = ", "),
                ".\n"))
  
  ### checks for tte.dist
  if (any(duplicated(tte.dist))) {
    warning("Duplicate entries removed from tte.dist.\n")
    tte.dist <- unique(tte.dist)
  }
  allowed_dists <- c("w","dw","pgw")
  if (any(is.na(match(tte.dist, allowed_dists))))
    stop(paste0("Argument tte.dist must be out of: ",
                paste(allowed_dists, collapse = ", "),
                ".\n"))
  
  # special case: only frequentist -> overwrite prior info input args 
  if (identical(est.approach, "f")) {
    prior.dist <- NA_character_
    fitpars.list <- list(w = data.frame(), dw = data.frame(), pgw = data.frame())
    post.ci.type = "ETI"
    warning("Arguments prior.dist, fitpars.list, post.ci.type and sensitivity.option are ignored. Reason: only frequentist estimation approach considered in simulation study. \n")
  }
  
  if (!identical(est.approach, "f")) { #if "b" considered, check for prior.dist and fitpars.list input
    ### checks for prior.dist
    if (any(duplicated(prior.dist))) {
      warning("Duplicate entries removed from prior.dist.\n")
      prior.dist <- unique(prior.dist)
    }
    allowed_priors <- c("fg","fl","gg","ll")
    if (any(is.na(match(prior.dist, allowed_priors))))
      stop(paste0("Argument prior.dist must be out of: ", paste(allowed_priors, collapse = ", "),".\n"))
    
    ### checks regarding fitpars.list
    if(!setequal(names(fitpars.list), c("w", "dw", "pgw", "prior.sds")))
      stop("Argument fitpars.list must be a list with elements $w, $dw, $pgw and $prior.sds.\nSetup a fitpars.list template with sim.priors_template.\n")
    # check whether tte.dist matches tte.dist defined fitpars.list
    tte.dist.in.fitpars.list <- names(Filter(function(x) is.data.frame(x) && nrow(x) > 0, fitpars.list))
    # Check matching distributions
    if (!setequal(tte.dist, tte.dist.in.fitpars.list))
      stop("Argument tte.dist must match non-empty fitpars.list data.frames.\nSetup a fitpars.list template with sim.priors_template.\n")
    # Validate each dataframe
    for (dist in tte.dist) {
      df <- fitpars.list[[dist]]
      # numeric check except first col
      if (!all(sapply(df[-1], is.numeric))) 
        stop(paste0("All parameter columns in fitpars.list$", dist," must be numeric (no characters or factors).\nSetup a fitpars.list template with sim.priors_template.\n"))
      # no remaining NA
      if (any(is.na(df[-1]))) 
        stop(paste0("All priors in fitpars.list$", dist, " must be filled - no NA allowed.\nSetup a fitpars.list template with sim.priors_template.\n"))
    }
  }
  
  
  ## for test parameters ---
  
  if (any(duplicated(post.ci.type))) {
    warning("Duplicate entries removed from post.ci.type.\n")
    post.ci.type <- unique(post.ci.type)
  }
  allowed_ci <- c("ETI","HDI")
  if (any(is.na(match(post.ci.type, allowed_ci))))
    stop(paste0("Argument post.ci.type must be out of: ",paste(allowed_ci, collapse = ", "),".\n"))
  
  if (!is.numeric(cred.level) || any(cred.level <= 0 | cred.level >= 1))
    stop("Argument cred.level must contain numeric values between 0 and 1.\n")
  # Credibility interval grid
  if (any(duplicated(cred.level))) {
    warning("Duplicate entries removed from cred.level.\n")
    cred.level <- unique(cred.level)
  }
  
  if (any(duplicated(sensitivity.option))) {
    warning("Duplicate entries removed from sensitivity.option.\n")
    sensitivity.option <- unique(sensitivity.option)
  }
  if (any(is.na(match(sensitivity.option, 1:3))))
    stop("Argument sensitivity.option must be out of 1, 2, 3.\n")
  
  ## for additional parameters ---
  
  if (!is.numeric(reps) || length(reps) != 1)
    stop("Argument reps must be a single numeric value.\n")
  if (reps <= 0)
    stop("Argument reps must be positive.\n")
  if (reps %% 1 != 0) {
    warning("Argument reps contains non-integer value and will be rounded.\n")
    reps <- round(reps)
  }
  
  if (!is.numeric(batch.size) || length(batch.size) != 1)
    stop("Argument batch.size must be a single numeric value.\n")
  if (batch.size <= 0)
    stop("Argument batch.size must be positive.\n")
  if (reps %% batch.size != 0) {
    stop("Argument batch.size must be an integer divisor of reps such that batch.size times an integer number of batches is exactly the specified number of reps.\n")
  }
  
  if (!is.character(resultpath) || length(resultpath) != 1)
    stop("Argument resultpath must be a single character string.\n")
  
  if (!is.numeric(stanmod.chains) || length(stanmod.chains) != 1)
    stop("Argument stanmod.chains must be a single numeric value.\n")
  if (stanmod.chains <= 0)
    stop("Argument stanmod.chains must be positive.\n")
  if (stanmod.chains %% 1 != 0) {
    warning("Argument stanmod.chains contains non-integer value and will be rounded.\n")
    stanmod.chains <- round(stanmod.chains)
  }
  
  if (!is.numeric(stanmod.iter) || length(stanmod.iter) != 1)
    stop("Argument stanmod.iter must be a single numeric value.\n")
  if (stanmod.iter <= 0)
    stop("Argument stanmod.iter must be positive.\n")
  if (stanmod.iter %% 1 != 0) {
    warning("Argument stanmod.iter contains non-integer value and will be rounded.\n")
    stanmod.iter <- round(stanmod.iter)
  }
  
  if (!is.numeric(stanmod.warmup) || length(stanmod.warmup) != 1)
    stop("Argument stanmod.warmup must be a single numeric value.\n")
  if (stanmod.warmup < 0)
    stop("Argument stanmod.warmup must be non-negative.\n")
  if (stanmod.warmup %% 1 != 0) {
    warning("Argument stanmod.warmup contains non-integer value and will be rounded.\n")
    stanmod.warmup <- round(stanmod.warmup)
  }
  
  # Set up list components -----------------------------------------------------
  dgp_pars = sim.setup_dgp_pars(N = N,
                                br = br,
                                adr.when = c(0, 0.25, 0.5, 0.75), # fixed (to reduce complexity)
                                adr.rate = adr.rate,
                                adr.relsd = adr.relsd,
                                study.period = study.period)
  
  test_pars = sim.setup_test_pars(post.ci.type = post.ci.type,
                                  cred.level = cred.level,
                                  sensitivity.option = sensitivity.option)
  
  add_pars = list(reps = reps,
                  batch.size = batch.size,
                  batch.nr = reps/batch.size,
                  resultpath = resultpath,
                  stanmod.chains = stanmod.chains,
                  stanmod.iter = stanmod.iter,
                  stanmod.warmup = stanmod.warmup
  )
  
  input_args = list(N = N, 
                    br = br, 
                    adr.rate = adr.rate, 
                    adr.when = c(0, 0.25, 0.5, 0.75), # fixed (reduce complexity)
                    adr.relsd = adr.relsd, 
                    study.period = study.period,
                    est.approach = est.approach,
                    tte.dist = tte.dist, 
                    prior.belief = c("none", "beginning", "middle", "end"), # fixed (matching adr.when)
                    prior.sds = fitpars.list$prior.sds,
                    prior.dist = prior.dist, 
                    post.ci.type = post.ci.type, 
                    cred.level = cred.level, 
                    sensitivity.option = sensitivity.option,
                    reps = reps, 
                    batch.size = batch.size, 
                    resultpath = resultpath, 
                    stanmod.chains = stanmod.chains, 
                    stanmod.iter = stanmod.iter, 
                    stanmod.warmup = stanmod.warmup
  )
  
  if (!identical(est.approach, "f")) {
    fit_pars = sim.setup_fit_pars(
      tte.dist = tte.dist,
      prior.belief = c("none", "beginning", "middle", "end"),
      prior.dist = prior.dist,
      fit_pars_list = fitpars.list
    )
    
    pc_table = list()
    if(nrow(fit_pars$w) > 0){
      pc_table$w = dplyr::cross_join(dgp_pars, fit_pars$w[c("tte.dist","prior.dist","prior.belief", "prior.sd")])
    } else pc_table$w = NULL
    
    if(nrow(fit_pars$dw) > 0){
      pc_table$dw = dplyr::cross_join(dgp_pars, fit_pars$dw[c("tte.dist","prior.dist","prior.belief", "prior.sd")])
    } else pc_table$dw = NULL
    
    if(nrow(fit_pars$pgw) > 0){
      pc_table$pgw = dplyr::cross_join(dgp_pars, fit_pars$pgw[c("tte.dist","prior.dist","prior.belief", "prior.sd")])
    } else pc_table$pgw = NULL
    
  } else {
    fit_pars = list(w = data.frame(), dw = data.frame(), pgw = data.frame())
    
    pc_table = dplyr::cross_join(
      dgp_pars,
      data.frame(
        tte.dist = tte.dist,
        prior.dist = NA_character_,
        prior.belief = NA_character_,
        prior.sd = NA_character_
      )
    )
  }
  
  pc_table = dplyr::bind_rows(pc_table)
  
  pc_table_freq = dplyr::select(pc_table, -dplyr::contains("prior"))
  pc_table_freq = dplyr::distinct(pc_table_freq)
  
  sim_pars = list(dgp = dgp_pars, fit = fit_pars, test = test_pars, add = add_pars, input = input_args, pc_table = pc_table)
  
  
  
  ## message regarding nr of parameter combinations and reps
  if (identical(est.approach, "b")) {
    message(
      "Each combination of sample scenario and prior specification leads to a total of ",
      nrow(pc_table),
      " Bayesian simulation settings."
    )
  }
  
  if (identical(est.approach, "f")) {
    message(
      "Each combination of sample scenario leads to a total of ",
      nrow(pc_table_freq),
      " frequentist simulation settings."
    )
  }
  
  if (setequal(est.approach, c("b","f"))) {
    message(
      "Each combination of sample scenario (and prior specification) leads to a total of ",
      nrow(pc_table),
      " Bayesian and ",
      nrow(pc_table_freq),
      " frequentist simulation settings."
    )
  }
  
  message(
    "Each simulation scenario's data generation and model estimation will be repeated ",
    reps,
    " times."
  )
  
  
  return(sim_pars)
}


## END OF DOC
