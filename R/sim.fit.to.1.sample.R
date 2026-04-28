#' Simstudy: apply all models to one specific sample
#'
#' Conducts all steps from data generation to simulation output generation including:
#' 
#' - data generation based on provided parameter combination
#'
#' - data preparation for stanmodel for all prior assumptions
#'
#' - fitting all model alternatives
#'
#' - extracting and returning relevant statistics from stan output:
#'   - general information "info" (noch meta), eg. prior mean and sds, fitting specifications etc.
#'   - information about the posterior distributions of nu and gamma "post"
#'   - credibility intervals for the tests (will be conducted afterwards)
#'   - 0 to 100 percentiles of the posterior distribution (to roughly estimate the
#'     probability mass within the region of practical equivalence for the null value)
#'
#' @param pc parameter combination 
#'
#' @return a data frame with 16 rows. Each row contains statistics for one of the 
#' 4x4 prior (fl, ll, fg, gg) and model (no ADR expected, ADR expected around 
#' 1st, 2nd, 3rd quartile of observation period) alternatives.
#'
#' @noRd



sim.fit.to.1.sample = function(pc, pc_list){
  est.approach = pc_list$input$est.approach
  
  ### Data simulation
  ttedat = sim.datagen_tte(genpar = c(pc$N, pc$br, pc$adr.rate, pc$adr.when, pc$adr.relsd, pc$study.period))

  
  ### Bayesian model fitting
  if("b" %in% est.approach){
    # tte and prior data preparation
    datstan = sim.fit.prep(ttedat = ttedat, pc = pc, pc_list = pc_list)
    
    mod = tryCatch(
      bwsp_model(datstan = datstan,
                 chains = pc_list$add$stanmod.chains,
                 iter = pc_list$add$stanmod.iter,
                 warmup = pc_list$add$stanmod.warmup),
      error = function(e) {
        return(NULL)
      }
    )
      
    btestres = tryCatch(
      bwsp_test(mod, 
                cred.level = pc_list$input$cred.level, 
                ci.type = pc_list$input$post.ci.type, 
                sensitivity.option = pc_list$input$sensitivity.option),
      error = function(e) {
        return(e)
        return(NULL)
      }
    )
    
    ### extracting Bayesian posterior statistics
    bstats = tryCatch(
      sim.stanfit.to.poststats(pc, 
                               stanfit.object = mod$fit
      ),
      error = function(e) {
        return(NULL)
      }
    )
    btests = c(bstats, btestres)
  }
  else{btests = NULL}
  
  ### Frequentist model fitting (MLE)
  if("f" %in% est.approach){
    mod.w = fwsp_model(dat = ttedat, tte.dist = "w")
    test.w = fwsp_test(mod.w, cred.level = pc_list$input$cred.level)
    mod.dw = fwsp_model(dat = ttedat, tte.dist = "dw")
    test.dw = fwsp_test(mod.dw, cred.level = pc_list$input$cred.level)
    
    mod.pgw = tryCatch(
      fwsp_model(dat = ttedat, tte.dist = "pgw"),
      error = function(e) {
        return(NULL)
      }
    )
    test.pgw = tryCatch(
      fwsp_test(mod.pgw, cred.level = pc_list$input$cred.level),
      error = function(e) {
        return(rep(NA, length(pc_list$input$cred.level)))
      }
    )
    
    ### formatting frequentist results
    ftests = c(pc, test.w, test.dw, test.pgw)
  }
  else{ftests = NULL}
  
  
  return(list(btests = btests, ftests = ftests)) # return both Bayesian and frequentist test results

}



## END OF DOC