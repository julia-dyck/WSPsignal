#' Load simulated scenario
#' 
#' Load all simulations of a simulated scenario (obtained using 
#' \code{\link{sim.repeat.1.scenario}}).
#'
#' @param wd working directory where ALL the simulation files lie
#' 
#' @param pc parameter combination (numeric version)
#' 
#' @param batchnr nr of files (batches) containing simulation output based on the 
#' same pc
#' 
#' @param bayes logical argument to specify whether to load the simulation results 
#' obtained with Bayesian (, if \code{bayes = TRUE}; default) or frequentist estimation 
#' approach
#' 
#' @return A dataframe \code{res.batch}.
#' 
#' @noRd


sim.load.scenario = function(wd, pc, batchnr = 1, bayes = TRUE){
  
  if(bayes == TRUE){
    filename = paste(c(pc, "bADR_sim", batchnr, ".RData") ,collapse="_")
  }
  if(bayes == FALSE){
    filename = paste(c(pc, "fADR_sim", batchnr, ".RData") ,collapse="_")
  }
  
  path = paste0(wd, "/", filename)
  
  # return(path)
  load(file = path)
  res = res.batch
  return(res)
}





## END OF DOC