#' Monitor the progress of simulationstudy
#' 
#' The function monitors the progress of the simulation study by going through 
#' raw result files.
#' This can be done to find the restarting point when the running of simulations 
#' was interrupted or to investigate which simulations could not be run
#' successfully (due to convergence issues).
#' 
#' @param pc_list simulation study parameters; mainly used here are pc_list$pc_table,
#' a data frame with all simulation scenarios and pc_list$add$reps stating the planned 
#' number of repetitions per scenario.
#' 
#' @return a vector with the number of successfully run simulations for each pc
#' 
#' @noRd


sim.monitor.progress = function(pc_list){
  progress = c()
  for(i in 1:nrow(pc_list$pc_table)){
    merged.res.one.pc = data.frame()
    for(ind.batch in 1:pc_list$add$batch.nr){
      
      tryCatch({
        # try merge as new row to existing part
        merged.res.one.pc = dplyr::bind_rows(merged.res.one.pc,
                                          sim.load.scenario(pc = pc_list$pc_table[i,], 
                                                            wd= pc_list$add$resultpath, 
                                                            batchnr = ind.batch,
                                                            bayes = T))
      },
      error=function(cond) {
        # if not, return a warning 
        warning(paste0(paste(c(pc_list$pc_table[i,], "bADR_sim", ind.batch, ".RData") ,collapse="_"), " not found."))
      }
      )
    }
    
    progress = c(progress, nrow(merged.res.one.pc))
    
  }
  
  # cbind(pc_list$pc_table, progress)
  
  return(progress)
}

