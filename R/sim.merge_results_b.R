#' Merge result table batches from simulation study
#' 
#' Merges result table batches from simulation study obtained from using 
#' \code{\link{sim.run}}.
#'
#' @param pc_list list of parameter combinations obtained from \link{sim.setup_simpars}
#' 
#' @return Dataframe containing all simulation results (one repetition of one 
#' simulation scenario per row). The 
#' simulation parameters are stored in the first 9 columns. The remaining columns 
#' contain 
#' \enumerate{
#'       \item posterior summary statistics and percentiles (ie information on the 
#'       posterior distribution) for each shape parameter and 
#'       \item posterior credibility 
#'       intervals (as specified in \code{$test} list element obtained from \link{sim.setup_simpars}).
#'       }
#'
#' @noRd

sim.merge_results_b = function(pc_list){
  
  # setup empty dfs
  merged.res.w =  data.frame()
  merged.res.dw =  data.frame()
  merged.res.pgw = data.frame()
  
  # go through all pc combinations in pc_list:
  
  for(ind.dgp in 1:nrow(pc_list$dgp)){      # go through dgp scenarios (per row)
    if(nrow(pc_list$fit$w)>0){
      for(ind.fitw in 1:nrow(pc_list$fit$w)){ # go through weibull fitting parameter combis
        # set up one dgp+fit combination
        pc_vect = sim.gather_pc_vect(pc_list$dgp[ind.dgp,], pc_list$fit$w[ind.fitw,c("tte.dist", "prior.dist", "prior.belief")])
        # go through batches for one parcombi
        for(ind.batch in 1:pc_list$add$batch.nr){
          tryCatch({
            # try merge as new row to existing part
            merged.res.w = dplyr::bind_rows(merged.res.w,
                                 sim.load.scenario(pc = pc_vect, wd= pc_list$add$resultpath, batchnr = ind.batch))
          },
          error=function(cond) {
            # if not, return a warning 
            warning(paste0(paste(c(pc_vect, "bADR_sim", ind.batch, ".RData") ,collapse="_"), " not found."))
          }
          )
          
        }
      }
    }
    
    if(nrow(pc_list$fit$dw)>0){
      for(ind.fitdw in 1:nrow(pc_list$fit$dw)){
        # set up one dgp+fit combination
        pc_vect = sim.gather_pc_vect(pc_list$dgp[ind.dgp,], pc_list$fit$dw[ind.fitdw,c("tte.dist", "prior.dist", "prior.belief")])
        # go through batches for one parcombi
        for(ind.batch in 1:pc_list$add$batch.nr){
          tryCatch({
            # try merge as new row to existing part
            merged.res.dw = dplyr::bind_rows(merged.res.dw,
                                 sim.load.scenario(pc = pc_vect, wd= pc_list$add$resultpath, batchnr = ind.batch))
          },
          error=function(cond) {
            # if not, return a warning 
            warning(paste0(paste(c(pc_vect, "bADR_sim", ind.batch, ".RData") ,collapse="_"), " not found."))
          }
          )
        }
      }
    }
    
    if(nrow(pc_list$fit$pgw)>0){
      for(ind.fitpgw in 1:nrow(pc_list$fit$pgw)){
        # set up one dgp+fit combination
        pc_vect = sim.gather_pc_vect(pc_list$dgp[ind.dgp,], pc_list$fit$pgw[ind.fitpgw,c("tte.dist", "prior.dist", "prior.belief")])
        # go through batches for one parcombi
        for(ind.batch in 1:pc_list$add$batch.nr){
          tryCatch({
            # try merge as new row to existing part
            merged.res.pgw = dplyr::bind_rows(merged.res.pgw,
                                 sim.load.scenario(pc = pc_vect, wd= pc_list$add$resultpath, batchnr = ind.batch))
          },
          error=function(cond) {
            # if not, return a warning 
            warning(paste0(paste(c(pc_vect, "bADR_sim", ind.batch, ".RData") ,collapse="_"), " not found."))
          }
          )
        }
      }
    }
  }
  
  # merge output for all tte.dists
  res_b = dplyr::bind_rows(merged.res.w, merged.res.dw, merged.res.pgw)
  return(res_b)
}

## END OF DOC