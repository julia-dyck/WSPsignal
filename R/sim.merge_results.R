#' Merge result table batches from simulation study
#' 
#' Merges result table batches from simulation study obtained from using 
#' \code{\link{sim.run}} or \code{\link{sim.run_parallel}} as preparation for evaluation..
#'
#' @param pc_list list of parameter combinations generated with \code{\link{sim.setup_sim_pars}}
#' @param save if \code{TRUE} (default), merged table is saved as res_b.RData or 
#' res_f.RData, respectively, in same path where batches are stored; 
#' else, result table(s) is/are returned to global environment
#' 
#' @return Dataframe or list of two dataframes containing all simulation results (one repetition of one 
#' simulation scenario per row). 
#'       
#' @export

sim.merge_results = function(pc_list, save = T){
  
  ## argument checks -----------------------------------------------------------
  # argument check for pc_list
  pc_list_is_valid <-
    is.list(pc_list) &&
    # pc_list$dgp must be a data.frame
    !is.null(pc_list$dgp) && 
    is.data.frame(pc_list$dgp) &&
    # pc_list$fit must be a list whose elements are all data.frames
    !is.null(pc_list$fit) && 
    is.list(pc_list$fit) &&
    length(pc_list$fit) > 0 &&
    all(vapply(pc_list$fit, is.data.frame, logical(1))) &&
    # pc_list$test must be a list
    !is.null(pc_list$test) && 
    is.list(pc_list$test) &&
    # pc_list$add must be a list with required numeric/character elements
    !is.null(pc_list$add) &&
    is.list(pc_list$add) &&
    is.numeric(pc_list$add$reps) &&
    is.numeric(pc_list$add$batch.size) &&
    is.numeric(pc_list$add$batch.nr) &&
    is.character(pc_list$add$resultpath) &&
    is.numeric(pc_list$add$stanmod.chains) &&
    is.numeric(pc_list$add$stanmod.iter) &&
    is.numeric(pc_list$add$stanmod.warmup) &&
    # pc_list$pc_table must be a non-empty data.frame
    !is.null(pc_list$pc_table) &&
    is.data.frame(pc_list$pc_table)
  
  if (!pc_list_is_valid) {
    stop("Argument pc_list has wrong format. It must be a list produced by sim.setup_sim_pars().\n")
  }
  
  # argument check for save
  if (!is.logical(save) || length(save) != 1 || is.na(save)) {
    stop("Argument 'save' must be a single logical value (TRUE or FALSE).\n")
  }
  
  ## fct body ------------------------------------------------------------------
  est.approach = pc_list$input$est.approach
  
  if("b" %in% est.approach){
    res_b = sim.merge_results_b(pc_list = pc_list)
    if(save == T){
      # save result
      path = pc_list$add$resultpath
      filename = "res_b.RData"
      save(res_b, file=paste0(path, "/", filename))
      message(sprintf("res_b saved to: %s", file.path(path, filename)))
    }
  }
  if("f" %in% est.approach){
    res_f = sim.merge_results_f(pc_list = pc_list)
    if(save == T){
      # save result
      path = pc_list$add$resultpath
      filename = "res_f.RData"
      save(res_f, file=paste0(path, "/", filename))
      message(sprintf("res_f saved to: %s", file.path(path, filename)))
    }
  }
  if(save == F){
    if("b" %in% est.approach && "f" %in% est.approach){
      return(list(res_b = res_b, res_f = res_f))
    } else if("b" %in% est.approach){
      return(res_b)
    } else if("f" %in% est.approach){
      return(res_f)
    } else {
      return(NULL)
    }
  }
}

## END OF DOC