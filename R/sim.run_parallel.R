#' Run simulation (in parallel)
#' 
#' Runs simulations in parallel for all specified data-generating processes, model 
#' and test configurations defined in \code{pc_list}. If simulations for some scenarios 
#' have already been completed, the function resumes from the remaining scenarios.
#' 
#' 
#' @param pc_list list of parameter combinations generated with \code{\link{sim.setup_sim_pars}}
#' @param subset_ind vector of integers specifying which rows of \code{pc_list$pc_table} 
#' to be considered in simulation runs; defaults to all rows
#' 
#' @details Parallelization needs to be initialized using the \code{\link[future]{plan}} 
#' command (see example).
#' 
#' 
#' @examples 
#' \dontrun{
#'  # install.packages(future)
#'  future::plan(multisession, workers = availableCores()) # or another plan strategy
#'  sim.run_parallel(pc_list) # run all simulations
#'  
#'  
#'  # To run only a subset of simulation scenarios, specify the indices of the rows
#'  # in pc_list$pc_table that you want to run. For example, to run the first 10 scenarios:
#'  
#'  # install.packages(future)
#'  future::plan(multisession, workers = availableCores()) # or another plan strategy
#'  sim.run_parallel(pc_list, subset_ind = 1:10) # run first 10 simulation scenarios
#' }
#' 
#' 
#' 
#' 
#' @export
#' 
#' 


sim.run_parallel = function(pc_list, subset_ind = NULL) {
  
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
  
  # argument check for subset_ind
  if (!is.null(subset_ind)) {
    subset_ok <- 
      is.numeric(subset_ind) &&
      !any(is.na(subset_ind)) &&
      all(subset_ind %% 1 == 0) &&          # integer-valued
      all(subset_ind > 0) &&                # positive
      all(subset_ind <= nrow(pc_list$pc_table))  # within range
    if (!subset_ok) {
      stop("Argument subset_ind must be NULL or a numeric vector of valid row indices of pc_list$pc_table.\n")
    }
    # auto-deduplicate with warning (consistent with earlier design)
    if (any(duplicated(subset_ind))) {
      warning("Duplicate entries removed from subset_ind.\n")
      subset_ind <- unique(subset_ind)
    }
  }
  
  ## fct body ------------------------------------------------------------------
  if (!is.null(subset_ind)) {
    pc_table = pc_list$pc_table[subset_ind, ]
  } else {
    pc_table = pc_list$pc_table## run 
  }
  
  pc_table_ext = dplyr::cross_join(pc_table, data.frame(batch_nr = 1:pc_list$add$batch.nr))
  
  # split table rows into list elements
  pc_table_as_list = split(pc_table_ext, seq_len(nrow(pc_table_ext)))
  
  # Parallelized execution using furrr::future_walk
  furrr::future_walk(
    .x = pc_table_as_list,
    .f = function(row) {
      pc_vect = row[1:10]
      ind.batch = row$batch_nr
      
      tryCatch({
        sim.load.scenario(pc = pc_vect, wd = pc_list$add$resultpath, batchnr = ind.batch)
      },
      error = function(cond) {
        sim.repeat.1.scenario(
          pc = pc_vect,
          pc_list = pc_list,
          batch.ind = ind.batch
        )
      })
    },
    .options = furrr::furrr_options(seed = NULL)  # optional seed if reproducibility isn't critical
  )
}



## END OF DOC