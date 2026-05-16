#' Evaluate number of non-convergence cases
#'
#' Only for Bayesian estimation approach.
#' Summarizes the number of planned vs. not successfully run simulations optionally
#' grouped by one or multiple model specifications.
#' This helps assess which of the modelling choices is suitable 
#' for CI+ROPE testing (along with other diagnostics such as 
#' \code{\link{eval.execution_times}} and \code{\link{eval.eff_sample_sizes}}).
#'
#' @param pc_list list of simulation parameters generated with \code{\link{sim.setup_sim_pars}}
#' @param group.by character vector specifying grouping variables; must be a 
#' subset of \code{c("tte.dist", "prior.dist", "prior.sd")}
#' 
#' @return A data frame with the following columns:
#'   - `tte.dist`: The tte distribution as grouping factor (if selected).
#'   - `prior.dist`: The prior distribution as grouping factor (if selected).
#'   - `prior.sd`: The prior standard deviation as grouping factor (if selected).
#'   - `total.planned`: The total number of planned repetitions in this group.
#'   - `total.notrun`: The total number of repetitions that were not run in this group.
#'   - `prop.notrun`: The proportion of repetitions that were not run in this group.
#'
#' @details
#' Calculations are based on the stored result file obtained with
#' \code{\link{sim.merge_results}}. 
#'
#' @seealso  \code{\link{eval.execution_times}}, \code{\link{eval.eff_sample_sizes}}
#' 
#' @examples
#' \dontrun{
#' # summarize non-convergence cases by tte and prior distribution and prior sd
#' eval.non_conv_cases(pc_list)
#'
#' # summarize non-convergence cases only by tte distribution
#' eval.non_conv_cases(pc_list, group.by = "tte.dist")
#' }
#' 
#' @export



eval.non_conv_cases = function(pc_list, group.by = c("tte.dist", "prior.dist", "prior.sd")){
  
  ## argument checks -----------------------------------------------------------
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
  
  if(!("b" %in% pc_list$input$est.approach)){
    stop("Functions eval.non_conv_cases(), eval.execution_times() and eval.eff_sample_sizes() only applicable if Bayesian estimation is included in simulation study.\n")
  }
  
  # group.by argument
  allowed_group_by <- c("tte.dist", "prior.dist", "prior.sd")
  
  if (any(duplicated(group.by))) {
    warning("Duplicate entries removed from group.by.\n")
    group.by <- unique(group.by)
  }
  
  if (!all(group.by %in% allowed_group_by)) {
    stop(paste0(
      "Argument 'group.by' must be a subset of: ",
      paste(allowed_group_by, collapse = ", "),
      ".\n"
    ))
  }
  
  ## fct body ------------------------------------------------------------------
  
  # shorten object names
  pc_table = pc_list$pc_table
  reps = pc_list$add$reps
  batch.size = pc_list$add$batch.size
  batch.nr = pc_list$add$batch.nr
  
  # monitor progress: how many batches were successfully run/ not run in simulation?
  prog = sim.monitor.progress(pc_list)
  
  # add information about planned, run and not run repetitions to pc table
  pc_table = cbind(pc_table, reps.done = prog, reps.planned = reps) 
  pc_table = cbind(pc_table, reps.notrun = pc_table$reps.planned - pc_table$reps.done)
  
  
  # sum up numbers for groups depending on tte.dist, prior.dist, prior.sd 
  # (or subset specified in group.by argument)
  nonconv.tab = dplyr::group_by(pc_table, dplyr::across(dplyr::all_of(group.by)))
  
  nonconv.tab = dplyr::summarise(
    nonconv.tab,
    total.planned = sum(reps.planned),
    total.notrun = sum(reps.notrun),
    prop.notrun = total.notrun / total.planned,
    .groups = "drop"
  )
  
  print(nonconv.tab)
  return(nonconv.tab)
  
}



## END OF DOC