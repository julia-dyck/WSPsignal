#' Evaluate execution times
#'
#' Only for Bayesian estimation approach.
#' Summarizes execution times of the models fitted optionally
#' grouped by one or multiple model specifications. 
#' This helps assess which of the modelling choices is suitable 
#' for CI+ROPE testing (along with other diagnostics such as 
#' \code{\link{eval.non_conv_cases}} and \code{\link{eval.eff_sample_sizes}}).
#'
#' @param pc_list list of simulation parameters generated with \code{\link{sim.setup_sim_pars}}
#' @param group.by character vector specifying grouping variables; must be a 
#' subset of \code{c("tte.dist", "prior.dist", "prior.sd")}
#' @param verbose logical; if \code{TRUE} (default), summary statistics
#' and plot are printed to the console
#' 
#' @details
#' Calculations are based on the stored result file obtained with
#' \code{\link{sim.merge_results}}. 
#' 
#' @return A list with summary statistics (`$summary`), a ggplot2 object (`$plot`), 
#' and the data (`$df`) on which summary and plot are based.
#' 
#' @seealso \code{\link{eval.non_conv_cases}}, \code{\link{eval.eff_sample_sizes}}
#'
#' @examples
#' # The package ships with a small precomputed toy simulation study with small
#' # numbers of repetitions (rep) and posterior sample sizes (stanmod.iter, stanmod.warmup). 
#' # The simulation parameters below match the shipped example results.
#'
#' #### prep:
#' toy_path <- system.file("extdata", "toysim", package = "WSPsignal")
#'
#' # setup prior template
#' fp_list <- sim.priors_template(tte.dist = c("w", "pgw"), prior.sds = 10)
#'
#' # fill in prior template with prior means
#' fp_list$w[,2] <- c(1, 1, 180, 300)
#' fp_list$w[,3] <- c(1, 0.207, 1, 4)
#'
#' fp_list$pgw[,2] <- c(1, 1, 20, 300)
#' fp_list$pgw[,3] <- c(1, 0.207, 5.5, 4)
#' fp_list$pgw[,4] <- c(1, 1, 14, 1)
#'
#' # recreate simulation settings used for the toy example
#' pc_list <- sim.setup_sim_pars(
#'   N = 500,
#'   br = 0.1,
#'   adr.rate = c(0, 1),
#'   adr.relsd = 0.05,
#'   study.period = 365,
#'   est.approach = c("f", "b"),
#'   tte.dist = c("w", "pgw"),
#'   prior.dist = "ll",
#'   fitpars.list = fp_list,
#'   post.ci.type = c("ETI", "HDI"),
#'   cred.level = seq(0.5, 0.9, by = 0.05),
#'   sensitivity.option = 1:3,
#'   reps = 6,
#'   batch.size = 3,
#'   resultpath = toy_path,
#'   stanmod.iter = 1100,
#'   stanmod.warmup = 100
#' )
#'
#' #### summarize execution times only by tte distribution (only varying factor)
#' ex_times = eval.execution_times(pc_list, group.by = "tte.dist")
#' 
#' 
#' @export


eval.execution_times = function(pc_list, 
                                group.by = c("tte.dist", "prior.dist", "prior.sd"),
                                verbose = TRUE){
  
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
 
  if (!exists("res_b")) { 
  # obtain res_b table
  tryCatch({
    load(paste0(pc_list$add$resultpath, "/res_b.RData"))
    message("res_b.RData successfully loaded")
  }, error = function(cond) {
    sim.merge_results(pc_list, save = TRUE)
    load(paste0(pc_list$add$resultpath, "/res_b.RData"))
    print(" batches merged and loaded")
  })
  }
  else{
    message("Object `res_b` loaded in current environment is used to extract execution times.")
  }
  
  # select relevant variables
  time.df = res_b[, c("tte.dist", "prior.dist", "prior.sd", "run.min")]
  # adjust format
  time.df$tte.dist <- as.factor(unlist(time.df$tte.dist))
  time.df$prior.dist <- as.factor(unlist(time.df$prior.dist))
  time.df$prior.sd <- as.factor(unlist(time.df$prior.sd))
  time.df$run.min = as.numeric(unlist(time.df$run.min)) 
  
  # create grouping variable (same as for plot)
  time.df$group <- do.call(interaction, c(time.df[group.by], sep = " - "))
  
  # summarise using that grouping
  time.summaries <- dplyr::summarise(
    dplyr::group_by(time.df, group),
    min = min(run.min, na.rm = TRUE),
    first_qu = stats::quantile(run.min, 0.25, na.rm = TRUE),
    median = stats::median(run.min, na.rm = TRUE),
    mean = mean(run.min, na.rm = TRUE),
    third_qu = stats::quantile(run.min, 0.75, na.rm = TRUE),
    max = max(run.min, na.rm = TRUE),
    .groups = "drop"
  )
  
  # group variable for plot
  time.df$group <- do.call(interaction, c(time.df[group.by], sep = " - "))
  
  # plot
  p = ggplot2::ggplot(time.df, ggplot2::aes(x = group, y = run.min)) +
    ggplot2::geom_boxplot(width = 0.5, fill = "lightgrey") +
    ggplot2::labs(
      x = paste(group.by, collapse = " - "),
      y = "Execution time (min)",
      title = "Execution time in minutes"
    ) +
    ggplot2::theme_minimal()
  
  if(verbose == TRUE){
    message("\nSummary:")
    print(time.summaries)
    print(p) 
  }
  
  
  return(list(summary = time.summaries,     # for overview
              plot = p,                     # for option to take plot as is and manipulate if further
              df = time.df  # for option to plot manually in preferred format
              ))
}



## END OF DOC
