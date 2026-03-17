#' Evaluate effective sample sizes by tte and prior distribution
#'
#' Only for Bayesian estimation approach.
#' Summarizes and visualizes effective sample sizes of the stan models fitted during 
#' the simulation study, grouped by time-to-event (tte) and prior distribution types. 
#' This helps assess which of the tte and prior distribution choices is suitable 
#' for HDI+ROPE testing (along with other diagnostics such as 
#' \code{\link{eval.non_conv_cases}} and \code{\link{eval.execution_times}}).
#'
#' @param pc_list a list containing simulation parameters (see \code{\link{sim.setup_sim_pars}})
#' @param threshold numeric threshold for effective sample size acceptable for 
#' HDI+ROPE testing (10000 by default as recommended by \insertCite{kruschke2015;textual}{WSPsignal})
#'
#' @return A list with summary statistics (`$summary`), a ggplot2 object (`$plot`), 
#' and the data (`$df`) on which summary and plot are based.
#' 
#' @seealso \code{\link{eval.non_conv_cases}}, \code{\link{eval.execution_times}}
#' 
#' @references
#'   \insertAllCited{}
#'
#' @export



eval.eff_sample_sizes = function(pc_list, threshold = 10000){
  
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
  
  if(!("b" %in% pc_list$input$est.approach)){
    stop("Functions eval.non_conv_cases(), eval.execution_times() and eval.eff_sample_sizes() only applicable if Bayesian estimation is included in simulation study.\n")
  }
  
  # argument check for threshold
  if (!is.numeric(threshold) || length(threshold) != 1 || is.na(threshold) || threshold < 0) {
    stop("Argument 'threshold' must be a single non-negative numeric value.\n")
  }
  
  
  ## fct body ------------------------------------------------------------------
  
  if (!exists("res_b")) { 
    # obtain res table
    tryCatch({
      load(paste0(pc_list$add$resultpath, "/res_b.RData"))
      message("res_b.RData successfully loaded")
    }, error = function(cond) {
      sim.merge_results(pc_list, save = T)
      load(paste0(pc_list$add$resultpath, "/res_b.RData"))
      message(" batches merged and loaded")
    })
  }
  else{
    message("Object `res_b` loaded in current environment is used to extract effective sample sizes.")
  }
  
  # NULL values ~> NA
  res_b$nu.po.n_eff <- lapply(res_b$nu.po.n_eff, function(x) if (length(x) == 0) NA else x)
  res_b$ga.po.n_eff <- lapply(res_b$ga.po.n_eff, function(x) if (length(x) == 0) NA else x)
  
  # Select and rename relevant cols for nu
  n_eff.nu <- dplyr::select(res_b, tte.dist, prior.dist, nu.po.n_eff)
  n_eff.nu <- dplyr::rename(n_eff.nu, n_eff = nu.po.n_eff)
  #n_eff.nu$n_eff = ifelse(is.null(n_eff.nu$n_eff), NA, n_eff.nu$n_eff)
  n_eff.nu$parameter <- "shape1"
  
  # Select and rename relevant cols for ga
  n_eff.ga <- dplyr::select(res_b, tte.dist, prior.dist, ga.po.n_eff)
  n_eff.ga <- dplyr::rename(n_eff.ga, n_eff = ga.po.n_eff)
  # Replace NULLs with NAs in a list-column ### HIER WEITER
  res_b$ga.po.n_eff <- lapply(res_b$ga.po.n_eff, function(x) if (is.null(x)) NA else x)
  n_eff.ga$parameter <- "shape2"
  
  # Combine into long format
  n_eff.long <- rbind(n_eff.nu, n_eff.ga)

  # Reorder columns
  n_eff.long <- n_eff.long[, c("tte.dist", "prior.dist", "parameter", "n_eff")]
  
  # Unlist all columns to remove any list-columns
  for (col in names(n_eff.long)) {
    n_eff.long[[col]] <- unlist(n_eff.long[[col]])
  }
  n_eff.long$n_eff <- as.numeric(n_eff.long$n_eff)
  
  # execution time dist per tte.dist and prior.dist
  ## WHAT DO WE REALLY NEED; MEAN; NO OF SIMS WITH NEFF < 10000; LOOK AT MANUSCRIPT
  
  n_eff.summaries <- suppressWarnings(
    dplyr::summarise(
      dplyr::group_by(n_eff.long, tte.dist, prior.dist, parameter),
      # min = min(n_eff, na.rm = TRUE),
      # first_qu = quantile(n_eff, 0.25, na.rm = TRUE),
      # median = median(n_eff, na.rm = TRUE),
      # mean = mean(n_eff, na.rm = TRUE),
      # third_qu = quantile(n_eff, 0.75, na.rm = TRUE),
      # max = max(n_eff, na.rm = TRUE),
      prob.above.thr = sum(n_eff > threshold, na.rm = TRUE) / dplyr::n(),
      .groups = "drop"
    )
  )
  
  #return(n_eff.summaries)
  
  # p = ggplot2::ggplot(n_eff.long, ggplot2::aes(x = prior.dist, y = n_eff, fill = parameter)) +
  #   ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.8)) +
  #   ggplot2::facet_wrap(~ tte.dist, ncol = 1) +
  #   ggplot2::labs(x = "prior.dist", y = "n_eff grouped by tte.dist", fill = "parameter") +
  #   ggplot2::theme_minimal() +
  #   ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
  # 
  
  n_eff.long$group <- interaction(n_eff.long$tte.dist, n_eff.long$prior.dist, sep = " - ")
  
  p <- ggplot2::ggplot(n_eff.long, ggplot2::aes(x = group, y = n_eff, fill = parameter)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.8), width = 0.5) +
    ggplot2::geom_hline(yintercept = 10000, linetype = "dashed", color = "black", size = 1) +
    ggplot2::labs(
      x = "tte.dist - prior.dist combination",
      y = "n_eff",
      title = "Effective sample sizes",
      fill = "Parameter"
    ) +
    ggplot2::theme_minimal() +
    
    ggplot2::theme(legend.position = "top")
  
  return(list(summary = n_eff.summaries,     # for overview
              plot = p,                     # for option to take plot as is and manipulate if further
              df = n_eff.long  # for option to plot manually in preferred format
  ))
}





## END OF DOC
