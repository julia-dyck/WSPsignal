#' Compute performance metrics for WSP test configurations
#'
#' Computes performance metrics for all the specified Weibull Shape Parameter (WSP) 
#' test configurations across simulated scenarios.
#' The output provides the base for a ranking of tests (see \code{\link{eval.rank_auc}}).
#'
#' @param pc_list list of simulation parameters generated with \code{\link{sim.setup_sim_pars}}
#'
#' @return A data frame containing one row per ADR–positive scenario, WSP model and
#' test configuration, and corresponding performance measurements in additional columns, namely 
#' the \eqn{auc}, \eqn{fpr}, \eqn{tpr}, \eqn{fnr} and \eqn{tnr}.
#' 
#' Scenarios with incomplete number of repetitions return `NA` for performance metrics.
#' Frequentist WSP tests return `NA` for scenario/model characteristics that are only relevant for 
#' Bayesian WSP test specification.
#' 
#' 
#'
#' @details 
#' Based on the merged simulation results obtained with 
#' \code{\link{sim.merge_results}}, the function performs WSP tests for all 
#' specified model and test configurations. Bayesian WSP tests depend on the 
#' combination of time-to-event (tte) distribution, prior model specification,
#' posterior credibility interval (CI) type, credibility level and sensitivity 
#' option (see \code{\link{bwsp_test}}). Frequentist WSP tests depend on the 
#' tte distribution and confidence level (see \code{\link{fwsp_test}}).
#' 
#' Given binary test results the function computes the
#' following performance measures:
#'   
#' \itemize{
#' 
#' \item False positive rate:
#' \deqn{ fpr = \frac{FP}{FP + TN} }
#'
#' \item True positive rate (sensitivity, recall):
#' \deqn{ tpr = \frac{TP}{TP + FN} }
#'
#' \item False negative rate:
#' \deqn{ fnr = \frac{FN}{TP + FN} }
#'
#' \item True negative rate (specificity):
#' \deqn{ tnr = \frac{TN}{FP + TN} }
#' }
#' with \eqn{FP} being the number of false positive cases, \eqn{TN} the number of true negative
#' cases, \eqn{TP} the number of true positive cases and \eqn{FN} the number of false negative cases
#' among simulation repetitions, as well as
#'
#' \itemize{
#' \item Area under the ROC curve (AUC):
#'   
#' The AUC is the area under the receiver operating characteristic (ROC) 
#' graph \insertCite{fawcett2004}{WSPsignal}. 
#' Here, the ROC curve with one threshold based on equal numbers
#' of ADR-positive and control 
#' scenarios is computed using the \code{\link[ROCR]{performance}} function. 
#' 
#' }
#' 
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
#'
#' #### compute performance metrics based on merged simulation results
#' 
#' perf = eval.calc_perf(pc_list)
#' head(perf, 10)
#' 
#' @references 
#' \insertAllCited{}
#'
#' @export
#' 

eval.calc_perf = function(pc_list){
  
  est.approach = pc_list$input$est.approach
  
  # calc performance measures for all Bayesian test types
  if("b" %in% est.approach){
    out_b = eval.calc_perf_b(pc_list)
    out_b_ext = cbind(test.type = rep("bwsp", nrow(out_b)), out_b)
  }
  else{out_b_ext = NULL}
  
  # calc performance measures for all frequentist test types
  if("f" %in% est.approach){
    out_f = eval.calc_perf_f(pc_list)
    out_f_ext = cbind(test.type = rep("fwsp", nrow(out_f)), out_f)
  }
  else{out_f_ext = NULL}
  
  out = dplyr::bind_rows(out_b_ext, out_f_ext)
  
  return(out)
}


## END OF DOC