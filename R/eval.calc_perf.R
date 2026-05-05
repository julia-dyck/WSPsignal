#' Calculate performance metrics for WSP test configurations in simulation
#'
#' Computes performance metrics for all the specified Bayesian and  
#' frequentist Weibull Shape Parameter (BWSP, BWSP) test configurations across simulated scenarios.
#' The output provides the base for a ranking of tests (see \code{\link{eval.rank_auc}}).
#'
#' @param pc_list list of simulation parameters generated with \code{\link{sim.setup_sim_pars}}
#'
#' @return A data frame containing one row per ADR–positive scenario, WSP model and
#' test configuration, and corresponding performance measurements in additional columns, namely 
#' the \eqn{auc}, \eqn{fpr}, \eqn{tpr}, \eqn{fnr} and \eqn{tnr}.
#' 
#' Scenarios with incomplete amount of simulations return `NA` for performance metrics.
#' FWSP tests return `NA` for scenario/model characteristics only relevant for 
#' BWSP test specification.
#' 
#' 
#'
#' @details 
#' 
#' 
#' Based the merged simulation results (obtained with 
#' \code{\link{sim.merge_results}}), the function performs WSP tests for all 
#' specified test setups - for Bayesian tests depending on the combination of posterior 
#' credibility interval (CI) type, credibility level and the sensitivity option 
#' (see \code{\link{bwsp_test}} and \code{pc_list$test}); for frequentist tests 
#' depending on the credibility level (see \code{\link{fwsp_test}} and pc_list$input$cred.level).
#' 
#' Given binary test results the function calculates
#' for each subset representing one scenario and test combination the 
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
#' among simulation repetitions, as well, as
#'
#' \itemize{
#' \item Area under the ROC curve:
#'   
#' The AUC is the area under the receiver operating characteristic (ROC) 
#' graph \insertCite{fawcett2004}{WSPsignal}. 
#' Here, we use the ROC curve with one threshold based on equal numbers
#' of ADR-positive and control 
#' scenarios. It is computed using the \code{\link[ROCR]{performance}} function. 
#' 
#' }
#' 
#' 
#' @references 
#' \insertAllCited{}
#'
#' @export
#' 

eval.calc_perf = function(pc_list){
  require(dplyr)
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