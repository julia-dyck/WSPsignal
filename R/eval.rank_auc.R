#' Ranking of WSP test configurations
#'
#' Ranks all model and test specifications grouped by simulation scenarios in terms
#' of the corresponding area under the curve (AUC) value for Weibull Shape Parameter 
#' (WSP) tests.
#'
#' @param perf data frame containing performance results for WSP tests returned 
#' by \code{\link{eval.calc_perf}}
#' @param test.type.subset character to filter for Bayesian and frequentist
#' WSP test types to be considered in the ranking; 
#' must be a subset of \code{c("bwsp", "fwsp")}
#'
#' @param tte.dist.subset character to filter for the time-to-event (tte)
#'   distributions considered in the ranking, must be a subset of
#'   \code{c("w", "dw", "pgw")}
#'
#' @param prior.dist.subset character to filter for the prior distribution
#'   (relevant only for BWSP tests), must be a subset of
#'   \code{c("fg", "fl", "gg", "ll")}
#'   
#' @param prior.sd.subset numeric to filter for the prior standard deviation
#'   (relevant only for BWSP tests), must be a subset of prior.sds considered 
#'   in the simulation
#'
#' @return A list containing 
#' \itemize{
#' \item \code{$rank.tab}: Ranking of fit and WSP test specifications according 
#' to AUC averaged over all sample scenarios (for BWSP tests given a correct 
#' specification of prior belief)
#' \item \code{$effect.of.N}: Effect of sample size on AUC for the optimal fit 
#' and WSP test (for BWSP tests given a correct specification of prior belief)
#' \item \code{$effect.of.br}: Effect of background rate on AUC for the optimal fit
#' and WSP test (for BWSP tests given a correct specification of prior belief)
#' \item \code{$effect.of.adr.rate}: Effect of ADR rate on AUC for the optimal fit
#' and WSP test (for BWSP tests given a correct specification of prior belief)
#' \item \code{$effect.of.adr.when}: Effect of true expected event times on AUC 
#' for the optimal fit and WSP test (for BWSP tests given a correct 
#' specification of prior belief)
#' \item \code{$effect.of.adr.relsd}: Effect of relative standard deviation of event 
#' time on AUC for the optimal fit and WSP test (for BWSP tests given a correct 
#' specification of prior belief)
#' \item \code{$effect.of.dist.prior.to.truth}: Effect of distance of prior belief 
#' to true \code{adr.when} on AUC for the optimal fit and WSP test (only BWSP)
#' }
#' 
#' @details For definitions of the performance metrics \code{AUC, FPR, TPR, FNR} 
#' and \code{TNR} returned in output, 
#' see the details section of \code{\link{eval.calc_perf}}.
#' 
#' The filter mechanism enables filtering for a subset of test specifications. 
#' This is helpful for example when tte distributions, prior distributions or 
#' an estimation approach are no longer under consideration for instance after 
#' inspecting the model diagnostics with \code{\link{eval.execution_times}}, 
#' \code{\link{eval.non_conv_cases}} and \code{\link{eval.eff_sample_sizes}}.
#' 
#' 
#' 
#' 
#' 
#' 
#' @examples
#' \dontrun{
#' # loading of performance metrics returned by eval.calc_perf function,
#' load(paste0(pc_list$add$resultpath, "/perf.RData"))
#' 
#' # ranking of all WSP tests considered in simulation setup
#' rank = eval.rank_auc(perf)
#' 
#' # ranking of subset of all Bayesian WSP tests considered in the simulation setup
#' rank_b = eval.rank_auc(perf, test.type.subset = "bwsp")
#'
#' # ranking of subset of all frequentist Weibull tests considered
#' # in the simulation setup
#' rank_w = eval.rank_auc(perf, tte.dist.subset = "w")
#' 
#' # ranking of subset of all Bayesian WSP tests with "ll"-prior considered in 
#' # the simulation setup
#' rank_b_ll = eval.rank_auc(perf, test.type.subset = "bwsp",
#'                        prior.dist.subset = "ll")
#'                        
#' }
#'
#' @export
#'




eval.rank_auc = function(perf, 
                         test.type.subset = c("bwsp", "fwsp"), 
                         tte.dist.subset = c("w", "dw", "pgw"), 
                         prior.dist.subset = c("fg", "fl", "gg", "ll"),
                         prior.sd.subset = NULL){

  
  # 0. argument checks ---------------------------------------------------------
  # for perf
  required_cols = c(
    "test.type", "N","br","adr.rate","adr.when","adr.relsd","study.period",
    "tte.dist","prior.dist","prior.belief", "prior.sd","dist.prior.to.truth",
    "post.ci.type","cred.level","sensitivity.option",
    "auc","fpr","tpr","fnr","tnr"
  )
  num_cols = c(
    "N","br","adr.rate","adr.when","adr.relsd","study.period",
    "prior.sd", "cred.level", "sensitivity.option","auc","fpr","tpr","fnr","tnr"
  )
  chr_cols <- c(
    "test.type", "tte.dist","prior.dist","prior.belief",
    "dist.prior.to.truth","post.ci.type"
  )
  if (
    !is.data.frame(perf) ||
    !all(required_cols %in% names(perf)) ||
    !all(vapply(perf[num_cols], is.numeric, logical(1))) ||
    !all(vapply(perf[chr_cols], function(x) is.character(x), logical(1)))
  ) {
    stop("Argument `perf` has wrong format. It must be the output of eval.calc_perf().")
  }
  # check: test.type.subset
  if (!all(test.type.subset %in% c("bwsp", "fwsp"))) {
    stop("Argument `test.type.subset` must be a subset of c('bwsp','fwsp').")
  }
  # check: tte.dist.subset
  if (!all(tte.dist.subset %in% c("w","dw","pgw"))) {
    stop("Argument `tte.dist.subset` must be a subset of c('w','dw','pgw').")
  }
  # check: prior.dist.subset
  if (!all(prior.dist.subset %in% c("fg","fl","gg","ll"))) {
    stop("Argument `prior.dist.subset` must be a subset of c('fg','fl','gg','ll').")
  }
  # warning if bwsp is filtered out, but prior.dist subset specified
  if (!("bwsp" %in% test.type.subset) && !identical(prior.dist.subset, c("fg","fl","gg","ll"))) {
    warning("No 'bwsp' tests are considered. Argument `prior.dist.subset` has no effect.")
  }
  
  # check: prior.sd.subset
  if (!is.null(prior.sd.subset) && (
      !is.numeric(prior.sd.subset) ||
      !all(prior.sd.subset %in% unique(stats::na.omit(perf$prior.sd)))
    )
  ) {
    stop(
      paste0(
        "Argument `prior.sd.subset` must be NULL or a subset of prior.sd values ",
        paste(sort(unique(stats::na.omit(perf$prior.sd))), collapse = ", "),
        "."
      )
    )
  }
  
  ## fct body ------------------------------------------------------------------
  # apply filter options
  perf = perf %>% dplyr::filter(tte.dist %in% tte.dist.subset) %>%
    dplyr::filter(prior.dist %in% c(prior.dist.subset,NA))
  
  if (!is.null(prior.sd.subset)) {
    perf = perf %>% dplyr::filter(prior.sd %in% c(prior.sd.subset, NA))
  }
  
  # apply ranking fcts
  if(("bwsp" %in% test.type.subset) & (sum(perf$test.type == "bwsp") > 0)){
    perf_b = perf %>% dplyr::filter(test.type == "bwsp")
    rlist_b = eval.rank_auc_b(perf_b)
    rtab_b = rlist_b$ranking
    rtab_b_ext = cbind(test.type = rep("bwsp", nrow(rtab_b)), rtab_b)
  }
  else{ # generate table, reduce to zero rows, to keep column structure
    rtab_b_ext = data.frame(test.type = "bwsp", tte.dist = "w", prior.dist = "ll",
                           post.ci.type = "hdi", cred.level = 0.9, sensitivity.option = 1,
                           AUC = 0.999, FPR = 0.99, TPR = 0.99, FNR = 0.99, TNR = 0.99)
    rtab_b_ext = rtab_b_ext[0,]
  } 
  
  if(("fwsp" %in% test.type.subset) & (sum(perf$test.type == "fwsp") > 0)){ 
    perf_f = perf %>% dplyr::filter(test.type == "fwsp") 
    rlist_f = eval.rank_auc_f(perf_f)
    rtab_f = rlist_f$ranking
    rtab_f_ext = cbind(test.type = rep("fwsp", nrow(rtab_f)), rtab_f)
  }
  else{ # generate table, reduce to zero rows, to keep column structure
    rtab_f_ext = data.frame(test.type = "fwsp", tte.dist = "w", cred.level = 0.9,
                            AUC = 0.999, FPR = 0.99, TPR = 0.99, FNR = 0.99, TNR = 0.99)
    rtab_f_ext = rtab_f_ext[0,]
  } 
  
  # combine bwsp & fwsp ranking tables
  if (nrow(rtab_b_ext) > 0 & nrow(rtab_f_ext) > 0) {
    rtab = dplyr::bind_rows(rtab_b_ext, rtab_f_ext) %>%
      dplyr::arrange(dplyr::desc(AUC))
  } else if (nrow(rtab_b_ext) > 0) {
    rtab = rtab_b_ext %>%
      dplyr::arrange(dplyr::desc(AUC))
  } else if (nrow(rtab_f_ext) > 0) {
    rtab = rtab_f_ext %>%
      dplyr::arrange(dplyr::desc(AUC))
  } else {
    rtab = rtab_b_ext  # empty with correct structure
  }
  
  if(nrow(rtab) == 0){
    # ie. when filter leads to empty ranking table
    out = list(rank.tab = rtab, 
               effect.of.N = NULL,
               effect.of.br = NULL,
               effect.of.adr.rate = NULL,
               effect.of.adr.when = NULL,
               effect.of.adr.relsd = NULL,
               effect.of.dist.prior.to.truth = NULL)
    warning("Filter mechanisms removed all test candidates from ranking.") 
  }
  else if(rtab$test.type[1] == "bwsp"){
    out = list(rank.tab = rtab, 
               effect.of.N = rlist_b$effect.of.N,
               effect.of.br = rlist_b$effect.of.br,
               effect.of.adr.rate = rlist_b$effect.of.adr.rate,
               effect.of.adr.when = rlist_b$effect.of.adr.when,
               effect.of.adr.relsd = rlist_b$effect.of.adr.relsd,
               effect.of.dist.prior.to.truth = rlist_b$effect.of.dist.prior.to.truth)
  }
  else if(rtab$test.type[1] == "fwsp"){
    out = list(rank.tab = rtab, 
               effect.of.N = rlist_f$effect.of.N,
               effect.of.br = rlist_f$effect.of.br,
               effect.of.adr.rate = rlist_f$effect.of.adr.rate,
               effect.of.adr.when = rlist_f$effect.of.adr.when,
               effect.of.adr.relsd = rlist_f$effect.of.adr.relsd,
               effect.of.dist.prior.to.truth = NULL)
  }
  out$rank.tab = out$rank.tab %>%
    dplyr::mutate(
      dplyr::across(c(AUC, FPR, TPR, FNR, TNR), ~ round(.x, 3))
    )
  print(out$rank.tab)
  return(out)
}







## END OF DOC
