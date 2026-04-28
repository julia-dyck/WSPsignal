#'
#' @noRd

# for the Bayesian tests 

eval.rank_auc_b = function(perf_b){
  
  
  tab.df = perf_b
  
  # 1. summarize and rank fit & test specifications
  
  # ranking under correct prior belief
  tab.ranked <- tab.df %>%
    dplyr::filter(dist.prior.to.truth == "correct specification") %>% # only for correct prior spec for now
    dplyr::group_by(tte.dist, prior.dist, prior.sd, post.ci.type, cred.level, sensitivity.option) %>% # group by fit&test specifications
    dplyr::summarise(AUC = mean(auc), FPR = mean(fpr), TPR = mean(tpr), FNR = mean(fnr), TNR = mean(tnr), 
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(AUC)) # ranking by AUC
  
  
  # 2. filter for best fit & test specification and investigate effect of scenario parameters 
  
  opti.tte.dist = tab.ranked$tte.dist[1]
  opti.prior.dist = tab.ranked$prior.dist[1]
  opti.prior.sd = tab.ranked$prior.sd[1]
  opti.post.ci.type = tab.ranked$post.ci.type[1]
  opti.cred.level = tab.ranked$cred.level[1]
  opti.sensitivity.option = tab.ranked$sensitivity.option[1]

  
  # optimal fit & test specification according to ranking under correct prior belief
  tab.opti = tab.df %>% filter(tte.dist == opti.tte.dist, 
                               prior.dist == opti.prior.dist,
                               prior.sd == opti.prior.sd,
                               post.ci.type == opti.post.ci.type,
                               cred.level == opti.cred.level,
                               sensitivity.option == opti.sensitivity.option,
                               dist.prior.to.truth == "correct specification")
  
  # Effect of N
  tab.opti.N = tab.opti %>% 
    dplyr::filter(dist.prior.to.truth == "correct specification") %>%
    group_by(N) %>% 
    dplyr::summarise(AUC = mean(auc), FPR = mean(fpr), TPR = mean(tpr), FNR = mean(fnr), TNR = mean(tnr), 
                     .groups = "drop")
  
  # Effect of br
  tab.opti.br = tab.opti %>% 
    dplyr::filter(dist.prior.to.truth == "correct specification") %>%
    group_by(br) %>% 
    dplyr::summarise(AUC = mean(auc), FPR = mean(fpr), TPR = mean(tpr), FNR = mean(fnr), TNR = mean(tnr), 
                     .groups = "drop")
  
  # Effect of adr.rate
  tab.opti.adr.rate = tab.opti %>% 
    dplyr::filter(dist.prior.to.truth == "correct specification") %>%
    group_by(adr.rate) %>% 
    dplyr::summarise(AUC = mean(auc), FPR = mean(fpr), TPR = mean(tpr), FNR = mean(fnr), TNR = mean(tnr), 
                     .groups = "drop")
  
  # Effect of adr.when
  tab.opti.adr.when = tab.opti %>% 
    dplyr::filter(dist.prior.to.truth == "correct specification") %>%
    group_by(adr.when) %>% 
    dplyr::summarise(AUC = mean(auc), FPR = mean(fpr), TPR = mean(tpr), FNR = mean(fnr), TNR = mean(tnr), 
                     .groups = "drop")
  
  # Effect of adr.relsd
  tab.opti.adr.relsd = tab.opti %>% 
    dplyr::filter(dist.prior.to.truth == "correct specification") %>%
    group_by(adr.relsd) %>% 
    dplyr::summarise(AUC = mean(auc), FPR = mean(fpr), TPR = mean(tpr), FNR = mean(fnr), TNR = mean(tnr), 
                     .groups = "drop")
  
  # 3. Robustness wrt prior specification
  
  # Effect of distance of prior belief to true adr.when
  tab.robu.dist.prior.to.truth = tab.df %>% 
    filter(tte.dist == opti.tte.dist, 
           prior.dist == opti.prior.dist,
           post.ci.type == opti.post.ci.type,
           cred.level == opti.cred.level,
           sensitivity.option == opti.sensitivity.option) %>% 
    group_by(dist.prior.to.truth) %>% 
    dplyr::summarise(AUC = mean(auc), FPR = mean(fpr), TPR = mean(tpr), FNR = mean(fnr), TNR = mean(tnr), 
                     .groups = "drop")
  
  out = list(ranking = tab.ranked,
             effect.of.N = tab.opti.N,
             effect.of.br = tab.opti.br,
             effect.of.adr.rate = tab.opti.adr.rate,
             effect.of.adr.when = tab.opti.adr.when,
             effect.of.adr.relsd = tab.opti.adr.relsd,
             effect.of.dist.prior.to.truth = tab.robu.dist.prior.to.truth
  )
  return(out)
  
}

## END OF DOC
