#'
#' @noRd

## eval.calc_perf for frequentist tests (inner fct)



eval.calc_perf_f = function(pc_list) {
  # require(dplyr)
  
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
  # fct body -------------------------------------------------------------------
  
  # Load results
  if (!exists("res_f")) {
    tryCatch({
      load(paste0(pc_list$add$resultpath, "/res_f.RData"))
      message("res_f.RData successfully loaded")
    }, error = function(cond) {
      sim.merge_results(pc_list, save = TRUE, bayes = FALSE)
      load(paste0(pc_list$add$resultpath, "/res_f.RData"))
      message("Batches merged and loaded")
    })
  } else {
    message("Object `res_f` currently loaded in environment is used to calculate performance measures.")
  }
  
  # 0. -------------------------------------------------------------------------
  #### reshape and remove NA test results
  
  pc_cols = names(res_f)[1:10]
  
  fwsp_w_cols   = grep("^fwsp_w_", names(res_f), value = TRUE)
  fwsp_dw_cols  = grep("^fwsp_dw_", names(res_f), value = TRUE)
  fwsp_pgw_cols = grep("^fwsp_pgw_", names(res_f), value = TRUE)
  
  # build clean subtables
  res_f_w = res_f %>%
    dplyr::filter(tte.dist == "w") %>%
    dplyr::select(dplyr::all_of(c(pc_cols, fwsp_w_cols))) %>%
    rename_with(~ sub("^fwsp_w_", "fwsp_", .x), dplyr::starts_with("fwsp_w"))
  
  res_f_dw = res_f %>%
    dplyr::filter(tte.dist == "dw") %>%
    dplyr::select(dplyr::all_of(c(pc_cols, fwsp_dw_cols))) %>%
    rename_with(~ sub("^fwsp_dw_", "fwsp_", .x), dplyr::starts_with("fwsp_dw"))
  
  res_f_pgw = res_f %>%
    dplyr::filter(tte.dist == "pgw") %>%
    dplyr::select(dplyr::all_of(c(pc_cols, fwsp_pgw_cols))) %>%
    rename_with(~ sub("^fwsp_pgw_", "fwsp_", .x), dplyr::starts_with("fwsp_pgw"))
  
  # combine
  res_f_long = dplyr::bind_rows(res_f_w, res_f_dw, res_f_pgw)
  
  
  # 1. -------------------------------------------------------------------------
  #### add label for true adr status
  res_f_long$lab = ifelse(res_f_long$adr.rate > 0, 1, 0)
  res.ext = res_f_long
  
  
  # 2. ------------------------------------------------------------------------- 
  #### calculate AUC for each simulation scenario (= one row of pc_list$pc_table)
  pc.pos = dplyr::filter(pc_list$pc_table, adr.rate > 0)
  pc.pos = unique(pc.pos[, -c(8:10)])  # Drop prior info (frequentist)
  
  fwsp_cols = grep("^fwsp_", names(res.ext), value = TRUE)
  nr.combined.tests = length(fwsp_cols)
  
  aucs = fprs = tprs = fnrs = tnrs = 
    matrix(NA, nrow = nrow(pc.pos), ncol = nr.combined.tests)
  colnames(aucs) = sub("^fwsp", "auc", fwsp_cols)
  colnames(fprs) = sub("^fwsp", "fpr", fwsp_cols)
  colnames(tprs) = sub("^fwsp", "tpr", fwsp_cols)
  colnames(fnrs) = sub("^fwsp", "fnr", fwsp_cols)
  colnames(tnrs) = sub("^fwsp", "tnr", fwsp_cols)
  
  run.reps = c()
  
  for (i in 1:nrow(pc.pos)) {
    scenario <- pc.pos[i, ]
    
    res.test = res.ext %>%
      dplyr::filter((adr.rate == 0 | adr.rate == scenario$adr.rate),
             (is.na(adr.when) | adr.when == scenario$adr.when),
             N == scenario$N,
             br == scenario$br,
             (is.na(adr.relsd) | adr.relsd == scenario$adr.relsd),
             tte.dist == scenario$tte.dist)
    
    run.reps[i] = nrow(res.test)
    
    labels = matrix(res.test$lab, nrow = run.reps[i], ncol = nr.combined.tests, byrow = FALSE)
    predictions = as.matrix(res.test[, fwsp_cols])
    
    # Calculate AUCs
    pred.obj = ROCR::prediction(predictions, labels)
    aucs[i, ] = ROCR::performance(pred.obj, "auc")@y.values %>% unlist()
    
    # Calculate tp, fp, tn, fn
    
    tp = colSums(predictions == 1 & labels == 1)
    fp = colSums(predictions == 1 & labels == 0)
    tn = colSums(predictions == 0 & labels == 0)
    fn = colSums(predictions == 0 & labels == 1)
    
    tprs[i, ] = tp / (tp + fn)
    fprs[i, ] = fp / (fp + tn)
    fnrs[i, ] = fn / (tp + fn)
    tnrs[i, ] = tn / (fp + tn)
  }
  
  fprs = cbind(pc.pos, fprs)
  tprs = cbind(pc.pos, tprs)
  fnrs = cbind(pc.pos, fnrs)
  tnrs = cbind(pc.pos, tnrs)
  aucs = cbind(pc.pos, aucs)
  
  
  # reshape performance measure matrices to long format
  
  fpr_cols <- grep("^fpr_", names(fprs), value = TRUE)
  # Reshape
  fprs_long <- reshape(
    fprs[, c(setdiff(names(fprs), fpr_cols), fpr_cols)],
    varying = fpr_cols,
    v.names = "fpr",
    timevar = "test_spec",
    times = fpr_cols,
    direction = "long"
  )
  rownames(fprs_long) <- NULL
  # Extract cred.level from test_spec and save in its own column
  fprs_long$cred.level <- as.numeric(sub(".*_(\\d+\\.?\\d*)$", "\\1", fprs_long$test_spec))
  
  tpr_cols <- grep("^tpr_", names(tprs), value = TRUE)
  # Reshape
  tprs_long <- reshape(
    tprs[, c(setdiff(names(tprs), tpr_cols), tpr_cols)],
    varying = tpr_cols,
    v.names = "tpr",
    timevar = "test_spec",
    times = tpr_cols,
    direction = "long"
  )
  rownames(tprs_long) <- NULL
  # Extract cred.level from test_spec and save in its own column
  tprs_long$cred.level <- as.numeric(sub(".*_(\\d+\\.?\\d*)$", "\\1", tprs_long$test_spec))
  
  fnr_cols <- grep("^fnr_", names(fnrs), value = TRUE)
  # Reshape
  fnrs_long <- reshape(
    fnrs[, c(setdiff(names(fnrs), fnr_cols), fnr_cols)],
    varying = fnr_cols,
    v.names = "fnr",
    timevar = "test_spec",
    times = fnr_cols,
    direction = "long"
  )
  rownames(fnrs_long) <- NULL
  # Extract cred.level from test_spec and save in its own column
  fnrs_long$cred.level <- as.numeric(sub(".*_(\\d+\\.?\\d*)$", "\\1", fnrs_long$test_spec))
  
  tnr_cols <- grep("^tnr_", names(tnrs), value = TRUE)
  # Reshape
  tnrs_long <- reshape(
    tnrs[, c(setdiff(names(tnrs), tnr_cols), tnr_cols)],
    varying = tnr_cols,
    v.names = "tnr",
    timevar = "test_spec",
    times = tnr_cols,
    direction = "long"
  )
  rownames(tnrs_long) <- NULL
  # Extract cred.level from test_spec and save in its own column
  tnrs_long$cred.level <- as.numeric(sub(".*_(\\d+\\.?\\d*)$", "\\1", tnrs_long$test_spec))
  
  # reshape aucs to long format
  auc_cols <- grep("^auc_", names(aucs), value = TRUE)
  # Reshape
  aucs_long <- reshape(
    aucs[, c(setdiff(names(aucs), auc_cols), auc_cols)],
    varying = auc_cols,
    v.names = "auc",
    timevar = "test_spec",
    times = auc_cols,
    direction = "long"
  )
  rownames(aucs_long) <- NULL
  # Extract cred.level from test_spec and save in its own column
  aucs_long$cred.level <- as.numeric(sub(".*_(\\d+\\.?\\d*)$", "\\1", aucs_long$test_spec))
  
  # merge performance measure matrices fprs, tprs, auc to one
  # merge fprs and tprs
  pm_long = dplyr::left_join(fprs_long, tprs_long, by = c(
    "N", "br", "adr.rate", "adr.when", "adr.relsd", "study.period", 
    "tte.dist", "cred.level"
  ))
  # merge fprs, tprs and fnrs
  pm_long = dplyr::left_join(pm_long, fnrs_long, by = c(
    "N", "br", "adr.rate", "adr.when", "adr.relsd", "study.period", 
    "tte.dist", "cred.level"
  ))
  # merge fprs, tprs, fnrs and tnrs
  pm_long = dplyr::left_join(pm_long, tnrs_long, by = c(
    "N", "br", "adr.rate", "adr.when", "adr.relsd", "study.period", 
    "tte.dist", "cred.level"
  ))
  # merge fprs, tprs, fnrs, tnrs and aucs
  # (this is the final result)
  pm_long = dplyr::left_join(pm_long, aucs_long, by = c(
    "N", "br", "adr.rate", "adr.when", "adr.relsd", "study.period", 
    "tte.dist", "cred.level"
  ))
  
  # select relevant columns:
  pm_long <- pm_long[, c(
    "N", "br", "adr.rate", "adr.when", "adr.relsd", "study.period", 
    "tte.dist", "cred.level", 
    "auc", "fpr", "tpr", "fnr", "tnr"
  )]
  
  # add missing columns to match required format (relevant when "b" not in est.approach)
  pm_long <- pm_long %>%
    mutate(
      prior.dist = NA_character_,
      prior.belief = NA_character_,
      prior.sd = NA_real_,
      dist.prior.to.truth = NA_character_,
      post.ci.type = NA_character_,
      sensitivity.option = NA_real_
    )
  pm_long = pm_long[, c(
    "N", "br", "adr.rate", "adr.when", "adr.relsd", "study.period", "tte.dist",
    "prior.dist", "prior.belief", "prior.sd", "dist.prior.to.truth",
    "post.ci.type", "cred.level", "sensitivity.option", 
    "auc", "fpr", "tpr", "fnr", "tnr"
  )]
  
  return(pm_long)
  
}



## END OF DOC