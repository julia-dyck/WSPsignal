#'
#' @noRd


## eval.calc_perf for Bayesian tests (inner fct)



eval.calc_perf_b = function(pc_list){
  require(dplyr) # for the pipe operator
  
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
  
  # 0. -------------------------------------------------------------------------
  #### load res table
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
    message("Object `res_b` currently loaded in environment is used to calculate performance measures.")
  }
  
  # 1. -------------------------------------------------------------------------
  #### add label for true adr status
  res_b$lab = ifelse(res_b$adr.rate > 0, 1, 0) # 1 = ADR, 0 = no ADR)
  
  # 2. -------------------------------------------------------------------------
  #### add ropes for all bwsp tests to be performed
  
  rope.infos = dplyr::filter(do.call(dplyr::bind_rows, pc_list$fit), prior.belief == "none")
  
  # calculate ropes for each tte.dist & prior.dist combi
  ropes = do.call(rbind, lapply(1:nrow(rope.infos), function(i) {
    # extract row
    rope.infos.row = rope.infos[i, ]
    # apply eval.calc_rope
    ropes = eval.calc_rope(cred.levels = pc_list$input$cred.level, rope.infos.row = rope.infos.row)
    return(ropes)
  }))
  
  rope.infos$prior.belief <- NULL # to not mess up the merge
  # remove cols unnecessary? (at the same time, they do not hurt)
  rope.infos = cbind(rope.infos, ropes)
  
  res.ext = merge(res_b, rope.infos, by = c("tte.dist", "prior.dist"), all.x = TRUE) #merge
  
  # 3. -------------------------------------------------------------------------
  #### perform all bwsp tests and save binary test result per row and per test specification as new res cols
  
  # names of the lower and upper interval bounds you need
  cred.levels = pc_list$input$cred.level
  
  # posterior CI types
  types = tolower(pc_list$input$post.ci.type)
  
  # define the different sensitivity options
  options = pc_list$input$sensitivity.option
  
  
  # loop over sensitivity options
  for (opt in options) {
    
    # loop over types (hdi, eti)
    for (type in types) {
      
      # loop over credibility levels
      for (lev in cred.levels) {
        
        # create the column names dynamically for credregion
        cred_cols = c(
          paste0("nu.", type, lev, "l"),
          paste0("nu.", type, lev, "u"),
          paste0("ga.", type, lev, "l"),
          paste0("ga.", type, lev, "u")
        )
        
        # create the column names dynamically for nullregion
        null_cols = c(
          paste0("nu.", "rope", lev, "l"),
          paste0("nu.", "rope", lev, "u"),
          paste0("ga.", "rope", lev, "l"),
          paste0("ga.", "rope", lev, "u")
        )
        
        # define a new column name for the result
        bwsp_col = paste0("bwsp_", type, "_", lev, "_opt", opt)
        
        # apply the test rowwise
        res.ext[[bwsp_col]] = apply(
          res.ext[, c(cred_cols, null_cols, "tte.dist")], 
          MARGIN = 1,
          FUN = function(x) {
            # first 4 entries = credregion, next 4 entries = nullregion
            credregion = as.numeric(unlist(x[1:4]))
            nullregion = as.numeric(unlist(x[5:8]))
            tte.dist = x[9]
            
            if(tte.dist == "dw" || tte.dist == "pgw"){
              sim.bwsp_test(
                credregion = credregion,
                nullregion = nullregion,
                option = opt,
                tte.dist =  tte.dist
              )
            } else if(tte.dist == "w"){
              sim.bwsp_test(
                credregion = credregion[1:2],
                nullregion = nullregion[1:2],
                option = opt,
                tte.dist =  tte.dist
              )
            }
          }
        ) # end of apply
      } # end of loop over cred.levels
    } # end of loop over cred.levels
  } # end of loop over options
  
  
  # 4. -------------------------------------------------------------------------
  #### calculate AUC for each simulation scenario (= one row of pc_list$pc_table)
  
  ## control cases are matched to each ADR-positive scenario for AUC calc
  pc.pos = filter(pc_list$pc_table, adr.rate > 0) # only ADR-positive scenarios
  
  # number of tests
  nr.combined.tests = length(grep("^bwsp_", names(res.ext))) 
  
  # Identify all bwsp_test result columns
  bwsp_cols = grep("^bwsp_", names(res.ext), value = TRUE)
  
  # prep empty matrix for performance measure results
  # false positive rate
  fprs = matrix(rep(NA, nr.combined.tests*nrow(pc.pos)), ncol = nr.combined.tests)
  colnames(fprs) = sub("^bwsp", "fpr", bwsp_cols)
  # true positive rate
  tprs = matrix(rep(NA, nr.combined.tests*nrow(pc.pos)), ncol = nr.combined.tests)
  colnames(tprs) = sub("^bwsp", "tpr", bwsp_cols)
  # false negative rate
  fnrs = matrix(rep(NA, nr.combined.tests*nrow(pc.pos)), ncol = nr.combined.tests)
  colnames(fnrs) = sub("^bwsp", "fnr", bwsp_cols)
  # true negative rate
  tnrs = matrix(rep(NA, nr.combined.tests*nrow(pc.pos)), ncol = nr.combined.tests)
  colnames(tnrs) = sub("^bwsp", "tnr", bwsp_cols)
  # area under the receiver operating characteristic (ROC) curve
  aucs = matrix(rep(NA, nr.combined.tests*nrow(pc.pos)), ncol = nr.combined.tests)
  colnames(aucs) = sub("^bwsp", "auc", bwsp_cols)
  
  run.reps = c()
  # go through every ADR-positive scenario linked with control
  for(i in 1:nrow(pc.pos)){
    N_i = pc.pos$N[i]
    br_i = pc.pos$br[i]
    adr.rate_i = pc.pos$adr.rate[i]
    adr.when_i = pc.pos$adr.when[i]
    adr.relsd_i = pc.pos$adr.relsd[i]
    # no study.period, as supposed to be only one value
    tte.dist_i = pc.pos$tte.dist[i]
    prior.dist_i = pc.pos$prior.dist[i]
    prior.belief_i = pc.pos$prior.belief[i]
    
    res.test = res.ext %>%
      dplyr::filter((adr.rate == 0 | adr.rate == adr.rate_i),
                    (is.na(adr.when) | adr.when == adr.when_i),
                    N == N_i,
                    br == br_i,
                    (is.na(adr.relsd) | adr.relsd == adr.relsd_i),
                    tte.dist == tte.dist_i,
                    prior.dist == prior.dist_i,
                    prior.belief == prior.belief_i)
    
    run.reps[i] = nrow(res.test) # number of repetitions obtained for this scenario
    if(run.reps[i] == 2*pc_list$add$reps){
      
      # set up labels and predictions in a matrix
      labels = matrix(res.test$lab, nrow = run.reps[i], ncol = nr.combined.tests, byrow = F)
      predictions = data.frame(res.test)[,bwsp_cols] %>%
        as.matrix()
      
      # calculate tpr, fpr, tnr, fnr manually   
      
      # true positives
      tp = colSums(predictions == 1 & labels == 1)
      # false positives
      fp = colSums(predictions == 1 & labels == 0)
      # true negatives
      tn = colSums(predictions == 0 & labels == 0)
      # false negatives
      fn = colSums(predictions == 0 & labels == 1)
      
      fprs[i,] = fp / (fp + tn) # false positive rate
      tprs[i,] = tp / (tp + fn) # true positive rate
      fnrs[i,] = fn / (tp + fn) # false negative rate
      tnrs[i,] = tn / (fp + tn) # true negative rate
      
      
      # calculate AUCs
      pred.obj <- ROCR::prediction(predictions, labels) # creating prediction object
      
      aucs[i,] = ROCR::performance(pred.obj, "auc") %>%
        .@y.values %>%
        as.numeric()
    }
    else{
      fprs[i,] = rep(NA, nr.combined.tests)
      tprs[i,] = rep(NA, nr.combined.tests)
      fnrs[i,] = rep(NA, nr.combined.tests)
      tnrs[i,] = rep(NA, nr.combined.tests)
      aucs[i,] = rep(NA, nr.combined.tests)
    }
  }
  
  # 5. -------------------------------------------------------------------------
  #### add info and reshape table format
  
  # add scenario information to performance measure matrices
  
  ## inner fct
  # add description of deviance between prior belief and simulated truth
  prior.correctness = function(pc_row) {
    adr.when = as.numeric(as.character(pc_row[4]))
    prior.belief = as.character(pc_row[9])
    out = ifelse((adr.when == 0.25 && prior.belief == "beginning") ||
                   (adr.when == 0.5 && prior.belief == "middle") ||
                   (adr.when == 0.75 && prior.belief == "end"), 
                 "correct specification", 
                 ifelse((adr.when == 0.25 && prior.belief == "middle") ||
                          (adr.when == 0.5 && prior.belief == "beginning") ||
                          (adr.when == 0.5 && prior.belief == "end") ||
                          (adr.when == 0.75 && prior.belief == "middle"),
                        "one quarter off",
                        ifelse((adr.when == 0.25 && prior.belief == "end") ||
                                 (adr.when == 0.75 && prior.belief == "beginning"),
                               "two quarters off",
                               "no ADR assumed")))
    return(out)
  }
  
  # add deviance to pc.pos (for grouped mean calculation later)
  dist.prior.to.truth = apply(pc.pos, 1, prior.correctness) 
  pc.pos = cbind(pc.pos, dist.prior.to.truth)
  
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
  # Extract post.ci.type, cred.level, and sensitivity.option from test_spec
  fprs_long$post.ci.type = sub("^fpr_([^_]+)_.*", "\\1", fprs_long$test_spec)
  fprs_long$cred.level <- as.numeric(sub("^fpr_[^_]+_([0-9\\.]+)_.*", "\\1", fprs_long$test_spec))
  fprs_long$sensitivity.option <- as.integer(sub(".*opt([0-9]+)$", "\\1", fprs_long$test_spec))
  
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
  # Extract post.ci.type, cred.level, and sensitivity.option from test_spec
  tprs_long$post.ci.type = sub("^tpr_([^_]+)_.*", "\\1", tprs_long$test_spec)
  tprs_long$cred.level <- as.numeric(sub("^tpr_[^_]+_([0-9\\.]+)_.*", "\\1", tprs_long$test_spec))
  tprs_long$sensitivity.option <- as.integer(sub(".*opt([0-9]+)$", "\\1", tprs_long$test_spec))
  
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
  # Extract post.ci.type, cred.level, and sensitivity.option from test_spec
  fnrs_long$post.ci.type = sub("^fnr_([^_]+)_.*", "\\1", fnrs_long$test_spec)
  fnrs_long$cred.level <- as.numeric(sub("^fnr_[^_]+_([0-9\\.]+)_.*", "\\1", fnrs_long$test_spec))
  fnrs_long$sensitivity.option <- as.integer(sub(".*opt([0-9]+)$", "\\1", fnrs_long$test_spec))
  
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
  # Extract post.ci.type, cred.level, and sensitivity.option from test_spec
  tnrs_long$post.ci.type = sub("^tnr_([^_]+)_.*", "\\1", tnrs_long$test_spec)
  tnrs_long$cred.level <- as.numeric(sub("^tnr_[^_]+_([0-9\\.]+)_.*", "\\1", tnrs_long$test_spec))
  tnrs_long$sensitivity.option <- as.integer(sub(".*opt([0-9]+)$", "\\1", tnrs_long$test_spec))
  
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
  # Extract post.ci.type, cred.level, and sensitivity.option from test_spec
  aucs_long$post.ci.type = sub("^auc_([^_]+)_.*", "\\1", aucs_long$test_spec)
  aucs_long$cred.level <- as.numeric(sub("^auc_[^_]+_([0-9\\.]+)_.*", "\\1", aucs_long$test_spec))
  aucs_long$sensitivity.option <- as.integer(sub(".*opt([0-9]+)$", "\\1", aucs_long$test_spec))
  
  
  # merge performance measure matrices fprs, tprs, auc to one
  # merge fprs and tprs
  pm_long = dplyr::left_join(fprs_long, tprs_long, by = c(
    "N", "br", "adr.rate", "adr.when", "adr.relsd", "study.period", "tte.dist", 
    "prior.dist", "prior.belief", "dist.prior.to.truth",
    "post.ci.type", "cred.level", "sensitivity.option"
  ))
  # merge fprs, tprs and fnrs
  pm_long = dplyr::left_join(pm_long, fnrs_long, by = c(
    "N", "br", "adr.rate", "adr.when", "adr.relsd", "study.period", "tte.dist", 
    "prior.dist", "prior.belief", "dist.prior.to.truth",
    "post.ci.type", "cred.level", "sensitivity.option"
  ))
  # merge fprs, tprs, fnrs and tnrs
  pm_long = dplyr::left_join(pm_long, tnrs_long, by = c(
    "N", "br", "adr.rate", "adr.when", "adr.relsd", "study.period", "tte.dist", 
    "prior.dist", "prior.belief", "dist.prior.to.truth",
    "post.ci.type", "cred.level", "sensitivity.option"
  ))
  # merge fprs, tprs, fnrs, tnrs and aucs
  # (this is the final result)
  pm_long = dplyr::left_join(pm_long, aucs_long, by = c(
    "N", "br", "adr.rate", "adr.when", "adr.relsd", "study.period", "tte.dist", 
    "prior.dist", "prior.belief", "dist.prior.to.truth",
    "post.ci.type", "cred.level", "sensitivity.option"
  ))
  
  # select relevant columns:
  pm_long <- pm_long[, c(
    "N", "br", "adr.rate", "adr.when", "adr.relsd", "study.period", 
    "tte.dist", "prior.dist", "prior.belief", "dist.prior.to.truth",
    "post.ci.type", "cred.level", "sensitivity.option", 
    "auc", "fpr", "tpr", "fnr", "tnr"
  )]
  
  return(pm_long)
  
}


## END OF DOC