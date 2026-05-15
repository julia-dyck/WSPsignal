#' Plot ROC curves for top WSP test specifications
#'
#' Plots receiver operating characteristic (ROC) curves for the top-ranked test 
#' specifications.
#'
#' @param rank.tab data frame of ranked test specifications obtained from 
#' \code{\link{eval.rank_auc}} (\code{output$rank.tab})
#' @param n number of top-ranked test specifications to plot (10 by default)
#'
#' @return A ggplot object displaying ROC curves with shaded AUC regions.
#' 
#' @details The function returns the receiver ROC curves for the top n WSP test 
#' configurations based 
#' on the ranking returned by \code{\link{eval.rank_auc}} by plotting the
#' true positive rate (TPR) on the y-axis against the false positive rate (FPR)
#' on the x-axis \insertCite{fawcett2004}{WSPsignal}.
#' Here, we use the ROC curve with one threshold based on equal numbers
#' of ADR-positive and control scenarios.
#' 
#' For definitions of the performance metrics \code{AUC, FPR, TPR, FNR} 
#' and \code{TNR} returned in printed output, 
#' see the details section of \code{\link{eval.calc_perf}}.
#'
#' The ggplot output can be adjusted to individual needs by adding \code{ggplot2}
#' layers to the output.
#' 
#'
#' @examples
#' \dontrun{
#' # Given object rank as output of eval.rank_auc:
#' eval.roc_curve(rank$rank.tab, n = 1) # top 1 test
#' eval.roc_curve(rank$rank.tab, n = 10) # top 10 tests
#' 
#' # exemplary further processing of the plot
#' library(patchwork) # again, formatted next to each other
#' # Left plot: n=1, no legend, no title
#' p1 <- eval.roc_curve(rank$rank.tab, n = 1) +
#'   ggplot2::guides(color = "none", fill = "none") +
#'   ggplot2::labs(title = NULL)
#' Right plot: n=10, legend, no title
#' p2 <- eval.roc_curve(rank$rank.tab, n = 10) +
#'   ggplot2::labs(title = NULL)
#' # Combine with patchwork and add global title
#' (p1 | p2) + patchwork::plot_annotation(title = "ROC curves")
#' }
#'
#' @references 
#' \insertAllCited{}
#'
#' @export

eval.roc_curve = function(rank.tab, n = 10) {

  # argument checks ------------------------------------------------------------
  
  # arg check rank.tab
  required_cols = c("test.type", "tte.dist", "prior.dist", "prior.sd", "post.ci.type", 
                    "cred.level","sensitivity.option", "AUC", "FPR", "TPR", "FNR", "TNR")
  
  valid_ranktab = is.data.frame(rank.tab) && all(required_cols %in% names(rank.tab))
  
  if (!valid_ranktab){
    stop("Argument rank.tab must be a data frame returned by eval.rank_auc (output$rank.tab).\n")
  }
  
  # arg check n
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 1) {
    stop("Argument n must be an integer value >= 1.\n")
  }
  
  # warn + clamp n if too large
  if (n > nrow(rank.tab)) {
    warning("Argument n exceeds number of rows in rank.tab; using n = nrow(rank.tab).\n")
  }
  n = min(n, nrow(rank.tab))
  
  
  # arrange command to make sure, even merged tables are ranked top down!
  rank.tab = dplyr::arrange(rank.tab, dplyr::desc(AUC))
  # print top n tests in console as well
  print(rank.tab[1:n,])
  
  # Add combined label
  rank.tab = rank.tab %>%
    dplyr::slice(1:n) %>%
    dplyr::mutate(
      spec = paste(tte.dist, prior.dist, prior.sd, post.ci.type, cred.level, sensitivity.option, sep = " | "),
      spec = factor(spec, levels = unique(spec))  # preserve legend order
    )
  
  # Build long-format ROC + polygon bottom closure
  roc_df = rank.tab %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      FPR_full = list(c(0, FPR, 1, 1, 0)),  # extended polygon
      TPR_full = list(c(0, TPR, 1, 0, 0))
    ) %>%
    tidyr::unnest(cols = c(FPR_full, TPR_full)) %>%
    dplyr::group_by(spec) %>%
    dplyr::mutate(point = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  # Base ROC curve points for line/points (without bottom closure)
  line_df = rank.tab %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      FPR_vec = list(c(0, FPR, 1)),
      TPR_vec = list(c(0, TPR, 1))
    ) %>%
    tidyr::unnest(cols = c(FPR_vec, TPR_vec)) %>%
    dplyr::mutate(spec = factor(paste(tte.dist, prior.dist, prior.sd, post.ci.type, cred.level, sensitivity.option, sep = " | "),
                                levels = unique(rank.tab$spec))) %>%
    dplyr::group_by(spec) %>%
    dplyr::mutate(point = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  ggplot2::ggplot() +
    # ROC area (polygon to bottom)
    ggplot2::geom_polygon(data = roc_df,
                          ggplot2::aes(x = FPR_full, y = TPR_full, group = spec, fill = spec),
                          alpha = 0.1, color = NA) +
    
    # ROC line
    ggplot2::geom_line(data = line_df,
                       ggplot2::aes(x = FPR_vec, y = TPR_vec, color = spec, group = spec),
                       alpha = 0.6, size = 1.2) +
    
    # ROC points
    ggplot2::geom_point(data = line_df,
                        ggplot2::aes(x = FPR_vec, y = TPR_vec, color = spec, group = spec),
                        alpha = 0.6, size = 2) +
    
    # Diagonal line
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", size = 0.8) +
    
    ggplot2::labs(
      x = "False Positive Rate (FPR)",
      y = "True Positive Rate (TPR)",
      color = "Test specification \n(descending ranking order)",
      fill = "Test specification \n(descending ranking order)",
      title = paste0("ROC curves of top ", n, " test specifications")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")
  
}


## END OF DOC