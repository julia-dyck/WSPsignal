#' The 'WSPsignal' package.
#'
#' @description An R package to perform Bayesian and frequentist Weibull Shape Parameter (WSP)
#' tests for signal detection.
#'
#' @docType package
#' @name WSPsignal-package
#' @aliases WSPsignal
#' @useDynLib WSPsignal, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertRef{sauzet2022}{WSPsignal}
#' 
#' \insertRef{sauzet2024}{WSPsignal}
#' 
#' \insertRef{dyck2024bpgwsppreprint}{WSPsignal}
#'
#' Stan Development Team (NA). RStan: the R interface to Stan. 
#' R package version 2.32.6. https://mc-stan.org
NULL
