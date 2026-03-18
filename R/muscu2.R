#' Simulated musculoskeletal pain time-to-event dataset
#'
#' A simulated time-to-event (tte) dataset ispired by tte data about musculoskeletal pain
#' following bisphosphonate intake.
#'
#' @format
#' A data frame with 1 208 rows and 2 variables:
#' \describe{
#'   \item{time}{event time (in days) or censoring time (365 days),}
#'   \item{status}{event indicator; 1 = event observed, 0 = censored.}
#' }
#'
#' @details
#' The data was generated using 
#' \code{\link{sim.datagen_tte(genpar = c(1208, 0.01, 0.5, 100/365, 0.1, 365))}}.
#' Data generation was guided by the data presented in the case study in
#' \insertCite{dyck2024bpgwsppreprint;textual}{WSPsignal}.
#' 
#'
#' @seealso \code{\link{sim.datagen_tte}}
#' @references \insertAllCited{}
#'
"muscu2"
