#' Simulated musculoskeletal pain time-to-event dataset
#'
#' A simulated time-to-event (tte) dataset representing musculoskeletal pain
#' following bisphosphonate intake.
#'
#' @format
#' A data frame with 19 777 rows and 2 variables:
#' \describe{
#'   \item{time}{event time (in days) or censoring time (365 days),}
#'   \item{status}{event indicator; 1 = event observed, 0 = censored.}
#' }
#'
#' @details
#' The data was generated using 
#' \code{\link{sim.datagen_tte(genpar = c(19777, 0.01, 0.89, 160/365, 0.1, 365))}} 
#' with parameters derived from the case study presented in 
#' \insertCite{dyck2024bpgwsppreprint;textual}{WSPsignal}.
#' 
#'
#' @seealso \code{\link{sim.datagen_tte}}
#' @references \insertAllCited{}
#'
"muscu1"
