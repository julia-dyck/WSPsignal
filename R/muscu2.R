#' Simulated musculoskeletal pain time-to-event dataset
#'
#' A simulated time-to-event (tte) dataset on time to musculoskeletal pain
#' from first bisphosphonate intake.
#'
#' @format
#' A data frame with 1 208 rows and 2 variables:
#' \describe{
#'   \item{time}{event time (in days) if an event was observed or censoring time (365 days) if no event was observed,}
#'   \item{status}{event status; 1 if an event was observed, 0 if no event was observed.}
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
