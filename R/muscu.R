#' Simulated musculoskeletal pain time-to-event dataset
#'
#' A simulated time-to-event (tte) dataset on time to musculoskeletal pain
#' from first bisphosphonate intake.
#'
#' @format
#' A data frame with 19 777 rows and 2 variables:
#' \describe{
#'   \item{time}{event time (in days) if an event was observed or censoring time (365 days) if no event was observed,}
#'   \item{status}{event status; 1 if an event was observed, 0 if no event was observed.}
#' }
#'
#' @details
#' The data was generated using 
#' \code{sim.datagen_tte(genpar = c(19777, 0.01, 0.89, 160/365, 0.1, 365))}
#' with parameters derived from the case study presented in 
#' \insertCite{dyck2024bpgwsppreprint;textual}{WSPsignal}.
#' 
#'
#' @seealso \code{\link{sim.datagen_tte}}
#' @references \insertAllCited{}
#'
"muscu"
