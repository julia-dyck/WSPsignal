#' test time-to-event dataset
#'
#' Simulated time-to-event (tte) data generated with \code{sim.datagen_tte(c(100, 0.1, 1, 0.5, 0.05, 365))}.
#' 
#' @format  A data frame with 100 rows and 2 variables:
#' \describe{
#'   \item{time}{event time if an event was observed, or censoring time if no event was observed,}
#'   \item{status}{event status; 1 if an event was observed, 0 if no event was observed.}
#'}
#' 
#' @seealso \code{\link{sim.datagen_tte}}
#'
"tte"