#' Fit frequentist model to time-to-event data
#' 
#' Fits a frequentist model to time-to-event (tte) data via maximum
#' likelihood (ML) estimation.
#' 
#' @param dat data frame or matrix with time information in first column and 
#' event information (binary status) in second column
#' @param tte.dist character specifying the distribution for the
#' model out of \code{"w", "dw", "pgw"} (see details)
#' 
#' @return
#' A list with components depending on \code{tte.dist}:
#' \itemize{
#'   \item \code{fit}: fitted model object
#'     (\code{summary.survreg} for \code{"w"}, a list of two \code{summary.survreg} objects for \code{"dw"},
#'     and an \code{nlm} output list for \code{"pgw"})
#'   \item \code{estimates}: data frame of parameter estimates in standard parametrization (scale, shape, powershape)
#'   \item \code{tte.dist}: character indicating the fitted tte distribution
#' }
#' 
#' @details The model can be a Weibull \code{("w")},
#' a double Weibull (\code{"dw"}, estimating two Weibull models - one to the data as is and 
#' one to the data censored at mid of observation period), 
#' or a power generalized Weibull \code{("pgw")} model.
#' 
#' The likelihood used in ML estimation is 
#' \deqn{\mathcal{L}(t) = \prod_{i=1}^N S(t_i)^{1-d_i}\cdot f(t_i)^{d_i}} 
#' with \eqn{S(t)}  being the survival function of the chosen distribution, \eqn{f(t)} the
#' density \insertCite{nikulin2016}{WSPsignal}, and \eqn{(t_i, d_i)} 
#' the (right-censored) tte observations.
#' 
#' For the estimation of the Weibull models\code{("w", "dw")}, the \code{\link[survival]{survreg}} 
#' function (with no covariates) is used.
#' 
#' The \code{"pgw"} model is estimated by numerically minimizing the corresponding
#'  negative log-likelihood function with \code{\link[stats]{nlm}}.
#' 
#' 
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' head(tte)
#' fwsp_model(tte, tte.dist = "w") # Weibull model
#' fwsp_model(tte, tte.dist = "dw") # double Weibull model
#' fwsp_model(tte, tte.dist = "pgw") # power generalized Weibull model
#' 
#' 
#' @export
#' 


fwsp_model = function(dat, 
                      tte.dist = c("dw")
                      ) {
  
  ### checks for tte.dist
  if (any(duplicated(tte.dist))) {
    warning("Duplicate entries removed from tte.dist.\n")
    tte.dist <- unique(tte.dist)
  }
  allowed_dists <- c("w","dw","pgw")
  if (any(is.na(match(tte.dist, allowed_dists))))
    stop(paste0("Argument tte.dist must be out of: ",
                paste(allowed_dists, collapse = ", "),
                ".\n"))
  
  # model fitting dep. on tte.dist
  if(tte.dist == "w"){
    # fit model
    res.w = summary(survival::survreg(survival::Surv(time = dat[,1], event = dat[,2])~1, dist = "weibull"))
    # print estimates in dweibull parametrization
    #   survreg's scale  =    1/(rweibull shape)
    #   survreg's intercept = log(rweibull scale)
    outprint = data.frame(parameter = c("scale", "shape"),
                          estimate = c(exp(res.w$coefficients[1]), 1/res.w$scale),
                          row.names = NULL)
    message("Weibull parameter estimates in rweibull parametrization:")
    print(outprint)
    return(list(fit = res.w, estimates = outprint, tte.dist = tte.dist))
  }
  
  if(tte.dist == "dw"){
    # extract data censored a middle of observation period
    dat.c = dat
    half_op = ceiling(max(dat[,1])/2)
    dat.c[dat$time > half_op,1] = half_op
    dat.c[dat$time > half_op,2] = 0
    # fit models
    res.w = summary(survival::survreg(survival::Surv(time = dat[,1], event = dat[,2])~1, dist = "weibull"))
    res.c.w = summary(survival::survreg(survival::Surv(time = dat.c[,1], event = dat.c[,2])~1, dist = "weibull"))
    outprint = data.frame(parameter = c("scale", "shape", "scale.c", "shape.c"),
                          estimate = c(exp(res.w$coefficients[1]), 1/res.w$scale, 
                                       exp(res.c.w$coefficients[1]), 1/res.c.w$scale),
                          row.names = NULL)
    message("Double Weibull parameter estimates in rweibull parametrization:")
    print(outprint)
    return(list(fit = list(uncens = res.w, cens = res.c.w), estimates = outprint, tte.dist = tte.dist))
  }
  
  if(tte.dist == "pgw"){
    # fit model
    res.pgw = try(nlm(mllk_pgw, p = c(0,0,0), dat = dat, hessian = T))
    outprint = data.frame(parameter = c("scale", "shape", "powershape"),
                          estimate = c(exp(res.pgw$estimate[1]), exp(res.pgw$estimate[2]), exp(res.pgw$estimate[3])),
                          row.names = NULL)
    message("Power generalized Weibull parameter estimates in rpgw parametrization:")
    print(outprint)
    return(list(fit = res.pgw, estimates = outprint, tte.dist = tte.dist))
  }
  
  return(mod)
}
