#' Fit frequentist model to time-to-event data
#' 
#' Fits a frequentist model to time-to-event (tte) data via maximum
#' likelihood (ML) estimation.
#' 
#' @param dat data frame or matrix with time information in first column and 
#' event information (binary status) in second column
#' @param tte.dist character specifying the distribution for the
#' model; options are \code{"w", "dw", "pgw"} (see details)
#' 
#' @return Output of the fitted model. For \code{"w"}, a \code{summary.survreg} object;
#' for \code{"dw"}, a list of two \code{summary.survreg} objects
#' with \code{$uncens} containing the ML estimate
#' obtained from the model to the tte data as is, and \code{$cens} containing the 
#' ML estimate obtained from the model fitted to the tte data censored at mid of observation period; 
#' for \code{pgw}, a list containing the maximum likelihood estimate 
#' (log of all parameters) obtained from \code{\link[stats]{nlm}}.
#' 
#' @details The model can be a Weibull \code{("w")},
#' a double Weibull (\code{"dw"}, estimating two Weibull models - one to the data as is and 
#' one to the data censored at mid of observation period), 
#' or a power generalized Weibull \code{("pgw")} model.
#' 
#' The likelihood used in ML estimation is 
#' \deqn{\mathcal{L}(t) = \prod_{i=1}^N S(t_i)^{1-d_i}\cdot f(t_i)^{d_i}} 
#' with \eqn{S(t)} the survival function of the chosen distribution and \eqn{f(t)} the
#' density function \insertCite{nikulin2016}{WSPsignal}. The pair \eqn{(t_i, d_i)} 
#' are the tte observations.
#' 
#' Since the \code{survreg} function from the \code{survival} package uses a different parametrization,
#' the parameters transformed to the parametrization used in \code{\link[stats]{rweibull}} 
#' and \code{\link{rpgw}} are printed after function call.
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
  tte.dist <- match.arg(tte.dist)
  
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
