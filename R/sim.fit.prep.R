#' Prepare Bayesian pgW fitting for simulation study
#' 
#' 
#' Formatting of simulated data generated with \code{\link{sim.datagen_tte}} as 
#' building block for tuning the Bayesian WSP test via simulation study.
#'
#' @param ttedat time-to-event data set
#' @param pc vector representing one parameter combination used in simulation study
#' @param pc_list list containing all parameter combinations and additional
#'        information necessary for simulation study
#'
#' @return A list containing all information to be inserted as \code{datstan} 
#' argument into the \code{\link{bwsp_model}} function.
#'
#' @noRd


## CONCRETE SPECIFICATION DEPENDS ON tte.dist and prior.belief 

sim.fit.prep = function(ttedat, pc, pc_list){
  # ttedat = data set in time-event-format generated from pc
  # pc contains N, br, adr.rate, adr.when, adr.relsd, censor, tte.dist, prior.dist, prior.belief, prior.sd
  # pc_list contains all additional parameters necessary to specify the simulation study 
  # (eg. prior means depending on prior belief, 
  # additional parameters constant for all simulations)
  
  
  if(pc$tte.dist == "w"){
    # extract prior pars for prior belief
    belief.ind = which(pc_list$fit$w$prior.belief == pc$prior.belief)[1]
    pars = pc_list$fit$w[belief.ind,
                         c("scale.mean_w", "shape.mean_w")]
    
    # format data and prior pars accordingly
    datstan = tte2priordat(dat = ttedat, 
                           tte.dist = "w",
                           prior.dist = pc$prior.dist,
                           scale.mean = pars$scale.mean_w, 
                           scale.sd = pc$prior.sd,
                           shape.mean = pars$shape.mean_w, 
                           shape.sd = pc$prior.sd
    )
  }
  
  if(pc$tte.dist == "dw"){
    # extract prior pars for prior belief
    belief.ind = which(pc_list$fit$dw$prior.belief == pc$prior.belief)[1]
    pars = pc_list$fit$dw[belief.ind,
                          c("scale.mean_dw", "shape.mean_dw", 
                            "scale_c.mean_dw", "shape_c.mean_dw")]
    # format data and prior pars accordingly
    datstan = tte2priordat(dat = ttedat, 
                           tte.dist = "dw",
                           prior.dist = pc$prior.dist,
                           scale.mean = pars$scale.mean_dw, 
                           scale.sd = pc$prior.sd,
                           shape.mean = pars$shape.mean_dw, 
                           shape.sd = pc$prior.sd,
                           scale_c.mean = pars$scale_c.mean_dw,
                           scale_c.sd = pc$prior.sd,
                           shape_c.mean = pars$shape_c.mean_dw,
                           shape_c.sd = pc$prior.sd
    )
  }
  
  if(pc$tte.dist == "pgw"){
    # extract prior pars for prior belief
    belief.ind = which(pc_list$fit$pgw$prior.belief == pc$prior.belief)[1]
    pars = pc_list$fit$pgw[belief.ind,
                           c("scale.mean_pgw", "shape.mean_pgw",
                             "powershape.mean_pgw")]
    # format data and prior pars accordingly
    datstan = tte2priordat(dat = ttedat,
                           tte.dist = "pgw",
                           prior.dist = pc$prior.dist,
                           scale.mean = pars$scale.mean_pgw, 
                           scale.sd = pc$prior.sd,
                           shape.mean = pars$shape.mean_pgw, 
                           shape.sd = pc$prior.sd,
                           powershape.mean = pars$powershape.mean_pgw,
                           powershape.sd = pc$prior.sd
    )
  }
  
  return(datstan)
}

