#' extract posterior summary statistics and credibility intervals
#'
#' Given a stanfit object returned by \code{\link{bwsp_model}}, the function extracts statistics about the posterior
#' distribution of the shape parameters.
#' These shall give an overview and provide the basis for the signal detection
#' testing in one vector.
#'
#' @param stanfit.object the estimated stan model output
#' 
#' @return Information about the posterior distributions stored in a vector
#'         consisting of the following entry parts first for parameter nu followed by the same statistics regarding
#'         parameter gamma:
#'
#' \item{nu/ga.post.stats}{Summary statistics from the stanfit object about the posterior distribution, namely \code{mean, se_mean, sd},
#'                  and bayesian convergence diagnostic measures, namely \code{n_eff, Rhat}}
#'
#'  @details The storing in one long vector for each parameter is provided to run a
#'  simulation study, where all relevant statistics per run are stored column-wise
#'  and the repeated runs for one data scenario are stored row-wise.
#'
#' @noRd




# function --------------------------------------------------------------------
sim.stanfit.to.poststats = function(pc, stanfit.object){
  
  if(pc$tte.dist == "w"){
    obj = stanfit.object
    
    # extract running time
    run.min = sum(rstan::get_elapsed_time(obj))/60 # in minutes
    names(run.min) = "run.min"
    
    # summary statistics for shape parameters
    post_summary = rstan::summary(obj, pars = c("nu"), probs = c())$summary
    poststats = list(nu = post_summary["nu",])
    
    nu.post.stats = poststats$nu
    names(nu.post.stats) = paste0("nu.po.", names(nu.post.stats))
    
    # vector to be returned
    
    ret.vect = data.frame(t(c(run.min,
                              nu.post.stats))) #, nu.eti, nu.hdi, nu.per))) 
    
    ret.vect = cbind(pc, ret.vect)
    
  }
  
  if(pc$tte.dist == "dw"){
    obj = stanfit.object$uncens
    obj_c = stanfit.object$cens
    
    # extract running time
    
    run.min = sum(rstan::get_elapsed_time(obj) + rstan::get_elapsed_time(obj_c))/60 # in minutes
    names(run.min) = "run.min"
    
    # uncensored: 
    # summary statistics for shape parameter
    post_summary = rstan::summary(obj, pars = c("nu"), probs = c())$summary
    poststats = list(nu = post_summary["nu",])
    
    nu.post.stats = poststats$nu
    names(nu.post.stats) = paste0("nu.po.", names(nu.post.stats))
    
    
    # censored: 
    # summary statistics for shape parameter
    post_summary = rstan::summary(obj_c, pars = c("nu"), probs = c())$summary
    poststats = list(ga = post_summary["nu",])
    
    ga.post.stats = poststats$ga
    names(ga.post.stats) = paste0("ga.po.", names(ga.post.stats))
    
    # vector to be returned
    
    ret.vect = data.frame(t(c(run.min,
                              nu.post.stats, # nu.eti, nu.hdi, nu.per,
                              ga.post.stats  # , ga.eti, ga.hdi, ga.per)))
                              )))
    ret.vect = cbind(pc, ret.vect)
    
  }
  
  if(pc$tte.dist == "pgw"){
    obj = stanfit.object
    
    # extract running time
    
    run.min = sum(rstan::get_elapsed_time(obj))/60 # in minutes
    names(run.min) = "run.min"
    
    # summary statistics for shape parameters
    
    post_summary = rstan::summary(obj, pars = c("nu", "gamma"), probs = c())$summary
    poststats = list(nu = post_summary["nu",], ga = post_summary["gamma",])
    
    nu.post.stats = poststats$nu
    names(nu.post.stats) = paste0("nu.po.", names(nu.post.stats))
    
    ga.post.stats = poststats$ga
    names(ga.post.stats) = paste0("ga.po.", names(ga.post.stats))
    
    # vector to be returned
    
    ret.vect = data.frame(t(c(run.min,
                              nu.post.stats, # nu.eti, nu.hdi, nu.per,
                              ga.post.stats  #, ga.eti, ga.hdi, ga.per)))
                              )))
    ret.vect = cbind(pc, ret.vect)
  }

  rownames(ret.vect) = NULL
  return(ret.vect)
}



# ## test
# 
# # mod = inner mod object in sim.fit.to.1.sample fct
# post = sim.stanfit.to.poststats(pc = pc, stanfit.object = mod)
# View(post)

## END of Doc
