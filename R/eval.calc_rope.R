#' Derive ROPEs from setup parameter combinations in simulation study
#' 
#' 
#' @noRd


eval.calc_rope = function(rope.infos.row, cred.levels){
  # rope.infos.row is a row vector containing tte.dist, prior.dist, and prior means & sds for all shape pars
  # cred.levels is a vector of credibility levels (e.g. c(0.8, 0.9, 0.95))
  
  tte.dist = rope.infos.row$tte.dist     # for if question
  prior.dist = rope.infos.row$prior.dist # for if question
   
  rope.percentile.l = (1-cred.levels)/2
  rope.percentile.u = cred.levels + (1-cred.levels)/2
  rope.perc = cbind(rope.percentile.l, rope.percentile.u) # calc rope for each row
  
  
  if(tte.dist == "w"){
    # w distribution
    m = rope.infos.row$shape.mean_w
    sd = rope.infos.row$shape.sd_w
    
    if(prior.dist == "ll" || prior.dist == "fl"){
      # lognormal pars from mean and sd
      pars = logprior_repar(mean = m, sd = sd)
      # prepare storage
      ropes.list = vector("list", length(cred.levels))
      
      for (i in 1:length(cred.levels)) {
        cred = cred.levels[i]
        lower = stats::qlnorm((1 - cred) / 2, meanlog = pars[1], sdlog = pars[2])
        upper = stats::qlnorm(cred + (1 - cred) / 2, meanlog = pars[1], sdlog = pars[2])
        ropes.list[[i]] = c(lower, upper)
      }
    } 
    else if(prior.dist == "fg" || prior.dist == "gg"){
      # gamma pars from mean and sd
      pars = gamprior_repar(mean = m, sd = sd)
      # prepare storage
      ropes.list = vector("list", length(cred.levels))
      
      for (i in 1:length(cred.levels)) {
        cred = cred.levels[i]
        lower = stats::qgamma((1 - cred) / 2, shape = pars[1], rate = pars[2])
        upper = stats::qgamma(cred + (1 - cred) / 2, shape = pars[1], rate = pars[2])
        ropes.list[[i]] = c(lower, upper)
      }
    }
    # flatten and name
    ropes.vect = unlist(ropes.list)
    rope.names = c(rbind(
      paste0("nu.rope", cred.levels, "l"),
      paste0("nu.rope", cred.levels, "u")
    ))
    names(ropes.vect) = rope.names
    ropes.vect.nu = ropes.vect
    
    # add empty ga cols
    ropes.vect.ga = rep(NA, length(cred.levels) * 2)
    rope.names = c(rbind(
      paste0("ga.rope", cred.levels, "l"),
      paste0("ga.rope", cred.levels, "u")
    ))
    names(ropes.vect.ga) = rope.names
    
    ropes.vect = c(ropes.vect.nu, ropes.vect.ga)
    return(ropes.vect)
    
    
    return(ropes.vect)
    
  } else if (tte.dist == "dw"){
    # dw distribution
    
    # for shape aka first shape parameter aka nu
    m = rope.infos.row$shape.mean_dw
    sd = rope.infos.row$shape.sd_dw
    
    if(prior.dist == "ll" || prior.dist == "fl"){
      # lognormal pars from mean and sd
      pars = logprior_repar(mean = m, sd = sd)
      # prepare storage
      ropes.list = vector("list", length(cred.levels))
      
      for (i in 1:length(cred.levels)) {
        cred = cred.levels[i]
        lower = stats::qlnorm((1 - cred) / 2, meanlog = pars[1], sdlog = pars[2])
        upper = stats::qlnorm(cred + (1 - cred) / 2, meanlog = pars[1], sdlog = pars[2])
        ropes.list[[i]] = c(lower, upper)
      }
    } 
    else if(prior.dist == "fg" || prior.dist == "gg"){
      # gamma pars from mean and sd
      pars = gamprior_repar(mean = m, sd = sd)
      # prepare storage
      ropes.list = vector("list", length(cred.levels))
      
      for (i in 1:length(cred.levels)) {
        cred = cred.levels[i]
        lower = stats::qgamma((1 - cred) / 2, shape = pars[1], rate = pars[2])
        upper = stats::qgamma(cred + (1 - cred) / 2, shape = pars[1], rate = pars[2])
        ropes.list[[i]] = c(lower, upper)
      }
    }
    # flatten and name
    ropes.vect = unlist(ropes.list)
    rope.names = c(rbind(
      paste0("nu.rope", cred.levels, "l"),
      paste0("nu.rope", cred.levels, "u")
    ))
    names(ropes.vect) = rope.names
    ropes.vect.nu = ropes.vect
    
    
    # for shape_c aka second shape parameter aka ga
    m = rope.infos.row$shape_c.mean_dw
    sd = rope.infos.row$shape_c.sd_dw
    
    if(prior.dist == "ll" || prior.dist == "fl"){
      # lognormal pars from mean and sd
      pars = logprior_repar(mean = m, sd = sd)
      # prepare storage
      ropes.list = vector("list", length(cred.levels))
      
      for (i in 1:length(cred.levels)) {
        cred = cred.levels[i]
        lower = stats::qlnorm((1 - cred) / 2, meanlog = pars[1], sdlog = pars[2])
        upper = stats::qlnorm(cred + (1 - cred) / 2, meanlog = pars[1], sdlog = pars[2])
        ropes.list[[i]] = c(lower, upper)
      }
    } 
    else if(prior.dist == "fg" || prior.dist == "gg"){
      # gamma pars from mean and sd
      pars = gamprior_repar(mean = m, sd = sd)
      # prepare storage
      ropes.list = vector("list", length(cred.levels))
      
      for (i in 1:length(cred.levels)) {
        cred = cred.levels[i]
        lower = stats::qgamma((1 - cred) / 2, shape = pars[1], rate = pars[2])
        upper = stats::qgamma(cred + (1 - cred) / 2, shape = pars[1], rate = pars[2])
        ropes.list[[i]] = c(lower, upper)
      }
    }
    # flatten and name
    ropes.vect = unlist(ropes.list)
    rope.names = c(rbind(
      paste0("ga.rope", cred.levels, "l"),
      paste0("ga.rope", cred.levels, "u")
    ))
    names(ropes.vect) = rope.names
    ropes.vect.ga = ropes.vect
    
    # combine
    ropes.vect = c(ropes.vect.nu, ropes.vect.ga)
    return(ropes.vect)
    
  } else if (tte.dist == "pgw"){
    # pgw distribution
    
    # for shape aka first shape parameter aka nu
    m = rope.infos.row$shape.mean_pgw
    sd = rope.infos.row$shape.sd_pgw
    
    if(prior.dist == "ll" || prior.dist == "fl"){
      # lognormal pars from mean and sd
      pars = logprior_repar(mean = m, sd = sd)
      # prepare storage
      ropes.list = vector("list", length(cred.levels))
      
      for (i in 1:length(cred.levels)) {
        cred = cred.levels[i]
        lower = stats::qlnorm((1 - cred) / 2, meanlog = pars[1], sdlog = pars[2])
        upper = stats::qlnorm(cred + (1 - cred) / 2, meanlog = pars[1], sdlog = pars[2])
        ropes.list[[i]] = c(lower, upper)
      }
    } 
    else if(prior.dist == "fg" || prior.dist == "gg"){
      # gamma pars from mean and sd
      pars = gamprior_repar(mean = m, sd = sd)
      # prepare storage
      ropes.list = vector("list", length(cred.levels))
      
      for (i in 1:length(cred.levels)) {
        cred = cred.levels[i]
        lower = stats::qgamma((1 - cred) / 2, shape = pars[1], rate = pars[2])
        upper = stats::qgamma(cred + (1 - cred) / 2, shape = pars[1], rate = pars[2])
        ropes.list[[i]] = c(lower, upper)
      }
    }
    # flatten and name
    ropes.vect = unlist(ropes.list)
    rope.names = c(rbind(
      paste0("nu.rope", cred.levels, "l"),
      paste0("nu.rope", cred.levels, "u")
    ))
    names(ropes.vect) = rope.names
    ropes.vect.nu = ropes.vect
    
    
    # for powershape aka second shape parameter aka ga
    m = rope.infos.row$powershape.mean_pgw
    sd = rope.infos.row$powershape.sd_pgw
    
    if(prior.dist == "ll" || prior.dist == "fl"){
      # lognormal pars from mean and sd
      pars = logprior_repar(mean = m, sd = sd)
      # prepare storage
      ropes.list = vector("list", length(cred.levels))
      
      for (i in 1:length(cred.levels)) {
        cred = cred.levels[i]
        lower = stats::qlnorm((1 - cred) / 2, meanlog = pars[1], sdlog = pars[2])
        upper = stats::qlnorm(cred + (1 - cred) / 2, meanlog = pars[1], sdlog = pars[2])
        ropes.list[[i]] = c(lower, upper)
      }
    } 
    else if(prior.dist == "fg" || prior.dist == "gg"){
      # gamma pars from mean and sd
      pars = gamprior_repar(mean = m, sd = sd)
      # prepare storage
      ropes.list = vector("list", length(cred.levels))
      
      for (i in 1:length(cred.levels)) {
        cred = cred.levels[i]
        lower = stats::qgamma((1 - cred) / 2, shape = pars[1], rate = pars[2])
        upper = stats::qgamma(cred + (1 - cred) / 2, shape = pars[1], rate = pars[2])
        ropes.list[[i]] = c(lower, upper)
      }
    }
    # flatten and name
    ropes.vect = unlist(ropes.list)
    rope.names = c(rbind(
      paste0("ga.rope", cred.levels, "l"),
      paste0("ga.rope", cred.levels, "u")
    ))
    names(ropes.vect) = rope.names
    ropes.vect.ga = ropes.vect
    
    # combine
    ropes.vect = c(ropes.vect.nu, ropes.vect.ga)
    return(ropes.vect)
    
    
    
  }
  
}

# # test
# 
# rbind(
# eval.calc_rope(cred.levels = c(0.8, 0.9, 0.95), 
#                rope.infos.row = rope.infos[1,]) ,
# 
# eval.calc_rope(cred.levels = c(0.8, 0.9, 0.95), 
#               rope.infos.row = rope.infos[2,]) ,
# 
# eval.calc_rope(cred.levels = c(0.8, 0.9, 0.95), 
#               rope.infos.row = rope.infos[3,]) ,
# 
# eval.calc_rope(cred.levels = c(0.8, 0.9, 0.95),
#               rope.infos.row = rope.infos[4,]) 
# )

