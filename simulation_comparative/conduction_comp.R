#### simulationstudy
#### as R script

#### comment on version:
# This R script is executable only with the tagged R package version
# https://github.com/julia-dyck/WSPsignal/releases/tag/corresponding_to_BPgWSP_preprint
# NOT with the actual version

#### prep ----------------------------------------------------------------------
# packages

library(BWSPsignal) # signal detection test & simulation fcts.
library(dplyr)      # for handling results and 

####  prior elicitation --------------------------------------------------------
# try a few prior parameter combinations and see whether the resulting hazard 
# roughly matches the prior belief about the hazard form 

# Expected event time can also be taken into account for some guidance, but 
# should not be prioritized.
# The reason is that we do not expect the model to accurately fit the hazard of the 
# data, but only catch the rough form by distinguishing the cases
# constant vs decreasing vs unimodal vs increasing hazard.

# set prior means for Weibull parameters (ie powershape = 1 always):
365/4 * c(1,2,3) # expected event times under prior belief beginning, middle, end of study period
plot_pgw(scale = 1, shape = 1, powershape = 1)     # under prior belief "none"
plot_pgw(scale = 1, shape = 0.207, powershape = 1) # under prior belief "beginning"
plot_pgw(scale = 180, shape = 1, powershape = 1)   # under prior belief "middle" (cannot create a unimodal form)
plot_pgw(scale = 300, shape = 4, powershape = 1)   # under prior belief "end"

# set prior means for double Weibull setting:
# uncensored Weibull parameters are set equal to Weibull parameters
plot_pgw(scale = 1, shape = 1, powershape = 1)     # under prior belief "none"
plot_pgw(scale = 1, shape = 0.207, powershape = 1) # under prior belief "beginning"
plot_pgw(scale = 180, shape = 1, powershape = 1)   # under prior belief "middle" (cannot create a unimodal form)
plot_pgw(scale = 300, shape = 4, powershape = 1)   # under prior belief "end"
# for the censored Weibull parameters, imagine a priori expected hazard forms up to half of the study period
plot_pgw(scale = 1, shape = 1, powershape = 1)     # under prior belief "none"
plot_pgw(scale = 1, shape = 0.207, powershape = 1) # under prior belief "beginning"
plot_pgw(scale = 100, shape = 4, powershape = 1)   # under prior belief "middle" 
plot_pgw(scale = 1, shape = 1, powershape = 1)     # under prior belief "end"

# set prior means for Power generalized Weibull parameters:
plot_pgw(scale = 1, shape = 1, powershape = 1)     # under prior belief "none"
plot_pgw(scale = 20, shape = 5.5, powershape = 14) # under prior belief "beginning"
plot_pgw(scale = 180, shape = 1, powershape = 1)   # under prior belief "middle"
plot_pgw(scale = 300, shape = 4, powershape = 1)   # under prior belief "end"


#### specify parameter combinations for simulation study -----------------------

fp_list = sim.priors_template(tte.dist = c("pgw"),
                              prior.sds = 10) # setup prior template
# fill in prior template with values chosen in prior elicitation

# fp_list$w[,2] = c(1, 1, 180, 300) # scale prior means
# fp_list$w[,4] = c(1, 0.207, 1, 4) # shape prior means
# 
# fp_list$dw[,2] = c(1, 1, 180, 300) # scale prior means
# fp_list$dw[,4] = c(1, 0.207, 1, 4) # shape prior means
# fp_list$dw[,6] = c(1, 1, 100, 1)   # scale_c prior means
# fp_list$dw[,8] = c(1, 0.207, 4, 1) # shape_c prior means

fp_list$pgw[,2] = c(1, 1, 20, 300)   # scale prior means
fp_list$pgw[,4] = c(1, 0.207, 5.5, 4)# shape prior means
fp_list$pgw[,6] = c(1, 1, 14, 1)     # powershape prior means

fp_list # filled fitpars.list

pc_list = sim.setup_sim_pars(N = c(500, 3000, 5000),
                             br = 0.1,
                             adr.rate = c(0, 0.5, 1),
                             adr.relsd = 0.05,
                             study.period = 365,
                             
                             tte.dist = c("pgw"),
                             prior.dist = c("ll", "gg"),
                             fitpars.list = fp_list,
                             
                             post.ci.type = c("ETI", "HDI"),
                             cred.level = c(seq(0.5,0.9, by = 0.05), seq(0.91,0.99, by = 0.01), seq(0.991, 0.999, by = 0.001)), 
                             sensitivity.option = 1:3,
                             
                             reps = 100, # additional parameters
                             batch.size = 10,
                             
                             resultpath = "C:/Users/jdyck/sciebo/bADR/simulation_comparative",
                             stanmod.chains = 4,
                             stanmod.iter = 11000,
                             stanmod.warmup = 1000
)


save(pc_list, file = paste0(pc_list$add$resultpath, "/pc_list.RData"))


#### test for issues in prior choice -------------------------------------------

# test rope calculation given all prior distributions under none (exemplary for
# credibility level 80%):
qlnorm(p = c(0.1,0.9), meanlog = logprior_repar(1, 10)[1], sdlog = logprior_repar(1, 10)[2])
qgamma(p = c(0.1,0.9), shape = gamprior_repar(1, 10)[1], rate = gamprior_repar(1, 10)[2])


#### run simulation study ------------------------------------------------------

sim.run(pc_list = pc_list)


#### merge results -------------------------------------------------------------

# merge bayesian result batches and save
sim.merge_results(pc_list, save = T)
# merge frequentist result batches and save
sim.merge_results(pc_list, save = T, bayes = F)
# load merged results into environment
load(paste0(pc_list$add$resultpath, "/res_b.RData"))
load(paste0(pc_list$add$resultpath, "/res_f.RData"))


#### convergence issues, execution times and effective sample sizes ------------
# for different prior distributional choices

## no of simulation runs with convergence issues
eval.non_conv_cases(pc_list)

## execution times 
eval.execution_times(pc_list)

## effective sample sizes
eval.eff_sample_sizes(pc_list, threshold = 10000)


#### Performance measure calculation -------------------------------------------

# calculate performance measures per scenario with
perf = eval.calc_perf(pc_list = pc_list)
save(perf, file = paste0(pc_list$add$resultpath, "/perf.RData"))

load(paste0(pc_list$add$resultpath, "/perf.RData"))
# perf contains performance metrics for B & FWSP test types

#### Ranking results -----------------------------------------------------------

# AUC ranking of all PgWSP tests
rank = eval.rank_auc(perf, tte.dist.subset = "pgw")


# Top tests

## roc curves
eval.roc_curve(rank$rank.tab, n = 1)
eval.roc_curve(rank$rank.tab, n = 10)

library(patchwork) # again, formatted next to each other
# Left plot: n=1, no legend, no title
p1 <- eval.roc_curve(rank$rank.tab, n = 1) +
  ggplot2::guides(color = "none", fill = "none") +
  ggplot2::labs(title = NULL)
# Right plot: n=10, legend, no title
p2 <- eval.roc_curve(rank$rank.tab, n = 10) +
  ggplot2::labs(title = NULL)
# Combine with patchwork and add global title
(p1 | p2) + patchwork::plot_annotation(title = "ROC curves")


## Effects of dgp parameters on AUC of top test
rank


## END OF DOC
