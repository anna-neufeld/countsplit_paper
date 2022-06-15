#!/usr/local/bin/Rscript


library("argparse")
library("tidyr")
source("sim_functions_overdisp.R")

parser <- ArgumentParser()
parser$add_argument("--simname", default = "robust_ses",
                    help = "name of simulation")
parser$add_argument("--nreps", type = "double", default = 200,
                    help = "number of replicates for each set of params")
args <- parser$parse_args()



propLowMedHighs <- 2
regCoeffs <- c(0)

intercepts <- c(0.1,0.5,1,3,5,10,100)
ratios <- c(0.5,0.1,0.25,0.5,1,5,10)
ns <- c(200)
ps <- c(10)
#bs <- intercepts # LOL no idea the scale



## number of monte-carlo iterations per job
nreps_per_combo <- args$nreps
## set up grid of parameters
param_grid <- expand.grid(regCoeff=regCoeffs,
                          intercepts=intercepts,
                          propLowMedHigh = propLowMedHighs,
                          n=ns,
                          p=ps,
                          ratios=ratios)

param_grid$b <- param_grid$intercepts/param_grid$ratios

jobid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

current_dynamic_args <- param_grid[jobid, ]


## -----------------------------------------
## run the simulation nreps_per_job times
## -----------------------------------------
current_seed <- jobid
set.seed(current_seed)

filename <- paste("overdispres_cluster/",args$simname, jobid, ".txt", sep="")

eps=c(0.5)

system.time(replicate(args$nreps, 
                      one_trial_count_split_cluster(n=current_dynamic_args$n,
                                            p=current_dynamic_args$p,
                                            filename,
                                            k=1,
                                            intercepts=log(current_dynamic_args$intercepts), 
                                            eps=eps, 
                                            sig_strength=current_dynamic_args$regCoeff,
                                 propLowMedHigh = probMatrix[current_dynamic_args$propLowMedHigh,],
                                 overdisp=TRUE,b=current_dynamic_args$b)))

filename <- paste("overdispres_pt/",args$simname, jobid, ".txt", sep="")

#eps=c(0.1,0.25,0.5,0.75,0.9)

system.time(replicate(args$nreps, 
                      one_trial_count_split(n=current_dynamic_args$n,
                                                    p=current_dynamic_args$p,
                                                    filename,
                                                    k=1,
                                                    intercepts=log(current_dynamic_args$intercepts), 
                                                    eps=eps, 
                                                    sig_strength=current_dynamic_args$regCoeff,
                                                    propLowMedHigh = probMatrix[current_dynamic_args$propLowMedHigh,],
                                                    overdisp=TRUE,b=current_dynamic_args$b)))
