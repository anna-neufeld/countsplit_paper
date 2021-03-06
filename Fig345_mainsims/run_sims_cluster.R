#!/usr/local/bin/Rscript

## Run the simulation!

## -----------------------------------------
## load user-defined functions, packages
## -----------------------------------------

## get command line arguments
library("argparse")
## tidy stuff
library("tidyr")
## Needed functions.
source("sim_functions_cluster.R")

## -----------------------------------------
## load any command line arguments
## -----------------------------------------
parser <- ArgumentParser()
parser$add_argument("--simname", default = "robust_ses",
                    help = "name of simulation")
parser$add_argument("--nreps", type = "double", default = 200,
                    help = "number of replicates for each set of params")
#parser$add_argument("--nreps-per-job", type = "double", default = 200,
#                    help = "number of replicates per job")
args <- parser$parse_args()


## -----------------------------------------
## set up a grid of parameters to cycle over
## -----------------------------------------

propLowMedHighs <- 1:3
regCoeffs <- c(log(1.2), log(1.35), log(1.5), log(1.6), log(1.7), log(1.85), log(2), log(2.5), log(3), log(4), log(5), log(7), log(10), log(15), log(20))
#propImps <- c(0.05,0.1,0.2,0.3)


propImps <- c(0.1)
ns <- c(200)
#ps <- c(20, 100,200,400)
ps <- c(100)

## number of monte-carlo iterations per job
nreps_per_combo <- args$nreps
## set up grid of parameters
param_grid <- expand.grid(regCoeff=regCoeffs,
                          propImp=propImps,
                          propLowMedHigh = propLowMedHighs,
                          n=ns,
                          p=ps)

probMatrix <- rbind(c(0.5,0.5),
                    c(0,1),
                    c(1,0))

## -----------------------------------------
## get current dynamic arguments
## -----------------------------------------
## get job id from scheduler

jobid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

                                   #job_id <- 1

## current args
current_dynamic_args <- param_grid[jobid, ]


## -----------------------------------------
## run the simulation nreps_per_job times
## -----------------------------------------
current_seed <- jobid
#current_seed <- 1
set.seed(current_seed)

filename <- paste("res/",args$simname, jobid, ".txt", sep="")

eps=c(0.1,0.25,0.5,0.75,0.9)

system.time(replicate(args$nreps, 
                      one_trial_count_split(n=current_dynamic_args$n,
                                            p=current_dynamic_args$p,
                                            filename,
                                            k=1,
                                            propImp=current_dynamic_args$propImp, 
                                            eps=eps, 
                                            sig_strength=current_dynamic_args$regCoeff,
                                 propLowMedHigh = probMatrix[current_dynamic_args$propLowMedHigh,])))


filename <- paste("clusterres/",args$simname, jobid, ".txt", sep="")

#eps=c(0.1,0.25,0.5,0.75,0.9)

system.time(replicate(args$nreps,
                      one_trial_count_split_cluster(n=current_dynamic_args$n,
                                            p=current_dynamic_args$p,
                                            filename,
                                            k=1,
                                            propImp=current_dynamic_args$propImp,
                                            eps=eps,
                                            sig_strength=current_dynamic_args$regCoeff,
                                 propLowMedHigh = probMatrix[current_dynamic_args$propLowMedHigh,])))


