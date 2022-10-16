#!/bin/bash
#$ -cwd

Rscript run_sims_cluster.R --simname $1 --nreps $2
