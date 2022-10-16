#!/bin/bash
#$ -cwd

Rscript run_sims_overdisp.R --simname $1 --nreps $2
