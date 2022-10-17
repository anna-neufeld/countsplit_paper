#!/bin/bash

njobs=60

qsub -q w-bigmem.q -e iotrash/ -o iotrash/ -l h='biostat-b15|biostat-b16|biostat-b17|biostat-b18' -l h_vmem=10G -t 1-$njobs -tc 60 ./call_sim.sh $1 $2
