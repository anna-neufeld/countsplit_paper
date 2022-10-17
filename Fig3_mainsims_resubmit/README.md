# Fig345_mainsims

The contents of this folder carry out the main simulations described in Section 5 of our paper. 

The file ``sim_functions_cluster.R`` contains all of the important functions for running the simulations, and ``run_sims_cluster.R`` calls these functions for a single set of parameters. The files ``call_sim.sh`` and ``submit_sim.sh`` are used for running these simulations in parallel for many combinations of parameters. 

The results from running these functions across all parameter settings used in the paper should be stored in folders called ``res`` (for pseudotime simulations) and ``clusterres`` (for clustering simulations). Each file in these folders stores the results from a different set of parameters. These files were too large to include in the github repository, so should be reproduced using "run_sims_cluster.R"

The function ``make_figs_non_null.R`` loads in the results files and makes Figures 3,4, and 5 for the paper. 

