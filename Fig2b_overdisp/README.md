The file run_sims_overdisp.R runs the simulations that led to Figure 1b. As it was run on the cluster, call_sim.sh and submit_sim.sh were used to run the file itself. The file run_sims_overdisp.R calls functions that are define in sim_functions_overdisp.R.

The output of running run_sims_overdisp.R for many different parameter values are stored in the folder overdispres_pt. The file make_fig_overdisp.R reads in these results and creates Figure 1b. 

Note that the results are stored for several values of Lambda. The make_fig_overdisp.R subsets the results to only show cases where Lambda=5, as specified in the paper.