# Fig1_globalnull

This folder contains code used to create Figure 1 in our paper.

The file ``global_null_run_others.R`` runs the code for the simulations for all methods besides PseudotimeDE. The code ``global_null_fun_pseudotimeDE.R`` runs the code only for pseudotimeDE-- fewer iterations are run for this method as it is the slowest by far. The results from runs of these methods are stored in ``global_null_res_may_23.RData`` (all methods besides pseudotimeDE) and ``pvals_DE_pt.RData`` (pseudotimeDE only). 

The file ``myjackstraw.R`` contains our implementation of Jackstraw, and is called by ``global_null_run_others.R``.

The file ``global_null_plot.R`` reads in the saved results and creates Figure from our paper. 

