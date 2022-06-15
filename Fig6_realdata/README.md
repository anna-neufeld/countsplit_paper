# Fig6_realdata

The file myMonocle.R contains the pseudotime estimation strategy that is outlined in detail in Appendix F of our paper. The file run_myMonocle.R calls this function repeatedly to conduct all analyses used in Section 6 of the paper, and to create Figure 6.

The file sce_subset.rds contains a SingleCellExperiment object which contains the raw data and metadata provided by Elborany et al. for our real data analysis. The file metadata_subset.tsv contains come additional metaata for this dataset.

The function Real_Data_Overdisp.R creates Figure 8, which is discussed in Appendix E of the paper. 