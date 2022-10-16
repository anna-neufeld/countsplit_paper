# Fig6_realdata

The file myMonocle.R contains the pseudotime estimation strategy that is outlined in detail in Appendix F of our paper. The file run_myMonocle.R calls this function repeatedly to conduct all analyses used in Section 6 of the paper, and to create Figure 4.

The raw data is stored as separate train and test sets in Xtest.RData and Xtrain.RData. Count splitting was performed by our collaborators on a raw dataset with over 900,000 cells and over 30,000 genes, prior to the selection of a single lineage for further analysis. The file metadata_tsv contains come additional metaata for this dataset. 

The function Real_Data_Overdisp.R creates Figure 7, which is discussed in Appendix E of the paper. 