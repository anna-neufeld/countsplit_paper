### Start with the Seurat simulations
library(Seurat)

n=500
p=200
n_trials <- 500

pvals <- matrix(NA, nrow=n_trials, ncol=p)
pvals_CS <-  matrix(NA, nrow= n_trials, ncol=p)

for (i in 1:n_trials) {
  print(i)
  set.seed(i)
  X <- matrix(rpois(n*p, lambda=5), nrow=p)
  colnames(X) <- paste("cell", 1:n, sep="")
  rownames(X) <- paste("gene", 1:p, sep="")
  
  ep=0.5
  Xtrain <- apply(X, 2, function(u) rbinom(length(u), u, ep))
  Xtest <- X - Xtrain
  colnames(Xtrain) <- paste("cell", 1:n, sep="")
  rownames(Xtrain) <- paste("gene", 1:p, sep="")
  colnames(Xtest) <- paste("cell", 1:n, sep="")
  rownames(Xtest) <- paste("gene", 1:p, sep="")
  
  
  ### CREATE SEURAT OBJECTS FOR BOTH. 
  cds <- CreateSeuratObject(X)
  cds_train <- CreateSeuratObject(Xtrain)
  
  cds <- ScaleData(cds)
  cds <- RunPCA(cds, features=rownames(X))
  cds  <- FindNeighbors(cds , dims = 1:10)
  cds <- FindClusters(cds, resolution = 1)
  
  cds_train <- ScaleData(cds_train)
  cds_train <- RunPCA(cds_train, features=rownames(X))
  cds_train  <- FindNeighbors(cds_train , dims = 1:10)
  cds_train <- FindClusters(cds_train, resolution = 1)
  
  cds_train[['test']] <- CreateAssayObject(counts=Xtest)
  cds_train <- NormalizeData( cds_train , assay="test")
  cds_train <- ScaleData( cds_train , assay="test")
  
  if (length(unique(Idents(cds)))>1) {
    pvals[i,]  <- FindMarkers(cds, ident.1 = 1, ident.2 = 2, min.pct=0, min.cells.feature=0,
                              min.cells.group=0,
                              logfc.threshold=0)$p_val
  }
  
  if (length(unique(Idents(cds_train)))>1) {
    pvals_CS[i,]  <- FindMarkers(cds_train, ident.1 = 1, ident.2 = 2, min.pct=0, assay="test", min.cells.feature=0,
                                 min.cells.group=0,
                                 logfc.threshold=0)$p_val
  }
}

save(pvals, file="pvals_seurat.RData")
save(pvals_CS, file="pvals_seurat_CS.RData")

### Now do the monocle3 simulations
library(monocle3)

n=500
p=200
n_trials <- 500

pvals_graph_test <- matrix(NA, nrow=n_trials, ncol=p)
pvals_graph_test_CS <-  matrix(NA, nrow= n_trials, ncol=p)

for (i in 1:n_trials) {
  print(i)
  set.seed(i)
  X <- matrix(rpois(n*p, lambda=5), nrow=p)
  colnames(X) <- paste("cell", 1:n, sep="")
  rownames(X) <- paste("gene", 1:p, sep="")
  
  ep=0.5
  Xtrain <- apply(X, 2, function(u) rbinom(length(u), u, ep))
  Xtest <- X - Xtrain
  colnames(Xtrain) <- paste("cell", 1:n, sep="")
  rownames(Xtrain) <- paste("gene", 1:p, sep="")
  colnames(Xtest) <- paste("cell", 1:n, sep="")
  rownames(Xtest) <- paste("gene", 1:p, sep="")
  
  
  cds <- new_cell_data_set(X)
  cds_train <- new_cell_data_set(Xtrain)
  cds_test <- new_cell_data_set(Xtest)

  #### THESE ARE REQUIRED before order_cells, which the computes pseudotime
  cds <- preprocess_cds(cds)
  cds <- reduce_dimension(cds)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  
  cds_train <- preprocess_cds(cds_train)
  cds_train <- reduce_dimension(cds_train)
  cds_train <- cluster_cells(cds_train)
  cds_train <- learn_graph(cds_train)
  
  #### We got some numbers out for pseudotime so I guess that is good!
  cds <- order_cells(cds, root_cells="cell1")
  cds_train <- order_cells(cds_train, root_cells="cell1")
  
  ###
  #gene_fits_one <- fit_models(cds, model_formula_str = "~pseudotime")
  graph_test_fits <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
  
  
  ### This object now stores test set counts but training set pseudotime. 
  temp <- counts(cds_test)
  cds_test <- cds_train
  counts(cds_test) <- temp
  
  graph_test_fits_CS <-  graph_test(cds_test, neighbor_graph="principal_graph", cores=4)
  
  pvals_graph_test[i,] <- graph_test_fits$p_value
  pvals_graph_test_CS[i,] <- graph_test_fits_CS$p_value
  
  #pvalsCS[i,] <- ((coefficient_table(gene_fits_two)%>% filter(term == "pseudotime"))$p_value)
}
save(pvals_graph_test, file="pvals_graph_test.RData")
save(pvals_graph_test_CS, file="pvals_graph_test_CS.RData")



###### Now make the plot. 
load("pvals_seurat.RData")
load("pvals_seurat_CS.RData")
load("pvals_graph_test.RData")
load("pvals_graph_test_CS.RData")


library(tidyverse)
p1 <- ggplot(data=NULL)+
  geom_qq(aes(sample=as.numeric(pvals_graph_test[1:314,]), col="Double dipping"), distribution="qunif")+
  geom_qq(distribution="qunif", aes(sample=as.numeric(pvals_graph_test_CS[1:314,]), col="Count splitting"))+
  geom_abline()+
  coord_fixed()+
  labs(col="")+theme_bw()+
  xlab("Unif(0,1) Quantiles")+
  ylab("Sample Quantiles")+ggtitle("Monocle3")
p2 <- ggplot(data=NULL)+
  geom_qq(aes(sample=as.numeric(pvals), col="Double dipping"), distribution="qunif")+
  geom_qq(distribution="qunif", aes(sample=as.numeric(pvals_CS), col="Count splitting"))+
  geom_abline()+
  coord_fixed()+
  labs(col="")+theme_bw()+
  xlab("Unif(0,1) Quantiles")+
  ggtitle("Seurat")+
  ylab("Sample Quantiles")


library(patchwork)
p1+p2+plot_layout(guides="collect")
ggsave("~/Dropbox/Pseudotime : PCA NEW/Paper/Biostat Revision August 2022/Figures/Seurat_Monocle_appendix.png")





