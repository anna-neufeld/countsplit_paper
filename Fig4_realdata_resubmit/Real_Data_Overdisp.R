library(scran)
library(tidyverse)
library(scater)
library(patchwork)
library(tricycle)
library(monocle3)
library(Seurat)
setwd("~/countsplit_paper/Fig4_realdata_resubmit/")
source("myMonocle.R")

load("Xtrain.rda")
load("Xtest.rda")
metadata <- read.delim("metadata.tsv")
rownames(metadata) <- metadata$cell
#### Genes that were never expressed (only contained 0s) in either the training
#### or test matrices got automatically removed from the data by Seurat.
### Thus, we need to make sure the dimensions match before we add Xtrain
### and Xtest together. 
cells.to.keep <- intersect(colnames(Xtrain ), colnames(Xtest))
genes.to.keep <- intersect(rownames(Xtrain), rownames(Xtest))

metadata <- metadata[cells.to.keep,]
cm.train <- Xtrain[genes.to.keep, cells.to.keep]
cm.test <- Xtest[genes.to.keep, cells.to.keep]

X <- Xtrain+Xtest

cm <- SingleCellExperiment(assays=list(counts=X))



cds_full <- myMonocle(cm, metadata, numGenes=2500)
pt_full <- as.numeric(pseudotime(cds_full))

lib.sf.cm <- librarySizeFactors(cm)
sizeFactors(cm) <- lib.sf.cm
cm <- logNormCounts(cm)
dec.cm <- modelGeneVar(cm)
hvg.cm.var <- getTopHVGs(dec.cm, fdr.threshold=0.05)[1:2500]
X <- t(counts(cm))[, hvg.cm.var]
sfs <-lib.sf.cm
p <- ncol(X)


overdisps <- data.frame(matrix(NA, nrow=p, ncol=2))
names(overdisps) <- c("means",  "nb1")
pred_means <- X

for (i in 1:p) {
  print(i)
  overdisps[i,1] <- mean(X[,i])
  ### Based on this stack overflow thread, we see that the parameter theta is b. 
  ### https://stats.stackexchange.com/questions/10419/what-is-theta-in-a-negative-binomial-regression-fitted-with-r
  try1 <- try(mod <- MASS::glm.nb(X[,i]~pt_full+offset(log(sfs))))
  if (class(try1) != "try-error") {
    overdisps[i,2] <- mod$theta
    pred_means[,i] <- predict(mod, type="response")
  }
}


pred_means_2 <- apply(pred_means, 1, function(u) u/overdisps$nb1)
#ggplot(data=NULL, aes(x=as.numeric(pred_means_2)+1),after_stat(density))+geom_histogram()+scale_x_log10()+
#theme_bw()+xlab("Overdispersion")

save.image("real_data_overdisp.RData")


set.seed(1)
### Use random indices to avoid crashing the computer. 
randomindices <- sample(1:length(as.numeric(pred_means)), size=5e6)
ggplot(data=NULL, aes(x=as.numeric(t(pred_means_2))[randomindices], y=..density..))+
  geom_histogram()+scale_x_log10()+theme_bw()+
  geom_vline(xintercept=1, col="red")+
  xlab("Overdispersion (log scale)")+ylab("")+ggtitle("Histogram of Estimated Overdispersion Values")
ggsave("~/Dropbox/Pseudotime : PCA NEW/Paper/Biostat_Resubmit_October/overdisps_cm.eps")
