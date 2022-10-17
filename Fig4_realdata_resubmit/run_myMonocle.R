library(scran)
library(tidyverse)
library(scater)
library(patchwork)
library(monocle3)
library(tricycle)
library(Seurat)
setwd("~/countsplit_paper/Fig4_realdata_resubmit/")
source("myMonocle.R")


####  Read in the training and test counts for the 10,000 cells that we wish to work with.
#### The countsplitting on the FULL dataset (~900,000 cells) was done externally by our collaborators,
#### As was the lineage subsetting.
#### Details on this external process can be found in Elorbany et al. 
load("Xtrain.rda")
load("Xtest.rda")
metadata <- read.delim("metadata.tsv")
cells.to.keep <- intersect(colnames(Xtrain ), colnames(Xtest))
genes.to.keep <- intersect(rownames(Xtrain ), rownames(Xtest))
rownames(metadata) <- metadata$cell
metadata <- metadata[cells.to.keep,]

X <- Xtrain+Xtest
Xtrain0 <- Xtrain[, metadata$diffday==0]
Xtest0 <- Xtest[, metadata$diffday==0]
X0 <- X[, metadata$diffday==0]

### Store everything in the correct type of object. 
trainSCE <- SingleCellExperiment(assays=list(counts=Xtrain))
testSCE <- SingleCellExperiment(assays=list(counts=Xtest))
fullSCE <- SingleCellExperiment(assays=list(counts=X))
trainSCE_0 <- trainSCE[, metadata$diffday==0]
testSCE_0 <-  testSCE[, metadata$diffday==0]
fullSCE_0 <-  fullSCE[, metadata$diffday==0]
metadata0 <-  metadata[metadata$diffday==0,]

#### Compute pseudotime for each of the three objects!!!
### This is a little bit slow. 
numGenes <- 2500
cds_full <- myMonocle(fullSCE, metadata, numGenes)
cds_train <- myMonocle(trainSCE, metadata, numGenes)
cds_test <- myMonocle(testSCE, metadata, numGenes)
pt_full <- as.numeric(pseudotime(cds_full), numGenes)
pt_train <- as.numeric(pseudotime(cds_train), numGenes)
pt_test <- as.numeric(pseudotime(cds_test), numGenes)
print(cor(pt_train, pt_test))
print(cor(pt_full, pt_train))
print(cor(pt_full, pt_test))


ggplot(data=NULL, aes(x=as.factor(colData(cds_full)$diffday), y=pt_full))+geom_boxplot()
ggplot(data=NULL, aes(x=as.factor(colData(cds_full)$diffday), y=pt_train))+geom_boxplot()
ggplot(data=NULL, aes(x=as.factor(colData(cds_full)$diffday), y=pt_test))+geom_boxplot()


#### Now we should subset to the day0 cells only.


cds_full0 <- myMonocle(fullSCE_0 , metadata0, numGenes)
cds_train0 <- myMonocle(trainSCE_0, metadata0, numGenes)
cds_test0 <- myMonocle(testSCE_0, metadata0, numGenes)

pt_full0 <- as.numeric(pseudotime(cds_full0))
pt_train0 <- as.numeric(pseudotime(cds_train0))
pt_test0 <- as.numeric(pseudotime(cds_test0))

cor(pt_full0, pt_train0)
cor(pt_full0, pt_test0)
cor(pt_train0, pt_test0)



#### I think the idea is we "assume size factors are fixed" when justifying a
#### Poisson model, but of course in practice we need to esitmate them, and estimating
#### on training set *is* the most consistent with general philosophy of "do things
#### downstream". 
logSFs <- log(librarySizeFactors(fullSCE))
logSFs0 <- log(librarySizeFactors(fullSCE_0))
logSFs_train <- log(librarySizeFactors(trainSCE))
logSFs_train_0 <- log(librarySizeFactors(trainSCE_0))
logSFs_test <- log(librarySizeFactors(testSCE))
logSFs_test_0 <- log(librarySizeFactors(testSCE_0))
library(speedglm)

### Notes! Even though for the purposes of *computing pseudotime* I let 
### Training and Testing set compute their own high variance genes, now here I need to compute
### p-values for the SAME set of genes for all methods to make a plot. 
### And for computational reasons I don't want to do 21,000 genes. Also-- some of the 21,000 genes are so sparse
### that the GLMs you fit are kind of dumb.
### So I ended up only fitting GLMs for the high variance genes. 
### Basically, if computation was not an issue I would just fit the GLMs to all 21000 genes (and can do this eventually).
### But for our plots (especialy in null example) it is easiest to demonstrate *bad double dipping* if we subset
### to less sparse and/or high variance genes. So that's why there is some subsetting here!!!
fullSCE <- logNormCounts(fullSCE)
dec.cm<- modelGeneVar(fullSCE)
hvg.cm.var <- getTopHVGs(dec.cm,fdr.threshold = 0.001)[1:2500]
cm0 <- logNormCounts(fullSCE_0)
dec.cm0 <- modelGeneVar(cm0 )
hvg.cm.var0 <- getTopHVGs(dec.cm0, fdr.threshold = 0.001)[1:2500]


pvals_countsplit_full_hand <- apply(Xtest[hvg.cm.var, ], 1, function(u) summary(speedglm(u~pt_train, offset=logSFs_train,  maxit=20, family=poisson("log")))$coefficients[2,4])
pvals_full_naive_hand <- apply(X[hvg.cm.var,], 1, function(u) summary(speedglm(u~pt_full, offset=logSFs, maxit=20,  family=poisson("log")))$coefficients[2,4])
pvals_test_naive_hand <- apply(Xtest[hvg.cm.var, ], 1, function(u) summary(speedglm(u~pt_test, offset=logSFs_test, maxit=20, family=poisson("log")))$coefficients[2,4])
pvals_train_naive_hand <- apply(Xtrain[hvg.cm.var, ], 1, function(u) summary(speedglm(u~pt_train, offset=logSFs_train, maxit=20, family=poisson("log")))$coefficients[2,4])

pvals_countsplit_0 <- apply(Xtest0[hvg.cm.var0, ], 1, function(u) summary(speedglm(u~pt_train0, offset=logSFs_train_0, family=poisson("log")))$coefficients[2,4])
pvals_full_naive0 <- apply(X0[hvg.cm.var0, ], 1, function(u) summary(speedglm(u~pt_full0, offset=logSFs0, family=poisson("log")))$coefficients[2,4])
pvals_test_naive0 <- apply(Xtest0[hvg.cm.var0, ], 1, function(u) summary(speedglm(u~pt_test0, offset=logSFs_test_0, family=poisson("log")))$coefficients[2,4])
pvals_train_naive0 <- apply(Xtrain0[hvg.cm.var, ], 1, function(u) summary(speedglm(u~pt_train0, offset=logSFs_train_0, maxit=20, family=poisson("log")))$coefficients[2,4])


save.image("monocle_workspace_Oct17.RData")

library(tidyverse)
### Really affected by ridiculuously small numbers.
p1 <- ggplot(data=NULL)+
  geom_point(aes(x=pvals_countsplit_full_hand, y=pvals_full_naive_hand, col="Full Data"))+
  geom_point(aes(x=pvals_countsplit_full_hand, y=pvals_test_naive_hand, col="Test Data"))+
  geom_abline(a=0,b=1)+
  scale_x_log10()+scale_y_log10()+
  theme_bw()+coord_fixed()+xlab("Count Splitting p-values")+ylab("Naive p-values")+
  ggtitle("All Days")+
  coord_fixed()+
  labs(col="Naive gets to use")
p1

#### This ignores all p-values that are less than 2e-16.
p2 <- ggplot(data=NULL)+
  geom_point(aes(x=pvals_countsplit_full_hand, y=pvals_full_naive_hand, col="Full Data"))+
  geom_point(aes(x=pvals_countsplit_full_hand, y=pvals_test_naive_hand, col="Test Data"))+
  geom_abline(a=0,b=1)+
  xlim(2e-16,1)+ylim(2e-16,1)+
  scale_x_log10(limits=c(2e-16,1))+scale_y_log10(limits=c(2e-16,1))+
  theme_bw()+coord_fixed()+xlab("Count Splitting p-values")+ylab("Naive p-values")+
  ggtitle("All Days")+
  coord_fixed()+
  labs(col="Naive gets to use")
p2

#### This jitters all p-values that are less than numerical precision!!!
pvals_countsplit_full_hand2 <- pvals_countsplit_full_hand
pvals_full_naive_hand2 <- pvals_full_naive_hand
pvals_test_naive_hand2 <- pvals_test_naive_hand

pvals_countsplit_full_hand2[pvals_countsplit_full_hand < 2e-16] <- jitter(rep(2e-16, sum(pvals_countsplit_full_hand < 2e-16)), amount = 1e-30)
pvals_full_naive_hand2[pvals_full_naive_hand < 2e-16] <- jitter(rep(2e-16, sum(pvals_full_naive_hand < 2e-16)), amount = 1e-30)
pvals_test_naive_hand2[pvals_test_naive_hand < 2e-16] <- jitter(rep(2e-16, sum(pvals_test_naive_hand< 2e-16)), amount = 1e-30)

p3 <- ggplot(data=NULL)+
  geom_point(aes(x=pvals_countsplit_full_hand2, y=pvals_full_naive_hand2, col="Full Data"))+
  geom_point(aes(x=pvals_countsplit_full_hand2, y=pvals_test_naive_hand2, col="Test Data"))+
  geom_abline(a=0,b=1)+
  scale_x_log10()+scale_y_log10()+
  theme_bw()+coord_fixed()+xlab("Count splitting p-values")+ylab("Double dipping p-values")+
  ggtitle("All Days")+
  coord_fixed()+
  labs(col="Double dipping method uses:")
p3

#### NULL STORY!!!

### Make sure we don't include any genes that are just literally all 0s.
good_indices <- rowSums(X0[hvg.cm.var0,]) > 1
cols <- c("#00B81F", "#00A5FF", "#F8766D")
p4 <- ggplot(data=NULL)+
  geom_qq(aes(sample=pvals_countsplit_0[good_indices], col="Count spliting"), distribution="qunif")+
  geom_qq(aes(sample=pvals_full_naive0[good_indices], col="Full double dipping"), distribution="qunif")+
  geom_qq(aes(sample=pvals_test_naive0[good_indices], col="Test double dipping"), distribution="qunif")+
  geom_abline()+theme_bw()+coord_fixed()+xlab("Unif(0,1) Quantiles")+ylab("Sample Quantiles")+
  labs(col="Method")+
  scale_color_manual(values=cols[c(3,1,2,4)])+
  ggtitle("Day 0 Only")
p4


library(patchwork)
p3+p4
ggsave(filename="~/Dropbox/Pseudotime : PCA NEW/Paper/Biostat_Resubmit_October/monocle_both_NEW.eps")

