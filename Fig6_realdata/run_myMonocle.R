library(scran)
library(tidyverse)
library(scater)
library(patchwork)
library(monocle3)
library(tricycle)



#### Read in data !
setwd("~/countsplit_paper/Fig6_realdata/")
source("myMonocle.R")
cm <- readRDS("sce_subset.rds")
metadata <- read.delim("metadata_subset.tsv")
metadata$diffday <- ordered(metadata$diffday,levels=c(
  "day0", "day1", "day3", "day5", "day7", "day11", "day15")
)
rownames(metadata) <- metadata$cell

#### Count splitting on all 7 days data!
epsilon=0.5
set.seed(1)
X <- t(counts(cm))
Xtrain <-apply(X,2,function(u) rbinom(n=length(u), size=u, p=epsilon))
rownames(Xtrain) <- rownames(X)
colnames(Xtrain) <- colnames(X)
Xtest <- X-Xtrain
rownames(Xtest) <- rownames(X)
colnames(Xtest) <- colnames(X)
### Store new counts in SCE objects for convenience. 
cm_test <- cm
cm_train <- cm
counts(cm_test) <- t(Xtest)
counts(cm_train) <- t(Xtrain)

#### Compute pseudotime for each of the three objects!!!
### This is a little bit slow. 
cds_full <- myMonocle(cm, metadata)
cds_train <- myMonocle(cm_train, metadata)
cds_test <- myMonocle(cm_test, metadata)
pt_full <- as.numeric(pseudotime(cds_full))
pt_train <- as.numeric(pseudotime(cds_train))
pt_test <- as.numeric(pseudotime(cds_test))

### Sanity checks!! Since there is real signal, we should have correlation between full/train/test.
### Also if pseudotime doesn't track pretty well with day then we did something wrong.
### Here it is fine but not amazing. 
cor(pt_full, pt_train)
cor(pt_full, pt_test)
ggplot(data=NULL, aes(x=colData(cds_full)$diffday, y=pt_full))+geom_boxplot()
ggplot(data=NULL, aes(x=colData(cds_full)$diffday, y=pt_train))+geom_boxplot()
ggplot(data=NULL, aes(x=colData(cds_full)$diffday, y=pt_test))+geom_boxplot()


#### NOW DO COUNT SPLITTING on DAY 0 CELLS  
cm0 <- cm[,cm$diffday=="day0" & cm$type == "IPSC"]
metadata0 <- metadata[cm$diffday=="day0" & cm$type == "IPSC",]
epsilon=0.5
set.seed(1234)
X0 <- t(counts(cm0))
Xtrain0 <- apply(X0,2,function(u) rbinom(n=length(u), size=u, p=epsilon))
Xtest0 <- X0-Xtrain0
rownames(Xtrain0) <- rownames(X0)
colnames(Xtrain0) <- colnames(X0)
rownames(Xtest0) <- rownames(X0)
colnames(Xtest0) <- colnames(X0)
cm0_train <- cm0
cm0_test <- cm0
counts(cm0_train) <- t(Xtrain0)
counts(cm0_test) <- t(Xtest0)

#### Get the pseudotime for the Day0 Cells
cds_full0 <- myMonocle(cm0, metadata0)
cds_train0 <- myMonocle(cm0_train, metadata0)
cds_test0 <- myMonocle(cm0_test, metadata0)
pt_full0 <- as.numeric(pseudotime(cds_full0))
pt_train0 <- as.numeric(pseudotime(cds_train0))
pt_test0 <- as.numeric(pseudotime(cds_test0))



#### I think the idea is we "assume size factors are fixed" when justifying a
#### Poisson model, but of course in practice we need to esitmate them, and estimating
#### on training set *is* the most consistent with general philosophy of "do things
#### downstream". 
logSFs <- log(colData(cds_full)$Size_Factor)

#theirSFs <- colData(cds_full)$Size_Factor
#mySFs0 <- colSums(counts(cds_full))
#mySFs <- mySFs0/exp(mean(log(mySFs0)))
#mySFs2 <- mySFs0/mean(mySFs0)
#plot(theirSFs, mySFs)
#cor(theirSFs, mySFs)
#all.equal(mySFs2, librarySizeFactors(cm))


logSFs0 <- log(colData(cds_full0)$Size_Factor)
logSFs_train <- log(librarySizeFactors(cm_train))
logSFs_train_0 <- log(librarySizeFactors(cm0_train))
logSFs_test <- log(librarySizeFactors(cm_test))
logSFs_test_0 <- log(librarySizeFactors(cm0_test))
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
cm <- logNormCounts(cm)
dec.cm<- modelGeneVar(cm)
hvg.cm.var <- getTopHVGs(dec.cm, fdr.threshold=0.01)[1:2500]
cm0 <- logNormCounts(cm0)
dec.cm0<- modelGeneVar(cm0)
hvg.cm.var0 <- getTopHVGs(dec.cm0, fdr.threshold=0.005)[1:2500]

pvals_countsplit_full_hand <- apply(Xtest[, hvg.cm.var ], 2, function(u) summary(speedglm(u~pt_train, offset=logSFs_train,  maxit=20, family=poisson("log")))$coefficients[2,4])
pvals_full_naive_hand <- apply(X[, hvg.cm.var ], 2, function(u) summary(speedglm(u~pt_full, offset=logSFs, maxit=20,  family=poisson("log")))$coefficients[2,4])
pvals_test_naive_hand <- apply(Xtest[, hvg.cm.var ], 2, function(u) summary(speedglm(u~pt_test, offset=logSFs_test, maxit=20, family=poisson("log")))$coefficients[2,4])
pvals_countsplit_0<- apply(Xtest0[,hvg.cm.var0 ], 2, function(u) summary(speedglm(u~pt_train0, offset=logSFs_train_0, family=poisson("log")))$coefficients[2,4])
pvals_full_naive0 <- apply(X0[, hvg.cm.var0 ], 2, function(u) summary(speedglm(u~pt_full0, offset=logSFs0, family=poisson("log")))$coefficients[2,4])
pvals_test_naive0 <- apply(Xtest0[, hvg.cm.var0 ], 2, function(u) summary(speedglm(u~pt_test0, offset=logSFs_test_0, family=poisson("log")))$coefficients[2,4])


save.image("monocle_workspace.RData")

