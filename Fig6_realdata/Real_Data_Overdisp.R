library(scran)
library(tidyverse)
library(scater)
library(patchwork)
library(tricycle)
library(monocle3)

setwd("~/Dropbox/Pseudotime : PCA NEW/Paper/FinalCode/Monocle_CM_Work/")
source("myMonocle.R")
cm <- readRDS("sce_subset.rds")
metadata <- read.delim("metadata_subset.tsv")
metadata$diffday <- ordered(metadata$diffday,levels=c(
  "day0", "day1", "day3", "day5", "day7", "day11", "day15")
)
rownames(metadata) <- metadata$cell

X <- t(counts(cm))
cds_full <- myMonocle(cm, metadata)
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
  ### Based on this stack overflow thread, I am almost positive that theta is b, not
  ### 1/b. https://stats.stackexchange.com/questions/10419/what-is-theta-in-a-negative-binomial-regression-fitted-with-r
  try1 <- try(mod <- MASS::glm.nb(X[,i]~pt_full+offset(log(sfs))))
  if (class(try1) != "try-error") {
    overdisps[i,2] <- mod$theta
    pred_means[,i] <- predict(mod, type="response")
  }
}


pred_means_2 <- apply(pred_means, 1, function(u) u/overdisps$nb1)
ggplot(data=NULL, aes(x=as.numeric(pred_means_2)+1),after_stat(density))+geom_histogram()+scale_x_log10()+
theme_bw()+xlab("Overdispersion")


set.seed(1)
### Use random indices to avoid crashing the computer. 
randomindices <- sample(1:length(as.numeric(pred_means)), size=5e6)
ggplot(data=NULL, aes(x=as.numeric(t(pred_means_2))[randomindices], y=..density..))+
  geom_histogram()+scale_x_log10()+theme_bw()+
  geom_vline(xintercept=1, col="red")+
  xlab("Overdispersion (log scale)")+ylab("")+ggtitle("Histogram of Estimated Overdispersion Values")
