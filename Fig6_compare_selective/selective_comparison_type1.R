library(ggplot2)
library(patchwork)
library(clusterpval)
library(fastcluster)

n=200
p=10
Lambda <- matrix(5, nrow=n,ncol=p)
nTrials <- 1000
eps <- 0.5
pvals_countsplit10 <- rep(NA, nTrials)
pvals_lucy10 <- rep(NA, nTrials)
pvals_naive10 <-  rep(NA, nTrials)

for (i in 1:nTrials) {
  print(i)
  set.seed(i)
  X <- apply(Lambda,2,function(u) rpois(length(u), u))
  hX <- log(X+1)
  
  #### Lucy method on logged data
  hcl <- hclust(dist(hX, method="euclidean")^2, method="average") 
  clusters <- cutree(hcl, k=2)
  pvals_lucy10[i] <- as.numeric(test_hier_clusters_exact(hX, link="average", K=2, k1=1, k2=2, hcl=hcl)$pval)

  ### NAIVE
  diff_means <- colMeans(hX[clusters==1,,drop='F']) - colMeans(hX[clusters==2, , drop='F'])
  stat <- sqrt(sum(diff_means^2))
  n1 <- sum(clusters==1)
  n2 <- sum(clusters==2)
  squared_norm_nu <- 1/n1 + 1/n2
  geneSSEs_train <- apply(hX, 2, function(u) sum((u - mean(u))^2))
  sigma_hat_train <- sqrt(1/(n*p-p)*sum(geneSSEs_train))
  scale_factor <- squared_norm_nu*sigma_hat_train^2
  pvals_naive10[i] <-1 - pchisq(stat^2/scale_factor, df=p)
  
  
  #### Count split method, but with multivariate T-test
  Xtrain <- apply(X,2,function(u) rbinom(n=length(u), size=u, p=eps))
  Xtest <- X-Xtrain
  hXtrain <- log(Xtrain+1)
  hXtest <- log(Xtest+1)
  hcltrain <- hclust(dist(hXtrain, method="euclidean")^2, method="average")
  clusters <- cutree(hcltrain,k=2)
  
  diff_means <- colMeans(hXtest[clusters==1,,drop='F']) - colMeans(hXtest[clusters==2, , drop='F'])
  stat <- sqrt(sum(diff_means^2))
  n1 <- sum(clusters==1)
  n2 <- sum(clusters==2)
  squared_norm_nu <- 1/n1 + 1/n2
  geneSSEs_train <- apply(hXtest, 2, function(u) sum((u - mean(u))^2))
  sigma_hat_train <- sqrt(1/(n*p-p)*sum(geneSSEs_train))
  scale_factor <- squared_norm_nu*sigma_hat_train^2
  pvals_countsplit10[i]  <- 1 - pchisq(stat^2/scale_factor, df=p)
}


naivecol <- "#E763F3"
countsplitcol <- "#E7861B"
selectivecol <- "darkblue"

ggplot(data=NULL)+geom_qq(aes(sample=pvals_lucy10, col="bSelective inference"), distribution="qunif")+
  geom_qq(aes(sample=pvals_countsplit10, col="aCount splitting"), distribution="qunif")+
  geom_qq(aes(sample=pvals_naive10, col="cDouble dipping"), distribution="qunif")+
  ggtitle("Uniform QQ plot")+
  geom_abline(a=0, b=1, col="black")+theme_bw()+coord_fixed()+
  labs(col="Method")+xlab("Unif(0,1) Quantiles")+ylab("Sample Quantiles")+
 scale_color_manual(values=c(countsplitcol, selectivecol, naivecol), 
                    labels=c("Count splitting", "Selective inference", "Double dipping"))
ggsave("~/Dropbox/Pseudotime : PCA NEW/Paper/Biostat_Resubmit_October/selective.eps")

