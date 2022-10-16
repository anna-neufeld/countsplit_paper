suppressPackageStartupMessages(library(PseudotimeDE))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
setwd("~/countsplit_paper/Fig1a_globalnull/")


#### Simulate poisson data that is as boring as possible.
#### Every single element is just a poisson RV with mean 4.
n <- 200
p <- 10
c <- 1
library(tidyverse)
Lambda <- matrix(c(rep(1, n*(p/2)), rep(10, n*(p/2))), nrow=n, ncol=p)


nTrials <- 50

pvals_DE_pt <- data.frame(matrix(NA, nrow=nTrials*p, ncol=3))
names(pvals_DE_pt) <- c("DE", "Intercept", "Setting")
pvals_DE_cluster <- data.frame(matrix(NA, nrow=nTrials*p, ncol=3))
names(pvals_DE_cluster) <- c("DE", "Intercept", "Setting")

pvals_DE_pt$Intercept <- rep(Lambda[1,], nTrials)
pvals_DE_cluster$Intercept <- rep(Lambda[1,], nTrials)

counter <- 1

for (i in 1:nTrials) {
  print(i)
  set.seed(i)
  X <- apply(Lambda,2,function(u) rpois(length(u), u))
  pseudo_orig <- tibble(cell = 1:n, pseudotime = princomp(log(X+c))$scores[,1])
  cluster_orig <- tibble(cell = 1:n, pseudotime = kmeans(log(X+c), centers=2)$cluster)

  subsets <- lapply(1:100, function(x) {
    sample(x = c(1:dim(X)[1]), size = 0.8*dim(X)[1], replace = FALSE)
  })

  pseudo_sub <- lapply(subsets, function(indices, data) {
    data <- data[indices,]
    tbl <- tibble(cell = indices, pseudotime = princomp(log(data+1))$scores[,1])
  
    ## Make sure the direction of pseudotime is the same as the original pseudotime
    merge.tbl <- left_join(tbl, pseudo_orig, by = "cell")
    if(cor(merge.tbl$pseudotime.x, merge.tbl$pseudotime.y) < 0) {
      tbl <- dplyr::mutate(tbl, pseudotime = -1*pseudotime)
    }
    tbl
  }, data=X)
  
  cluster_sub <- lapply(subsets, function(indices, data) {
    data <- data[indices,]
    tbl <- tibble(cell = indices, pseudotime = kmeans(log(data+1), centers=2)$cluster)
    #merge.tbl <- left_join(tbl, cluster_orig, by = "cell")
    #if(cor(merge.tbl$pseudotime.x, merge.tbl$pseudotime.y) < 0) {
    #  tbl <- dplyr::mutate(tbl, pseudotime = -1*pseudotime)
    #}
    tbl
  }, data=X)
  
  rownames(X) <- paste("cell", 1:n, sep="")
  colnames(X) <- paste("gene", 1:p, sep="")
  X <- t(X)
  set.seed(5)
  res <- PseudotimeDE::runPseudotimeDE(gene.vec =  paste("gene", 1:p, sep=""),
                                       ori.tbl = pseudo_orig,
                                       sub.tbl = pseudo_sub, 
                                       mat = X, 
                                       model = "nb")
  res2 <- PseudotimeDE::runPseudotimeDE(gene.vec =  paste("gene", 1:p, sep=""),
                                       ori.tbl = cluster_orig,
                                       sub.tbl = cluster_sub, 
                                       mat = X, 
                                       model = "nb")
  
  
  pvals_DE_pt[counter:(counter+p-1),1] <- res$emp.pv
  pvals_DE_cluster[counter:(counter+p-1),1] <- res2$emp.pv
  counter <- counter+p
  
  save(pvals_DE_pt, file="pvals_DE_pt.RData")
  save(pvals_DE_cluster, file="pvals_DE_cluster.RData")
}



