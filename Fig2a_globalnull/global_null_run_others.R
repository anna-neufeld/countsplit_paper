setwd("~/countsplit_paper/Fig1a_globalnull/")
source("myjackstraw.R")

n <- 200
p <- 10
c <- 1

library(tidyverse)

Lambda <- matrix(c(rep(1, n*(p/2)), rep(10, n*(p/2))), nrow=n, ncol=p)
nTrials <- 1000

pvals_pseudotime <- data.frame(matrix(NA, nrow=nTrials*p, ncol=6))
names(pvals_pseudotime) <- c("Naive", "CellSplitting", "CountSplitting", "Jackstraw", "Intercept")

pvals_pseudotime$Intercept <- rep(Lambda[1,], nTrials)


counter <- 1

for (i in 1:nTrials) {
  set.seed(i)
  print(i)
  
  X <- apply(Lambda,2,function(u) rpois(length(u), u))
  hX <- log(X+c)
  
  #### NAIVE
  pseudotime <- princomp(hX)$scores[,1]
  pvals_pseudotime[counter:(counter+p-1),1] <- apply(X, 2, function(u) summary(glm(u~pseudotime, family="poisson"))$coefficients[2,4])

  #### CELL SPLITTING
  train <- sample(1:n, size=n/2)
  
  Xfullcenter <- apply(hX,2,function(u) u -mean(u))
  svdtrain <- svd(Xfullcenter[train,])
  testpt <- Xfullcenter[-train,]%*%svdtrain$v[,1]
  pvals_pseudotime[counter:(counter+p-1),2] <- apply(X[-train,], 2, function(u) summary(glm(u~testpt, family="poisson"))$coefficients[2,4])

  #### COUNT SPLITTING
  eps <- 0.5
  Xtrain <- apply(X,2,function(u) rbinom(n=length(u), size=u, p=eps))
  Xtest <- X-Xtrain
  
  hXtrain <- log(Xtrain+c)
  pseudotime <- princomp(hXtrain)$scores[,1]
  pvals_pseudotime[counter:(counter+p-1),3] <- apply(Xtest, 2, function(u) summary(glm(u~pseudotime, family="poisson"))$coefficients[2,4])

  
  #### JACKSTRAW
  pvals_pseudotime[counter:(counter+p-1),4] <- myjackstraw.z(X, Ltype="pseudotime",s=10,B=100)
  counter <- counter+p
}

save(pvals_pseudotime, file="global_null_res_may_23.RData")

