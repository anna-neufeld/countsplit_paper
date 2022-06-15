gen_Lambda <- function(n=200,p=100, k=1, Fs, intercepts) {
  if (k==1) {
    L <- rnorm(n)
    L <- L - mean(L)
  }
  
  LFT <- L%*%t(Fs)
  logLambda <- t(apply(LFT, 1, function(u) u+intercepts))
  Lambda <- exp(logLambda)
  
  return(Lambda)
  return(list(dat=X,mean=EX, coeffs=Fs, intercepts=intercepts))
}


gen_pois_data <- function(gammas, Lambda) {
  EX <- apply(Lambda, 2, function(u) u*gammas)
  X <- apply(EX,2,function(u) rpois(length(u),u))
  return(X)
}


one_trial_count_split <- function(n,p,filename, k=1,propImp=0.1, sig_strength=5, propLowMedHigh = c(1/2,1/2), eps=c(0.5), c=1) {
  
                                        
    
  #intercepts <- sampl, log(25), log(100)), prob=propLowMedHigh, replace=TRUE,size=p)

  intercepts <- sample(c(log(3), log(25)), prob=propLowMedHigh, replace=TRUE,size=p)
  num_non_null <- floor(propImp*p)
  numNull <- p-num_non_null
  
  Fs <- c(rep(sig_strength, num_non_null), rep(0, numNull))
  #Fs <- rep(c(-1,1), p/2)*Fs
  
  Lambda <- gen_Lambda(n,p,k,Fs,intercepts)
  gammas <- rgamma(n, 10, scale=1/10)
  
  X <- gen_pois_data(gammas,Lambda)
  
  for (ep in eps) {
    res <- cbind(1:p, Fs, intercepts, countsplit(X, Lambda, ep, c, gammas), ep, "known", propImp, n,p, propLowMedHigh[1])
    write(t(res), file=filename, append=TRUE, ncol=12)
  }
  
}


countsplit <- function(X,Lambda, ep,c=1, gammas=rep(1,NROW(X))) {
  Xtrain <- apply(X,2,function(u) rbinom(n=length(u), size=u, p=ep))
  Xtest <- X-Xtrain
  
  hXtrain <- log(diag(1/gammas)%*%(Xtrain+c))
  hXtraincenter <- apply(hXtrain,2,function(u) u-mean(u))
  

  #pseudotime <- princomp(hXtrain)$scores[,1]
  pseudotime <- svd(hXtraincenter)$u[,1]
  
  # Could it possibly matter here whether or not I put a 1-eps in the offset??? 
  pvals_pseudotime <- apply(Xtest, 2, function(u) summary(glm(u~pseudotime, offset=log(gammas), family="poisson"))$coefficients[2,4])

  
  true_coeffs <- suppressWarnings(apply(Lambda, 2, function(u) summary(glm(u~pseudotime, family="poisson"))$coefficients[2,1]))


  trueprincomp <- svd(apply(log(Lambda),2,function(u) u-mean(u)))$u[,1]

  true_cor <- cor(pseudotime, trueprincomp)
  
  return(cbind(pvals_pseudotime,true_coeffs, true_cor))
}


countsplit_cluster <- function(X,Lambda, ep,c=1, gammas=rep(1,NROW(X)), trueclusters=rep(1,NROW(X))) {
  Xtrain <- apply(X,2,function(u) rbinom(n=length(u), size=u, p=ep))
  Xtest <- X-Xtrain
  
  hXtrain <- log(diag(1/gammas)%*%(Xtrain+c))
  
  clusters <- as.factor(kmeans(apply(hXtrain,2,function(u) u-mean(u)), 2)$cluster)
  
  pvals_pseudotime <- apply(Xtest, 2, function(u) summary(glm(u~clusters, offset=log(gammas), family="poisson"))$coefficients[2,4])

  #### CONSIDER REMOVING 1-EPS FROM OFFSET
  true_coeffs <- suppressWarnings(apply(Lambda, 2, function(u) summary(glm(u~clusters, family="poisson"))$coefficients[2,1]))
  
  true_cor <- mclust::adjustedRandIndex(clusters, trueclusters)
  
  return(cbind(pvals_pseudotime,true_coeffs, true_cor))
}



one_trial_count_split_cluster <- function(n,p,filename, k=1,propImp=0.1, sig_strength=5, propLowMedHigh = c(1/3,1/3,1/3,1/3,1/3), eps=c(0.5), c=1, props=0.5) {
  
  #intercepts <- sample(c(log(1), log(3),log(25),log(100)), prob=propLowMedHigh, replace=TRUE, size=p)
  intercepts <- sample(c(log(3), log(25)), prob=propLowMedHigh, replace=TRUE,size=p)
  num_non_null <- floor(propImp*p)
  numNull <- p-num_non_null
  
  L <- sample(0:1, size=n, props)
  
  Fs <- c(rep(sig_strength, num_non_null), rep(0, numNull))
  #Fs <- rep(c(-1,1), p/2)*Fs
  
  LFT <- L%*%t(Fs)
  logLambda <- t(apply(LFT, 1, function(u) u+intercepts))
  Lambda <- exp(logLambda)
  
  gammas <- rgamma(n, 10, scale=1/10)
  
  #table(kmeans(apply(log(Lambda), 2, function(u) u-mean(u)), centers=2)$cluster, L)
  
  X <- gen_pois_data(gammas,Lambda)
  
  for (ep in eps) {
    res <- cbind(1:p, Fs, intercepts, countsplit_cluster(X, Lambda, ep, c, gammas, L), ep, "known", propImp, n,p, propLowMedHigh[1])
    write(t(res), file=filename, append=TRUE, ncol=12)
  }
  
}

