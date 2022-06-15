Ftest.glm<- function(u, L) {
  return(anova(glm(u~1, family="poisson"),glm(u~L, family="poisson"))$Deviance[2])
}

ztest.glm<- function(u, L) {
  return(summary(glm(u~L))$coefficients[2,3])
}

Ftest.lm<- function(u, L) {
  return(anova(lm(u~L))$"F value"[1])
}


myjackstraw.z <- function(X, s=10, B=100, c=1, Ltype="pseudotime", filename="test") {
  
  p <- NCOL(X)
  n <- NROW(X)
  
  
  Xtrans <- log(X+c)
  Xfullcenter <- apply(Xtrans, 2, function(u) (u-mean(u)))
  
  if (Ltype=="pseudotime") {
    Lhat <- svd(Xfullcenter)$u[,1]
  } else {
    Lhat = as.factor(kmeans(Xfullcenter,centers=2)$cluster)
  }
  
  origFs <- apply(X,2,ztest.glm, Lhat) 
  
  empiricalnull <- rep(0, s*B)
  for (i in 1:B) {
    permuteVars <- sample(1:p, size=s)
    Xpermute <- X
    Xpermute[,permuteVars] <- apply(Xpermute[,permuteVars], 2, function(u) sample(u, size=length(u)))
    Xtranspermute <-log(Xpermute+c)
    Xfullcenterpermute <- apply(Xtranspermute, 2, function(u) (u-mean(u)))
    
    if (Ltype=="pseudotime") {
      Lhatpermute <- svd(Xfullcenterpermute)$u[,1]
    } else {Lhatpermute <- kmeans(Xfullcenterpermute,centers=2)$cluster}
    
    empiricalnull[((i-1)*s+1):(i*s)]<- apply(Xpermute[,permuteVars],2,ztest.glm, Lhatpermute) 
  }
  
  pvals <- sapply(origFs, function(u) mean(abs(empiricalnull) > abs(u)))
  return(pvals)
}

myjackstraw <- function(X, s=10, B=100, c=1, Ltype="pseudotime", filename="test") {
  
  p <- NCOL(X)
  n <- NROW(X)
  

  Xtrans <- log(X+c)
  Xfullcenter <- apply(Xtrans, 2, function(u) (u-mean(u)))

  if (Ltype=="pseudotime") {
    Lhat <- svd(Xfullcenter)$u[,1]
  } else {
    Lhat = as.factor(kmeans(Xfullcenter,centers=2)$cluster)
  }
 
  origFs <- apply(X,2,Ftest.glm, Lhat) 
  
  empiricalnull <- rep(0, s*B)
  for (i in 1:B) {
    permuteVars <- sample(1:p, size=s)
    Xpermute <- X
    Xpermute[,permuteVars] <- apply(Xpermute[,permuteVars], 2, function(u) sample(u, size=length(u)))
    Xtranspermute <-log(Xpermute+c)
    Xfullcenterpermute <- apply(Xtranspermute, 2, function(u) (u-mean(u)))
    
    if (Ltype=="pseudotime") {
      Lhatpermute <- svd(Xfullcenterpermute)$u[,1]
    } else {Lhatpermute <- kmeans(Xfullcenterpermute,centers=2)$cluster}
    
    empiricalnull[((i-1)*s+1):(i*s)]<- apply(Xpermute[,permuteVars],2,Ftest.glm, Lhatpermute) 
  }
  
  pvals <- sapply(origFs, function(u) mean(empiricalnull > u))
  return(pvals)
}
