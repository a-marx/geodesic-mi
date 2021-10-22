#!/usr/bin/Rscript

source("volume_estimation.R")
source("marginale_forest_distance.R")

adjustDistance <- function(X, Xdist, epsilon=0.001){
  for(i in 1:nrow(Xdist)){
    X0.indicator = (Xdist[i,] == 0)
    min.d = min(Xdist[i,!X0.indicator])
    X0.indicator[i] = F # remove element itself
    nX0 = sum(X0.indicator)
    X0 = c(c(1:nrow(Xdist))[X0.indicator])
    if(nX0 > 0){
      X0.l2 <- sapply(X0, function (k)
      {
        l2dist(X[i,], X[k,])
      })
      max.l2 = max(X0.l2)
      X0.l2.adjusted = X0.l2 / max.l2 * (min.d * (1-epsilon))
      Xdist[i,X0] = X0.l2.adjusted
    }
  }
  return(Xdist)
}

mi_ksg_estimate_forest <- function(X, Y, k, geodesic=T, fm=NULL, Xdist=NULL, Ydist=NULL, trees=300, normalizeData = F, splitCrit="fastbic", pathSimilarity=T)
{
  N <- nrow(X)
  if (geodesic)
  {
    m <- "geodesic"
  } else {
    m <- "euclidean"
  }
  # compute distances
  Z <- cbind(X, Y)
  # distance in joint space
  if (geodesic)
  {
    if(is.null(fm))
    {
      fm <- Urerf(Z, trees=trees, normalizeData=normalizeData, splitCrit=splitCrit)
    }
    if(is.null(Xdist))
    {
      Xdist <- marginale_forest_distance(X=Z, forest=fm$forest, b=c( rep(1, dim(X)[2]), rep(0,dim(Y)[2]) ), normalizeData=normalizeData)
    }
    if(is.null(Ydist))
    {
      Ydist <- marginale_forest_distance(X=Z, forest=fm$forest, b=c( rep(0, dim(X)[2]), rep(1,dim(Y)[2]) ), normalizeData=normalizeData)
    }
  } else {
    Xdist <- as.matrix(dist(X, method="maximum"))
    Ydist <- as.matrix(dist(Y, method="maximum"))
  }
  if(geodesic){ ## adjust for larger set of nodes always in the same leaf
    Xdist = adjustDistance(X,Xdist)
    Ydist = adjustDistance(Y,Ydist)
  }
  Zdist <- do.call(rbind, lapply(1:nrow(Xdist), function (i)
  {
    pmax(Xdist[i,], Ydist[i,])
  }))
  Zdist <- t(apply(Zdist,1,sort))
  # sort distances for faster computation
  MI <- digamma(k) + digamma(N)
  MI <- MI - 1/N*sum(sapply(1:N, function(i)
  {
    # determine epsilon_i (distance of sample i to kth neighbour in joint space)
    eps_i <- Zdist[i,k+1]
    # count number of samples with less distance in each subspace
    n_xi <- sum(sapply(1:N, function (j) {Xdist[i,j] < eps_i})) - 1
    n_yi <- sum(sapply(1:N, function (j) {Ydist[i,j] < eps_i})) - 1
    return(digamma(n_xi + 1) + digamma(n_yi + 1))
  }
  ))
  return(MI)
}
