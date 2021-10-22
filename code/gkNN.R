#!/usr/bin/Rscript

source("volume_estimation.R")

## Directly compute the entropy
entropy_estimate_gknn <- function(X, k, Xdist)
{
  require("FastKNN")
  N <- nrow(X)
  d <- ncol(X)
  Xknn <- matrix(0,N,k)
  for (i in 1:N){
    Xknn[i,] <- k.nearest.neighbors(i, Xdist, k = k)
  }
  # first compute "constants"
  H_X <- log(N) + log(pi^(d/2)/gamma(1+d/2))
  # compute distances between each data point and its kth neighbour
  H_X <- H_X + d/N*sum(sapply(1:N, function (i)
  {
    log(l2dist(X[i,], X[Xknn[i,k],]))
  }))
  H_X <- H_X + 1/N*sum(sapply(1:N, function (i)
  {
    # as.matrix is necessary for 1D data
    kneighood <- NULL
    if(d == 1){
      kneighood <- rbind(as.matrix(X[i,]), as.matrix(X[Xknn[i,],]))
    }else{
      kneighood <- rbind(X[i,], X[Xknn[i,],])
    }
    centroid <- colMeans(kneighood)
    Y_i <- do.call(rbind, lapply(1:nrow(kneighood), function (rIdx)
    {kneighood[rIdx,] - centroid}))
    # compute the SVD and retrieve its singular values
    svd_Yi <- svd(Y_i)
    sing_Yi <- svd_Yi$d
    # compute the number of points that lie in ellipsoid
    #   first center the values at ellipsoid center (which is now x_i)
    Z_i <- do.call(rbind, lapply(1:k, function (j)
    {
        X[Xknn[i, j], ] - X[i, ]
    }))
    #   compute number of points, at least kth neighbour is contained, 
    #   we don't have to check for that one
    k_i <- sum(sapply(1:nrow(Z_i), function (j)
    {
      hyperellipsoid_check(svd_Yi, Z_i[j,])
    }))
    res <- -log(k_i)
    
    # compute the fractions of sv to determine volume
    res <- res + sum(sapply(1:d, function(l)
    {
      log(sing_Yi[l]/sing_Yi[1])
    }
    ))
    return(res)
  }))
  return(H_X)
}
mi_gKNN <- function(X, Y, k){
    Xdist = as.matrix(dist(X))
    Ydist = as.matrix(dist(Y))
    XYdist = as.matrix(dist(cbind(X,Y)))
    hx = entropy_estimate_gknn(X,k,Xdist)
    hy = entropy_estimate_gknn(Y,k,Ydist)
    hxy = entropy_estimate_gknn(cbind(X,Y),k,XYdist)
    return(hx + hy - hxy)
}
