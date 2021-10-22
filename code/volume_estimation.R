#!/usr/bin/Rscript

source("kNN_estimation.R")

hyperellipsoid_check <- function(svd, x)
{
  # project x into coordinate system given by svd, i.e. the axes of the ellipse
  x_proj <- NULL
  if(length(diag(svd$d)) == 1){
    x_proj <- x %*% (svd$v)
  }else{
    x_proj <- x %*% (svd$v)
  }
  # check if the projected x lies within boundaries, due to projection we can
  # use the "simple" inequality check
  return(sum(sapply(1:length(x_proj), function(i)
  {
    (x_proj[i] / svd$d[i])^2
  })) <= 1)
}

l2dist <- function (x, y)
{
  return(sqrt(sum((x-y)^2)))
}

# @data   full data matrix (rows are samples, columns are features)
# @idx    index of sample x_i we consider
# @knns   the k nearest neighbours of x_i
computeVolEllipsoid <- function(data, idx, knns)
{
  d <- ncol(data)
  k <- length(knns)
  # compute k-neighbourhood
  kneighood <- rbind(data[idx,], data[knns,])
  # center data
  centroid <- 1/(k+1) * colSums(kneighood)
  Y_i <- do.call(rbind, lapply(1:nrow(kneighood), function (rIdx)
    {kneighood[rIdx,] - centroid}))
  # compute the l2 distance between point and its k nearest neighbour
  eps_i <- l2dist(kneighood[1,], kneighood[nrow(kneighood),])
  # compute the SVD and retrieve its singular values in decreasing order
  sing_Yi <- svd(Y_i)$d
  # compute volume
  #       for the whole dataset
  V_i <- pi^(d/2)/gamma(1+d/2)
  V_i <- V_i * eps_i * prod(sapply(1:d, function (l) {sing_Yi[l]/sing_Yi[1]}))
  return (V_i)
}

require(mvtnorm)
test <- function()
{
  # generate 50 points
  X <- rmvnorm(50, c(0,0,1), matrix(c(1,0.3,0.8,   0.3,1,0.5,   0.8,0.5,1), nrow=3))
  knns <- getKNNMatrix(X, 5, geodesic=F)
  print(knns)
  Vol1 <- computeVolEllipsoid(X, 1, knns[1,])
  print(Vol1)
  Vol2 <- computeVolEllipsoid(X, 2, knns[2,])
  print(Vol2)
}

