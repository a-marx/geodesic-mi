#!/usr/bin/Rscript

# maybe require("RANN") ## finds kNN in N log N runtime; fast implementation
# maybe use class RcppAnnoy --- calculates approximate nn based on forests
# maybe "Seurat" package for nearest neighbour calculation with any distance matrix
require("FastKNN")
require("rerf")
# needed is also packages "rerf" for forests (see https://github.com/neurodata/SPORF/tree/staging/R-Project)
# installation instructions:
# get latest valid version of rerf at https://cran.r-project.org/src/contrib/Archive/rerf/
# install tar.gz via install_local (remotes package)

# @data data matrix expected to be of type matrix
getKNNMatrix <- function(data, k=10, geodesic=T){
  n <- dim(data)[1]
  dm <- NULL
  if(geodesic){
    ## train unsupervised forests
    forest.model <- Urerf(data, trees=500)
    ## get distance matrix
    dm <- (1 - forest.model$similarityMatrix)
  }else{
    dm <- as.matrix(dist(data))
  }
  ## extract knn
  nn <- matrix(0,n,k)
  for (i in 1:n){
    nn[i,] <- k.nearest.neighbors(i, dm, k = k)
  }
  ## for geodesic distance locally sort nn by euclidian distance
  if(geodesic){
    sortedNN <- sapply(1:n, function (rIdx)
    {
      res <- sapply(1:k, function (idNN)
      {
        l2dist(data[rIdx,], data[nn[rIdx,idNN],])
      })
      nn[rIdx,order(res)]
    })
    return(t(sortedNN))
  }else{
    return(nn)
  }
}

# computes the distance matrix and sorts each row
getDistances <- function(data, measure = "euclidean")
{
  if(measure == "geodesic"){
    ## train unsupervised forests
    forest.model <- Urerf(data, splitCrit = 'bicfast', trees=100)
    ## get distance matrix
    dm <- (1 - forest.model$similarityMatrix)
  }else{
    dm <- as.matrix(dist(data))
  }
  return(dm)
}
