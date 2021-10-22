#!/usr/bin/Rscript
source("get_marginalized_forest.R")
source("fast_gf_similarity.R")

## Input:
# forest
# data matrix p x n
# vector {0,1}^p, where a 0 indicates that a dimension is not regarded for the distances and a 1 indicates that the corresponding dimension is used
## Output: a similarity matrix nxn between 0 and 1

marginale_forest_distance <- function(X, forest, b, normalizeData=FALSE)
{
  if(length(b) != dim(X)[2]){
    error("Number of rows for X has to match the length of the vecotor b.")
  }
  forest.m <- get_marginalized_forest(X, forest, b, normalizeData)
  # get similarity matrix
  mat <- fast_gf_similarity(forest.m)
  return(1 - mat)
}
