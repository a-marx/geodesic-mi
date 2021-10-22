#!/usr/bin/Rscript
## Input:
# forest
## Output: a distance matrix nxn
## The distance is smaller, the more frequently two nodes i,j appear in the same leaf
fast_gf_similarity <- function(Forest)
{
  tS <- Forest[[1]]$TrainSize
  numTrees <- length(Forest)
  tree_memberships <- matrix(0, nrow=tS, ncol=numTrees)
  simMatrix <- matrix(0, nrow = tS, ncol = tS)
  for (i in 1:numTrees)
  {
    childNodes <- which(Forest[[i]]$Children[, 1] == 0)
    for(j in childNodes)
    {
      for(k in 1:length(Forest[[i]]$ALeaf[[j]]))
      {
        tree_memberships[Forest[[i]]$ALeaf[[j]][k],i] <- j
      }
    }
  }
  for (i in 1:tS) {
    simMatrix[i,i] <- numTrees
  }
  for (i in 2:tS) {
    for (j in 1:(i-1)) {
      simij <- sum(tree_memberships[i,] == tree_memberships[j,])
      simMatrix[i,j] <- simij
      simMatrix[j,i] <- simij
    }
  }
  simMatrix <- simMatrix / numTrees
  if(any(diag(simMatrix) != 1))
  {
    print("diag not zero")
    diag(simMatrix) <- 1
  }
  return(simMatrix)
}
