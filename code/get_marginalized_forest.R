#!/usr/bin/Rscript
normalizeTheData <- function(X, normData)
{
  if(normData)
  {
    X <- sweep(X, 2, apply(X, 2, min), "-")
    return(sweep(X, 2, apply(X, 2, max), "/"))
  }
  else
  {
    return(X)
  }
}

## Input:
# forest
# data matrix p x n
# vector {0,1}^p, where a 0 indicates that a dimension is not regarded for the distances and a 1 indicates that the corresponding dimension is used
## Output: a distance matrix nxn

get_marginalized_forest <- function(X, forest, b, normalizeData=TRUE)
{
  if(length(b) != dim(X)[2]){
    error("Number of rows for X has to match the length of the vecotor b.")
  }
  # normalize data if necessary
  X = normalizeTheData(X, normalizeData)
  # init forest for marginalization
  forest.m <- vector("list", length(forest))
  for(tree.i in 1:length(forest))
  {
    num.leaves <- dim(forest[[tree.i]]$Children)[1]
    # init emptry tree
    Assigned2Leaf <- vector("list", num.leaves)
    stack <- c(0,1)
    Assigned2Leaf[[1]] <- 1:(dim(X)[1])
    while(stack[length(stack)] != 0)
    {
      # get element from stack
      current <- stack[length(stack)]
      stack <- c(stack[1:length(stack)-1])
      # check if the current node has children
      if(forest[[tree.i]]$Children[current,][1] != 0){
        # prepare vector a' <- <a,b>
        matA <- forest[[tree.i]]$matA[[current]]
        a.dash <- rep(0,length(b))
        for(i.a in 1:(length(matA) / 2)){
          pos.col <- matA[i.a]
          val.A <- matA[i.a + (length(matA) / 2)]
          a.dash[pos.col] <- val.A * b[pos.col]
        }
        # project X based on a.dash
        x.proj <- X[Assigned2Leaf[[current]], ] %*% a.dash
        x.left <- c(x.proj) < forest[[tree.i]]$CutPoint[current]
        # assigen indices to left child
        l.child <- forest[[tree.i]]$Children[current,][1]
        stack <- c(stack, l.child)
        Assigned2Leaf[[l.child]] <- Assigned2Leaf[[current]][x.left]
        # assigen indices to right child
        r.child <- forest[[tree.i]]$Children[current,][2]
        stack <- c(stack, r.child)
        Assigned2Leaf[[r.child]] <- Assigned2Leaf[[current]][!x.left]
      }
    }
    forest.m[[tree.i]] <- list(ALeaf = Assigned2Leaf, Children = forest[[tree.i]]$Children, TrainSize = forest[[tree.i]]$TrainSize)
  }
  return(forest.m)
}
