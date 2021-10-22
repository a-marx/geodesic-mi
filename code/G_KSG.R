#!/usr/bin/Rscript

source("information_measures.R")

G_KSG <- function(X,Y,k,nTrees=300,splitCrit="bicfast",normalizeData=FALSE,minParent=round((2*n)^.5)){
  data = cbind(X,Y)
  fm <- Urerf(data, trees=nTrees, normalizeData=normalizeData, splitCrit=splitCrit, min.parent=minParent)
  Xdist <- marginale_forest_distance(X=data, forest=fm$forest, b=c( rep(1, dim(X)[2]), rep(0,dim(Y)[2]) ), normalizeData=normalizeData)
  Ydist <- marginale_forest_distance(X=data, forest=fm$forest, b=c( rep(0, dim(X)[2]), rep(1,dim(Y)[2]) ), normalizeData=normalizeData)
  estimate <- mi_ksg_estimate_forest(X=X,Y=Y, fm=fm, Xdist=Xdist, Ydist=Ydist, k=k, geodesic=T, normalizeData=normalizeData)
  return(estimate)
}
