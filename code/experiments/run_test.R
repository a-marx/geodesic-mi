#!/usr/bin/Rscript
source("experiments/data_generators.R")
source("information_measures.R")
source("gkNN.R")
require("rmi")
require("rerf")

## Note
# These settings apply for running G-KSG. If gKNN or LNN should be executed set geodesic in line 27 to F, uncomment the corresponding application in lines 104-120. Also, you need to change the name of the output file (oFile).

######################################
# Specify Parameters for Experiments #
######################################
# Which data?
generator <- "Sphere" ## Gaussian, Uniform, Sphere, Spiral
p.data <- list(alpha=0.01)  ## use cov=0.9 for Gaussian and alpha=0.01 for rest
# How much data?
n <- 500
dims <- 0:10 # 0:10 # c(0,50,100,150,200,250,300)
# How many repetitions?
loops <- 20 # 20 3
d <- rep(0,length(dims))
results <- data.frame(D=dims, KSG3M=d, KSG3MSE=d, KSG5M=d, KSG5MSE=d, KSG10M=d, KSG10MSE=d)
# What algorithm setting?
geodesic <- T
nTrees = 300

############################################################
# Generate automatic output file name based on input above #
############################################################
# default settings for forest
normD <- F
pathSimilarity <- F
splitCrit <- "bicfast"
minParent <- round((2*n)^.5) # if changed, change in line 41
unifNoise = F
paramVec <- c()
if(geodesic){
  paramVec <- c(paramVec, "G")
  if(normD){
    paramVec <- c(paramVec, "N")
  }
  if(hybridSimilarity){
    paramVec <- c(paramVec, "Hs")
  }else if(pathSimilarity){
    paramVec <- c(paramVec, "Ps")
  }
  if(splitCrit == "twomeans"){
    paramVec <- c(paramVec, "SCtm")
  }else if(splitCrit == "bicfast"){
    paramVec <- c(paramVec, "SCfb")
  }else{
    paramVec <- c(paramVec, "SCr")
  }
  paramVec <- c(paramVec, "MP2n")
  paramVec <- c(paramVec, paste("T", nTrees, sep=""))
}else{
  paramVec <- c(paramVec, "E")
}
oFile <- paste(c("experiments/results/", generator ,"GNoiseN", n, "D", max(dims), "L", loops, "_", paramVec, ".tab"), collapse="")

# Print file name to confirm settings
print(oFile)
set.seed(1)
for(di in 1:length(dims)){
  dimension = dims[di]
  KSG3 = rep(0,loops)
  KSG5 = rep(0,loops)
  KSG10 = rep(0,loops)
  for(l in 1:loops){
    # generate base data
    gen.data <- generate_data(n,p.data,generator)
    X <- gen.data$X
    Y <- gen.data$Y
    trueI <- gen.data$trueI
    # add noise dimensions
    if(dimension > 0){
      for(iter in 1:dimension){
        if(unifNoise){
          X = cbind(X, runif(n, min=0, max=1))
          Y = cbind(Y, runif(n, min=0, max=1))
        }
      }
    }
    data = cbind(X,Y)
    # precompute tree and marginalizations
    fm <- NULL
    Xdist <- NULL
    Ydist <- NULL
    if(geodesic){
      fm <- Urerf(data, trees=nTrees, normalizeData=normD, splitCrit=splitCrit, min.parent=minParent)
      Xdist <- marginale_forest_distance(X=data, forest=fm$forest, b=c( rep(1, dim(X)[2]), rep(0,dim(Y)[2]) ), normalizeData=normD)
      Ydist <- marginale_forest_distance(X=data, forest=fm$forest, b=c( rep(0, dim(X)[2]), rep(1,dim(Y)[2]) ), normalizeData=normD)
    }
    # compute MI estimates
    # G-KSG
    KSG3[l] = mi_ksg_estimate_forest(X=X,Y=Y, fm=fm, Xdist=Xdist, Ydist=Ydist, k=3, geodesic=geodesic, normalizeData=normD)
    # gKNN
    #KSG3[l] <- mi_gKNN(X=X, Y=Y, k=20)
    # LNN
    #KSG3[l] <- lnn_mi(data, c(ncol(X),ncol(Y)), k=3)
    # G-KSG
    KSG5[l] = mi_ksg_estimate_forest(X=X,Y=Y, fm=fm, Xdist=Xdist, Ydist=Ydist, k=5, geodesic=geodesic, normalizeData=normD)
     # gKNN
    #KSG5[l] <- mi_gKNN(X=X, Y=Y, k=30)
    # LNN
    #KSG5[l] <- lnn_mi(data, c(ncol(X),ncol(Y)), k=5)
    # G-KSG
    KSG10[l] = mi_ksg_estimate_forest(X=X,Y=Y, fm=fm, Xdist=Xdist, Ydist=Ydist, k=10, geodesic=geodesic, normalizeData=normD)
    # gKNN
    #KSG10[l] <- mi_gKNN(X=X, Y=Y, k=50)
    # LNN
    #KSG10[l] <- lnn_mi(data, c(ncol(X),ncol(Y)), k=10, tr=50)
  }
  results$KSG3M[di] = mean(KSG3)
  results$KSG3MSE[di] = mean((KSG3 - trueI)^2)
  results$KSG5M[di] = mean(KSG5)
  results$KSG5MSE[di] = mean((KSG5 - trueI)^2)
  results$KSG10M[di] = mean(KSG10)
  results$KSG10MSE[di] = mean((KSG10 - trueI)^2)
  
  write.table(results, file=oFile, sep="\t", row.names=F, col.names=T, quote=F)
}
