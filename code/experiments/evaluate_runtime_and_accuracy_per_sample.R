#!/usr/bin/Rscript
source("experiments/data_generators.R")
source("information_measures.R")
source("gkNN.R")
require("rmi")
require("rerf")

## Note
# These settings apply for running G-KSG. If gKNN or LNN should be executed set geodesic in line 27 to F, comment line 90 out and uncomment line 92 resp. 92. Also, you need to change the name of the output file (oFile).

######################################
# Specify Parameters for Experiments #
######################################
# Which data?
generator <- "Uniform" ## Gaussian, Uniform, Sphere, Spiral
p.data <- list(alpha=0.01)  ## use cov=0.9 for Gaussian and alpha=0.01 for rest
# How much data?
n_steps <- c(1:20)  #c(1:5) c(1:20)
n_factor <- 100  # 2000 100
unifNoise = F
# How many repetitions?
loops <- 10 # 5
d <- rep(0,length(n_steps))
results <- data.frame(N=n_steps*n_factor, Mean=d, MSE=d, Time=d)
# What algorithm setting?
geodesic <- T
nTrees <- 300

############################################################
# Generate automatic output file name based on input above #
############################################################
# default settings for forest
normD <- F
pathSimilarity <- F
splitCrit <- "bicfast"
paramVec <- c()
if(geodesic){
  paramVec <- c(paramVec, "G")
  if(normD){
    paramVec <- c(paramVec, "N")
  }
  if(pathSimilarity){
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
oFile <- paste(c("experiments/results/Samples", generator ,"GNoiseN", n_factor, "L", loops, "_", paramVec, ".tab"), collapse="")

# Print file name to confirm settings
print(oFile)
set.seed(1)
for(ni in 1:length(n_steps)){
  n = ni * n_factor
  MeanV = rep(0,loops)
  TimeV = rep(0,loops)
  for(l in 1:loops){
    # generate base data
    gen.data <- generate_data(n,p.data,generator)
    X <- gen.data$X
    Y <- gen.data$Y
    trueI <- gen.data$trueI
    data = cbind(X,Y)
    # precompute tree and marginalizations
    fm <- NULL
    Xdist <- NULL
    Ydist <- NULL
    start_time <- Sys.time()
    if(geodesic){
      minParent <- round((2*n)^.5) # if changed, change in line 41
      fm <- Urerf(data, trees=nTrees, normalizeData=normD, splitCrit=splitCrit, min.parent=minParent)
      Xdist <- marginale_forest_distance(X=data, forest=fm$forest, b=c( rep(1, dim(X)[2]), rep(0,dim(Y)[2]) ), normalizeData=normD)
      Ydist <- marginale_forest_distance(X=data, forest=fm$forest, b=c( rep(0, dim(X)[2]), rep(1,dim(Y)[2]) ), normalizeData=normD)
    }
    # compute MI estimates
    # G-KSG
    MeanV[l] <- mi_ksg_estimate_forest(X=X,Y=Y, fm=fm, Xdist=Xdist, Ydist=Ydist, k=3, geodesic=geodesic, normalizeData=normD)
    # gKNN
    #MeanV[l] <- mi_gKNN(X=X, Y=Y, k=20)
    # LNN
    #MeanV[l] <- lnn_mi(data, c(ncol(X),ncol(Y)), k=3)
    total_time <- as.numeric(Sys.time() - start_time, units = "secs")
    TimeV[l] <- total_time
  }
  results$Mean[ni] = mean(MeanV)
  results$MSE[ni] = mean((MeanV - trueI)^2)
  results$Time[ni] = mean(TimeV)

  write.table(results, file=oFile, sep="\t", row.names=F, col.names=T, quote=F)
}
