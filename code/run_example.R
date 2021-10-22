#!/usr/bin/Rscript

require("MASS")
require("rerf")
source("G_KSG.R")
source("experiments/data_generators.R")
# Generate Data
set.seed(123)
n = 500
# generates spiral (or helix) data
gen.data <- generate_data(n,list(alpha=0.01),"Spiral")
X <- gen.data$X
Y <- gen.data$Y
trueI <- gen.data$trueI
# Run G_KSG
G_KSG(X,Y,k=3)
# [1] 2.828043
plot_helix_data(Y,X[,1],X[,2],X[,3])
