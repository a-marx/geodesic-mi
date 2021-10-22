#!/usr/bin/Rscript
source("gen_helix_data.R")
source("gen_sphere_data.R")
require("MASS")

# A data gernerator has arguments n and p for parameter.
# @n number of samples to be generated
# @p a list of parameters (e.g. list(alhpa=0.1) for uniform, or list(tCov=0.9) for Gaussian)
# @return list(X of type matrix, Y of type matrix, trueI of type double)
get_gaussian <- function(n,p){
  if(is.null(p$cov)){
    stop("Generator -> get_gaussian: Covariance ('p$cov') not specified!")
  }
  trueI <- -0.5 * log(1-p$cov^2)
  covM <- matrix(c(1,p$cov,p$cov,1), nrow=2, ncol=2)
  data = as.matrix(mvrnorm(n=n, Sigma=covM, mu=c(0,0)))
  X = as.matrix(data[,1])
  Y = as.matrix(data[,2])
  return(list(X=X,Y=Y,trueI=trueI))
}
get_uniform <- function(n,p){
  if(is.null(p$alpha)){
    stop("Generator -> get_uniform: Noise level ('p$alpha') not specified!")
  }
  trueI = -log(p$alpha) + p$alpha/2.0
  X <- as.matrix(runif(n, min=0, max=1))
  N <- as.matrix(runif(n, min=0, max=1))
  Y <- as.matrix(X + (p$alpha)*N)
  return(list(X=X,Y=Y,trueI=trueI))
}
get_sphere <- function(n,p){
  if(is.null(p$alpha)){
    stop("Generator -> get_sphere: Noise level ('p$alpha') not specified!")
  }
  # Note: the true information is >= trueI
  trueI = -log(p$alpha) + p$alpha/2.0
  dat = gen_sphere_data_unif(n=n,alpha=p$alpha)
  X = cbind(dat$x,dat$y,dat$z)
  Y = cbind(dat$u1, dat$u2)
  return(list(X=X,Y=Y,trueI=trueI))
}
get_spiral <- function(n,p){
  if(is.null(p$alpha)){
    stop("Generator -> get_spiral: Noise level ('p$alpha') not specified!")
  }
  trueI = -log(p$alpha) + p$alpha/2.0
  dat = gen_helix_data_unif(n=n,min=0,max=1,alpha=p$alpha)
  X = cbind(dat$x,dat$y,dat$z)
  Y = as.matrix(dat$r)
  return(list(X=X,Y=Y,trueI=trueI))
}
generate_data <- function(n,p,generator){
  if(generator == "Gaussian"){
    return(get_gaussian(n,p))
  }else if(generator == "Uniform"){
    return(get_uniform(n,p))
  }else if(generator == "Sphere"){
    return(get_sphere(n,p))
  }else if(generator == "Spiral"){
    return(get_spiral(n,p))
  }else{
    stop(paste(c("Generator '", generator, "' not found!"), collapse=""))
  }
}
