require("plot3D")

# Add small dots on basal plane and on the depth plane
scatter3D_fancy <- function(x, y, z,..., colvar = z)
{
  panelfirst <- function(pmat) {
    XY <- trans3D(x, y, z = rep(min(z), length(z)), pmat = pmat)
    scatter2D(XY$x, XY$y, col = "black", pch = ".", 
              cex = 0.2, add = TRUE, colkey = FALSE)
    
    XY <- trans3D(x = rep(min(x), length(x)), y, z, pmat = pmat)
    scatter2D(XY$x, XY$y, col = "black", pch = ".", 
              cex = 0.2, add = TRUE, colkey = FALSE)
  }
  scatter3D(x, y, z, ..., colvar = colvar, panel.first=panelfirst) 
}


## Function to generate data of a 3d helix
## Data is computed by a value r (the "point" along the helix) and the x,y,z coordinates of the helix
## Optionally, Gaussian noise with zero mean can be added to the coordinates
##
## Arguments:
## n        number of samples to generate
## min      minimum value for r
## max      maximum value of r
## n_sigma  standard deviation of noise (0 means no noise)
##
## ground truth mutual information between r and xyz (without noise) is ln(max-min)
gen_helix_data <- function(n, min=0, max=2, n_sigma=0)
{
  r <- runif(n, min=min, max=max)
  x <- sapply(r, function(r_1) {
    res = 5*pi + 3*pi*r_1
    res = res*cos(res)/(8*pi)
    return(res)
  })
  y <- sapply(r, function(r_1) {
    res = 5*pi + 3*pi*r_1
    res = res*sin(res)/(8*pi)
    return(res)
  })
  ## set the denominator to 2*pi for more space in Z direction
  z <- sapply(r, function(r_1) {
    res = 5*pi + 3*pi*r_1/(8*pi)
    return(res)
  })
  
  if (n_sigma > 0)
  {
    x <- x + rnorm(n, sd=n_sigma)
    y <- y + rnorm(n, sd=n_sigma)
    z <- z + rnorm(n, sd=n_sigma)
  }
  return(list(r=r, x=x, y=y, z=z))
}

gen_helix_data_unif <- function(n, min=0, max=1, alpha=0.01)
{
  r <- runif(n, min=min, max=max)
  x <- sapply(r, function(r_1) {
    res = 5*pi + 3*pi*r_1
    res = res*cos(res)/(8*pi)
    return(res)
  })
  y <- sapply(r, function(r_1) {
    res = 5*pi + 3*pi*r_1
    res = res*sin(res)/(8*pi)
    return(res)
  })
  z <- sapply(r, function(r_1) {
    res = 5*pi + 3*pi*r_1/(8*pi)
    return(res)
  })
  
  if (alpha > 0)
  {
    x <- x + alpha*runif(n, min=0,max=1)
    y <- y + alpha*runif(n, min=0,max=1)
    z <- z + alpha*runif(n, min=0,max=1)
  }
  return(list(r=r, x=x, y=y, z=z))
}

plot_helix_data <- function(r, x, y, z)
{
  scatter3D_fancy(x,y,z, colvar=r, pch=16, cex=0.5,
            ticktype="detailed", bty="g", theta=30, phi=30
  )
}

