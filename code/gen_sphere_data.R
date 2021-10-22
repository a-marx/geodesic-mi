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


## Function to generate data of a 3d sphere
## Data is computed by two values u1 and u2, out of which longitude and latitute is derived, from which we can obtain xyz coords.
## see https://math.stackexchange.com/questions/1585975/how-to-generate-random-points-on-a-sphere for more info.
## Optionally, Gaussian noise with zero mean can be added to the coordinates.
##
## Arguments:
## n        number of samples to generate
## r        radius (this parameter should matter if you get few samples and smaller balls, then it gets
##                  more difficult to distinguish manifold from ambient space, use plots to find a good r)
## n_sigma  standard deviation of noise (0 means no noise)
##
gen_sphere_data <- function(n, r=1, n_sigma=0)
{
  u1 <- runif(n, min=0, max=1)
  u2 <- runif(n, min=0, max=1)
  
  ## get long and lat
  lat <- acos(2*u1 - 1) - pi/2
  long <- 2*pi*u2
  
  ## generate x,y,z coordinates
  x <- cos(lat)*cos(long)*r
  y <- cos(lat)*sin(long)*r
  z <- sin(lat)*r
  
  if (n_sigma > 0)
  {
    x <- x + rnorm(n, sd=n_sigma)
    y <- y + rnorm(n, sd=n_sigma)
    z <- z + rnorm(n, sd=n_sigma)
  }
  return(list(u1=u1, u2=u2, x=x, y=y, z=z))
}

gen_sphere_data_unif <- function(n, r=1, alpha=0.1)
{
  u1 <- runif(n, min=0, max=1)
  u2 <- runif(n, min=0, max=1)
  
  ## get long and lat
  lat <- acos(2*u1 - 1) - pi/2
  long <- 2*pi*u2
  
  ## generate x,y,z coordinates
  x <- cos(lat)*cos(long)*r
  y <- cos(lat)*sin(long)*r
  z <- sin(lat)*r
  
  if (alpha > 0)
  {
    x <- x + alpha * runif(n, min=-0.5,max=0.5)
    y <- y + alpha * runif(n, min=-0.5,max=0.5)
    z <- z + alpha * runif(n, min=-0.5,max=0.5)
  }
  return(list(u1=u1, u2=u2, x=x, y=y, z=z))
}

plot_sphere_data <- function(u1, u2, x, y, z)
{
  scatter3D_fancy(x,y,z, colvar=u1, pch=16, cex=0.5,
                  ticktype="detailed", bty="g", theta=30, phi=30
  )
}
