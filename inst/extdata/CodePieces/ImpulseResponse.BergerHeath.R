
ImpulseResponse.BergerHeath <- function(z, d, z0 = 0) {
  ### Deep solution for Berger and Heath; e.g. page 65 in Officer and Lynch
  #z =  depth
  #d =  mixing depth
  #z0 = depth of the pulse without bioturbation
  x <- z0 + d - z
  epsilon <- x / d
  result <- 1 / d * exp(-epsilon)
  result[z > (z0 + d)] <- 0
  return(result)
}

z <- seq(-100, 100, length.out = 1000)
plot(z, ImpulseResponse.BergerHeath(z, 50, 0), type = "l")
abline(v = 50)
integrate(ImpulseResponse.BergerHeath, lower = -100, upper = 100, d = 10)
