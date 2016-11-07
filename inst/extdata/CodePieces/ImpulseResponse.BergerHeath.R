
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

z <- seq(0, 15, length.out = 100)
plot(z, ImpulseResponse.BergerHeath(z, 10), type = "l")
abline(v = 10)

integrate(ImpulseResponse.BergerHeath, lower = 0, upper = 20, d = 10)
