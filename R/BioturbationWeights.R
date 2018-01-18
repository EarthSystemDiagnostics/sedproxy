#' Bioturbation weights
#' @description For a given focal depth (or time), this function returns the probability
#' that material collected from that depth was orignially deposited at depth(s)
#' z. In other words, that the material would have been found at depth z if there
#' had been no bioturbation. It is the convolution of the depth solution from
#' Berger and Heath (1968) with a uniform distribution to account for the width
#' of the sediment layer from which samples
#' were picked/extracted. It is a probability density function.
#' @inheritParams ClimToProxyClim
#' @param z A vector of times or depths at which to evaluate the bioturbation weights
#' @param focal.depth The depth (or time) for which source dates are wanted
#' @param scale whether to scale depths by sediment accumulation rate to give
#' positions in terms of time
#' @return a vector of weights
#' @export
#' @references Berger, W. H., & Heath, G. R. (1968).
#' Vertical mixing in pelagic sediments.
#' Journal of Marine Research, 26(2), 134–143.
#' @examples
#' z <- 0:10000
#' w <- BioturbationWeights(z, focal.depth = 2000, layer.width = 1, sed.acc.rate = 5 / 1000, bio.depth = 10)
#' plot(z, w, "l")
BioturbationWeights <- function(z, focal.depth, layer.width=1, sed.acc.rate, bio.depth, scale = c("time", "depth")){

  sed.acc.rate <- sed.acc.rate / 1000

  scale <- match.arg(scale)

  if (scale == "time"){
    lwy <- (layer.width / sed.acc.rate)
    mdy <- (bio.depth / sed.acc.rate)
  }else{
    lwy <- (layer.width)
    mdy <- (bio.depth)
  }

  fd <- focal.depth

  C <- lwy/2
  lam <- 1/mdy

  z <- z - fd + mdy

  if (mdy <= 1){
    fz <- dunif(z, -C, C)
  }else if (lwy == 0){
    fz <- dexp(z, 1/mdy)
  }else{
    fz <- (z < -C) * 0 +
      (z >= -C & z <= C) * (lam*(1/lam-exp(-lam*C-lam*z)/lam))/(2*C)  +
      (z > C) * (lam*(exp(lam*C-lam*z)/lam-exp(-lam*C-lam*z)/lam))/(2*C)
  }

  return(fz)
}

z <- seq(0, 500, length.out = 1000)
w <- BioturbationWeights(z, focal.depth = 50, layer.width = 1, sed.acc.rate = 50, bio.depth = 10, scale = "depth")
plot(z, w, "l")
