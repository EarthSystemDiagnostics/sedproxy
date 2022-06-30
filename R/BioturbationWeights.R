#' Bioturbation weights
#' @description For a given focal depth (or time), this function returns the probability
#' that material collected from that depth was originally deposited at depth(s)
#' z. In other words, that the material would have been found at depth z if there
#' had been no bioturbation. It is the convolution of the depth solution from
#' Berger and Heath (1968) with a uniform distribution to account for the width
#' of the sediment layer from which samples
#' were picked/extracted. It is a probability density function.
#' @inheritParams ClimToProxyClim
#' @param z A vector of times or depths at which to evaluate the bioturbation weights
#' @param focal.z The depth (or time) for which source dates are wanted
#' @param scale Whether to scale depths by sediment accumulation rate to give
#' positions in terms of time. Defaults to time.
#' @return a numerical vector of weights.
#' @export
#' @references Berger, W. H., & Heath, G. R. (1968).
#' Vertical mixing in pelagic sediments.
#' Journal of Marine Research, 26(2), 134–143.
#' @examples
#' z <- 0:10000
#' w <- BioturbationWeights(z, focal.z = 4000, layer.width = 1, sed.acc.rate = 5, bio.depth = 10)
#' plot(z, w, "l")
BioturbationWeights <- function(z, focal.z, layer.width=1, sed.acc.rate, bio.depth, scale = c("time", "depth")){

  sed.acc.rate <- sed.acc.rate / 1000

  scale <- match.arg(scale)

  if (scale == "depth"){
    z <- z / sed.acc.rate
    focal.z <- focal.z / sed.acc.rate
  }


  lwy <- ceiling(layer.width / sed.acc.rate)
  #lwy[lwy < 1] <- 1 #lwy can be 0
  mdy <- ceiling(bio.depth / sed.acc.rate)

  if (lwy == 0 & mdy == 0) lwy <- 1

  C <- lwy/2
  lam <- 1/mdy

  z <- z - focal.z + mdy

  if (mdy <= 1){
    fz <- stats::dunif(z, -C, C)
  }else if (lwy == 0){
    fz <- stats::dexp(z, 1/mdy)
  }else{
    fz <- (z < -C) * 0 +
      (z >= -C & z <= C) * (lam*(1/lam-exp(-lam*C-lam*z)/lam))/(2*C)  +
      (z > C) * (lam*(exp(lam*C-lam*z)/lam-exp(-lam*C-lam*z)/lam))/(2*C)
  }
  if (sum(fz) == 0){fz}else{
    fz <- fz / sum(fz, na.rm = T)
  }

  return(fz)
}
