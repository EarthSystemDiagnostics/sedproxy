#' Title
#'
#' @param z
#' @param focal.depth
#' @param layer.width
#' @param sed.acc.rate
#' @param mix.depth
#'
#' @return
#' @export
#'
#' @examples
#' z <- 0:10000
#' w <- BioturbationWeights(z, focal.depth = 2000, layer.width = 1, sed.acc.rate = 5 / 1000, mix.depth = 10)
#' plot(z, w, "l")
BioturbationWeights <- function(z = NULL, focal.depth, layer.width, sed.acc.rate, mix.depth){

  lwy <- (layer.width / sed.acc.rate)
  mdy <- (mix.depth / sed.acc.rate)

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



