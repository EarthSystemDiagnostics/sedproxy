#' ImpulseResponse
#'
#' @param z depth relative to focus horizon
#' @param d mixing depth
#' @param z0 depth of the pulse without bioturbation
#'
#' @description Depth solution for Berger and Heath; e.g. page 65 in Officer and
#'   Lynch (1983).
#'
#'   For a given mixing depth \code{d}, it gives the probability
#'   that a particle has been moved from a depth horizon \code{z} units away
#'   from the focus horizon.
#'
#' @return Given a vector of depths, \code{z}, relative to a given horizon, the
#'   function returns a vector of probabilities that a given particle at that
#'   horizon has been moved there by bioturbation from depth \code{z}
#' @export
#'
#' @examples
#' z <- seq(-100, 100, length.out = 1000)
#' d <- 10
#'
#' plot(z, ImpulseResponse(z, d), type = "l")
#' abline(v = d)
#'
ImpulseResponse <- function(z, d, z0 = 0) {
  x <- z0 + d - z
  epsilon <- x / d
  result <- 1 / d * exp(-epsilon)
  result[z > (z0 + d)] <- 0
  return(result)
}
