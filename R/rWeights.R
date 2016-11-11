#' rWeights
#'
#' @param n
#' @param weights
#'
#' @return
#' @export
#'
#' @examples
rWeights<-function(n, weights)
  ### Draw random number between 1 and length(weights) with the distribution
  ### from weights
  #Input
  #n Number of random numbers to be drawn
  #weights= vector of numbers; defines the distribution
  #uses runuran package Josef Leydold and Wolfgang H\{}"ormann
{
  dpv <- new(className("unuran.discr", "Runuran"), pv=weights, lb=1)
  gen <- unuran.new(dpv)
  return(ur(gen,n))
}
