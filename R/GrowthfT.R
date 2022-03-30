#' @title Foraminifer Growth Rate Function from Lombard et al. (2009)
#' @description Implements the function for foraminifer growth rate as a function of
#' temperature from Lombard et al. (2009) with parametrization from FAME 1.0
#' (Roche et al, 2018).
#'
#' @param foram Name of foram.
#' @param temperature_K Temperature in Kelvin
#' @param norm Optional normalizing factor
#' @param min.growth.thresh Sets a lower cutoff for growth as a proportion of
#'   the maximum growth rate for that taxon. For example in Roche et al (2018) a cutoff of
#'   0.1 was used, meaning all growth rates less than 0.1*max were set to zero.
#' @references
#'
#' Lombard, F., Labeyrie, L., Michel, E., Spero, H. J., and Lea, D. W.:
#' Modelling the temperature dependent growth rates of planktic foraminifera,
#' Marine Micropaleontology, 70, 1–7,
#' https://doi.org/10.1016/j.marmicro.2008.09.004, 2009.
#'
#' Roche, D. M., Waelbroeck, C., Metcalfe, B. and Caley, T.: FAME
#'   (v1.0): a simple module to simulate the effect of planktonic foraminifer
#'   species-specific habitat on their oxygen isotopic content, Geosci. Model
#'   Dev. Discuss., 2017, 1–22, doi:10.5194/gmd-2017-251, 2017.
#'
#' @return A numerical vector or matrix with the same dimensions as the object
#'   passed to temperature_K. Units are daily growth rate, unless norm == TRUE.
#' @export
#'
#' @examples
#' ForamGrowthfT(foram = 'ruber', temperature_K = (c(280, 290)), norm = 1)
ForamGrowthfT <- function(foram = c("sacculifer", "bulloides", "pachy_d",
                                           "siphonifera", "universa", "pachy_s",
                                           "dutertrei", "ruber"),
                              temperature_K, norm = FALSE, min.growth.thresh = 0) {

  foram <- match.arg(foram)

  pars <- l09_cnsts_dic[[foram]]
  names(pars) <- c("muT1", "TA", "TL", "TH", "TAL", "TAH")

  pars <- as.list(pars)

  gT <- function(muT1, TA, TL, TH, TAL, TAH) {
    muT1 * exp(TA / 293 - TA / temperature_K) /
      (1 + exp(TAL / temperature_K - TAL / TL) + exp(TAH / TH - TAH / temperature_K))
  }


  out <- do.call(gT, pars)

  max_func = l09_maxgrowth_dic[[foram]]

  out[out < min.growth.thresh * max_func] <- 0

  if (norm == TRUE)  {
    out <- out / max_func
  }

  out[temperature_K < -2 + 273.15] <- 0

  colnames(out) <- colnames(temperature_K)
  return(out)
}













