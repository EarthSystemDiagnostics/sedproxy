#' @title Foraminifer Growth Rate Function from FAME 1.0
#' @description The function for foraminifer growth rate as a function of
#' temperature \code{growth_rate_l09} and tables of parameters
#' \code{l09_maxgrowth_dic} and \code{l09_cnsts_dic}
#' from FAME 1.0 (Roche et al, 2018). Function and data objects
#' translated into R from the original Python code supplied with the below cited
#' discussion paper.
#'
#' The code is made available under the GNU General Public License
#' https://www.gnu.org/licenses/gpl.html
#'
#' @param nm_foram Name of foram.
#' @param input_T Input temperature in Kelvin
#' @param norm Optional normalizing factor
#' @param min.growth.thresh Sets a lower cutoff for growth as a proportion of
#'   the maximum growth rate for that taxon. In Roche et al (2018) a cutoff of
#'   0.1 was used, meaning all growth rates less than 0.1*max were set to zero.
#' @references Roche, D. M., Waelbroeck, C., Metcalfe, B. and Caley, T.: FAME
#'   (v1.0): a simple module to simulate the effect of planktonic foraminifer
#'   species-specific habitat on their oxygen isotopic content, Geosci. Model
#'   Dev. Discuss., 2017, 1–22, doi:10.5194/gmd-2017-251, 2017.
#'
#' @return Growth rate/day
#' @export
#'
#' @examples
#' growth_rate_l09_R(nm_foram = 'ruber', input_T = (c(280, 290)), norm = 1)
growth_rate_l09_R <- function(nm_foram = c(
  "sacculifer", "bulloides", "pachy_d", "siphonifera", "universa", "pachy_s",
  "dutertrei", "ruber"),
  input_T, norm = FALSE, min.growth.thresh = 0) {
  # Growth rate in day^-1 ...  Inputs: foram name as string
  # input_T temperature in °K

  # Optional input norm is a normalization factor

  # mut, TA, TL, TH, TAL, TAH = l09_cnsts_dic[nm_foram]
  pars <- l09_cnsts_dic[[nm_foram]]
  names(pars) <- c("mut", "TA", "TL", "TH", "TAL", "TAH")
  pars <- as.list(pars)

  out <- with(pars, {
    {
      if (TAL > 0) {
        (mut * exp(TA/293 - TA/input_T)/(1 + exp(TAL/input_T -
                                                   TAL/TL) + exp(TAH/TH - TAH/input_T)))
      } else {
        (mut * exp(TA/293 - TA/input_T)/(2 + exp(TAH/TH -
                                                   TAH/input_T)))
      }
    }
  })

  max_func = l09_maxgrowth_dic[[nm_foram]]

  out[out < min.growth.thresh * max_func] <- 0

  if (norm == TRUE)  {
    out <- out / max_func
  }

  out[input_T < -2 + 273.15] <- 0


  colnames(out) <- colnames(input_T)
  return(out)
}













