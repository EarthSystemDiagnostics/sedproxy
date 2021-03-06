% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FAME.R
\name{growth_rate_l09_R}
\alias{growth_rate_l09_R}
\title{Foraminifer Growth Rate Function from FAME 1.0}
\usage{
growth_rate_l09_R(
  nm_foram = c("sacculifer", "bulloides", "pachy_d", "siphonifera", "universa",
    "pachy_s", "dutertrei", "ruber"),
  input_T,
  norm = FALSE,
  min.growth.thresh = 0
)
}
\arguments{
\item{nm_foram}{Name of foram.}

\item{input_T}{Input temperature in Kelvin}

\item{norm}{Optional normalizing factor}

\item{min.growth.thresh}{Sets a lower cutoff for growth as a proportion of
the maximum growth rate for that taxon. In Roche et al (2018) a cutoff of
0.1 was used, meaning all growth rates less than 0.1*max were set to zero.}
}
\value{
Growth rate/day
}
\description{
The function for foraminifer growth rate as a function of
temperature \code{growth_rate_l09} and tables of parameters
\code{l09_maxgrowth_dic} and \code{l09_cnsts_dic}
from FAME 1.0 (Roche et al, 2018). Function and data objects
translated into R from the original Python code supplied with the below cited
discussion paper.

The code is made available under the GNU General Public License
https://www.gnu.org/licenses/gpl.html
}
\examples{
growth_rate_l09_R(nm_foram = 'ruber', input_T = (c(280, 290)), norm = 1)
}
\references{
Roche, D. M., Waelbroeck, C., Metcalfe, B. and Caley, T.: FAME
  (v1.0): a simple module to simulate the effect of planktonic foraminifer
  species-specific habitat on their oxygen isotopic content, Geosci. Model
  Dev. Discuss., 2017, 1–22, doi:10.5194/gmd-2017-251, 2017.
}
