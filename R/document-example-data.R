#' Example data
#'
#' @name example.data
#' @title Data for running the examples in package \code{sedproxy}
#' @description A set of data objects for running the examples in \code{sedproxy}
#'
#' \code{\link{N41.proxy}}
#'
#' \code{\link{N41.proxy.details}}
#'
#' \code{\link{N41.G.ruber.seasonality}}
#'
#' \code{\link{N41.t21k.climate}}
NULL

#' @title Mg/Ca proxy based temperature reconstruction for core MD97-2141
#' @description Mg/Ca proxy based temperature reconstruction for core MD97-2141, originally published in
#' Rosenthal et al. (2003), extracted from Shakun et al. (2012).
#' @format A data frame with 216 rows and 4 variables:
#' \tabular{ll}{
#'   \cr \code{Published.age} \tab Age in yr BP
#'   \cr \code{Published.temperature} \tab Reconstructed temperature in degrees C
#'   \cr \code{Sed.acc.rate.cm.ka} \tab Sediment accumulation rate in cm per kyr
#'   }
#' @details Published age and published temperature were extracted from
#' Shakun, J. D., Clark, P. U., He, F., Marcott, S. A.,
#' Mix, A. C., Liu, Z., … Bard, E. (2012). Global warming preceded by increasing
#' carbon dioxide concentrations during the last deglaciation. Nature,
#' 484(7392), 49–54. \doi{10.1038/nature10915}.
#'
#' Sediment accumulation rates were estimated by fitting splines to the published
#' age and depth estimates.
#'
#' @source Original reference: Rosenthal, Y., Oppo, D. W., & Linsley, B. K.
#' (2003). The amplitude and phasing of climate change during the last
#' deglaciation in the Sulu Sea, western equatorial Pacific. Geophysical
#' Research Letters, 30(8), 1428. \doi{10.1029/2002GL016612}.
#'
"N41.proxy"

#' @title Metadata for datset \code{N41.proxy}
#' @description Metadata for core MD97-2141 taken from Shakun et al. (2012)
#' @format A data frame with 1 rows and 17 variables:
#' \tabular{ll}{
#'    \cr \code{Number} \tab Proxy ID number from Shakun et al. (2012)
#'    \cr \code{ID.no} \tab Proxy ID number from Shakun et al. (2012) with prefix "N"
#'    \cr \code{Core} \tab ID code of sediment core
#'    \cr \code{Location} \tab Location of core
#'    \cr \code{Proxy} \tab Proxy type
#'    \cr \code{Lat} \tab Latutide of core
#'    \cr \code{Lon} \tab Longitude of core
#'    \cr \code{Elevation} \tab Depth of core-top below sea level in metres
#'    \cr \code{Reference} \tab Original reference for the proxy record
#'    \cr \code{Resolution} \tab Average time resolution in years
#'    \cr \code{Calibration.ref} \tab Reference for Mg/Ca calibration
#'    \cr \code{Calibration} \tab Mg/Ca calibration formula
#'    \cr \code{Foram.sp} \tab Foram species analyses for proxy
#'    \cr \code{Ref.14C} \tab Reference for carbon dating
#'    \cr \code{Notes} \tab
#'    \cr \code{Geo.cluster} \tab Coarse geographic location
#'    \cr \code{Archive.type} \tab Type of proxy archive
#'}
#' @source Shakun, J. D., Clark, P. U., He, F., Marcott, S. A.,
#' Mix, A. C., Liu, Z., … Bard, E. (2012). Global warming preceded by increasing
#' carbon dioxide concentrations during the last deglaciation. Nature,
#' 484(7392), 49–54. \doi{10.1038/nature10915}
"N41.proxy.details"

#' @title Seasonality of Globigerinoides ruber at core MD97-2141
#' @description Seasonality of Globigerinoides ruber at core MD97-2141 predicted
#' by the PLAFOM model (Fraile et al. 2008).
#' @format A vector of 12 values
#' @source Fraile, I., Schulz, M., Mulitza, S., & Kucera, M. (2008).
#' Predicting the global distribution of planktonic foraminifera using a
#' dynamic ecosystem model. Biogeosciences, 5(3), 891–911.
#'
"N41.G.ruber.seasonality"

#' @title Climate (surface temperature) at core MD97-2141 from TraCE-21ka
#' @description Modelled surface temperature at core MD97-2141. Model output
#' from TraCE-21ka simulations.
#' @format A matrix with 22040 rows and 12 columns
#' @source Liu, Z., Otto-Bliesner, B. L., He, F., Brady, E. C., Tomas, R.,
#' Clark, P. U., … Cheng, J. (2009).
#' Transient Simulation of Last Deglaciation with a New Mechanism for
#' Bølling-Allerød Warming. Science, 325(5938), 310–314.
#' https://doi.org/10.1126/science.1171041
"N41.t21k.climate"

#' @title Scussolini et al. (2013) Table 1
#' @description Data from table 1 in Scussolini et al. (2013)
#' @format A dataframe with 22 rows and 6 columns
#' @source Scussolini, P., van Sebille, E., & Durgadoo, J. V. (2013).
#' Paleo Agulhas rings enter the subtropical gyre during the penultimate
#' deglaciation. Climate of the Past, 9(6), 2631–2639.
#' https://doi.org/10.5194/cp-9-2631-2013
"scussolini.tab1"

#' @title sedproxy parameters
#' @description Parameters and variables required to generate a pseudo-proxy with \code{ClimToProxyClim}
#' @format A data frame with 12 rows and 4 variables:
#' \describe{
#'   \item{\code{Function argument}}{character Argument name in \code{ClimToProxyClim}}
#'   \item{\code{Description}}{character Description of argument and corresponding variable/parameter}
#'   \item{\code{Possible sources}}{character Source or possible sources of values for arguments}
#'   \item{\code{Default}}{character Default values for arguments}
#'}
"param.tab"


#' @title gisp2 ice core data at annual resolution 
#' @description gisp2 ice core data interpolated to annual resolution 
#' @format A data frame with 49885 rows and 3 variables:
#' \describe{
#'   \item{\code{age.yr.bp}}{Age in years BP}
#'   \item{\code{temperature}}{temperature in deg C}
#'   \item{\code{temperature.rescaled}}{temperature rescaled to resemble 
#'   a d18O record as in Löwemark et al 2008} 
#'}
#' @details GISP2 data from Alley, R.B.. 2004. GISP2 Ice Core Temperature and
#' Accumulation Data. IGBP PAGES/World Data Center for Paleoclimatology Data
#' Contribution Series 2004-013. NOAA/NGDC Paleoclimatology Program, Boulder CO,
#' USA. Interpolated to annual resolution and additionally rescaled to resemble
#' a d18O record. As in Fig.2 of Löwemark, L., Konstantinou, K. I. and Steinke,
#' S.: Bias in foraminiferal multispecies reconstructions of paleohydrographic
#' conditions caused by foraminiferal abundance variations and bioturbational
#' mixing: A model approach, Marine Geology, 256(1–4), 101–106,
#' doi:10.1016/j.margeo.2008.10.005, 2008.
"gisp2.ann"
