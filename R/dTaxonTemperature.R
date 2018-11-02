#' Foram production as a function of temperature, using published temperature preference parameters.
#'
#' @param temperature temperature in degrees celsius
#' @param taxon taxon ID in table taxa.temperature.prefs,
#'
#' @return
#' @export
#'
#' @examples
dTaxonTemperature <- function(temperature, taxon = taxa.temperature.prefs$ID){

  tax <- match.arg(taxon)

  args <- sedproxy::taxa.temperature.prefs[sedproxy::taxa.temperature.prefs$ID == tax, ]
  d <- dnorm(temperature, mean = args$Topt, sd = args$Tsigma)
  d[temperature >= args$Tmax] <- 0
  d[temperature <= args$Tmin] <- 0
  return(d)
}

