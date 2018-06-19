#' Convert between Temperature in Degrees C and Proxy Units
#'
#' @description A wrapper function for accessing proxy - temperature conversion functions
#'
#' @param temperature Temperature in degrees C
#' @param proxy.value Temperature in proxy units
#' @param proxy.calibration.type Type of proxy, e.g. UK37 or MgCa
#' @param point.or.sample Use the "best estimate" calibration parameters, or
#'   parameters sampled from the fitted calibration model
#' @param n the number of replicate conversions to make in the case of sampled
#'   calibration parameters
#' @param taxon The name of a specific taxon for which calibration parameters
#'   are provided by sedproxy. Currently applies only to proxy.calibration.type MgCa.
#' @param slp.int.means Optional user supplied vector of values for the slope
#'   and intercept of the calibration function. Overides the defaults.
#' @param slp.int.vcov Optional user supplied variance covariance matrix
#'   calibration parameters. Overides the defaults.
#' @details Valid entries for taxon are: "10 Foram Taxa", "G. aequilateralis_350-500", "G.
#'   aequilateralis_500-1000", "G. conglobatus_350-500", "G. hirsuta_350-500",
#'   "G. inflata_350-500", "G. ruber pink_250-350", "G. ruber pink_350-500", "G.
#'   ruber white_250-350", "G. ruber white_350-500", "G. sacculifer with
#'   sac_350-500", "G. sacculifer without sac_350-500", "G.
#'   truncatulinoides_350-500", "G. truncatulinoides_500-1000", "N.
#'   dutertrei_350-500", "O. univesa_350-500", "P. obliquiloculata_350-500"
#' @return a vector of temperatures or proxy values
#' @export
#' @importFrom mvtnorm rmvnorm
#' @family calib
#'
#' @examples
#' # From temperature to UK'37
#' ## With fixed calibration
#' ProxyConversion(temperature = c(1, 2), point.or.sample = "point", proxy.calibration.type = "UK37")
#'
#' ## With random calibration, 5 replicates
#' ProxyConversion(temperature = c(1, 2), n = 5, point.or.sample = "sample", proxy.calibration.type = "UK37")
#'
#'
#' ## Back-transformation with same calibration
#' ProxyConversion(proxy.value = as.vector(CalibUK37(temperature = c(21, 22), point.or.sample = "point"))
#'            , point.or.sample = "point", proxy.calibration.type = "UK37")
#'
#' ## Back-transformation with random calibration
#' ProxyConversion(proxy.value = as.vector(CalibUK37(temperature = c(21, 22), point.or.sample = "point"))
#'            , n = 5, point.or.sample = "sample", proxy.calibration.type = "UK37")
#'
#' ## Incompatible arguments
#' \dontrun{
#' CalibUK37(temperature = 1, proxy.value = 1)
#' }
ProxyConversion <- function(temperature = NULL, proxy.value = NULL,
                            proxy.calibration.type = c("MgCa", "UK37"),
                            slp.int.means = NULL, slp.int.vcov = NULL,
                            taxon = NULL,
                            point.or.sample = c("point", "sample"), n = 1){

  if (is.null(temperature) & is.null(proxy.value) |
      is.null(temperature) == FALSE & is.null(proxy.value) == FALSE){
    stop("One and only one of temperature or proxy.value must be supplied")
  }

  proxy.calibration.type <- match.arg(proxy.calibration.type)

  out <- switch(proxy.calibration.type,
                MgCa = CalibMgCa(temperature = temperature, proxy.value = proxy.value,
                                  point.or.sample = point.or.sample, n = n,
                                 taxon = taxon,
                                 slp.int.means = slp.int.means, slp.int.vcov = slp.int.vcov),
                UK37 = CalibUK37(temperature = temperature, proxy.value = proxy.value,
                                 point.or.sample = point.or.sample, n = n,
                                 slp.int.means = slp.int.means, slp.int.vcov = slp.int.vcov)
                )
  return(out)
}


#' Convert between UK'37 and Temperature in degrees C
#'
#' @inheritParams ProxyConversion
#'
#' @return a vector of temperatures or proxy values
#' @export
#' @family calib
#'
#' @examples
#' # From temperature to UK'37
#' ## With fixed calibration
#' CalibUK37(temperature = c(1, 2), point.or.sample = "point")
#'
#' ## With random calibration, 5 replicates
#' CalibUK37(temperature = c(20, 25), n = 5, point.or.sample = "sample")
#'
#'
#' ## Back-transformation with same calibration
#' CalibUK37(proxy.value = as.vector(CalibUK37(temperature = c(21, 22), point.or.sample = "point"))
#'            , point.or.sample = "point")
#'
#' ## Back-transformation with random calibration
#' CalibUK37(proxy.value = as.vector(CalibUK37(temperature = c(21, 22), point.or.sample = "point"))
#'            , n = 5, point.or.sample = "sample")
#'
#' ## Incompatible arguments
#' \dontrun{
#' CalibUK37(temperature = 1, proxy.value = 1)
#' }
#' @details To get the UK'37 calibration parameters
#' \preformatted{
#' uk37.dat <- ecusdata::mueller.uk37.sst
#' lm.uk37 <- lm(`UK'37`~ 1 + `SST (1-12) [°C]`, data = uk37.dat)
#' means.mueller <- coef(lm.uk37)[2:1]
#' vcov.mueller <- vcov(lm.uk37)[2:1, 2:1]
#' }
#'
CalibUK37 <- function(temperature = NULL, proxy.value = NULL,
                      point.or.sample = c("point", "sample"), n = 1,
                      slp.int.means = NULL, slp.int.vcov = NULL
                      ){

  if (is.null(temperature) & is.null(proxy.value) |
      is.null(temperature) == FALSE & is.null(proxy.value) == FALSE){
    stop("One and only one of temperature or proxy.value must be supplied")
  }

  type <- match.arg(point.or.sample)

  if (is.null(slp.int.means)){
    cfs <- UK37.pars$mueller.uk37$means
    cfs <-  matrix(cfs, ncol = 2, byrow = TRUE)
    }else{cfs <- matrix(slp.int.means, nrow = 1)}

  if (is.null(slp.int.vcov)){
    vcov <- UK37.pars$mueller.uk37$vcov
  }else{
      vcov <- slp.int.vcov
    }

  if (type == "sample"){

    if (is.null(slp.int.means) == FALSE & is.null(slp.int.vcov))
      warning("Sampling calibration parameters using user supplied values
              for the mean slope and intercept but the variance covariance matrix for the
              default or named taxon.")

    cfs <- mvtnorm::rmvnorm(n=n, mean=cfs, sigma=vcov)
  }

  # convert from temperature to UK'37
  if (is.null(proxy.value)){
    out <- t(cfs[, 2] + (outer(cfs[, 1], temperature, FUN = "*")))
  }

  # convert from UK'37 to temperature
  if (is.null(temperature)){
    out <- t(t(outer(proxy.value, cfs[, 2], FUN = "-")) / cfs[, 1])
  }
  return(out)
}

#' Convert between MgCa and Temperature in degrees C
#'
#' @inheritParams ProxyConversion
#'
#' @return a vector of temperatures or proxy values
#' @export
#' @family calib
#'
#' @examples
#' # From temperature to MgCa
#' ## With fixed calibration
#' CalibMgCa(temperature = c(20, 25), point.or.sample = "point")
#'
#' ## With random calibration, 5 replicates
#' CalibMgCa(temperature = c(20, 25), n = 5, point.or.sample = "sample")
#'
#'
#' ## Back-transformation with same calibration
#' CalibMgCa(proxy.value = as.vector(CalibMgCa(temperature = c(21, 22), point.or.sample = "point"))
#'            , point.or.sample = "point")
#'
#' ## Back-transformation with random calibration
#' CalibMgCa(proxy.value = as.vector(CalibMgCa(temperature = c(21, 22), point.or.sample = "point"))
#'            , n = 5, point.or.sample = "sample")
#'
#' ## Incompatible arguments
#' \dontrun{
#' CalibMgCa(temperature = 1, proxy.value = 1)
#' }
CalibMgCa <- function(temperature = NULL, proxy.value = NULL,
                      point.or.sample = c("point", "sample"), n = 1,
                      slp.int.means = NULL, slp.int.vcov = NULL,
                      taxon = c("10 Foram Taxa", "G. aequilateralis_350-500", "G. aequilateralis_500-1000",
                                "G. conglobatus_350-500", "G. hirsuta_350-500", "G. inflata_350-500",
                                "G. ruber pink_250-350", "G. ruber pink_350-500", "G. ruber white_250-350",
                                "G. ruber white_350-500", "G. sacculifer with sac_350-500", "G. sacculifer without sac_350-500",
                                "G. truncatulinoides_350-500", "G. truncatulinoides_500-1000",
                                "N. dutertrei_350-500", "O. univesa_350-500", "P. obliquiloculata_350-500")
                      ){

  if (is.null(temperature) & is.null(proxy.value) |
      is.null(temperature) == FALSE & is.null(proxy.value) == FALSE){
    stop("One and only one of temperature or proxy.value must be supplied")
  }

  type <- match.arg(point.or.sample)
  taxon <- if (is.null(taxon)) {"10 Foram Taxa"} else {match.arg(taxon)}

  if (is.null(slp.int.means)){
    cfs <- MgCa.foram.pars[[taxon]]$means[c("slope", "intercept")]
    cfs <-  matrix(cfs, ncol = 2, byrow = TRUE)
  }else{cfs <- matrix(slp.int.means, nrow = 1)}

  if (is.null(slp.int.vcov)){
    vcov <- MgCa.foram.pars[[taxon]]$vcov[c("slope", "intercept"), c("slope", "intercept")]
    }else{
      vcov <- slp.int.vcov
      }

  if (type == "sample"){

    if (is.null(slp.int.means) == FALSE & is.null(slp.int.vcov))
      warning("Sampling calibration parameters using user supplied values
              for the mean slope and intercept but the variance covariance matrix for the
              default or named taxon.")

    cfs <- mvtnorm::rmvnorm(n=n, mean=cfs, sigma=vcov)
  }

  cfs[,2] <- exp(cfs[,2])

  # convert from temperature to MgCa
  if (is.null(proxy.value)){
    out <- t(cfs[, 2] * exp(outer(cfs[, 1], temperature, FUN = "*")))
  }

  # convert from MgCa to temperature
  if (is.null(temperature)){
    out <- t(t(log(outer(proxy.value, cfs[, 2], FUN = "/"))) / cfs[, 1])
  }

  return(out)
}



#' Calibration Uncertainty
#'
#' @description Calculates the calibration uncertainty at a given temperature
#'   from the variance-covariance matrix of the coefficients of a fitted
#'   calibration regression equation. Also returns the mean proxy value at the a
#'   given temperature.
#' @param temperature temperature at which to evaluate proxy uncertainty in
#'   degrees celsius
#' @param means intercept and slope of fitted calibration regression
#' @param vcov variance-covariance matrix of slope and intercept from fitted
#'   calibration regression
#'
#' @return a tibble (dataframe) with temperature, mean proxy unit and 1 x sigma
#'   proxy unit
#' @export
#'
#' @examples
CalibrationUncertainty <- function(temperature, means, vcov){
  mm <- cbind(1, temperature)
  vars <- mm %*% vcov %*% t(mm)
  sds <- sqrt(diag(vars))
  mu <- as.vector(mm %*% means)
  return(tibble::tibble(temperature, mu=mu, sigma=sds))
}


