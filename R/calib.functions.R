#' Convert between Temperature in Degrees C and Proxy Units
#'
#' @description A wrapper function for accessing proxy - temperature conversion functions
#'
#' @param temperature Temperature in degrees C
#' @param proxy.value Temperature in proxy units
#' @param proxy.calibration.type Type of proxy, e.g. UK37 or MgCa
#' @param point.or.sample Use the "best estimate" calibration parameters,
#' or parameters sampled from the fitted calibration model
#' @param n the number of replicate conversions to make in the case of sampled calibration parameters
#'
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
                            point.or.sample = c("point", "sample"), n = 1){

  if (is.null(temperature) & is.null(proxy.value) |
      is.null(temperature) == FALSE & is.null(proxy.value) == FALSE){
    stop("One and only one of temperature or proxy.value must be supplied")
  }

  proxy.calibration.type <- match.arg(proxy.calibration.type)

  out <- switch(proxy.calibration.type,
                MgCa = CalibMgCa(temperature = temperature, proxy.value = proxy.value,
                                  point.or.sample = point.or.sample, n = n),
                UK37 = CalibUK37(temperature = temperature, proxy.value = proxy.value,
                                  point.or.sample = point.or.sample, n = n)
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
#' CalibUK37(temperature = c(1, 2), n = 5, point.or.sample = "sample")
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
#'
#' lm.uk37 <- lm(`UK'37`~`SST (1-12) [°C]`, data = uk37.dat)
#'
#' vcov(lm.uk37)
#'
#' }
CalibUK37 <- function(temperature = NULL, proxy.value = NULL,
                      point.or.sample = c("point", "sample"), n = 1){

  if (is.null(temperature) & is.null(proxy.value) |
      is.null(temperature) == FALSE & is.null(proxy.value) == FALSE){
    stop("One and only one of temperature or proxy.value must be supplied")
  }

  type <- match.arg(point.or.sample)

  cfs.mueller <- matrix(c(0.0686612340110185, 0.0328750614815548),
                        ncol = 2, byrow = TRUE)

  vcov.mueller <-
    structure(c(6.06536807458765e-05, -2.80815422746781e-06,
                -2.80815422746781e-06, 1.46053818728255e-07),
              .Dim = c(2L, 2L))

  if (type == "sample"){
    cfs.mueller <- mvtnorm::rmvnorm(n=n, mean=cfs.mueller, sigma=vcov.mueller)
  }

  # convert from temperature to UK'37
  if (is.null(proxy.value)){
    out <- t(cfs.mueller[, 1] + (outer(cfs.mueller[, 2], temperature, FUN = "*")))
  }

  # convert from UK'37 to temperature
  if (is.null(temperature)){
    out <- t(t(outer(proxy.value, cfs.mueller[, 1], FUN = "-")) / cfs.mueller[, 2])
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
#' CalibMgCa(temperature = c(1, 2), point.or.sample = "point")
#'
#' ## With random calibration, 5 replicates
#' CalibMgCa(temperature = c(1, 2), n = 5, point.or.sample = "sample")
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
                      point.or.sample = c("point", "sample"), n = 1){

  if (is.null(temperature) & is.null(proxy.value) |
      is.null(temperature) == FALSE & is.null(proxy.value) == FALSE){
    stop("One and only one of temperature or proxy.value must be supplied")
  }

  type <- match.arg(point.or.sample)

  cfs.anand <- matrix(c(0.0958964702727252, -1.08117546459042),
                        ncol = 2, byrow = TRUE)

  # Need the covariance between A and B
  # > dput(colMeans(sma.boot[,c("slope", "elevation")]))
  # structure(c(0.0958964702727252, -1.08117546459042), .Names = c("slope",
  # "elevation")) > dput(var(sma.boot[,c("slope", "elevation")]))
  # structure(c(9.87197167306858e-06, -0.000219581347799779,
  # -0.000219581347799779, 0.00501539679249436), .Dim = c(2L, 2L), .Dimnames =
  # list(c("slope", "elevation"), c("slope", "elevation")))
  # "slope"), c("elevation", "slope")))

   vcov.anand <- structure(c(9.87197167306858e-06, -0.000219581347799779,
                             -0.000219581347799779, 0.00501539679249436),
                           .Dim = c(2L, 2L),
                           .Dimnames = list(c("A", "B"), c("A", "B")))

  if (type == "sample"){
    cfs.anand <- mvtnorm::rmvnorm(n=n, mean=cfs.anand, sigma=vcov.anand)
  }

  cfs.anand[,2] <- exp(cfs.anand[,2])

  # convert from temperature to MgCa
  if (is.null(proxy.value)){
    out <- t(cfs.anand[, 2] * exp(outer(cfs.anand[, 1], temperature, FUN = "*")))
  }

  # convert from MgCa to temperature
  if (is.null(temperature)){
    out <- t(t(log(outer(proxy.value, cfs.anand[, 2], FUN = "/"))) / cfs.anand[, 1])
  }

  return(out)
}
