#' Convert between Temperature in Degrees C and Proxy Units
#'
#' @description A wrapper function for accessing proxy - temperature conversion functions
#'
#' @param temperature Temperature in degrees C
#' @param proxy.value Temperature in proxy units
#' @param calibration.type Type of proxy, e.g. Uk37 or MgCa
#' @param point.or.sample Use the "best estimate" calibration parameters, or
#'   parameters sampled from the fitted calibration model
#' @param n the number of replicate conversions to make in the case of sampled
#'   calibration parameters
#' @param calibration The name of a specific calibration for which calibration parameters
#'   are provided by sedproxy. Currently applies only to calibration.type MgCa.
#' @param slp.int.means Optional user supplied vector of values for the slope
#'   and intercept of the calibration function. Overides the defaults.
#' @param slp.int.vcov Optional user supplied variance covariance matrix
#'   calibration parameters. Overides the defaults.
#' @details Valid entries for calibration are: "Ten planktonic species_350-500", "G. aequilateralis_350-500", "G.
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
#' @examples
#' # From temperature to UK'37
#' ## With fixed calibration
#' ProxyConversion(temperature = c(10, 20), point.or.sample = "point",
#'                 calibration.type = "Uk37")
#' 
#' ## With random calibration, 5 replicates
#' ProxyConversion(temperature = c(1, 2), n = 5, point.or.sample = "sample",
#'                 calibration.type = "Uk37")
#' 
#' 
#' ## Back-transformation with same calibration
#' ProxyConversion(
#'   proxy.value = as.vector(
#'     ProxyConversion(
#'       temperature = c(21, 22),
#'       calibration.type = "Uk37",
#'       point.or.sample = "point"
#'     )
#'   ),
#'   point.or.sample = "point",
#'   calibration.type = "Uk37"
#' )
#' 
#' ## Back-transformation with random calibration
#' ProxyConversion(
#'   proxy.value = as.vector(
#'     ProxyConversion(
#'       temperature = c(21, 22),
#'       calibration.type = "Uk37",
#'      point.or.sample = "point"
#'     )
#'   )
#'   ,
#'   n = 5,
#'   point.or.sample = "sample",
#'   calibration.type = "Uk37"
#' )
#' 
#' ## Incompatible arguments
#' \dontrun{
#' ProxyConversion(temperature = 1, proxy.value = 1)
#' }
ProxyConversion <- function(temperature = NULL, proxy.value = NULL,
                            calibration.type = "identity",
                            slp.int.means = NULL, slp.int.vcov = NULL,
                            calibration = switch(calibration.type,
                                           identity = NA,
                                           Uk37 = "Mueller global",
                                           MgCa = "Ten planktonic species_350-500"),
                            point.or.sample = c("point", "sample"), n = 1){
 
  pct <- match.arg(calibration.type,
                   choices = c("identity", "MgCa", "Uk37"))
  
  point.or.sample <- match.arg(point.or.sample)
  

  if (is.null(temperature) & is.null(proxy.value) |
      is.null(temperature) == FALSE & is.null(proxy.value) == FALSE){
    stop("One and only one of temperature or proxy.value must be supplied")
  }

  ## Check dimensions if matrix
  if (any(is.matrix(temperature), is.matrix(proxy.value))){
    if (point.or.sample == "sample" & max(ncol(temperature), ncol(proxy.value)) != n) {
      stop("If input is matrix and point.or.sample == 'sample', n must equal ncol(input)")
    }
  }

  if (point.or.sample == "point" & n > 1){
    stop("Multiple replicates only returned if point.or.sample == 'sample'")}


  ## Get calibration parameters
  if (pct != "identity"){

    prs <- calibration.parameters
    cfs.vcov <- prs[prs$calibration.type == pct & prs$calibration == calibration, ]
    #if (is.null(calibration)) calibration <- 1
    if (is.null(slp.int.means)){
      cfs <-  matrix(c(cfs.vcov$slope, cfs.vcov$intercept), ncol = 2, byrow = TRUE)
      # cfs <-
      #   sedproxy::CalibrationParameters[[pct]][[calibration]][[
      #     "means"]][c("slope", "intercept")]
      # cfs <-  matrix(cfs, ncol = 2, byrow = TRUE)
    }else{cfs <- matrix(slp.int.means, nrow = 1)}

    if (is.null(slp.int.vcov)){
      vcov <- cfs.vcov$vcov[[1]]
      # vcov <-
      #   sedproxy::CalibrationParameters[[pct]][[calibration]][[
      #     "vcov"]][c("slope", "intercept"), c("slope", "intercept")]

    }else{
      vcov <- slp.int.vcov
    }

    if (point.or.sample == "sample"){

      if (is.null(slp.int.means) == FALSE & is.null(slp.int.vcov))
        warning("Sampling calibration parameters using user supplied values
                for the mean slope and intercept but the variance covariance matrix for the
                default or named calibration.")

      cfs <- mvtnorm::rmvnorm(n=n, mean=cfs, sigma=vcov)
      }
  }

  # Do conversion

  ## check if vector input and convert to 1 column matrix
  is.vec <- any(is.vector(temperature), is.vector(proxy.value))

  if (is.vec){
    if (is.null(temperature)){
      proxy.value <- matrix(rep(proxy.value, n), ncol = n)
    } else if (is.null(proxy.value)){
      temperature <- matrix(rep(temperature, n), ncol = n)
    }
  }

  switch(pct,
         identity = {
           out <- if (is.null(temperature)){
             proxy.value
           } else {
             temperature
           }},
         #
         MgCa = {
           cfs[,2] <- exp(cfs[,2])

           # convert from temperature to MgCa
           out <-
             if (is.null(proxy.value)){
               t(cfs[, 2] * exp(t(temperature) * cfs[, 1]))
             } else if (is.null(temperature)){
               # convert from MgCa to temperature
               t(log(t(proxy.value) / cfs[,2]) / cfs[,1])
             }
         },

         Uk37 = {
           # convert from temperature to UK'37
           out <- if (is.null(proxy.value)){
             t(cfs[, 2] + t(temperature) * cfs[, 1])
           } else if (is.null(temperature)){
             # convert from UK'37 to temperature
             t((t(proxy.value) - cfs[, 2]) / cfs[, 1])
           }
         }
  )

  # convert back to vector if vector input and n == 1
  if (is.vec & n == 1){
    out <- as.vector(out)
  }

  return(out)
}

# # Visual checks of calibration data and parameters
#
# ## Uk37
#
# df <- tibble(x = 1:30,
#              y = ProxyConversion(temperature = x, point.or.sample = "point",
#                                  calibration.type = "Uk37"))
#
# dat <- climproxycalibration::mueller.uk37.sst
#
# dat %>%
#   ggplot(aes(x = `SST (1-12) [°C]`, y = `UK'37`)) +
#   geom_point( ) +
#   geom_line(data = df, aes(x=x, y =y))
#
# t1 <- 1:30
#
# p1 <- ProxyConversion(temperature = t1, point.or.sample = "point",
#                       calibration.type = "Uk37")
#
# t2 <- ProxyConversion(proxy.value = p1, point.or.sample = "point",
#                       calibration.type = "Uk37")
#
# all.equal(t1, t2)
#
#
# ## MgCa
#
# df <- tibble(x = 1:30,
#              y = ProxyConversion(temperature = x, point.or.sample = "point",
#                                  calibration.type = "MgCa"))
#
# dat <- climproxycalibration::anand.2003.hydro
#
# dat %>%
#   ggplot(aes(x = calc.t, y = `Mg/Ca (mmol/mol)`)) +
#   geom_point( ) +
#   geom_line(data = df, aes(x=x, y =y))
#
# t1 <- 1:30
#
# p1 <- ProxyConversion(temperature = t1, point.or.sample = "point",
#                       calibration.type = "MgCa")
#
# t2 <- ProxyConversion(proxy.value = p1, point.or.sample = "point",
#                       calibration.type = "MgCa")
#
# all.equal(t1, t2)