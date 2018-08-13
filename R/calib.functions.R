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
#' ProxyConversion(temperature = c(10, 20), point.or.sample = "point", proxy.calibration.type = "UK37")
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
                            proxy.calibration.type = "identity",
                            slp.int.means = NULL, slp.int.vcov = NULL,
                            taxon = NULL,
                            point.or.sample = c("point", "sample"), n = 1){
  
  if (is.null(temperature) & is.null(proxy.value) |
      is.null(temperature) == FALSE & is.null(proxy.value) == FALSE){
    stop("One and only one of temperature or proxy.value must be supplied")
  }
  
  pct <- match.arg(proxy.calibration.type,
                   choices = c("identity", "MgCa", "UK37"))
  
  point.or.sample <- match.arg(point.or.sample)
  
  if (point.or.sample == "point" & n > 1){
    stop("Multiple replicates only returned if point.or.sample == 'sample'")}
  
  
  ## Get calibration parameters
  if (pct != "identity"){
    if (is.null(taxon)) taxon <- 1
    if (is.null(slp.int.means)){
      cfs <-
        sedproxy::CalibrationParameters[[pct]][[taxon]][[
          "means"]][c("slope", "intercept")]
      cfs <-  matrix(cfs, ncol = 2, byrow = TRUE)
    }else{cfs <- matrix(slp.int.means, nrow = 1)}
    
    if (is.null(slp.int.vcov)){
      vcov <- 
        sedproxy::CalibrationParameters[[pct]][[taxon]][[
          "vcov"]][c("slope", "intercept"), c("slope", "intercept")]
      
    }else{
      vcov <- slp.int.vcov
    }
    
    if (point.or.sample == "sample"){
      
      if (is.null(slp.int.means) == FALSE & is.null(slp.int.vcov))
        warning("Sampling calibration parameters using user supplied values
              for the mean slope and intercept but the variance covariance matrix for the
              default or named taxon.")
      
      cfs <- mvtnorm::rmvnorm(n=n, mean=cfs, sigma=vcov)
    }
  }
  
  # Do conversion
  
  ## check if matrix input
  if (any(is.matrix(temperature), is.matrix(proxy.value))){
    if (point.or.sample == "sample" & max(ncol(temperature), ncol(proxy.value)) != n) {
      stop("If input is matrix and point.or.sample == 'sample', n must equal ncol(input)")
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
           
           UK37 = {
             # convert from temperature to UK'37
             out <- if (is.null(proxy.value)){
               t(cfs[, 2] + t(temperature) * cfs[, 1])
             } else if (is.null(temperature)){
               # convert from UK'37 to temperature
               t((t(proxy.value) - cfs[, 2]) / cfs[, 1])
             }
           }
           
    )
    ## Vector input
  }else if (any(is.vector(temperature), is.vector(proxy.value))){
    switch(pct,
           
           identity = {
             out <- if (is.null(temperature)){
               matrix(rep(proxy.value, n), ncol = n)
             } else {
               matrix(rep(temperature, n), ncol = n)
             }},
           #  
           MgCa = {
             cfs[,2] <- exp(cfs[,2])
             
             # convert from temperature to MgCa
             out <- 
               if (is.null(proxy.value)){
                 
                 t(cfs[, 2] * exp(outer(cfs[, 1], temperature, FUN = "*")))
               } else if (is.null(temperature)){
                 # convert from MgCa to temperature
                 t(t(log(outer(proxy.value, cfs[, 2], FUN = "/"))) / cfs[, 1])
               }
           },
           
           UK37 = {
             # convert from temperature to UK'37
             out <- if (is.null(proxy.value)){
               t(cfs[, 2] + (outer(cfs[, 1], temperature, FUN = "*")))
             } else if (is.null(temperature)){
               # convert from UK'37 to temperature
               t(t(outer(proxy.value, cfs[, 2], FUN = "-")) / cfs[, 1])
             }
           }
    )
  }
  

  return(out)
}


ScaleError <- function(mean.temperature = NULL,
                       sd.temperature = NULL,
                       mean.proxy.value = NULL,
                       sd.proxy.value = NULL,
                       proxy.calibration.type,
                       slp.int.means = NULL,
                       slp.int.vcov = NULL,
                       taxon = NULL,
                       point.or.sample = c("point", "sample")) {
  
  if (is.null(sd.temperature) & is.null(sd.proxy.value) |
      is.null(sd.temperature) == FALSE &
      is.null(sd.proxy.value) == FALSE) {
    stop("One and only one of sd.temperature or sd.proxy.value must be supplied 
         with corresponding mean value")
  }
  

  if ((class(sd.temperature) == class(mean.temperature)) == FALSE |
      (class(sd.proxy.value) == class(mean.proxy.value)) == FALSE) {
    stop("Both mean and SD of either temperature or proxy value must be provided.")
  }
  
  
  proxy.calibration.type <- match.arg(proxy.calibration.type,
                                      choices = c("MgCa", "UK37"))
  
  t.mean.plus.sd <- if(is.null(mean.temperature)){NULL} else {
    mean.temperature + sd.temperature}
  p.mean.plus.sd <- if(is.null(mean.proxy.value)){NULL} else {
    mean.proxy.value + sd.proxy.value}
  
  out <- ProxyConversion(
      temperature = t.mean.plus.sd,
      proxy.value = p.mean.plus.sd,
      point.or.sample = point.or.sample,
      proxy.calibration.type = proxy.calibration.type,
      taxon = taxon,
      slp.int.means = slp.int.means,
      slp.int.vcov = slp.int.vcov
    ) -
    ProxyConversion(
      temperature = mean.temperature,
      proxy.value = mean.proxy.value,
      point.or.sample = point.or.sample,
      proxy.calibration.type = proxy.calibration.type,
      taxon = taxon,
      slp.int.means = slp.int.means,
      slp.int.vcov = slp.int.vcov
    ) 
  
  return(as.vector(out))
}

# ScaleError(mean.temperature = c(20, 10), sd.temperature = 2,
#            proxy.calibration.type = "MgCa")
# 
# ScaleError(mean.proxy.value = c(2, 4), sd.proxy.value = 0.1,
#            proxy.calibration.type = "MgCa")
# 
# ScaleError(mean.temperature = c(20, 20), sd.temperature = c(2, 1),
#            proxy.calibration.type = "UK37")
# 
# ScaleError(mean.temperature = c(10, 20), sd.temperature = c(2, 2),
#            proxy.calibration.type = "UK37")
# 
# ScaleError(sd.temperature = c(2, 1),
#            proxy.calibration.type = "UK37")

#' Calibration Uncertainty
#'
#' @description Calculates the calibration uncertainty at a given temperature
#'   from the variance-covariance matrix of the coefficients of a fitted
#'   calibration regression equation. Also returns the mean proxy value at the
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

# Temp <- 1:30
#
# C.U <- CalibrationUncertainty(Temp, means = UK37.pars$mueller.uk37$means,
#                        vcov = UK37.pars$mueller.uk37$vcov)
#
# C.U %>%
#   ggplot(aes(x = temperature, y = mu)) +
#   geom_ribbon(aes(ymax = mu + sigma, ymin = mu - sigma,
#                   colour = "1 SD", fill = "1 SD")) +
#   geom_line(aes(colour = "Mean", fill = "Mean")) +
#   scale_fill_discrete("") +
#   scale_color_discrete("")
#
#
# # For Mg/Ca the proxy units need exonentiating
# C.U <- CalibrationUncertainty(10:28, means = MgCa.foram.pars$`G. aequilateralis_350-500`$means,
#                               vcov = MgCa.foram.pars$`G. aequilateralis_350-500`$vcov) %>%
#   mutate(Temp2 = as.vector(ProxyConversion(proxy.value = exp(mu), proxy.calibration.type = "MgCa")),
#          T.upr = as.vector(ProxyConversion(proxy.value = exp(mu+sigma), proxy.calibration.type = "MgCa")),
#          T.lwr = as.vector(ProxyConversion(proxy.value = exp(mu-sigma), proxy.calibration.type = "MgCa")))
#
#
# C.U %>%
#   ggplot(aes(x = temperature, y = exp(mu))) +
#   geom_ribbon(aes(ymax = exp(mu + sigma), ymin = exp(mu - sigma),
#                   colour = "1 SD", fill = "1 SD")) +
#   geom_line(aes(colour = "Mean", fill = "Mean")) +
#   scale_fill_discrete("") +
#   scale_color_discrete("")
#
# C.U %>%
#   ggplot(aes(x = temperature, y = Temp2)) +
#   geom_ribbon(aes(ymax = T.upr, ymin = T.lwr,
#                   colour = "1 SD", fill = "1 SD")) +
#   geom_line(aes(colour = "Mean", fill = "Mean")) +
#   scale_fill_discrete("") +
#   scale_color_discrete("")
#
# C.U %>%
#   ggplot(aes(x = exp(mu), y = Temp2)) +
#   geom_ribbon(aes(ymax = T.upr, ymin = T.lwr,
#                   colour = "1 SD", fill = "1 SD")) +
#   geom_line(aes(colour = "Mean", fill = "Mean")) +
#   scale_fill_discrete("") +
#   scale_color_discrete("")

#
# C.U %>%
#   ggplot(aes(x = temperature, y = Temp2)) +
#   geom_line(aes(colour = "Mean", fill = "Mean")) +
#   scale_fill_discrete("") +
#   scale_color_discrete("")

# ProxyConversion(proxy.value = as.vector(CalibMgCa(temperature = c(21, 22), point.or.sample = "point"))
#                             , point.or.sample = "point", proxy.calibration.type = "MgCa")


CalibrationUncertainty <- function(temperature, means, vcov){
  mm <- cbind(temperature, 1)
  vars <- mm %*% vcov %*% t(mm)
  sds <- sqrt(diag(vars))
  mu <- as.vector(mm %*% means)
  return(tibble::tibble(temperature, mu=mu, sigma=sds))
}


# CalibrationUncertainty.2 <- function(temperature = NULL, proxy.value = NULL,
#                                      proxy.calibration.type,
#                                      slp.int.means = NULL, slp.int.vcov = NULL,
#                                      taxon = NULL){
#
#   if (is.null(temperature) & is.null(proxy.value) |
#       is.null(temperature) == FALSE & is.null(proxy.value) == FALSE){
#     stop("One and only one of temperature or proxy.value must be supplied")
#   }
#
#   proxy.calibration.type <- match.arg(proxy.calibration.type,
#                                       choices = c("MgCa", "UK37"))
#
#   if (proxy.calibration.type == "UK37"){
#     if (is.null(slp.int.means)){
#       slp.int.means <- sedproxy::UK37.pars$mueller.uk37$means[c("slope", "intercept")]
#       #slp.int.means <- matrix(slp.int.means, ncol = 2, byrow = TRUE)
#     }else{slp.int.means <- slp.int.means}
#
#     if (is.null(slp.int.vcov)){
#       slp.int.vcov <- sedproxy::UK37.pars$mueller.uk37$vcov[c("slope", "intercept"), c("slope", "intercept")]
#     }else{
#       slp.int.vcov <- slp.int.vcov
#     }
#   }else if (proxy.calibration.type == "MgCa"){
#     taxon <- if (is.null(taxon)) {"10 Foram Taxa"} else {match.arg(taxon)}
#
#     if (is.null(slp.int.means)){
#       slp.int.means <- sedproxy::MgCa.foram.pars[[taxon]]$means[c("slope", "intercept")]
#
#     }else{slp.int.means <- matrix(slp.int.means, nrow = 1)}
#
#     if (is.null(slp.int.vcov)){
#       slp.int.vcov <- sedproxy::MgCa.foram.pars[[taxon]]$vcov[c("slope", "intercept"), c("slope", "intercept")]
#     }else{
#       slp.int.vcov <- slp.int.vcov
#     }
#   }
#
#
#   mm <- cbind(temperature, 1)
#   vars <- mm %*% slp.int.vcov %*% t(mm)
#   sd.proxy <- sqrt(diag(vars))
#   mu.proxy <- as.vector(mm %*% slp.int.means)
#   if (proxy.calibration.type == "MgCa") {
#     sd.proxy = exp(mu.proxy + sd.proxy) - exp(mu.proxy)
#     mu.proxy = exp(mu.proxy)
#     }
#   mu.temperature <- as.vector(ProxyConversion(proxy.value = mu.proxy,
#                                     proxy.calibration.type = proxy.calibration.type,
#                                     slp.int.means = slp.int.means))
#
#   sd.temperature <- as.vector(ProxyConversion(proxy.value = mu.proxy + sd.proxy,
#                                               proxy.calibration.type = proxy.calibration.type,
#                                               slp.int.means = slp.int.means)) -
#     as.vector(ProxyConversion(proxy.value = mu.proxy,
#                               proxy.calibration.type = proxy.calibration.type,
#                               slp.int.means = slp.int.means))
#
#   if (proxy.calibration.type == "MgCa"){
#     # noise SD needs to be divided by the mean temperature in proxy units in
#     # order to maintain a consistent SD in temperature units.
#     sd.temperature <- sd.temperature
#   }
#
#
#   return(tibble::tibble(temperature, mu.proxy,
#                         sigma.proxy = sd.proxy,
#                         mu.temperature, sigma.temperature = sd.temperature))
# }
#


