ScaleError <- function(mean.temperature = NULL,
                       sd.temperature = NULL,
                       mean.proxy.value = NULL,
                       sd.proxy.value = NULL,
                       calibration.type,
                       slp.int.means = NULL,
                       slp.int.vcov = NULL,
                       calibration = NULL,
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


  calibration.type <- match.arg(calibration.type,
                                choices = c("MgCa", "Uk37"))

  t.mean.plus.sd <- if(is.null(mean.temperature)){NULL} else {
    mean.temperature + sd.temperature}
  p.mean.plus.sd <- if(is.null(mean.proxy.value)){NULL} else {
    mean.proxy.value + sd.proxy.value}

  out <- ProxyConversion(
    temperature = t.mean.plus.sd,
    proxy.value = p.mean.plus.sd,
    point.or.sample = point.or.sample,
    calibration.type = calibration.type,
    calibration = calibration,
    slp.int.means = slp.int.means,
    slp.int.vcov = slp.int.vcov
  ) -
    ProxyConversion(
      temperature = mean.temperature,
      proxy.value = mean.proxy.value,
      point.or.sample = point.or.sample,
      calibration.type = calibration.type,
      calibration = calibration,
      slp.int.means = slp.int.means,
      slp.int.vcov = slp.int.vcov
    )

  return(as.vector(out))
  }

# ScaleError(mean.temperature = c(20, 10), sd.temperature = 2,
#            calibration.type = "MgCa")
#
# ScaleError(mean.proxy.value = c(2, 4), sd.proxy.value = 0.1,
#            calibration.type = "MgCa")
#
# ScaleError(mean.temperature = c(20, 20), sd.temperature = c(2, 1),
#            calibration.type = "Uk37")
#
# ScaleError(mean.temperature = c(10, 20), sd.temperature = c(2, 2),
#            calibration.type = "Uk37")
#
# ScaleError(sd.temperature = c(2, 1),
#            calibration.type = "Uk37")

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
# C.U <- CalibrationUncertainty(Temp, means = Uk37.pars$mueller.uk37$means,
#                        vcov = Uk37.pars$mueller.uk37$vcov)
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
#   mutate(Temp2 = as.vector(ProxyConversion(proxy.value = exp(mu), calibration.type = "MgCa")),
#          T.upr = as.vector(ProxyConversion(proxy.value = exp(mu+sigma), calibration.type = "MgCa")),
#          T.lwr = as.vector(ProxyConversion(proxy.value = exp(mu-sigma), calibration.type = "MgCa")))
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
CalibrationUncertainty <- function(temperature, means, vcov){
  mm <- cbind(temperature, 1)
  vars <- mm %*% vcov %*% t(mm)
  sds <- sqrt(diag(vars))
  mu <- as.vector(mm %*% means)
  return(tibble::tibble(temperature, mu=mu, sigma=sds))
}


# CalibrationUncertainty.2 <- function(temperature = NULL, proxy.value = NULL,
#                                      calibration.type,
#                                      slp.int.means = NULL, slp.int.vcov = NULL,
#                                      calibration = NULL){
#
#   if (is.null(temperature) & is.null(proxy.value) |
#       is.null(temperature) == FALSE & is.null(proxy.value) == FALSE){
#     stop("One and only one of temperature or proxy.value must be supplied")
#   }
#
#   calibration.type <- match.arg(calibration.type,
#                                       choices = c("MgCa", "Uk37"))
#
#   if (calibration.type == "Uk37"){
#     if (is.null(slp.int.means)){
#       slp.int.means <- sedproxy::Uk37.pars$mueller.uk37$means[c("slope", "intercept")]
#       #slp.int.means <- matrix(slp.int.means, ncol = 2, byrow = TRUE)
#     }else{slp.int.means <- slp.int.means}
#
#     if (is.null(slp.int.vcov)){
#       slp.int.vcov <- sedproxy::Uk37.pars$mueller.uk37$vcov[c("slope", "intercept"), c("slope", "intercept")]
#     }else{
#       slp.int.vcov <- slp.int.vcov
#     }
#   }else if (calibration.type == "MgCa"){
#     calibration <- if (is.null(calibration)) {"Ten planktonic species_350-500"} else {match.arg(calibration)}
#
#     if (is.null(slp.int.means)){
#       slp.int.means <- sedproxy::MgCa.foram.pars[[calibration]]$means[c("slope", "intercept")]
#
#     }else{slp.int.means <- matrix(slp.int.means, nrow = 1)}
#
#     if (is.null(slp.int.vcov)){
#       slp.int.vcov <- sedproxy::MgCa.foram.pars[[calibration]]$vcov[c("slope", "intercept"), c("slope", "intercept")]
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
#   if (calibration.type == "MgCa") {
#     sd.proxy = exp(mu.proxy + sd.proxy) - exp(mu.proxy)
#     mu.proxy = exp(mu.proxy)
#     }
#   mu.temperature <- as.vector(ProxyConversion(proxy.value = mu.proxy,
#                                     calibration.type = calibration.type,
#                                     slp.int.means = slp.int.means))
#
#   sd.temperature <- as.vector(ProxyConversion(proxy.value = mu.proxy + sd.proxy,
#                                               calibration.type = calibration.type,
#                                               slp.int.means = slp.int.means)) -
#     as.vector(ProxyConversion(proxy.value = mu.proxy,
#                               calibration.type = calibration.type,
#                               slp.int.means = slp.int.means))
#
#   if (calibration.type == "MgCa"){
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
