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
#'
#'   The uncertainty in proxy units can be scaled to temperature units by
#'   dividing by the calibration slope.
#'
#'   For linear calibrations on transformed variabes such as Mg/Ca, this is the
#'   uncertainty on the log(Mg/Ca) values.
#'
#' @param temperature Temperature at which to evaluate proxy uncertainty in
#'   degrees celsius
#' @param means slope and intercept of fitted calibration regression
#' @param vcov variance-covariance matrix of slope and intercept from fitted
#'   calibration regression
#'
#' @return a tibble (dataframe) with temperature, mean proxy unit and 1 x sigma
#'   proxy unit
#' @export
#'
#' @examples
#' Temp <- 10:20
#'
#' cp <- with(calibration.parameters, calibration.parameters[calibration.type == "Uk37", ])
#'
#' C.U <- CalibrationUncertainty(Temp, means = c(cp$slope, cp$intercept),
#'                               vcov = cp$vcov[[1]])
#'
#' C.U %>%
#'   ggplot(aes(x = temperature, y = mu.proxy)) +
#'   geom_ribbon(aes(ymax = mu.proxy + sigma.proxy, ymin = mu.proxy - sigma.proxy,
#'                   colour = "1 SD", fill = "1 SD")) +
#'   geom_line() +
#'   scale_fill_discrete("") +
#'   scale_color_discrete("")
#'
#'
#' # For Mg/Ca the proxy units need exonentiating
#' cp <- with(calibration.parameters, calibration.parameters[calibration == "G. ruber pink_250-350", ])
#'
#' Temp <- 1:30
#'
#'
#' C.U <- CalibrationUncertainty(Temp, means = c(cp$slope, cp$intercept),
#'                               vcov = cp$vcov[[1]]) %>%
#'   mutate(MgCa = exp(mu.proxy),
#'          MgCa.upr = exp(mu.proxy + sigma.proxy),
#'          MgCa.lwr = exp(mu.proxy - sigma.proxy))
#'
#' C.U %>%
#'   ggplot(aes(x = temperature, y = MgCa)) +
#'   geom_ribbon(aes(ymax = MgCa.upr, ymin = MgCa.lwr,
#'                   colour = "1 SD", fill = "1 SD")) +
#'   geom_line() +
#'   scale_fill_discrete("") +
#'   scale_color_discrete("") +
#'   scale_x_continuous(limits = c(-1, 31)) +
#'   scale_y_continuous(limits = c(0, 7))
#'
#'
#' C.U %>%
#'   ggplot(aes(x = MgCa, y = temperature)) +
#'   geom_ribbon(aes(ymax = temperature + sigma.temperature, ymin = temperature - sigma.temperature,
#'                   colour = "1 SD", fill = "1 SD")) +
#'   geom_line() +
#'   scale_fill_discrete("") +
#'   scale_color_discrete("") +
#'   scale_y_continuous(limits = c(-1, 31)) +
#'   scale_x_continuous(limits = c(0, 7)) +
#'   coord_flip()
#'
CalibrationUncertainty <- function(temperature, means, vcov){
  mm <- cbind(temperature, 1)
  vars <- mm %*% vcov %*% t(mm)
  sds <- sqrt(diag(vars))
  mu <- as.vector(mm %*% means)

  return(tibble::tibble(temperature, mu.proxy=mu, sigma.proxy=sds, sigma.temperature = sds / means[1]))
}



# CalibrationUncertainty.2 <- function(temperature = NULL, proxy.value = NULL,
#                                      calibration.type,
#                                      slp.int.means = NULL, slp.int.vcov = NULL,
#                                      calibration = NULL){
#
#   pct <- match.arg(calibration.type,
#                    choices = c("identity", "MgCa", "Uk37"))
#
#
#     if (is.null(temperature) & is.null(proxy.value) |
#       is.null(temperature) == FALSE & is.null(proxy.value) == FALSE){
#     stop("One and only one of temperature or proxy.value must be supplied")
#   }
#
#
#   ## Get calibration parameters
#   if (pct != "identity"){
#     prs <- calibration.parameters
#     cfs.vcov <- prs[prs$calibration.type == pct & prs$calibration == calibration, ]
#     if (is.null(slp.int.means)){
#       cfs <-  matrix(c(cfs.vcov$slope, cfs.vcov$intercept), ncol = 2, byrow = TRUE)
#      }else{cfs <- matrix(slp.int.means, nrow = 1)}
#
#     if (is.null(slp.int.vcov)){
#       vcov <- cfs.vcov$vcov[[1]]
#       }else{
#       vcov <- slp.int.vcov
#       }
#     }
#
#   mm <- cbind(temperature, 1)
#   vars <- mm %*% slp.int.vcov %*% t(mm)
#   sd.proxy <- sqrt(diag(vars))
#   mu.proxy <- as.vector(mm %*% slp.int.means)
#
#   if (calibration.type == "MgCa") {
#     sd.proxy = exp(mu.proxy + sd.proxy) - exp(mu.proxy)
#     mu.proxy = exp(mu.proxy)
#   }
#
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




#' Calibration (regression) error from RSE and x values of fitted regression
#'
#' @param x_star New x values for which error in y is desired
#' @param rse Residual Standard Error of fitted regression
#' @param x x values from fitted regression
#'
#' @return vector
#' @export
#'
#' @examples
#' n <- 20
#' df <- data.frame(x = runif(n, 50, 100))
#' df$y <- 1 + 0.5 * df$x + rnorm(n, 0, 10)
#'
#' plot(df)
#' lm1 <- lm(y~x, data = df)
#' abline(lm1)
#'
#' ####
#'
#' library(dplyr)
#' library(ggplot2)
#'
#' n <- 20
#' df <- data.frame(x = runif(n, 50, 100))
#' df$y <- 1 + 0.5 * df$x + rnorm(n, 0, 10)
#'
#' plot(df)
#' lm1 <- lm(y~x, data = df)
#' abline(lm1)
#'
#' CalErr(c(23, 56, 79), summary(lm1)$sigma, df$x)
#'
#' df2 <- tibble(x = 0:200,
#'               y = predict(lm1, newdata = data.frame(x=x))) %>%
#'   mutate(y.err = CalErr(x_star = x, rse = summary(lm1)$sigma, x = df$x))
#'
#' df %>%
#'   ggplot(aes(x , y)) +
#'   geom_point() +
#'   geom_smooth(method = "lm") +
#'   geom_line(data = df2, aes(x, y = y + 2*y.err))+
#'   geom_line(data = df2, aes(x, y = y - 2*y.err))
#'
CalErr <- function(x_star, rse, x){
  n <- length(x)
  mu_x <- mean(x)
  rse * sqrt(1/n + (x_star - mu_x)^2 / sum((x - mu_x)^2))
}


# From Gebregiorgis, D., Hathorne, E. C., Giosan, L., Clemens, S., Nürnberg, D.
# and Frank, M.: Southern Hemisphere forcing of South Asian monsoon
# precipitation over the past ~1 million years, Nature Communications, 9(1),
# 4702, doi:10.1038/s41467-018-07076-2, 2018.

# sig.sq_T_Mg_Ca <- function(Mg_Ca, a, sigma_a, b, sigma_b, sigma_Mg_Ca){
#
#   dT_da <- - 1/a^2 * log(Mg_Ca / b)
#
#   dT_db <- -1/(a*b)
#
#   dT_dMg_Ca <- 1/a * 1/Mg_Ca
#
#
#   sig.sq_T <- (dT_da * sigma_a)^2 + (dT_db * sigma_b)^2 + (dT_dMg_Ca * sigma_Mg_Ca)^2
#
#   return(sig.sq_T)
# }



