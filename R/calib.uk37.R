# UK'37 calibration

# uk37.dat <- ecusdata::mueller.uk37.sst
#
# lm.uk37 <- lm(`UK'37`~`SST (1-12) [°C]`, data = uk37.dat)
# summary(lm.uk37, correlation = TRUE)

#' Convert between UK'37 and Temperature in degrees C
#'
#' @param temperature
#' @param UK37
#' @param point.or.sample
#' @param n
#'
#' @return
#' @export
#'
#' @examples
calib.uk37 <- function(temperature = NULL, UK37 = NULL, point.or.sample = c("point", "sample"), n = 1){

  if (is.null(temperature) & is.null(UK37) |
      is.null(temperature) == FALSE & is.null(UK37) == FALSE){
    stop("One and only one of temperature or UK37 must be supplied")
  }

  type <- match.arg(point.or.sample)

  cfs.mueller <- matrix(c(0.0686612340110185, 0.0328750614815548), ncol = 2, byrow = TRUE)

  vcov.mueller <- structure(c(6.06536807458765e-05, -2.80815422746781e-06, -2.80815422746781e-06,
                              1.46053818728255e-07), .Dim = c(2L, 2L),
                            .Dimnames = list(c("(Intercept)",
                                               "`SST (1-12) [°C]`"), c("(Intercept)", "`SST (1-12) [°C]`")))

  if (type == "sample"){
  cfs.mueller <- mvtnorm::rmvnorm(n=n, mean=cfs.mueller, sigma=vcov.mueller)
  }

  # convert from temperature to UK'37
  if (is.null(UK37)){
    ret <- t(cfs.mueller[, 1] + (outer(cfs.mueller[, 2], temperature, FUN = "*")))
  }

  # convert from UK'37 to temperature
  if (is.null(temperature)){
    ret <- t(t(outer(UK37, cfs.mueller[, 1], FUN = "-")) / cfs.mueller[, 2])
  }

  return(ret)
}


# calib.uk37(temperature = c(1, 2), point.or.sample = "point")
# calib.uk37(temperature = c(1, 2), n = 5, point.or.sample = "sample")
#
#
# calib.uk37(UK37 = as.vector(calib.uk37(temperature = c(1, 2), point.or.sample = "point"))
#            , point.or.sample = "point")
#
#
# calib.uk37(UK37 = as.vector(calib.uk37(temperature = c(1, 2), point.or.sample = "point"))
#            , n = 5, point.or.sample = "sample")
#
#
# calib.uk37(temperature = 1, UK37 = 1)
#
# Tmps <- runif(100, 0, 20)
# UKs <- calib.uk37(temperature = Tmps
#            , n = 3, point.or.sample = "sample")
#
# plot(Tmps, UKs[, 1])
#
#
# Tmps <- runif(100, 0, 20)
# UKs <- calib.uk37(UK37 = Tmps
#                   , n = 3, point.or.sample = "sample")
#
# plot(Tmps, UKs[, 1])

# mod_boot <- function(data, model){
#   model <- model[[1]]
#   data$resp <- simulate(model)[,1]
#   lm1 <-  update(model, resp ~ ., data=data)
#   coef(lm1)
# }
