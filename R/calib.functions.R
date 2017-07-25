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
#' # From temperature to UK'37
#' ## With fixed calibration
#' calib.uk37(temperature = c(1, 2), point.or.sample = "point")
#'
#' ## With random calibration, 5 replicates
#' calib.uk37(temperature = c(1, 2), n = 5, point.or.sample = "sample")
#'
#'
#' ## Back-transformation with same calibration
#' calib.uk37(UK37 = as.vector(calib.uk37(temperature = c(21, 22), point.or.sample = "point"))
#'            , point.or.sample = "point")
#'
#' ## Back-transformation with random calibration
#' calib.uk37(UK37 = as.vector(calib.uk37(temperature = c(21, 22), point.or.sample = "point"))
#'            , n = 5, point.or.sample = "sample")
#'
#' ## Incompatible arguments
#' calib.uk37(temperature = 1, UK37 = 1)
calib.uk37 <- function(temperature = NULL, UK37 = NULL, point.or.sample = c("point", "sample"), n = 1){

  if (is.null(temperature) & is.null(UK37) |
      is.null(temperature) == FALSE & is.null(UK37) == FALSE){
    stop("One and only one of temperature or UK37 must be supplied")
  }

  type <- match.arg(point.or.sample)

  cfs.mueller <- matrix(c(0.0686612340110185, 0.0328750614815548),
                        ncol = 2, byrow = TRUE)

  vcov.mueller <-
    structure(c(6.06536807458765e-05, -2.80815422746781e-06,
                -2.80815422746781e-06, 1.46053818728255e-07),
              .Dim = c(2L, 2L),
              .Dimnames = list(
                c("(Intercept)", "`SST (1-12) [°C]`"),
                c("(Intercept)", "`SST (1-12) [°C]`")))

  if (type == "sample"){
    cfs.mueller <- mvtnorm::rmvnorm(n=n, mean=cfs.mueller, sigma=vcov.mueller)
  }

  # convert from temperature to UK'37
  if (is.null(MgCa)){
    out <- t(cfs.mueller[, 1] + (outer(cfs.mueller[, 2], temperature, FUN = "*")))
  }

  # convert from UK'37 to temperature
  if (is.null(temperature)){
    out <- t(t(outer(MgCa, cfs.mueller[, 1], FUN = "-")) / cfs.mueller[, 2])
  }
  return(out)
}


calib.mgca <- function(temperature = NULL, MgCa = NULL, point.or.sample = c("point", "sample"), n = 1){

  if (is.null(temperature) & is.null(MgCa) |
      is.null(temperature) == FALSE & is.null(MgCa) == FALSE){
    stop("One and only one of temperature or MgCa must be supplied")
  }

  type <- match.arg(point.or.sample)

  cfs.anand <- matrix(c(0.09, 0.38),
                        ncol = 2, byrow = TRUE, dimnames = list(NULL, c("A", "B")))

  # Need the covariance between A and B
  if (type == "sample"){
    cfs.anand <- matrix(c(rnorm(n, 0.09, 0.003), rnorm(n, 0.38, 0.02)),
                        ncol = 2, byrow = FALSE, dimnames = list(NULL, c("A", "B")))
  }

  # convert from temperature to UK'37
  if (is.null(MgCa)){
    out <- t(cfs.anand[, "B"] * exp(outer(cfs.anand[, "A"], temperature, FUN = "*")))
  }

  # convert from UK'37 to temperature
  if (is.null(temperature)){
    out <- t(t(log(outer(MgCa, cfs.anand[, "B"], FUN = "/"))) / cfs.anand[, "A"])
  }

  return(out)
}

#
# calib.mgca(temperature = c(15, 25), point.or.sample = "point")
# calib.mgca(MgCa = c(2, 4), point.or.sample = "point")
#
# calib.mgca(temperature = c(15, 25), point.or.sample = "sample", n = 2)
#
# calib.mgca(MgCa = c(2, 4), point.or.sample = "sample", n = 2)
#
# calib.mgca(MgCa = as.vector(
#   calib.mgca(temperature = c(15, 25),point.or.sample = "point"))
#   , point.or.sample = "point")
#
#
# Tmps <- runif(100, 12, 28)
# MgCas <- calib.mgca(temperature = Tmps
#                   , n = 15, point.or.sample = "sample") %>%
#   tbl_df() %>%
#   mutate(Temperature = Tmps) %>%
#   gather(Rep, MgCa, -Temperature)
#
# MgCas %>%
#   ggplot(aes(x = Temperature, y = log(MgCa), group = Rep)) +
#   geom_point()
#
#
#
# MgCas <- runif(100, 1, 6)
#
# Tmps <- calib.mgca(MgCa = MgCas
#                     , n = 105, point.or.sample = "sample") %>%
#   tbl_df() %>%
#   mutate(MgCa = MgCas) %>%
#   gather(Rep, Temperature, -MgCa)
#
# Tmps %>%
#   ggplot(aes(x = MgCa, y = (Temperature), group = Rep)) +
#   geom_line(alpha = 0.05)
#
#
# mean(calib.mgca(temperature = mean(1:10), point.or.sample = "point"))
# mean(calib.mgca(temperature = 1:10, point.or.sample = "point"))

