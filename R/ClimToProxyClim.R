#' Title
#'
#' @param clim.signal The "assumed true" climate signal, e.g. climate model output or
#'   instrumental record. A years x 12 (months) matrix of temperatures.
#' @param timepoints The timepoints for which the proxy record is to be modelled
#' @param seas.prod The seasonal pattern of productivity for the organism(s)
#'   archived in the proxy. Either a vector of 12 values or a matrix of the same
#'   dimensions as clim.signal. Defaults to a uniform seasonal distribution.
#' @param bio.depth Depth of the bioturbated layer in metres, defaults to 0.1 m
#' @param sed.acc.rate Sediment accumulation rate in metres per year. Defaults
#'   to 5e-04 m per year (0.5 m / kyr). Either a single value, or vector of same
#'   length as "timepoints"
#' @param meas.noise The amount of noise to add to each simulated proxy value.
#'   Defined as the standard deviation of a normal distribution with mean = 0
#' @param meas.bias The amount of bias to add to each simulated proxy
#'   time-series. Each replicate proxy time-series has a constant bias added,
#'   drawn from a normal distribution with mean = 0, sd = meas.bias. Bias
#'   defaults to zero.
#' @param n.samples Number of e.g. foraminifera sampled per timepoint
#' @param n.replicates Number of replicate proxy time-series to simulate from
#'   the climate signal
#'
#' @return a list with three elements: a dataframe \code{simulated.proxy}, a
#'   dataframe \code{smoothed.signal}, and a list \code{everything}
#'
#'   The dataframe \code{simulated.proxy} contains a single realisation of the
#'   final forward modelled proxy as well as the intermediate stages and the
#'   original climate signal at the requested timepoints.
#'
#' \tabular{ll}{
#' \bold{Variable} \tab \bold{Description} \cr
#' timepoints                 \tab requested timepoints                                                                                                               \cr
#' clim.signal.timepoints.100 \tab 100 year means of climate signal evaluated at the requested timepoints                                                             \cr
#' biot.sig.inf               \tab proxy after bioturbation                                                                                                           \cr
#' proxy.sig.inf              \tab proxy after bioturbation + seasonal bias                                                                                           \cr
#' sed.acc.rate               \tab sediment accumulation rates for each timepoint                                                                                     \cr
#' smoothing.width            \tab weighted mean time span represented in a sample after bioturbation                                                                 \cr
#' proxy.sig.samp             \tab proxy after bioturbation, seasonal bias and finite sampling                                                                        \cr
#' proxy.sig.inf.b            \tab proxy after bioturbation, seasonal bias and calibration bias                                                                       \cr
#' proxy.sig.inf.b.n          \tab proxy after bioturbation, seasonal bias, calibration bias and measurement noise                                                    \cr
#' proxy.sig.samp.b           \tab proxy after bioturbation, seasonal bias, finite sampling and calibration bias                                                      \cr
#' proxy.sig.samp.b.n         \tab proxy after bioturbation, seasonal bias, finite sampling, calibration bias and measurement noise                                   \cr
#' simulated.proxy            \tab final simulated proxy, this will be same as proxy.sig.inf.b.n when n.samples = Inf, and proxy.sig.samp.b.n when n.samples is finite
#'}
#'
#'   The dataframe \code{smoothed.signal} contains the original climate signal
#'   at 100 year resolution.
#'
#'   The list \code{everything} contains all of the above, but if n.replicates >
#'   1, stages of the forward model containing stochastic elements are returned
#'   as matrices with \bold{n.timepoints} rows and \bold{n.replicates} columns.
#'
#' @importFrom dplyr tbl_df
#' @export
#'
#' @examples
ClimToProxyClim <- function(clim.signal,
                            timepoints,
                            seas.prod = rep(1, 12),
                            bio.depth = 0.1,
                            sed.acc.rate = 5e-04,
                            meas.noise = 0,
                            meas.bias = 0,
                            n.samples = Inf,
                            n.replicates = 1) {
  # Check inputs --------
  stopifnot(is.matrix(clim.signal))
  stopifnot(length(sed.acc.rate) == length(timepoints) |
              length(sed.acc.rate) == 1)
  if (is.matrix(seas.prod))
    stop("Matrix form of seasonality not yet supported")

  # Calculate timepoint invariant values ------
  max.clim.signal.i <- nrow(clim.signal)
  sig.years.i <- 1:max.clim.signal.i

  n.timepoints <- length(timepoints)

  if (length(sed.acc.rate) == 1) {
    sed.acc.rate <- rep(sed.acc.rate, n.timepoints)
  }

  # At 100 yr intervals
  timepoints.100 <- seq(1, nrow(clim.signal), by = 100)
  clim.signal.100 <- sapply(timepoints.100, function(tp){

    avg.window.i.1 <- (-50:49) + tp

    if (max(avg.window.i.1) > nrow(clim.signal)) {
      warning("Climate 100 yr average window extends below end of clim.signal")
    }

    avg.window.i <-
      avg.window.i.1[avg.window.i.1 > 0 &
                       avg.window.i.1 < nrow(clim.signal)]

    stopifnot(avg.window.i > 0)
    stopifnot(nrow(clim.signal) > max(avg.window.i))

    mean(clim.signal[avg.window.i, , drop = FALSE])
  })

  # For each timepoint ------
  out <- sapply(1:n.timepoints, function(tp) {
    # Get bioturbation window ----------
    bio.depth.timesteps <- round(bio.depth / sed.acc.rate[tp])
    bioturb.window <- (-1*bio.depth.timesteps):(3*bio.depth.timesteps)

    # Get bioturbation weights --------
    bioturb.weights <-
      ImpulseResponse(-bioturb.window, bio.depth.timesteps, z0 = 0)
    bioturb.weights <- bioturb.weights / sum(bioturb.weights)

    smoothing.width = sum(bioturb.weights*abs(bioturb.window))
    # Check depth and time order match
    # plot(bioturb.weights, (bioturb.window), type = "l", ylim = rev(range(bioturb.window)))

    # get portion of clim.signal corresponding to bioturbation window -------
    sig.window.i.1 <- bioturb.window + timepoints[tp]

    if (max(sig.window.i.1) > nrow(clim.signal)) {
      warning("Bioturbation window extends below end of clim.signal")
    }

    sig.window.i <-
      sig.window.i.1[sig.window.i.1 > 0 &
                       sig.window.i.1 < nrow(clim.signal)]

    stopifnot(sig.window.i > 0)
    stopifnot(nrow(clim.signal) > max(sig.window.i))

    clim.sig.window <- clim.signal[sig.window.i, , drop = FALSE]


    # Get bioturbation X no-seasonality weights matrix ---------
    biot.sig.weights <- bioturb.weights %o% rep(1, 12)
    biot.sig.weights <-
      biot.sig.weights[sig.window.i.1 > 0 &
                         sig.window.i.1 < nrow(clim.signal), , drop = FALSE]
    biot.sig.weights <- biot.sig.weights / sum(biot.sig.weights)


    # Get bioturbation X seasonality weights matrix ---------
    clim.sig.weights <- bioturb.weights %o% seas.prod
    clim.sig.weights <-
      clim.sig.weights[sig.window.i.1 > 0 &
                         sig.window.i.1 < nrow(clim.signal), , drop = FALSE]
    clim.sig.weights <- clim.sig.weights / sum(clim.sig.weights)

    # Check weights sum to 1, within tolerance
    weight.err <- abs(sum(clim.sig.weights) - 1)
    if ((weight.err < 1e-10) == FALSE) stop(paste0("weight.err = ", weight.err))


    # Calculate mean clim.signal -------
    biot.sig.inf <- sum(biot.sig.weights * clim.sig.window)

    proxy.sig.inf <- sum(clim.sig.weights * clim.sig.window)

    if (is.infinite(n.samples)) {
      proxy.sig.samp <- rep(NA, n.replicates)
    } else if (is.finite(n.samples)) {
      # call sample once for all replicates together, then take means of
      # groups of n.samples
      samp <-  sample(clim.sig.window,
                      n.samples * n.replicates,
                      prob = clim.sig.weights,
                      replace = TRUE)

      samp <- matrix(samp, nrow = n.samples)
      proxy.sig.samp <- plyr::aaply(samp, 2, mean)
    }

    # get 100 year clim.average at timepoints -------
    avg.window.i.1 <- (-50:49) + timepoints[tp]

    if (max(avg.window.i.1) > nrow(clim.signal)) {
      warning("Climate average window extends below end of clim.signal")
    }

    avg.window.i <-
      avg.window.i.1[avg.window.i.1 > 0 &
                       avg.window.i.1 < nrow(clim.signal)]

    stopifnot(avg.window.i > 0)
    stopifnot(nrow(clim.signal) > max(avg.window.i))

    clim.signal.timepoints.100 <- mean(clim.signal[avg.window.i, , drop = FALSE])


    # Gather output ----------
    list(
      smoothing.width = smoothing.width,
      clim.signal.timepoints.100 = clim.signal.timepoints.100,
      biot.sig.inf = biot.sig.inf,
      proxy.sig.inf = proxy.sig.inf,
      proxy.sig.samp = proxy.sig.samp)
  })

  #out <- apply(out, 1, function(x) simplify2array(x))
  # use plyr::alply to always return a list
  out <- plyr::alply(out, 1, function(x) simplify2array(x), .dims = TRUE)

  # remove extra attributes added by alply
  attr(out, "split_type") <- NULL
  attr(out, "split_labels") <- NULL

  out$proxy.sig.samp <- t(out$proxy.sig.samp)

  # Add bias and noise to infinite sample --------
  if (meas.bias != 0) {
    bias <- rnorm(n = n.replicates, mean = 0, sd = meas.bias)
  } else{
    bias <- rep(0, n.replicates)
  }
  if (meas.noise != 0) {
    noise <- rnorm(n = n.replicates * n.timepoints, mean = 0, sd = meas.noise)
  }else{
    noise <- rep(0, n.replicates)
  }

  out$proxy.sig.inf.b <- outer(out$proxy.sig.inf, bias, FUN = "+")
  out$proxy.sig.inf.b.n <- out$proxy.sig.inf.b + noise

  # Add bias and noise to finite sample --------
  out$proxy.sig.samp.b <- out$proxy.sig.samp + bias
  out$proxy.sig.samp.b.n <- out$proxy.sig.samp.b + noise

  # Add items to output list -----------
  out$timepoints = timepoints
  out$clim.signal.ann = rowSums(clim.signal[timepoints,  , drop = FALSE]) / ncol(clim.signal)
  out$sed.acc.rate = sed.acc.rate
  out$timepoints.100 = timepoints.100
  out$clim.signal.100 = clim.signal.100

  # Organise output -------
  simulated.proxy <-
    dplyr::tbl_df(out[c(
      "timepoints",
      "clim.signal.timepoints.100",
      "biot.sig.inf",
      "proxy.sig.inf",
      "sed.acc.rate",
      "smoothing.width"
    )])

  simulated.proxy$proxy.sig.samp <- out$proxy.sig.samp[1, ]
  simulated.proxy$proxy.sig.inf.b <- out$proxy.sig.inf.b[, 1]
  simulated.proxy$proxy.sig.inf.b.n <- out$proxy.sig.inf.b.n[, 1]
  simulated.proxy$proxy.sig.samp.b <- out$proxy.sig.samp.b[1, ]
  simulated.proxy$proxy.sig.samp.b.n <- out$proxy.sig.samp.b.n[1, ]

  if (is.finite(n.samples)) {simulated.proxy$simulated.proxy <- simulated.proxy$proxy.sig.samp.b.n}else{
    simulated.proxy$simulated.proxy <- simulated.proxy$proxy.sig.inf.b.n
  }


  smoothed.signal <- dplyr::tbl_df(out[c(
    "timepoints.100",
    "clim.signal.100"
    )])

  return(list(simulated.proxy=simulated.proxy,
              smoothed.signal=smoothed.signal,
              everything = out))
}



#' Convert "everything" part of output from ClimToProxyClim to dataframe
#'
#' @param PFM output from ClimToProxyClim
#'
#' @return
#' @export
#'
#' @examples
MakePFMDataframe <- function(PFM){
  df <- data.frame(
    proxy.sig.samp = as.vector(PFM$proxy.sig.samp),
    proxy.sig.inf.b = as.vector(PFM$proxy.sig.inf.b),
    proxy.sig.samp.b = as.vector(PFM$proxy.sig.samp.b),
    proxy.sig.inf.b.n = as.vector(PFM$proxy.sig.inf.b.n),
    proxy.sig.samp.b.n = as.vector(PFM$proxy.sig.samp.b.n),
    stringsAsFactors = FALSE)

  df$Age <- PFM$timepoints
  df$replicate <- rep(1:ncol(PFM$proxy.sig.inf.b), each = length(PFM$timepoints))
  df <- tbl_df(df) %>%
    gather(Stage, value, -Age, -replicate)

  #df$proxy.sig.inf = as.vector(PFM$proxy.sig.inf)

  biot.inf <- data.frame(
    replicate = 1,
    Age = PFM$timepoints,
    Stage = "biot.sig.inf",
    value = PFM$biot.sig.inf,
    stringsAsFactors = FALSE) %>%
    tbl_df()

  sig.inf <- data.frame(
    replicate = 1,
    Age = PFM$timepoints,
    Stage = "proxy.sig.inf",
    value = PFM$proxy.sig.inf,
    stringsAsFactors = FALSE) %>%
    tbl_df()

  clim <- data.frame(
    replicate = 1,
    Age = PFM$timepoints,
    Stage = "clim.signal.ann",
    value = PFM$clim.signal.ann,
    stringsAsFactors = FALSE) %>%
    tbl_df()

  clim2 <- data.frame(
    replicate = 1,
    Age = PFM$timepoints,
    Stage = "clim.signal.timepoints.100",
    value = PFM$clim.signal.timepoints.100,
    stringsAsFactors = FALSE) %>%
    tbl_df()

  clim3 <- data.frame(
    replicate = 1,
    Age = PFM$timepoints.100,
    Stage = "clim.signal.100",
    value = PFM$clim.signal.100,
    stringsAsFactors = FALSE) %>%
    tbl_df()

  rtn <- bind_rows(df, biot.inf,  sig.inf, clim, clim2, clim3)

  return(rtn)
}
