#' Title
#'
#' @param clim.signal The "true" climate signal, i.e. model output or
#'   instrumental record. A years x 12 (months) matrix of temperatures.
#' @param timepoints The timepoints for which the proxy record is to be modelled
#' @param seas.prod The seasonal pattern of productivity for the organism(s)
#'   archived in the proxy. Either a vector of 12 values or a matrix of the same
#'   dimensions as clim.signal. Defaults to uniform seasonal distribution.
#' @param bio.depth Depth of the bioturbated layer in metres, defaults to 0.1 m
#' @param acc.rate Sediment accumulation rate in metres per year. Defaults to
#'   5e-04 m per year (0.5 m / kyr). Either a single value, or vector of same
#'   length as "timepoints"
#' @param meas.noise The amount of noise to add to each simulated proxy value.
#'   Defined as the standard deviation of a normal distribution with mean = 0
#' @param meas.bias The amount of bias to add to each simulated proxy time-series.
#'   Each replicate proxy time-series has a constant bias added, drawn from a normal distribution with
#'   mean = 0, sd = meas.bias. Bias defaults to zero.
#' @param n.samples Number of e.g. foraminifera sampled per timepoint
#' @param n.replicates Number of replicate proxy time-series to simulate from
#'   the climate signal
#'
#' @return a list
#' @export
#'
#' @examples
ClimToProxyClim <- function(clim.signal,
                            timepoints,
                            seas.prod = rep(1, 12),
                            bio.depth = 0.1,
                            acc.rate = 5e-04,
                            meas.noise = 0,
                            meas.bias = 0,
                            n.samples = Inf,
                            n.replicates = 1) {
  # Check inputs --------
  stopifnot(is.matrix(clim.signal))
  stopifnot(length(acc.rate) == length(timepoints) |
              length(acc.rate) == 1)
  if (is.matrix(seas.prod))
    stop("Matrix form of seasonality not yet supported")

  # Calculate timepoint invariant values ------
  max.clim.signal.i <- nrow(clim.signal)
  sig.years.i <- 1:max.clim.signal.i

  n.timepoints <- length(timepoints)

  if (length(acc.rate) == 1) {
    acc.rate <- rep(acc.rate, n.timepoints)
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

    clim.100.avg <- mean(clim.signal[avg.window.i, , drop = FALSE])
  })

  # For each timepoint ------
  out <- sapply(1:n.timepoints, function(tp) {
    # Get bioturbation window ----------
    bio.depth.timesteps <- round(bio.depth / acc.rate[tp])
    bioturb.window <- (-1*bio.depth.timesteps):(3*bio.depth.timesteps)

    # Get bioturbation weights --------
    bioturb.weights <-
      ImpulseResponse(-bioturb.window, bio.depth.timesteps, z0 = 0)
    bioturb.weights <- bioturb.weights / sum(bioturb.weights)

    bioturb.window.width = sum(bioturb.weights*abs(bioturb.window))
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
      proxy.sig.samp <- apply(samp, 2, mean)
    }

    # get 100 year clim.average -------
    avg.window.i.1 <- (-50:49) + timepoints[tp]

    if (max(avg.window.i.1) > nrow(clim.signal)) {
      warning("Climate average window extends below end of clim.signal")
    }

    avg.window.i <-
      avg.window.i.1[avg.window.i.1 > 0 &
                       avg.window.i.1 < nrow(clim.signal)]

    stopifnot(avg.window.i > 0)
    stopifnot(nrow(clim.signal) > max(avg.window.i))

    clim.100.avg <- mean(clim.signal[avg.window.i, , drop = FALSE])


    # Gather output ----------
    list(
      bioturb.window.width = bioturb.window.width,
      clim.100.avg = clim.100.avg,
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
  if (is.finite(n.samples)){
    out$proxy.sig.samp.b <- t(t(out$proxy.sig.samp) + bias)
    out$proxy.sig.samp.b.n <- out$proxy.sig.samp.b + noise
  }

  if(is.infinite(n.samples)){
    out$proxy.sig.samp.b <- out$proxy.sig.samp.b.n <- NA
  }

  # Add items to output list -----------
  out$timepoints = timepoints
  out$clim.timepoints = rowSums(clim.signal[timepoints,  , drop = FALSE]) / ncol(clim.signal)
  out$sed.acc.rate = acc.rate
  out$timepoints.100 = timepoints.100
  out$clim.signal.100 = clim.signal.100

  return(out)
}
