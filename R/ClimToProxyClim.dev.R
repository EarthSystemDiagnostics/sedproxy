#' Simulate sediment archived proxy records from an "assumed true" climate signal.
#'
#' @md
#' @inheritParams ClimToProxyClim
#' @param seas.prod Either the seasonal pattern of productivity for the organism(s)
#'   recording / producing the proxy as a vector of 12 values, or a function that
#'   produces an index of productivity as a function of temperature.
#'   Defaults to a uniform seasonal distribution.
#' @param seas.prod.args A named list of parameters to be passed to a function named in seas.prod
#' @inherit ClimToProxyClim return
#' @inherit ClimToProxyClim description
#' @importFrom dplyr tbl_df
#' @export
#'
#' @examples
ClimToProxyClim.dev <- function(clim.signal,
                            timepoints,
                            proxy.calibration.type = c("identity", "UK37", "MgCa"),
                            smoothed.signal.res = 100,
                            seas.prod = rep(1, 12),
                            seas.prod.args = NULL,
                            bio.depth = 0.1,
                            sed.acc.rate = 5e-04,
                            meas.noise = 0,
                            meas.bias = 0,
                            n.samples = Inf,
                            n.replicates = 1) {
  # Check inputs --------
  n.timepoints <- length(timepoints)

  if((length(n.samples) == 1 | length(n.samples)==n.timepoints)==FALSE)
    stop("n.sample must be either a single value, or a vector the same
         length as timepoints")

  if (all(is.finite(n.samples))==FALSE & all(is.infinite(n.samples))==FALSE)
    stop("n.samples cannot be a mix of finite and infinite")

  stopifnot(is.matrix(clim.signal))
  stopifnot(length(sed.acc.rate) == n.timepoints |
              length(sed.acc.rate) == 1)

  if (any(is.function(seas.prod), (is.numeric(seas.prod)&length(seas.prod) == ncol(clim.signal))) == FALSE)
    stop("seas.prod must be either a vector of weights with length = ncol(clim.signal), or a function.
         Function names should be given unquoted, e.g. dnorm, not \"dnorm\"")

  # Calculate timepoint invariant values ------
  max.clim.signal.i <- nrow(clim.signal)
  sig.years.i <- 1:max.clim.signal.i

  if (length(sed.acc.rate) == 1) {
    sed.acc.rate <- rep(sed.acc.rate, n.timepoints)
  }

  # Replicate n.samples if not vector
  if (length(n.samples) == 1) {
    n.samples <- rep(n.samples, n.timepoints)
  }

  # Check whether bioturbation window will extend beyond climate signal for any of the timepoints

  max.min.windows <- t(sapply(1:length(timepoints), function(tp){
    bio.depth.timesteps <- round(bio.depth / sed.acc.rate[tp])
    if (bio.depth.timesteps == 0){
      bio.depth.timesteps <- 0
      bioturb.window <- 1
    }else{
      bioturb.window <- (-1*bio.depth.timesteps):(3*bio.depth.timesteps)
    }
    #print(range(bioturb.window))
    return(c(max = max(bioturb.window + timepoints[tp]),
             min = min(bioturb.window + timepoints[tp])))
  }))

  #print(max.min.windows)
  max.ind <- max.min.windows[,"max"] >= max.clim.signal.i
  min.ind <- max.min.windows[,"min"] <  1

  if (any(max.ind))
    warning(paste0("One or more requested timepoints is too old. Bioturbation window(s) for timepoint(s) ",
                   paste(timepoints[max.ind], collapse = ", "),
                   " extend(s) beyond end of input climate signal. Returning pseudo-proxy for valid timepoints."))

  if (any(max.min.windows[,"min"] < 1))
    warning(paste0("One or more requested timepoints is too recent. Bioturbation window(s) for timepoint(s) ",
                timepoints[max.min.windows[, "min"] < 1],
                " extend(s) above start of input climate signal. Returning pseudo-proxy for valid timepoints."))


  timepoints <- timepoints[max.ind == FALSE & min.ind == FALSE]
  n.timepoints <- length(timepoints)

  # Trim timepoint invariant values ------
  sed.acc.rate <- sed.acc.rate[max.ind == FALSE & min.ind == FALSE]
  n.samples <- n.samples[max.ind == FALSE & min.ind == FALSE]



  # Generate productivity weights from function if supplied
  if (is.function(seas.prod)){
    FUN <- match.fun(seas.prod)
    seas.prod.weights <- do.call(FUN, args = c(list(x = clim.signal), seas.prod.args))
    seas.prod.weights <- seas.prod.weights / sum(seas.prod.weights)
  }

  # Ensure seasonal productivities are weights and matrix
  if (is.numeric(seas.prod)){
  seas.prod.weights <- seas.prod / sum(seas.prod)
  seas.prod.weights <- matrix(rep(seas.prod.weights, nrow(clim.signal)),
                              nrow = nrow(clim.signal), byrow = TRUE)
  }


  proxy.calibration.type <- match.arg(proxy.calibration.type)

  if (proxy.calibration.type != "identity") {
    mean.temperature <-  mean(as.vector(clim.signal))
    proxy.clim.signal <-
      matrix(
        ProxyConversion(
          temperature = as.vector(clim.signal),
          proxy.calibration.type = proxy.calibration.type,
          point.or.sample = "point",
          n = 1
        )[, 1],
        ncol = ncol(clim.signal),
        byrow = FALSE
      )
    meas.noise <- as.vector(ProxyConversion(temperature = mean.temperature + meas.noise,
                                            proxy.calibration.type = proxy.calibration.type) -
                              ProxyConversion(temperature = mean.temperature,
                                              proxy.calibration.type = proxy.calibration.type))
  } else{
    proxy.clim.signal <- clim.signal
  }

  # Create smoothed climate signal --------
  if (is.na(smoothed.signal.res)) {
    timepoints.smoothed <- NA
    clim.signal.smoothed <- NA
  } else{
    timepoints.smoothed <- seq(1, max.clim.signal.i, by = smoothed.signal.res)
    clim.signal.smoothed <- ChunkMatrix(timepoints.smoothed, smoothed.signal.res,
                                        proxy.clim.signal)
  }

  # For each timepoint ------
  out <- sapply(1:n.timepoints, function(tp) {
    # Get bioturbation window ----------
    bio.depth.timesteps <- round(bio.depth / sed.acc.rate[tp])
    if (bio.depth.timesteps == 0){
      bio.depth.timesteps <- 0
      bioturb.window <- 1
    }else{
      bioturb.window <- (-1*bio.depth.timesteps):(3*bio.depth.timesteps)
      #print(paste0("loop = ", range(bioturb.window)))
    }

    # Get bioturbation weights --------
    bioturb.weights <-
      ImpulseResponse(-bioturb.window, bio.depth.timesteps, z0 = 0)

    # Check depth and time order match
    # plot(bioturb.weights, (bioturb.window), type = "l", ylim = rev(range(bioturb.window)))

    # Get portion of clim.signal corresponding to bioturbation window -------


    # Correction 2017.09.29
    # Do not shift by Tau: bioturbation does not cause a time shift (assuming constant
    # sedimentation rate) once out of the bioturbated layer!

    # [my speculative reasoning AMD]
    # This will be generally true for any signal with a constant first derivative,
    # like time if sedimentation rate is constant, or on average for stationary white noise.

    ## shift by bio.depth.timesteps (tau in Torben's notation)
    ## to remove timeshift due to bioturbation, which would effect dating in the same way
    sig.window.i.1 <- bioturb.window + timepoints[tp] #+ bio.depth.timesteps


    if (max(sig.window.i.1) >= max.clim.signal.i) {
      warning("Bioturbation window extends below end of clim.signal")
    }

    valid.window.logical <- sig.window.i.1 > 0 &
      sig.window.i.1 <= max.clim.signal.i

    bioturb.weights <- bioturb.weights[valid.window.logical]

    sig.window.i <-
      sig.window.i.1[valid.window.logical]


    stopifnot(sig.window.i > 0)
    stopifnot(max.clim.signal.i >= max(sig.window.i))

    clim.sig.window <- proxy.clim.signal[sig.window.i, , drop = FALSE]

    # this is estimating mean deviation MD, (not MAD or SD)
    # no need to estimate this from the psuedo data
    # MD = 2/(exp(1)/std) for exponential, where std = lambda = bio.depth.timesteps
    smoothing.width = sum(bioturb.weights*abs(bioturb.window))

    # Get bioturbation X no-seasonality weights matrix ---------
    biot.sig.weights <- bioturb.weights %o% rep(1, ncol(clim.signal))
    biot.sig.weights <- biot.sig.weights / sum(biot.sig.weights)


    # Get bioturbation X seasonality weights matrix ---------
    seas.prod.weights <- seas.prod.weights[sig.window.i.1, , drop = FALSE]
    seas.prod.weights <- seas.prod.weights / sum(seas.prod.weights)
    clim.sig.weights <- bioturb.weights * seas.prod.weights
    clim.sig.weights <- clim.sig.weights / sum(clim.sig.weights)

    # Check weights sum to 1, within tolerance
    weight.err <- abs(sum(clim.sig.weights) - 1)
    if ((weight.err < 1e-10) == FALSE) stop(paste0("weight.err = ", weight.err))


    # Calculate mean clim.signal -------

    # Just bioturbation
    proxy.bt <- sum(biot.sig.weights * clim.sig.window)

    # Bioturbation + seasonal bias
    proxy.bt.sb <- sum(clim.sig.weights * clim.sig.window)

    # Bioturbation + seasonal bias + aliasing
    if (is.infinite(n.samples[tp])) {
      proxy.bt.sb.sampY <- rep(NA, n.replicates)
      proxy.bt.sb.sampYM <- rep(NA, n.replicates)
    } else if (is.finite(n.samples[tp])) {
      # call sample once for all replicates together, then take means of
      # groups of n.samples
      # Get indices not values
      samp.indices <-  sample(length(clim.sig.window),
                              n.samples[tp] * n.replicates,
                              prob = clim.sig.weights,
                              replace = TRUE)

      # convert vector to matrix (cheap only attributes changed), then means
      # can be taken across columns to get per replicate means
      samp <- matrix(clim.sig.window[samp.indices], nrow = n.samples[tp])
      #proxy.bt.sb.sampYM <- apply(samp, 2, mean)
      proxy.bt.sb.sampYM <- colMeans(samp)

      # Get without seasonal aliasing (bioturbation aliasing only)

      clim.sig.window.ann <- rowSums(clim.sig.window * seas.prod.weights)
      row.indices <- (samp.indices-1) %% nrow(clim.sig.window) + 1

      samp.bt <- matrix(clim.sig.window.ann[row.indices], nrow = n.samples[tp])
      proxy.bt.sb.sampY <- colMeans(samp.bt)

    }


    # Gather output ----------
    list(
      smoothing.width = smoothing.width,
      proxy.bt = proxy.bt,
      proxy.bt.sb = proxy.bt.sb,
      proxy.bt.sb.sampY = proxy.bt.sb.sampY,
      proxy.bt.sb.sampYM = proxy.bt.sb.sampYM)
  })

  #out <- apply(out, 1, function(x) simplify2array(x))
  out <- plyr::alply(out, 1, function(x) simplify2array(x), .dims = TRUE)

  # remove extra attributes added by alply
  attr(out, "split_type") <- NULL
  attr(out, "split_labels") <- NULL

  #print(out$proxy.bt.sb.sampYM)
  if (n.replicates == 1) out$proxy.bt.sb.sampYM <- matrix(out$proxy.bt.sb.sampYM, nrow = 1)
  out$proxy.bt.sb.sampYM <- t(out$proxy.bt.sb.sampYM)

  if (n.replicates == 1) out$proxy.bt.sb.sampY <- matrix(out$proxy.bt.sb.sampY, nrow = 1)
  out$proxy.bt.sb.sampY <- t(out$proxy.bt.sb.sampY)
  #print(out$proxy.bt.sb.sampYM)

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

  out$proxy.bt.sb.inf.b <- outer(out$proxy.bt.sb, bias, FUN = "+")
  out$proxy.bt.sb.inf.b.n <- out$proxy.bt.sb.inf.b + noise

  if (all(is.finite(n.samples))){
    out$proxy.bt.sb.inf.b[,] <- NA
    out$proxy.bt.sb.inf.b.n[,] <- NA
  }

  # Add bias and noise to finite sample --------
  out$proxy.bt.sb.sampYM.b <- out$proxy.bt.sb.sampYM + bias
  out$proxy.bt.sb.sampYM.b.n <- out$proxy.bt.sb.sampYM.b + noise

  # set intermediate bias stages to NA if no bias modelled
  if (meas.bias == 0) {
    out$proxy.bt.sb.inf.b[,] <- NA
    out$proxy.bt.sb.sampYM.b[,] <- NA
  }

  # Calculate chunked climate at timepoints
  
  # Create smoothed climate signal
  if (is.na(smoothed.signal.res)) {
    out$clim.timepoints.ssr <- NA

  } else{
    out$clim.timepoints.ssr <- ChunkMatrix(timepoints, smoothed.signal.res, proxy.clim.signal)
  }

  # Add items to output list -----------
  out$timepoints = timepoints
  out$clim.signal.ann = rowSums(proxy.clim.signal[timepoints,  , drop = FALSE]) / ncol(proxy.clim.signal)
  out$sed.acc.rate = sed.acc.rate
  out$timepoints.smoothed = timepoints.smoothed
  out$clim.signal.smoothed = clim.signal.smoothed

  # Organise output -------
  simulated.proxy <-
    dplyr::tbl_df(out[c(
      "timepoints",
      "clim.signal.ann",
      "clim.timepoints.ssr",
      "proxy.bt",
      "proxy.bt.sb",
      "sed.acc.rate",
      "smoothing.width"
    )])

  simulated.proxy$proxy.bt.sb.sampY <- out$proxy.bt.sb.sampY[, 1, drop = TRUE]
  simulated.proxy$proxy.bt.sb.sampYM <- out$proxy.bt.sb.sampYM[, 1, drop = TRUE]
  simulated.proxy$proxy.bt.sb.inf.b <- out$proxy.bt.sb.inf.b[, 1, drop = TRUE]
  simulated.proxy$proxy.bt.sb.inf.b.n <- out$proxy.bt.sb.inf.b.n[, 1, drop = TRUE]
  simulated.proxy$proxy.bt.sb.sampYM.b <- out$proxy.bt.sb.sampYM.b[, 1, drop = TRUE]
  simulated.proxy$proxy.bt.sb.sampYM.b.n <- out$proxy.bt.sb.sampYM.b.n[, 1, drop = TRUE]

  if (all(is.finite(n.samples))) {
    simulated.proxy$simulated.proxy <- simulated.proxy$proxy.bt.sb.sampYM.b.n
    out$simulated.proxy <- out$proxy.bt.sb.sampYM.b.n
  } else{
    simulated.proxy$simulated.proxy <- simulated.proxy$proxy.bt.sb.inf.b.n
    out$simulated.proxy <- out$proxy.bt.sb.inf.b.n
  }


  smoothed.signal <- dplyr::tbl_df(out[c(
    "timepoints.smoothed",
    "clim.signal.smoothed"
  )])

  smoothed.signal <- dplyr::rename(smoothed.signal,
                                   timepoints = timepoints.smoothed,
                                   value = clim.signal.smoothed)

  smoothed.signal$Stage <- "clim.signal.smoothed"

  everything <- MakePFMDataframe(out)

  return(list(simulated.proxy=simulated.proxy,
              smoothed.signal=smoothed.signal,
              everything = everything))
  #return(everything)
}

## Some tests

# library(sedproxy)
# library(tidyverse)
#
# clim.in <- N41.t21k.climate[nrow(N41.t21k.climate):1,] - 273.15
# PFM <- ClimToProxyClim(clim.signal = clim.in,
#                        timepoints = c(round(N41.proxy$Published.age), 22030),
#                        proxy.calibration.type = "identity",
#                        sed.acc.rate = c(N41.proxy$Sed.acc.rate.m.yr, 3.047077e-04),
#                        seas.prod = N41.G.ruber.seasonality,
#                        meas.noise = 0.46, n.samples = 30,
#                        n.replicates = 10)
#
# set.seed(20171023)
# PFM.dev <- ClimToProxyClim.dev(clim.signal = clim.in,
#                                timepoints = c(round(N41.proxy$Published.age), 22030),
#                                proxy.calibration.type = "identity",
#                                seas.prod = dnorm,
#                                seas.prod.args = list(mean = 26, sd = 2),
#                                sed.acc.rate = c(N41.proxy$Sed.acc.rate.m.yr, 3.047077e-04),
#                                meas.noise = 0.46, n.samples = 30,
#                                n.replicates = 10)
#
#
# PFM.dev <- ClimToProxyClim.dev(clim.signal = clim.in,
#                                timepoints = c(1000, round(N41.proxy$Published.age)),
#                                proxy.calibration.type = "identity",
#                                seas.prod = dnorm,
#                                seas.prod.args = list(mean = 26, sd = 2),
#                                sed.acc.rate = c(3.047077e-04, N41.proxy$Sed.acc.rate.m.yr),
#                                meas.noise = 0.46, n.samples = 30,
#                                n.replicates = 10)
#
#
# PFM.dev <- ClimToProxyClim.dev(clim.signal = clim.in,
#                                timepoints = c(round(N41.proxy$Published.age)),
#                                proxy.calibration.type = "identity",
#                                seas.prod = dnorm,
#                                seas.prod.args = list(mean = 26, sd = 2),
#                                sed.acc.rate = c(N41.proxy$Sed.acc.rate.m.yr),
#                                meas.noise = 0.46, n.samples = 30,
#                                n.replicates = 10)
#
#
# PFM$everything %>%
#   PlotPFMs(max.replicates = 1)
#
# PFM.dev$everything %>%
#   PlotPFMs(max.replicates = 1)

