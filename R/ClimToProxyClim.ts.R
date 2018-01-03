#' Simulate sediment archived proxy records from an "assumed true" climate signal.
#'
#' @md
#' @inheritParams ClimToProxyClim
#' @param layer.width the width of the sediment layer from which samples were taken
#' e.g. foraminifera were picked or alkenones were extracted, in cm.
#' @inherit ClimToProxyClim return
#' @inherit ClimToProxyClim description
#' @importFrom dplyr tbl_df
#' @export
#'
#' @examples
#' library(ggplot2)
#' set.seed(26052017)
#' clim.in <- N41.t21k.climate[(nrow(N41.t21k.climate)-40):1, ] - 273.15
#' clim.in.ts <- ts(N41.t21k.climate[(nrow(N41.t21k.climate)):1, ] - 273.15, start = -39)
#'
#' tpts <- c(1, round(N41.proxy$Published.age), 21*1000)
#' #'tpts <- 10000
#' n.samps <- round(runif(length(tpts), 10, 40))
#' sed.acc.rates <- c(25, N41.proxy$Sed.acc.rate.cm.ka, 25)
#' #'sed.acc.rates <- 505+pi
#' system.time({
#'   set.seed(1)
#'   PFM <- ClimToProxyClim(clim.signal = clim.in,
#'                          timepoints = tpts,
#'                          proxy.calibration.type = "identity",
#'                          seas.prod = N41.G.ruber.seasonality,
#'                          sed.acc.rate = sed.acc.rates,
#'                          meas.noise = 0.46, n.samples = n.samps,
#'                          smoothed.signal.res = 100, meas.bias = 1,
#'                          n.replicates = 10)
#' })
#'
#' system.time({
#'   set.seed(1)
#'   PFM.ts <- ClimToProxyClim.ts(clim.signal = clim.in.ts,
#'                                timepoints = tpts,
#'                                proxy.calibration.type = "identity",
#'                                seas.prod = N41.G.ruber.seasonality,
#'                                sed.acc.rate = sed.acc.rates,
#'                                layer.width = 0,
#'                                meas.noise = 0.46, n.samples = n.samps,
#'                                smoothed.signal.res = 100, meas.bias = 1,
#'                                n.replicates = 10)
#' })
#'
#' PlotPFMs(PFM$everything, max.replicates = 1, stage.order = "seq") +
#'   facet_wrap(~stage)
#'
#' PlotPFMs(PFM.ts$everything, max.replicates = 1, stage.order = "seq") +
#'   facet_wrap(~stage)
#'
#' table(PFM$simulated.proxy$simulated.proxy - PFM.ts$simulated.proxy$simulated.proxy)
#' table(PFM$simulated.proxy$proxy.bt - PFM.ts$simulated.proxy$proxy.bt)
#'
#' plot(PFM$simulated.proxy$simulated.proxy, PFM.ts$simulated.proxy$simulated.proxy)
#' abline(0,1)
#'
#' sub <- PFM.ts$smoothed.signal %>%
#'   filter(timepoints > 0)
#'
#' plot(PFM$smoothed.signal$value, sub$value)
#' table(PFM$smoothed.signal$value - sub$value)
#'
#' plot(PFM$simulated.proxy$clim.timepoints.ssr, PFM.ts$simulated.proxy$clim.timepoints.ssr)
#' table((PFM$simulated.proxy$clim.timepoints.ssr - PFM.ts$simulated.proxy$clim.timepoints.ssr) == 0)
ClimToProxyClim.ts <- function(clim.signal,
                            timepoints,
                            proxy.calibration.type = c("identity", "UK37", "MgCa"),
                            smoothed.signal.res = 100,
                            seas.prod = rep(1, 12),
                            bio.depth = 10,
                            sed.acc.rate = 50,
                            layer.width = 1,
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
  stopifnot(is.ts(clim.signal))
  stopifnot(length(sed.acc.rate) == n.timepoints |
              length(sed.acc.rate) == 1)

  if (is.matrix(seas.prod))
    stop("Matrix form of seasonality not yet supported")

  # check no production weights match dimensions of climate
  #print(paste0("seas.prod = ", seas.prod))
  stopifnot(ncol(clim.signal) == length(seas.prod))

  # Ensure seasonal productivities are weights
  seas.prod <- seas.prod / sum(seas.prod)

  # Calculate timepoint invariant values ------
  max.clim.signal.i <- end(clim.signal)[1]
  min.clim.signal.i <- start(clim.signal)[1]

  # Rescale sed.acc.rate to per year
  sed.acc.rate <- sed.acc.rate / 1000

  if (length(sed.acc.rate) == 1) {
    sed.acc.rate <- rep(sed.acc.rate, n.timepoints)
  }

  # Replicate n.samples if not vector
  if (length(n.samples) == 1) {
    n.samples <- rep(n.samples, n.timepoints)
  }

  # Check whether bioturbation window will extend beyond climate signal for any of the timepoints

  # bioturbation window will be focal.timepoint - bio.depth.timesteps - layer.width.years / 2 to
  # focal.timepoint + 3*bio.depth.timesteps

  max.min.windows <- matrix(t(sapply(1:length(timepoints), function(tp) {
    bio.depth.timesteps <- round(bio.depth / sed.acc.rate[tp])
    layer.width.years <- ceiling(layer.width / sed.acc.rate[tp])
    return(c(max = timepoints[tp] + 3 * bio.depth.timesteps,
             min = timepoints[tp] - bio.depth.timesteps - layer.width.years / 2))
  })), ncol = 2)

  colnames(max.min.windows) <- c("max", "min")

  #print(max.min.windows)
  max.ind <- max.min.windows[,"max"] >= max.clim.signal.i
  min.ind <- max.min.windows[,"min"] <  min.clim.signal.i

  if (any(max.ind))
    warning(paste0("One or more requested timepoints is too old. Bioturbation window(s) for timepoint(s) ",
                   paste(timepoints[max.ind], collapse = ", "),
                   " extend(s) beyond end of input climate signal. Returning pseudo-proxy for valid timepoints."))

  if (any(max.min.windows[,"min"] < min.clim.signal.i))
    warning(paste0("One or more requested timepoints is too recent. Bioturbation window(s) for timepoint(s) ",
                   timepoints[max.min.windows[, "min"] < min.clim.signal.i],
                   " extend(s) above start of input climate signal. Returning pseudo-proxy for valid timepoints."))


  timepoints <- timepoints[max.ind == FALSE & min.ind == FALSE]
  n.timepoints <- length(timepoints)

  # Remove timepoints that exceed clim.signal ------
  sed.acc.rate <- sed.acc.rate[max.ind == FALSE & min.ind == FALSE]
  n.samples <- n.samples[max.ind == FALSE & min.ind == FALSE]
  max.min.windows <- max.min.windows[max.ind == FALSE & min.ind == FALSE, , drop = FALSE]


  # Convert to proxy units if requested --------
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
    timepoints.smoothed <- seq(min.clim.signal.i, max.clim.signal.i, by = smoothed.signal.res)
    clim.signal.smoothed <- ChunkMatrix(timepoints.smoothed, smoothed.signal.res,
                                        proxy.clim.signal)
  }


  # For each timepoint ------
  out <- sapply(1:n.timepoints, function(tp) {
    # Get bioturbation window ----------
    first.tp <- max.min.windows[tp, "min"]
    last.tp <- max.min.windows[tp, "max"]
    bioturb.window <- first.tp:last.tp

    # Get bioturbation weights --------
    bioturb.weights <- BioturbationWeights(z = bioturb.window, focal.depth = timepoints[tp],
                                           layer.width = layer.width, sed.acc.rate = sed.acc.rate[tp],
                                           mix.depth = bio.depth)

    bioturb.weights <- bioturb.weights / sum(bioturb.weights)

    # Get portion of clim.signal corresponding to bioturbation window -------
    ## window is slow
    # clim.sig.window <- stats:::window.ts(proxy.clim.signal,
    #                           start = first.tp,
    #                           end = last.tp)
    clim.sig.window <-  clim.signal[first.tp:last.tp - min.clim.signal.i+1, ]

    # this is estimating mean deviation MD, (not MAD or SD)
    # no need to estimate this from the psuedo data
    # MD = 2/(exp(1)/std) for exponential, where std = lambda = bio.depth.timesteps
    smoothing.width = sum(bioturb.weights*abs(bioturb.window))

    # Get bioturbation X no-seasonality weights matrix ---------
    biot.sig.weights <- bioturb.weights %o% rep(1, ncol(clim.signal))
    biot.sig.weights <- biot.sig.weights / sum(biot.sig.weights)

    # Get bioturbation X seasonality weights matrix ---------
    clim.sig.weights <- bioturb.weights %o% seas.prod
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
      clim.sig.window.ann <- rowSums(clim.sig.window %*% diag(seas.prod))

      # weights passed as a matrix are applied columnwise, so
      # modulo on nrows is need here
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
    bias <- stats::rnorm(n = n.replicates, mean = 0, sd = meas.bias)
  } else{
    bias <- rep(0, n.replicates)
  }
  if (meas.noise != 0) {
    noise <- stats::rnorm(n = n.replicates * n.timepoints, mean = 0, sd = meas.noise)
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

  smoothed.signal$stage <- "clim.signal.smoothed"

  everything <- MakePFMDataframe(out)

  return(list(simulated.proxy=simulated.proxy,
              smoothed.signal=smoothed.signal,
              everything = everything))
  #return(everything)
}
