#' Simulate sediment archived proxy records from an "assumed true" climate signal.
#'
#' @md
#' @description \code{ClimToProxyClim} simulates the creation of a proxy climate record
#'   from a climate signal that is assumed to be true.
#'
#'   The following aspects of proxy creation are currently modelled.
#'
#'   1. Seasonal bias in the encoding of a proxy due to the interaction between
#'   climate seasonality and any seasonality in the life cycle of the organism
#'   encoding the climate signal (e.g. Foraminifera for Mg/Ca ratios, or
#'   phytoplankton for Alkenone unsaturation indices).
#'
#'   2. Bioturbation of the sediment archived proxy. For each requested
#'   timepoint, the simulated proxy consists of a weighted mean of the climate
#'   signal over a time window that is determined by the sediment accumulation
#'   rate \{sed.acc.rate} and the bioturbation depth \{bio.depth} which defaults
#'   to 0.1 m. The weights are given by the depth solution to an impulse
#'   response function (Berger and Heath, 1968).
#'
#'   3. Aliasing of seasonal and inter-annual climate variation onto to
#'   bioturbated (smoothed) signal.
#'   For proxies measured on a small number of discrete particles both seasonal
#'   and inter-annual climate variation is aliased into the proxy record. For
#'   example, Foraminifera have a life-cycle of approximately 1 month, so they
#'   record something like the mean temperature from a single month. If Mg/Ca is
#'   measured on e.g. \code{n.samples} = 30 individuals, the measured proxy
#'   signal is a mean of 30 distinct monthly mean temperatures and will thus be
#'   a stochastic sample of the true mean climate.
#'
#'   4. Measurement noise/error is added as a pure Gaussian white noise
#'   process with mean = 0, standard deviation = \code{meas.noise}.
#'
#'   5. Additionally, a random *bias* can be added to each realisation of a
#'   proxy record. Bias is simulated as a Gaussian random variable with mean =
#'   0, standard deviation = \code{meas.bias}. The same randomly generated bias
#'   value is applied to all timepoints in a simulated proxy record, when
#'   multiple replicate proxies are generated (\{n.replicates} > 1) each
#'   replicate has a different bias applied.
#'
#'   \code{ClimToProxyClim} returns one or more replicates of the final simulated proxy
#'   as well as several intermediate stages (see section **Value** below).
#'
#'
#' @param clim.signal The "assumed true" climate signal, e.g. climate model output or
#'   instrumental record. A years x 12 (months) matrix of temperatures.
#' @param timepoints The timepoints for which the proxy record is to be modelled
#' @param smoothed.signal.res The resolution, in years, of the smoothed
#' (block averaged) version of the input climate signal returned for plotting.
#' This does not affect what the proxy model uses as input.
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
#' @param n.samples Number of e.g. Foraminifera sampled per timepoint
#' @param n.replicates Number of replicate proxy time-series to simulate from
#'   the climate signal
#'
#' @return \code{ClimToProxyClim} returns a list with three elements:
#'
#'  1. a dataframe \code{simulated.proxy}
#'  2. a dataframe \code{smoothed.signal}
#'  3. a list \code{everything}
#'
#'
#' The dataframe \code{simulated.proxy} contains a single realisation of the
#' final forward modelled proxy, as well as the intermediate stages and the
#' original climate signal at the requested timepoints.
#'
#' The dataframe \code{smoothed.signal} contains a block averaged version the
#' input climate signal, defaults to 100 year means but this is set by the parameter
#' smoothed.signal.res. This is useful for plotting against the resulting simulated
#' proxy.
#'
#' The list \code{everything} contains all of the above, but where a stage contains stochastically
#' generated noise, rather than a vector, a \code{timepoints} **by** \code{n.replicates} matrix is returned.
#'
#' **Named elements of the returned proxy record:**
#'
#' \tabular{ll}{
#' \bold{Variable} \tab \bold{Description} \cr
#' timepoints                 \tab requested timepoints                                                                                                               \cr
#' clim.timepoints.1000 \tab 1000 year means of climate signal evaluated at the requested timepoints                                                             \cr
#' clim.timepoints.100 \tab 100 year means of climate signal evaluated at the requested timepoints                                                             \cr
#' clim.timepoints.50 \tab 50 year means of climate signal evaluated at the requested timepoints                                                             \cr
#' proxy.bt               \tab proxy after bioturbation                                                                                                           \cr
#' proxy.bt.sb              \tab proxy after bioturbation + seasonal bias                                                                                           \cr
#' sed.acc.rate               \tab sediment accumulation rates for each timepoint                                                                                     \cr
#' smoothing.width            \tab weighted mean time span represented in a sample after bioturbation                                                                 \cr
#' proxy.bt.sb.sampY             \tab proxy after bioturbation, seasonal bias and finite sampling of years but not seasonality                                                                        \cr
#' proxy.bt.sb.sampYM             \tab proxy after bioturbation, seasonal bias and finite sampling of years and months                                                                        \cr
#' proxy.bt.sb.inf.b            \tab proxy after bioturbation, seasonal bias and calibration bias                                                                       \cr
#' proxy.bt.sb.inf.b.n          \tab proxy after bioturbation, seasonal bias, calibration bias and measurement noise                                                    \cr
#' proxy.bt.sb.sampYM.b           \tab proxy after bioturbation, seasonal bias, finite sampling and calibration bias                                                      \cr
#' proxy.bt.sb.sampYM.b.n         \tab proxy after bioturbation, seasonal bias, finite sampling, calibration bias and measurement noise                                   \cr
#' simulated.proxy            \tab final simulated proxy, this will be same as proxy.bt.sb.inf.b.n when n.samples = Inf, and proxy.bt.sb.sampYM.b.n when n.samples is finite
#'}
#'
#' @importFrom dplyr tbl_df
#' @export
#'
#' @examples
ClimToProxyClim <- function(clim.signal,
                            timepoints,
                            smoothed.signal.res = 100,
                            seas.prod = rep(1, 12),
                            bio.depth = 0.1,
                            sed.acc.rate = 5e-04,
                            meas.noise = 0,
                            meas.bias = 0,
                            n.samples = Inf,
                            n.replicates = 1) {
  # Check inputs --------
  n.timepoints <- length(timepoints)

  stopifnot(is.matrix(clim.signal))
  stopifnot(length(sed.acc.rate) == n.timepoints |
              length(sed.acc.rate) == 1)
  if (is.matrix(seas.prod))
    stop("Matrix form of seasonality not yet supported")

  # Ensure seasonal productivities are weights
  seas.prod <- seas.prod / sum(seas.prod)

  # Calculate timepoint invariant values ------
  max.clim.signal.i <- nrow(clim.signal)
  sig.years.i <- 1:max.clim.signal.i



  if (length(sed.acc.rate) == 1) {
    sed.acc.rate <- rep(sed.acc.rate, n.timepoints)
  }

  # Create smoothed climate signal
  timepoints.smoothed <- seq(1, max.clim.signal.i, by = smoothed.signal.res)
  clim.signal.smoothed <- ChunkMatrix(timepoints.smoothed, smoothed.signal.res, clim.signal)

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

    if (max(sig.window.i.1) > max.clim.signal.i) {
      warning("Bioturbation window extends below end of clim.signal")
    }

    sig.window.i <-
      sig.window.i.1[sig.window.i.1 > 0 &
                       sig.window.i.1 < max.clim.signal.i]

    stopifnot(sig.window.i > 0)
    stopifnot(max.clim.signal.i > max(sig.window.i))

    clim.sig.window <- clim.signal[sig.window.i, , drop = FALSE]


    # Get bioturbation X no-seasonality weights matrix ---------
    biot.sig.weights <- bioturb.weights %o% rep(1, 12)
    biot.sig.weights <-
      biot.sig.weights[sig.window.i.1 > 0 &
                         sig.window.i.1 < max.clim.signal.i, , drop = FALSE]
    biot.sig.weights <- biot.sig.weights / sum(biot.sig.weights)


    # Get bioturbation X seasonality weights matrix ---------
    clim.sig.weights <- bioturb.weights %o% seas.prod
    clim.sig.weights <-
      clim.sig.weights[sig.window.i.1 > 0 &
                         sig.window.i.1 < max.clim.signal.i, , drop = FALSE]
    clim.sig.weights <- clim.sig.weights / sum(clim.sig.weights)

    # Check weights sum to 1, within tolerance
    weight.err <- abs(sum(clim.sig.weights) - 1)
    if ((weight.err < 1e-10) == FALSE) stop(paste0("weight.err = ", weight.err))


    # Calculate mean clim.signal -------
    proxy.bt <- sum(biot.sig.weights * clim.sig.window)

    proxy.bt.sb <- sum(clim.sig.weights * clim.sig.window)

    if (is.infinite(n.samples)) {
      proxy.bt.sb.sampY <- rep(NA, n.replicates)
      proxy.bt.sb.sampYM <- rep(NA, n.replicates)
    } else if (is.finite(n.samples)) {
      # call sample once for all replicates together, then take means of
      # groups of n.samples
      # Get indices not values
      samp.indices <-  sample(length(clim.sig.window),
                              n.samples * n.replicates,
                              prob = clim.sig.weights,
                              replace = TRUE)

      # convert vector to matrix (cheap only attributes changed), then means
      # can be taken across columns to get per replicate means
      samp <- matrix(clim.sig.window[samp.indices], nrow = n.samples)
      #proxy.bt.sb.sampYM <- apply(samp, 2, mean)
      proxy.bt.sb.sampYM <- colMeans(samp)

      # Get without seasonal aliasing (bioturbation aliasing only)
      clim.sig.window.ann <- rowSums(clim.sig.window %*% diag(seas.prod))
      row.indices <- (samp.indices-1) %% nrow(clim.sig.window) + 1

      samp.bt <- matrix(clim.sig.window.ann[row.indices], nrow = n.samples)
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
  # use plyr::alply to always return a list
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

  # Add bias and noise to finite sample --------
  out$proxy.bt.sb.sampYM.b <- out$proxy.bt.sb.sampYM + bias
  out$proxy.bt.sb.sampYM.b.n <- out$proxy.bt.sb.sampYM.b + noise

  # Calculate chunked climate at timepoints
  # get 100 year clim.average at timepoints -------

  out$clim.timepoints.1000 <- ChunkMatrix(timepoints, 1000, clim.signal)
  out$clim.timepoints.100 <- ChunkMatrix(timepoints, 100, clim.signal)
  out$clim.timepoints.50 <- ChunkMatrix(timepoints, 50, clim.signal)



  # Add items to output list -----------
  out$timepoints = timepoints
  out$clim.signal.ann = rowSums(clim.signal[timepoints,  , drop = FALSE]) / ncol(clim.signal)
  out$sed.acc.rate = sed.acc.rate
  out$timepoints.smoothed = timepoints.smoothed
  out$clim.signal.smoothed = clim.signal.smoothed



  # Organise output -------
  simulated.proxy <-
    dplyr::tbl_df(out[c(
      "timepoints",
      "clim.timepoints.1000",
      "clim.timepoints.100",
      "clim.timepoints.50",
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

  if (is.finite(n.samples)) {simulated.proxy$simulated.proxy <- simulated.proxy$proxy.bt.sb.sampYM.b.n}else{
    simulated.proxy$simulated.proxy <- simulated.proxy$proxy.bt.sb.inf.b.n
  }


  smoothed.signal <- dplyr::tbl_df(out[c(
    "timepoints.smoothed",
    "clim.signal.smoothed"
    )])

  smoothed.signal <- dplyr::rename(smoothed.signal,
                            timepoints = timepoints.smoothed,
                            value = clim.signal.smoothed)

  smoothed.signal$Stage <- "clim.signal.smoothed"

  return(list(simulated.proxy=simulated.proxy,
              smoothed.signal=smoothed.signal,
              everything = out))
}



ChunkMatrix <- function(timepoints, width, climate.matrix){
  max.clim.signal.i <- nrow(climate.matrix)

  rel.wind <- 1:width -round(width/2)

  sapply(timepoints, function(tp){

    avg.window.i.1 <- (rel.wind) + tp

    if (max(avg.window.i.1) > max.clim.signal.i) {
      warning("Window extends below end of clim.signal")
    }

    avg.window.i <- avg.window.i.1[avg.window.i.1 > 0 &
                                     avg.window.i.1 < max.clim.signal.i]

    stopifnot(avg.window.i > 0)
    stopifnot(max.clim.signal.i > max(avg.window.i))

    mean(climate.matrix[avg.window.i, , drop = FALSE])
  })
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
    proxy.bt.sb.sampYM = as.vector(PFM$proxy.bt.sb.sampYM),
    proxy.bt.sb.inf.b = as.vector(PFM$proxy.bt.sb.inf.b),
    proxy.bt.sb.sampYM.b = as.vector(PFM$proxy.bt.sb.sampYM.b),
    proxy.bt.sb.inf.b.n = as.vector(PFM$proxy.bt.sb.inf.b.n),
    proxy.bt.sb.sampYM.b.n = as.vector(PFM$proxy.bt.sb.sampYM.b.n),
    stringsAsFactors = FALSE)

  df$Age <- PFM$timepoints
  df$replicate <- rep(1:ncol(PFM$proxy.bt.sb.inf.b), each = length(PFM$timepoints))
  df <- tbl_df(df) %>%
    dplyr::gather(Stage, value, -Age, -replicate)

  #df$proxy.bt.sb = as.vector(PFM$proxy.bt.sb)

  bt.inf <- data.frame(
    replicate = 1,
    Age = PFM$timepoints,
    Stage = "proxy.bt",
    value = PFM$proxy.bt,
    stringsAsFactors = FALSE) %>%
    tbl_df()

  sig.inf <- data.frame(
    replicate = 1,
    Age = PFM$timepoints,
    Stage = "proxy.bt.sb",
    value = PFM$proxy.bt.sb,
    stringsAsFactors = FALSE) %>%
    tbl_df()

  clim <- data.frame(
    replicate = 1,
    Age = PFM$timepoints,
    Stage = "clim.signal.ann",
    value = PFM$clim.signal.ann,
    stringsAsFactors = FALSE) %>%
    tbl_df()

  clim2a <- data.frame(
    replicate = 1,
    Age = PFM$timepoints,
    Stage = "clim.timepoints.1000",
    value = PFM$clim.timepoints.1000,
    stringsAsFactors = FALSE) %>%
    tbl_df()

  clim2 <- data.frame(
    replicate = 1,
    Age = PFM$timepoints,
    Stage = "clim.timepoints.100",
    value = PFM$clim.timepoints.100,
    stringsAsFactors = FALSE) %>%
    tbl_df()

  clim2b <- data.frame(
    replicate = 1,
    Age = PFM$timepoints,
    Stage = "clim.timepoints.50",
    value = PFM$clim.timepoints.50,
    stringsAsFactors = FALSE) %>%
    tbl_df()

  clim3 <- data.frame(
    replicate = 1,
    Age = PFM$timepoints.smoothed,
    Stage = "clim.signal.smoothed",
    value = PFM$clim.signal.smoothed,
    stringsAsFactors = FALSE) %>%
    tbl_df()

  rtn <- dplyr::bind_rows(df, bt.inf,  sig.inf, clim, clim2a, clim2, clim2b, clim3)

  return(rtn)
}
