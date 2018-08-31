#' @md
#' @title Simulate sediment archived proxy records from an "assumed true"
#'   climate signal.
#' @description \code{ClimToProxyClim} simulates the creation of a proxy climate
#'   record from a climate signal that is assumed to be true.
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
#'   to 10 cm. The weights are given by the depth solution to an impulse
#'   response function (Berger and Heath, 1968).
#'
#'   3. Aliasing of seasonal and inter-annual climate variation onto to
#'   bioturbated (smoothed) signal. For proxies measured on a small number of
#'   discrete particles both seasonal and inter-annual climate variation is
#'   aliased into the proxy record. For example, Foraminifera have a life-cycle
#'   of approximately 1 month, so they record something like the mean
#'   temperature from a single month. If Mg/Ca is measured on e.g.
#'   \code{n.samples} = 30 individuals, the measured proxy signal is a mean of
#'   30 distinct monthly mean temperatures and will thus be a stochastic sample
#'   of the true mean climate.
#'
#'   4. Measurement noise/error is added as a pure Gaussian white noise process
#'   with mean = 0, standard deviation =
#'    \code{sqrt(sigma.measurement^2 + sigma.individual^2/n.samples)}.
#'
#'   5. Additionally, a random *bias* can be added to each realisation of a
#'   proxy record. Bias is simulated as a Gaussian random variable with mean =
#'   0, standard deviation = \code{meas.bias}. The same randomly generated bias
#'   value is applied to all timepoints in a simulated proxy record, when
#'   multiple replicate proxies are generated (\{n.replicates} > 1) each
#'   replicate has a different bias applied.
#'
#'   \code{ClimToProxyClim} returns one or more replicates of the final
#'   simulated proxy as well as several intermediate stages (see section
#'   **Value** below).
#'
#' @param clim.signal The "assumed true" climate signal, e.g. climate model
#'   output or instrumental record. A \code{\link{ts}} object consisting of a
#'   years x 12 (months) x n habitats (e.g. depths) matrix of temperatures. The
#'   time series should be at annual resolution and in reverse, i.e. "most
#'   recent timepoint first" order.
#' @param timepoints The timepoints for which the proxy record is to be modelled
#' @param proxy.calibration.type Type of proxy, e.g. Uk'37 or MgCa, to which the
#'   clim.signal is converted before the archiving and measurement of the proxy
#'   is simulated. Defaults to "identity" which means no conversion takes place.
#' @param noise.type Determines whether additive or multiplicative measurement
#'   noise is added. The appropriate type depends on the units of the proxy.
#'   Defaults to multiplicative for MgCa, additive for Uk'37 and identity (none)
#'   calibration types. Can be overidden with a string, "additive" or
#'   "multiplicative" in the case that pre-converted climate signal and
#'   measurement noise values are used in combination with an "identity"
#'   calibration type.
#' @param smoothed.signal.res The resolution, in years, of the smoothed (block
#'   averaged) version of the input climate signal returned for plotting. This
#'   does not affect what the proxy model uses as input. If set to NA, no
#'   smoothed climate output is generated, this can speed up some simulations.
#' @param habitat.weights Production weights for the proxy / proxy-carrier
#' either as a vector of values with length = ncol(clim.signal), i.e. 1 weight
#' for each month x habitat combination, a matrix of the same dimensions as the
#' input climate signal matrix, or a function that produces an index of
#' productivity as a function of temperature.
#' Defaults to a vector of length = ncol(clim.signal) of equal weights.
#' @param habitat.wt.args A named list of parameter values to be passed to
#' a function named in habitat.weights.
#' @param bio.depth Depth of the bioturbated layer in cm, defaults to 10 cm.
#' @param layer.width the width of the sediment layer from which samples were
#'   taken, e.g. foraminifera were picked or alkenones were extracted, in cm.
#'   Defaults to 1 cm. If bio.depth and layer.width are both set to zero,
#'   each timepoint samples from a single year of the clim.signal, equivalent to
#'   sampling an annually laminated sediment core.
#' @param sed.acc.rate Sediment accumulation rate in cm per 1000 years. Defaults
#'   to 50 cm per ka. Either a single value, or vector of same length as
#'   "timepoints"
#' @param sigma.measurement The standard deviation of the measurement error
#' added to each simulated proxy value.
#' @param sigma.individual The standard deviation of error between individuals
#' (e.g. Forams) not otherwise modelled. This could included "vital effects" or
#' aliasing of depth habitat variation not modelled via a depth resolved input
#' climate signal and habitat weights. sigma.individual is scaled by n.samples
#' before being combined with sigma.measurement.
#' @param n.samples Number of e.g. Foraminifera sampled per timepoint, this can
#'   be either a single number, or a vector of length = timepoints
#' @param meas.bias The amount of bias to add to each simulated proxy
#'   time-series. Each replicate proxy time-series has a constant bias added,
#'   drawn from a normal distribution with mean = 0, sd = meas.bias. Bias
#'   defaults to zero.
#' @param n.replicates Number of replicate proxy time-series to simulate from
#'   the climate signal
#' @inheritParams ProxyConversion
#'
#' @return \code{ClimToProxyClim} returns a list with three elements:
#'
#'   1. a dataframe \code{simulated.proxy}
#'   2. a dataframe \code{smoothed.signal}
#'   3. a dataframe \code{everything}
#'
#'
#'   The dataframe \code{simulated.proxy} contains a single realisation of the
#'   final forward modelled proxy, as well as the intermediate stages and the
#'   original climate signal at the requested timepoints.
#'
#'   The dataframe \code{smoothed.signal} contains a block averaged version the
#'   input climate signal, defaults to 100 year means but this is set by the
#'   parameter smoothed.signal.res. This is useful for plotting against the
#'   resulting simulated proxy.
#'
#'   The dataframe \code{everything} contains all of the above but with multiple
#'   replicates of the pseudo-proxy records if requested. The data are in
#'   "long form", with the column "stage" inidcating the proxy stage or input
#'   climate resolution and column "value" giving the values.
#'
#' **Named elements of the returned proxy record:**
#'
#' \describe{
#'    \item{timepoints}{Requested timepoints}
#'    \item{clim.signal.ann}{Input climate signal at requested timepoints at annual resolution}
#'    \item{clim.signal.smoothed}{Input climate signal at regular time intervals and resolution = smoothed.signal.res}
#'    \item{clim.timepoints.ssr}{Input climate signal at requested timepoints, smoothed to resolution = smoothed.signal.res}
#'    \item{proxy.bt}{Climate signal after bioturbation}
#'    \item{proxy.bt.sb}{Climate signal after bioturbation and habitat bias}
#'    \item{proxy.bt.sb.inf.b}{Climate signal after bioturbation, habitat bias, and calibration bias}
#'    \item{proxy.bt.sb.inf.b.n}{Climate signal after bioturbation, habitat bias, and measurement error}
#'    \item{proxy.bt.sb.sampY}{Climate signal after bioturbation, habitat bias, and aliasing of inter-annual variation}
#'    \item{proxy.bt.sb.sampYM}{Climate signal after bioturbation, habitat bias, and aliasing of inter-annual and intra-annual variation such as monthly temperatures or depth habitats}
#'    \item{proxy.bt.sb.sampYM.b}{Climate signal after bioturbation, habitat bias, and aliasing of inter-annual and intra-annual variation such as monthly temperatures or depth habitats, and calibration bias}
#'    \item{proxy.bt.sb.sampYM.b.n}{Climate signal after bioturbation, habitat bias, aliasing, and measurement error}
#'    \item{simulated.proxy}{Final simulated pseudo-proxy, this will be same as proxy.bt.sb.inf.b.n when n.samples = Inf, and proxy.bt.sb.sampYM.b.n when n.samples is finite}
#'    \item{observed.proxy}{True observed proxy (when supplied)}
#' }
#'
#'@importFrom dplyr tbl_df rename
#'@importFrom plyr alply
#'@export
#'
#' @examples
#' library(ggplot2)
#' set.seed(26052017)
#' clim.in <- ts(N41.t21k.climate[nrow(N41.t21k.climate):1,] - 273.15, start = -39)
#'
#' PFM <- ClimToProxyClim(clim.signal = clim.in,
#'                        timepoints = round(N41.proxy$Published.age),
#'                        proxy.calibration.type = "identity",
#'                        habitat.weights = N41.G.ruber.seasonality,
#'                        sed.acc.rate = N41.proxy$Sed.acc.rate.cm.ka,
#'                        layer.width = 1,
#'                        sigma.measurement = 0.46,
#'                        sigma.individual = 0, n.samples = Inf,
#'                        smoothed.signal.res = 10, meas.bias = 1,
#'                        n.replicates = 10)
#'
#' PlotPFMs(PFM$everything, max.replicates = 1, stage.order = "seq") +
#'   facet_wrap(~stage)
#'
#' PlotPFMs(PFM$everything, max.replicates = 1, stage.order = "var")
#'
#' PlotPFMs(PFM$everything, stage.order = "var", plot.stages = "all")
#'
ClimToProxyClim <- function(clim.signal,
                            timepoints,
                            proxy.calibration.type = c("identity", "UK37", "MgCa"),
                            taxon = NULL,
                            slp.int.means = NULL, slp.int.vcov = NULL,
                            noise.type = switch(proxy.calibration.type,
                                                identity = "additive",
                                                UK37 = "additive",
                                                MgCa = "multiplicative"),
                            smoothed.signal.res = 100,
                            habitat.weights = rep(1/ncol(clim.signal),
                                                     ncol(clim.signal)),
                            habitat.wt.args = NULL,
                            bio.depth = 10,
                            sed.acc.rate = 50,
                            layer.width = 1,
                            sigma.measurement = 0,
                            sigma.individual = 0,
                            meas.bias = 0,
                            scale.noise = switch(proxy.calibration.type,
                                                 identity = FALSE,
                                                 UK37 = TRUE,
                                                 MgCa = TRUE),
                            n.samples = Inf,
                            n.replicates = 1) {
  # Check inputs --------
  n.timepoints <- length(timepoints)

  if((length(n.samples) == 1 | length(n.samples)==n.timepoints)==FALSE)
    stop("n.sample must be either a single value, or a vector the same
         length as timepoints")

  if((length(sigma.measurement) == 1 | length(sigma.measurement)==n.timepoints)==FALSE)
    stop("sigma.measurement must be either a single value, or a vector the same
         length as timepoints")

  if (all(is.finite(n.samples))==FALSE & all(is.infinite(n.samples))==FALSE)
    stop("n.samples cannot be a mix of finite and infinite")

  stopifnot(is.matrix(clim.signal))
  if (is.ts(clim.signal)==FALSE)
    stop("Since version 0.3.1 of sedproxy, ClimToProxyClim requires clim.signal to be a ts object")
  stopifnot(length(sed.acc.rate) == n.timepoints |
              length(sed.acc.rate) == 1)

  if (any(is.function(habitat.weights),
          (is.vector(habitat.weights)&length(habitat.weights) == ncol(clim.signal)),
          (is.matrix(habitat.weights)&dim(habitat.weights) == dim(clim.signal))) == FALSE)
    stop("habitat.weights must be either a vector of weights with length = ncol(clim.signal),
         a matrix of weights of the same dimensions as the input climate signal, or a function.
         Function names should be given unquoted, e.g. dnorm, not \"dnorm\"")



  # Ensure seasonal productivities are weights
  habitat.weights <- habitat.weights / sum(habitat.weights)

  # Calculate timepoint invariant values ------
  max.clim.signal.i <- end(clim.signal)[1]
  min.clim.signal.i <- start(clim.signal)[1]

  if (length(sed.acc.rate) == 1) {
    sed.acc.rate <- rep(sed.acc.rate, n.timepoints)
  }

  # Replicate n.samples if not vector
  if (length(n.samples) == 1) {
    n.samples <- rep(n.samples, n.timepoints)
  }

  if (length(sigma.measurement) == 1) {
    sigma.measurement <- rep(sigma.measurement, n.timepoints)
  }

  if (length(sigma.individual) == 1) {
    sigma.individual <- rep(sigma.individual, n.timepoints)
  }

  # Check whether bioturbation window will extend beyond climate signal for any of the timepoints

  # bioturbation window will be focal.timepoint - bio.depth.timesteps - layer.width.years / 2 to
  # focal.timepoint + 3*bio.depth.timesteps

  max.min.windows <- matrix(t(sapply(1:length(timepoints), function(tp) {
    bio.depth.timesteps <- round(1000 * bio.depth / sed.acc.rate[tp])
    layer.width.years <- ceiling(1000 * layer.width / sed.acc.rate[tp])
    return(c(max = timepoints[tp] + 3 * bio.depth.timesteps,
             min = timepoints[tp] - bio.depth.timesteps - layer.width.years / 2))
  })), ncol = 2)

  colnames(max.min.windows) <- c("max", "min")

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

  sigma.measurement <- sigma.measurement[max.ind == FALSE & min.ind == FALSE]
  sigma.individual <- sigma.individual[max.ind == FALSE & min.ind == FALSE]

  max.min.windows <- max.min.windows[max.ind == FALSE & min.ind == FALSE, , drop = FALSE]


  # Scale sigma.individual by n.samples and create combined error term
  sigma.ind.scl <- ifelse(is.finite(n.samples),
                          sigma.individual / sqrt(n.samples), 0)

  sigma.meas.ind <- sqrt(sigma.measurement^2 + sigma.ind.scl^2)


  # Generate productivity weights from function if supplied
  if (is.function(habitat.weights)){
    FUN <- match.fun(habitat.weights)
    habitat.weights <- do.call(FUN, args = c(list(x = clim.signal), habitat.wt.args))
    habitat.weights <- habitat.weights / sum(habitat.weights)
  }

  # If vector ensure seasonal productivities are weights and matrix
  if (is.vector(habitat.weights)){
    habitat.weights <- habitat.weights / sum(habitat.weights)
    habitat.weights <- matrix(rep(habitat.weights, nrow(clim.signal)),
                                 nrow = nrow(clim.signal), byrow = TRUE)
  }



  # Convert to proxy units if requested --------
  proxy.calibration.type <- match.arg(proxy.calibration.type)

  if (proxy.calibration.type != "identity") {
    proxy.clim.signal <-
        ProxyConversion(
          temperature = clim.signal,
          proxy.calibration.type = proxy.calibration.type,
          taxon = taxon,
          slp.int.means = slp.int.means, slp.int.vcov = slp.int.vcov,
          point.or.sample = "point",
          n = 1
        )
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
                                        clim.signal)
  }


  # For each timepoint ------
  out <- sapply(1:n.timepoints, function(tp) {
    # Get bioturbation window ----------
    first.tp <- max.min.windows[tp, "min"]
    last.tp <- max.min.windows[tp, "max"]
    bioturb.window <- first.tp:last.tp

    # Get bioturbation weights --------
    bioturb.weights <- BioturbationWeights(z = bioturb.window, focal.z = timepoints[tp],
                                           layer.width = layer.width, sed.acc.rate = sed.acc.rate[tp],
                                           bio.depth = bio.depth)

    bioturb.weights <- bioturb.weights / sum(bioturb.weights)

    # Get portion of clim.signal corresponding to bioturbation window -------
    ## window is slow
    # clim.sig.window <- stats:::window.ts(proxy.clim.signal,
    #                           start = first.tp,
    #                           end = last.tp)
    clim.sig.window <-  proxy.clim.signal[first.tp:last.tp - min.clim.signal.i+1, ]

    # this is estimating mean deviation MD, (not MAD or SD)
    # no need to estimate this from the psuedo data
    # MD = 2/(exp(1)/std) for exponential, where std = lambda = bio.depth.timesteps
    # smoothing.width = sum(bioturb.weights*abs(bioturb.window))

    # Get bioturbation X no-seasonality weights matrix ---------
    biot.sig.weights <- bioturb.weights %o% rep(1, ncol(proxy.clim.signal))
    biot.sig.weights <- biot.sig.weights / sum(biot.sig.weights)


    # Get bioturbation X seasonality weights matrix ---------
    habitat.weights <- habitat.weights[first.tp:last.tp - min.clim.signal.i+1, , drop = FALSE]
    habitat.weights <- habitat.weights / sum(habitat.weights)
    clim.sig.weights <- bioturb.weights * habitat.weights
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

      clim.sig.window.ann <- rowSums(clim.sig.window * habitat.weights)
      row.indices <- (samp.indices-1) %% nrow(clim.sig.window) + 1

      samp.bt <- matrix(clim.sig.window.ann[row.indices], nrow = n.samples[tp])
      proxy.bt.sb.sampY <- colMeans(samp.bt)

    }


    # Gather output ----------
    list(
      #smoothing.width = smoothing.width,
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

  if (n.replicates == 1) out$proxy.bt.sb.sampYM <- matrix(out$proxy.bt.sb.sampYM, nrow = 1)
  out$proxy.bt.sb.sampYM <- t(out$proxy.bt.sb.sampYM)

  if (n.replicates == 1) out$proxy.bt.sb.sampY <- matrix(out$proxy.bt.sb.sampY, nrow = 1)
  out$proxy.bt.sb.sampY <- t(out$proxy.bt.sb.sampY)


  # Add bias and noise --------
  # Rescale noise if using a calibration -----
  if (scale.noise != FALSE) {
    # scale.noise will either be TRUE because a non-identity calibration is being used
    # or it will be a string to identify the correct calibration if the input time-series
    # has already been converted.
    message("Rescaling noise")

    # If cal type is identity re-scaling still required
    pct <- if (proxy.calibration.type == "identity") {
      scale.noise
    } else{
      proxy.calibration.type
    }

    # mean temperature in temperature units at each timepoint - use bioturbated signal
    mean.temperature <-  as.vector(ProxyConversion(proxy.value = out$proxy.bt,
                                                   proxy.calibration.type = pct, taxon = taxon,
                                                   slp.int.means = slp.int.means, slp.int.vcov = slp.int.vcov))

     sigma.meas.ind <- ProxyConversion(temperature = mean.temperature + sigma.meas.ind,
                                      proxy.calibration.type = pct, taxon = taxon,
                                      slp.int.means = slp.int.means, slp.int.vcov = slp.int.vcov) -
      ProxyConversion(temperature = mean.temperature,
                      proxy.calibration.type = pct, taxon = taxon,
                      slp.int.means = slp.int.means, slp.int.vcov = slp.int.vcov)

    sigma.meas.ind <- as.vector(sigma.meas.ind)

    if (noise.type == "multiplicative"){
      # noise SD needs to be divided by the mean temperature in proxy units in
      # order to maintain a consistent SD in temperature units.
      sigma.meas.ind <- sigma.meas.ind / out$proxy.bt
    }
  }

  if (noise.type == "additive") {
    noise <- stats::rnorm(n = n.replicates * n.timepoints, mean = 0, sd = sigma.meas.ind)

    if (meas.bias != 0) {
      bias <- stats::rnorm(n = n.replicates, mean = 0, sd = meas.bias)
    } else{
      bias <- rep(0, n.replicates)
    }

    # Add bias and noise to infinite sample --------

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
  }else if (noise.type == "multiplicative"){
    noise <- exp(stats::rnorm(n.replicates * n.timepoints, 0, sigma.meas.ind))

    if (meas.bias != 0) {
      bias <- exp(stats::rnorm(n = n.replicates, mean = 0, sd = meas.bias))
    } else{
      bias <- rep(1, n.replicates)
    }

    # Add bias and noise to infinite sample --------
    out$proxy.bt.sb.inf.b <- outer(out$proxy.bt.sb, bias, FUN = "*")
    out$proxy.bt.sb.inf.b.n <- out$proxy.bt.sb.inf.b * noise

    if (all(is.finite(n.samples))){
      out$proxy.bt.sb.inf.b[,] <- NA
      out$proxy.bt.sb.inf.b.n[,] <- NA
    }

    # Add bias and noise to finite sample --------
    out$proxy.bt.sb.sampYM.b <- out$proxy.bt.sb.sampYM * bias
    out$proxy.bt.sb.sampYM.b.n <- out$proxy.bt.sb.sampYM.b * noise

    # set intermediate bias stages to NA if no bias modelled
    if (meas.bias == 0) {
      out$proxy.bt.sb.inf.b[,] <- NA
      out$proxy.bt.sb.sampYM.b[,] <- NA
    }
  }


  # Create smoothed climate signal -----
  if (is.na(smoothed.signal.res)) {
    out$clim.timepoints.ssr <- NA

  } else{
    out$clim.timepoints.ssr <- ChunkMatrix(timepoints, smoothed.signal.res, clim.signal)
  }

  # Add items to output list -----------
  out$timepoints = timepoints
  out$clim.signal.ann = rowSums(clim.signal[timepoints,  , drop = FALSE]) / ncol(clim.signal)
  #out$sed.acc.rate = sed.acc.rate
  out$timepoints.smoothed = timepoints.smoothed
  out$clim.signal.smoothed = clim.signal.smoothed

  # Organise output -------
  simulated.proxy <-
    dplyr::tbl_df(out[c(
      "timepoints",
      "clim.signal.ann",
      "clim.timepoints.ssr",
      "proxy.bt",
      "proxy.bt.sb"#,
      #"sed.acc.rate",
      #"smoothing.width"
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


  # Add calibration uncertainty -------
  # If n.replicates > 1
  # First convert back to temperature units with fixed parameters
  # Then re-convert to proxy units with random parameters
  if (proxy.calibration.type != "identity"){

    if (n.replicates > 1){
      out$simulated.proxy.cal.err <-
        ProxyConversion(proxy.value = out$simulated.proxy,
                        proxy.calibration.type = proxy.calibration.type,
                        taxon = taxon,
                        point.or.sample = "point", n = 1,
                        slp.int.means = slp.int.means,
                        slp.int.vcov = slp.int.vcov)

      out$simulated.proxy.cal.err <-
        ProxyConversion(temperature = out$simulated.proxy.cal.err,
                        proxy.calibration.type = proxy.calibration.type,
                        taxon = taxon,
                        point.or.sample = "sample", n = n.replicates,
                        slp.int.means = slp.int.means, slp.int.vcov = slp.int.vcov)
    }else{
      out$simulated.proxy.cal.err <- out$simulated.proxy
    }

    # Do this in all cases, not just if n.replicates == 1
    out$reconstructed.climate <-
      ProxyConversion(proxy.value = out$simulated.proxy,
                      proxy.calibration.type = proxy.calibration.type,
                      taxon = taxon,
                      point.or.sample = "point", n = 1,
                      slp.int.means = slp.int.means, slp.int.vcov = slp.int.vcov)

  }else{
    out$simulated.proxy.cal.err <- out$simulated.proxy
    out$reconstructed.climate <- out$simulated.proxy
  }

  simulated.proxy$simulated.proxy.cal.err <- out$simulated.proxy.cal.err[, 1, drop = TRUE]
  simulated.proxy$reconstructed.climate <- out$reconstructed.climate[, 1, drop = TRUE]

  everything <- MakePFMDataframe(out)

  slp.int.means <- if (is.null(slp.int.means)) {
    taxon <- if (proxy.calibration.type == "MgCa" & is.null(taxon)) {
      "10 Foram Taxa"
    } else {taxon}
    switch(proxy.calibration.type,
           MgCa = MgCa.foram.pars[[taxon]][["means"]],
           UK37 = UK37.pars[["mueller.uk37"]][["means"]])
  } else{
    slp.int.means
  }

  calibration.pars <- list(proxy.calibration.type = proxy.calibration.type,
                           taxon = taxon,
                           slp.int.means = slp.int.means)

  attr(simulated.proxy, "calibration.pars") <-  calibration.pars
  attr(everything, "calibration.pars") <-  calibration.pars

  out <- list(simulated.proxy=simulated.proxy,
              smoothed.signal=smoothed.signal,
              everything = everything,
              calibration.pars = calibration.pars)

  class(out) <- "sedproxy.pfm"

  return(out)

}

ChunkMatrix <- function(timepoints, width, climate.matrix){

  if (is.ts(climate.matrix)) {
     sapply(timepoints, function(tp){
      rel.wind <- 1:width -round(width/2)
      inds <- rel.wind + tp - start(climate.matrix)[1] + 1
      inds <- inds[inds > 0 & inds < nrow(climate.matrix)]
      m <- climate.matrix[inds, , drop = FALSE]
      mean(m)
    })}else{
    max.clim.signal.i <- nrow(climate.matrix)

    rel.wind <- 1:width -round(width/2)

    sapply(timepoints, function(tp){

      avg.window.i.1 <- (rel.wind) + tp

      if (max(avg.window.i.1) > max.clim.signal.i) {
        warning("In ChunkMatrix: window extends below end of clim.signal")
      }

      avg.window.i <- avg.window.i.1[avg.window.i.1 > 0 &
                                       avg.window.i.1 < max.clim.signal.i]

      stopifnot(avg.window.i > 0)
      stopifnot(max.clim.signal.i > max(avg.window.i))
      mean(climate.matrix[avg.window.i, , drop = FALSE])
    })}
}


#' Convert "everything" part of output from ClimToProxyClim to dataframe.
#' Used internally.
#'
#' @param PFM output from ClimToProxyClim
#' @return a dataframe
#' @importFrom dplyr bind_rows filter
#' @importFrom tidyr gather
MakePFMDataframe <- function(PFM){
  df <- data.frame(
    proxy.bt.sb.sampY = as.vector(PFM$proxy.bt.sb.sampY),
    proxy.bt.sb.sampYM = as.vector(PFM$proxy.bt.sb.sampYM),
    proxy.bt.sb.inf.b = as.vector(PFM$proxy.bt.sb.inf.b),
    proxy.bt.sb.sampYM.b = as.vector(PFM$proxy.bt.sb.sampYM.b),
    proxy.bt.sb.inf.b.n = as.vector(PFM$proxy.bt.sb.inf.b.n),
    proxy.bt.sb.sampYM.b.n = as.vector(PFM$proxy.bt.sb.sampYM.b.n),
    simulated.proxy = as.vector(PFM$simulated.proxy),
    simulated.proxy.cal.err = as.vector(PFM$simulated.proxy.cal.err),
    reconstructed.climate = as.vector(PFM$reconstructed.climate),
    stringsAsFactors = FALSE)

  df$timepoints <- PFM$timepoints
  df$replicate <- rep(1:ncol(PFM$proxy.bt.sb.inf.b), each = length(PFM$timepoints))
  df <- dplyr::tbl_df(df)
  df <- tidyr::gather(df, stage, value, -timepoints, -replicate)

  df2 <- data.frame(
    replicate = 1,
    timepoints = PFM$timepoints,
    proxy.bt = PFM$proxy.bt,
    proxy.bt.sb = PFM$proxy.bt.sb,
    clim.signal.ann = PFM$clim.signal.ann,
    clim.timepoints.ssr = PFM$clim.timepoints.ssr,
    stringsAsFactors = FALSE)
  df2 <- tidyr::gather(df2, stage, value, -timepoints, -replicate)

  df.smoothed <- data.frame(
    replicate = 1,
    timepoints = PFM$timepoints.smoothed,
    stage = "clim.signal.smoothed",
    value = PFM$clim.signal.smoothed,
    stringsAsFactors = FALSE)

  rtn <- dplyr::bind_rows(df, df2, df.smoothed)

  rtn <- droplevels(dplyr::filter(rtn, stats::complete.cases(value)))
  rtn <- dplyr::left_join(rtn, stages.key, by = "stage")

  return(rtn)
}

