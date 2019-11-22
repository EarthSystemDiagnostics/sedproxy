# Re-design of ClimToProxyClim

# Requirements
## Return IFA stats, summaries + optionally the individuals values

## Operate on a per-timepoint basis
## pre-allocate any storage
## All the statistics can be calculated from the matrix of weights + clim.matrix
## + samples
## What to do about individual noise - if keeping the individual values then
## this will actually need to be sampled.




#' @md
#' @title Simulate sediment archived proxy records from an input climate signal.
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
#'    \code{sqrt(sigma.meas^2 + sigma.ind^2/n.samples)}.
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
#' @param calibration.type Type of proxy, e.g. Uk'37 or MgCa, to which the
#'   clim.signal is converted before the archiving and measurement of the proxy
#'   is simulated. Defaults to "identity" which means no conversion takes place.
#' @param noise.type Determines whether additive or multiplicative measurement
#'   noise is added. The appropriate type depends on the units of the proxy.
#'   Defaults to multiplicative for MgCa, additive for Uk'37 and identity (none)
#'   calibration types. Can be overidden with a string, "additive" or
#'   "multiplicative" in the case that pre-converted climate signal and
#'   measurement noise values are used in combination with an "identity"
#'   calibration type.
#' @param plot.sig.res The resolution, in years, of the smoothed (block
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
#' @param sigma.meas The standard deviation of the measurement error
#' added to each simulated proxy value.
#' @param sigma.ind The standard deviation of error between individuals
#' (e.g. Forams) not otherwise modelled. This could included "vital effects" or
#' aliasing of depth habitat variation not modelled via a depth resolved input
#' climate signal and habitat weights. sigma.ind is scaled by n.samples
#' before being combined with sigma.meas.
#' @param n.samples Number of e.g. Foraminifera sampled per timepoint, this can
#'   be either a single number, or a vector of length = timepoints. Can be set
#'   to Inf for non-discrete proxies, e.g. for Uk’37.
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
#'   parameter plot.sig.res. This is useful for plotting against the
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
#'    \item{clim.signal.smoothed}{Input climate signal at regular time intervals and resolution = plot.sig.res}
#'    \item{clim.timepoints.ssr}{Input climate signal at requested timepoints, smoothed to resolution = plot.sig.res}
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
#'@export
#'
#'@examples
#' library(ggplot2)
#' set.seed(26052017)
#' clim.in <- ts(N41.t21k.climate[nrow(N41.t21k.climate):1,] - 273.15, start = -39)
#'
#' PFM <- ClimToProxyClim(clim.signal = clim.in,
#'                        timepoints = round(N41.proxy$Published.age),
#'                        calibration.type = "identity",
#'                        habitat.weights = N41.G.ruber.seasonality,
#'                        sed.acc.rate = N41.proxy$Sed.acc.rate.cm.ka,
#'                        layer.width = 1,
#'                        sigma.meas = 0.46,
#'                        sigma.ind = 0, n.samples = Inf,
#'                        plot.sig.res = 10, meas.bias = 1,
#'                        n.replicates = 10)
#'
#' PlotPFMs(PFM$everything, max.replicates = 1, stage.order = "seq") +
#'   facet_wrap(~stage)
#'
#' PlotPFMs(PFM$everything, max.replicates = 1, stage.order = "var")
#'
#' PlotPFMs(PFM$everything, stage.order = "var", plot.stages = "all")
#'
ClimToProxyClim_IFA <- function(clim.signal,
                                timepoints,
                                calibration.type = c("identity", "Uk37", "MgCa"),
                                calibration = switch(calibration.type,
                                                     identity = NA,
                                                     Uk37 = "Mueller global",
                                                     MgCa = "Ten planktonic species_350-500"),
                                slp.int.means = NULL, slp.int.vcov = NULL,
                                noise.type = switch(calibration.type,
                                                    identity = "additive",
                                                    Uk37 = "additive",
                                                    MgCa = "multiplicative"),
                                plot.sig.res = 100,
                                habitat.weights = rep(1/ncol(clim.signal),
                                                      ncol(clim.signal)),
                                habitat.wt.args = NULL,
                                bio.depth = 10,
                                sed.acc.rate = 50,
                                layer.width = 1,
                                sigma.meas = 0,
                                sigma.ind = 0,
                                meas.bias = 0,
                                scale.noise = switch(calibration.type,
                                                     identity = FALSE,
                                                     Uk37 = TRUE,
                                                     MgCa = TRUE),
                                n.samples = Inf,
                                n.replicates = 1,
                                n.bd = 3) {
  # Check inputs --------

  n.timepoints <- length(timepoints)

  CheckLength <- function(vec, correct.length){
    if((length(vec) == 1 | length(vec)==correct.length)==FALSE)
      stop(paste0(deparse(substitute(vec)), " must be either a single value,
                  or a vector the same length as ",
                  deparse(substitute(correct.length))
                  )
           )
  }

  CheckLength(n.samples, n.timepoints)
  CheckLength(sigma.meas, n.timepoints)
  CheckLength(sigma.ind, n.timepoints)
  CheckLength(sed.acc.rate, n.timepoints)

  if (all(is.finite(n.samples))==FALSE & all(is.infinite(n.samples))==FALSE)
    stop("n.samples cannot be a mix of finite and infinite")

  stopifnot(is.matrix(clim.signal))
  if (is.ts(clim.signal)==FALSE)
    stop("Since version 0.3.1 of sedproxy, ClimToProxyClim requires clim.signal
         to be a ts object")

  if (any(is.function(habitat.weights),
          (is.vector(habitat.weights)&length(habitat.weights) == ncol(clim.signal)),
          (is.matrix(habitat.weights)&dim(habitat.weights) == dim(clim.signal))) == FALSE)
    stop("habitat.weights must be either a vector of weights with length = ncol(clim.signal),
         a matrix of weights of the same dimensions as the input climate signal, or a function.
         Function names should be given unquoted, e.g. dnorm, not \"dnorm\"")

  # Convert to proxy units if requested --------
  calibration.type <- match.arg(calibration.type)

  if (calibration.type != "identity") {
    proxy.clim.signal <-
      ProxyConversion(
        temperature = clim.signal,
        calibration.type = calibration.type,
        calibration = calibration,
        slp.int.means = slp.int.means, slp.int.vcov = slp.int.vcov,
        point.or.sample = "point",
        n = 1
      )
  } else{
    proxy.clim.signal <- clim.signal
  }

  # Calculate timepoint invariant values ------
  max.clim.signal.i <- end(clim.signal)[1]
  min.clim.signal.i <- start(clim.signal)[1]

  # Create smoothed climate signal --------
  if (is.na(plot.sig.res)) {
    timepoints.smoothed <- NA
    clim.signal.smoothed <- NA
  } else{
    timepoints.smoothed <- seq(min.clim.signal.i-1, max.clim.signal.i-1, by = plot.sig.res)
    clim.signal.smoothed <- ChunkMatrix(timepoints.smoothed, plot.sig.res,
                                        clim.signal)
  }

  # Check whether bioturbation window will extend beyond climate signal for any of the timepoints

  # bioturbation window will be focal.timepoint - bio.depth.timesteps - layer.width.years / 2 to
  # focal.timepoint + 3*bio.depth.timesteps

  bio.depth.timesteps <- round(1000 * bio.depth / sed.acc.rate)
  layer.width.years <- ceiling(1000 * layer.width / sed.acc.rate)

  max.min.windows <- cbind(max = timepoints + n.bd * bio.depth.timesteps,
                           min = timepoints - bio.depth.timesteps - layer.width.years / 2)

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


  valid.inds <- max.ind == FALSE & min.ind == FALSE
  timepoints <- timepoints[valid.inds]
  n.timepoints <- length(timepoints)

  max.min.windows <- max.min.windows[valid.inds, , drop = FALSE]

  # Scale sigma.ind by n.samples and create combined error term
  sigma.ind.scl <- ifelse(is.finite(n.samples),
                          sigma.ind / sqrt(n.samples), 0)

  sigma.meas.ind <- sqrt(sigma.meas^2 + sigma.ind.scl^2)


  # Get bioturbation windows for all timepoints
  # Create vectors from "scalar" inputs
  if (length(sed.acc.rate) == 1) {
    sed.acc.rate <- rep(sed.acc.rate, n.timepoints)
  }

  if (length(layer.width) == 1) {
    layer.width <- rep(layer.width, n.timepoints)
  }

  if (length(n.samples) == 1) {
    n.samples <- rep(n.samples, n.timepoints)
  }

  if (length(sigma.meas) == 1) {
    sigma.meas <- rep(sigma.meas, n.timepoints)
  }

  if (length(sigma.ind) == 1) {
    sigma.ind <- rep(sigma.ind, n.timepoints)
  }

  # Remove timepoint specific parameters that exceed clim.signal ------

  if (length(sed.acc.rate) > n.timepoints) sed.acc.rate <- sed.acc.rate[valid.inds]
  if (length(layer.width) > n.timepoints) layer.width <- layer.width[valid.inds]
  if (length(n.samples) > n.timepoints) n.samples <-  n.samples[valid.inds]
  if (length(sigma.meas) > n.timepoints) sigma.meas <- sigma.meas[valid.inds]
  if (length(sigma.ind) > n.timepoints) sigma.ind <- sigma.ind[valid.inds]


 # Generate productivity weights from function if supplied
  if (is.function(habitat.weights)){
    FUN <- match.fun(habitat.weights)
    habitat.weights <- do.call(FUN, args = c(list(x = clim.signal), habitat.wt.args))
    habitat.weights <- habitat.weights / sum(habitat.weights)
  }

  # If vector ensure habitat.weights are weights and matrix
  if (is.vector(habitat.weights)){
    habitat.weights <- habitat.weights / sum(habitat.weights)
    habitat.weights <- matrix(rep(habitat.weights, nrow(clim.signal)),
                              nrow = nrow(clim.signal), byrow = TRUE)
  }

  bt.windows <- lapply(seq_along(timepoints), function(tp){
    max.min.windows[tp, "min"]:max.min.windows[tp, "max"]
  })

  bt.wts <- lapply(seq_along(timepoints), function(tp){
    bioturb.weights <- BioturbationWeights(z = bt.windows[[tp]], focal.z = timepoints[tp],
                                           layer.width = layer.width[tp],
                                           sed.acc.rate = sed.acc.rate[tp],
                                           bio.depth = bio.depth)

    biot.sig.weights <- bioturb.weights %o% rep(1, ncol(proxy.clim.signal))
    biot.sig.weights <- biot.sig.weights / sum(biot.sig.weights)
    return(biot.sig.weights)
  })

  hab.wts <- lapply(seq_along(timepoints), function(tp){
    habitat.weights <- habitat.weights[bt.windows[[tp]], , drop = FALSE]
    habitat.weights <- habitat.weights / sum(habitat.weights)
    })

  clim.sig.wts <- lapply(seq_along(timepoints), function(tp) {
    clim.sig.wts <- hab.wts[[tp]] * bt.wts[[tp]]
    clim.sig.wts / sum(clim.sig.wts)
    })

  sample.inds <- lapply(seq_along(timepoints), function(tp) {
    sample(length(clim.sig.wts[[tp]]),
           n.samples[tp] * n.replicates,
           prob = clim.sig.wts[[tp]],
           replace = TRUE)
  })

  # get sampled rows

  row.inds.rel <- lapply(seq_along(timepoints), function(tp) {
    inds <- (sample.inds[[tp]]-1) %% nrow(clim.sig.wts[[tp]]) + 1
    inds
  })

  row.inds.abs <- lapply(seq_along(timepoints), function(tp) {
    row.inds.rel[[tp]] + max.min.windows[tp, "min"] - min.clim.signal.i
  })

  col.inds <- lapply(seq_along(timepoints), function(tp) {
    ceiling(sample.inds[[tp]] / nrow(hab.wts[[tp]]))
    })

  clim.sig.windows <- lapply(seq_along(timepoints), function(tp) {
    proxy.clim.signal[bt.windows[[tp]], , drop=FALSE]
  })

  samples <- lapply(seq_along(timepoints), function(tp) {
   samp <- clim.sig.windows[[tp]][sample.inds[[tp]]]
   matrix(samp, nrow = n.samples[tp])
  })

  proxy.bt <- sapply(seq_along(timepoints), function(tp) {
    sum(clim.sig.windows[[tp]] * bt.wts[[tp]])
  })

  proxy.bt.sb <- sapply(seq_along(timepoints), function(tp) {
    sum(clim.sig.windows[[tp]] * clim.sig.wts[[tp]])
  })

  proxy.bt.sampY <- sapply(seq_along(timepoints), function(tp) {
    hw <- hab.wts[[tp]]
    hw <- hw / rowSums(hw)
    rows <- rowSums(clim.sig.windows[[tp]] * hw)[row.inds.rel[[tp]]]
    colMeans(matrix(rows, nrow = n.samples[tp]))
  })

  proxy.bt.sb.sampYM <- sapply(seq_along(timepoints), function(tp) {
    # can be taken across columns to get per replicate means
    samp <- matrix(samples[[tp]], nrow = n.samples[tp])
    colMeans(samp)
  })

  # ind.errs <- lapply(seq_along(timepoints), function(tp) {
  #   errs <- rnorm(n.samples[tp] * n.replicates, 0, sigma.ind)
  #   matrix(errs, nrow = n.samples[tp])
  # })

  # meas.errs <- rnorm(length(timepoints) * n.replicates, 0, sigma.meas)

  IFA <- list(row.inds = simplify2array(row.inds.abs),
              row.inds.rel = simplify2array(row.inds.rel),
              col.inds = simplify2array(col.inds))

  out <- list(
    timepoints = timepoints,
    proxy.bt = proxy.bt,
    proxy.bt.sb = proxy.bt.sb,
    proxy.bt.sampY = t(proxy.bt.sampY),
    proxy.bt.sb.sampYM = t(proxy.bt.sb.sampYM)
  )

  # # #return(out)
  # RestructOut <- function(out, n.replicates){
  #   if (n.replicates == 1){
  #     tmp <- apply(out, 1, function(x) simplify2array(x))
  #     as.list(data.frame(apply(tmp, 2, as.vector)))
  #   }else{
  #     apply(out, 1, function(x) simplify2array(x))
  #   }
  # }
  #
  # out <- RestructOut(out, n.replicates = n.replicates)
  # #out <- apply(out, 1, function(x) simplify2array(x))
  # #out <- plyr::alply(out, 1, function(x) simplify2array(x), .dims = TRUE)
  # #
  # #  # remove extra attributes added by alply
  # #  attr(out, "split_type") <- NULL
  # #  attr(out, "split_labels") <- NULL

  #print(out)

  # if (n.replicates == 1) out$proxy.bt.sb.sampYM <- matrix(out$proxy.bt.sb.sampYM, nrow = 1)
  # out$proxy.bt.sb.sampYM <- t(out$proxy.bt.sb.sampYM)
  #
  # if (n.replicates == 1) out$proxy.bt.sb.sampY <- matrix(out$proxy.bt.sb.sampY, nrow = 1)
  # out$proxy.bt.sb.sampY <- t(out$proxy.bt.sb.sampY)


  # Add bias and noise --------
  # Rescale noise if using a calibration -----
  if (scale.noise != FALSE) {
    # scale.noise will either be TRUE because a non-identity calibration is being used
    # or it will be a string to identify the correct calibration if the input time-series
    # has already been converted.
    message("Rescaling noise")

    # If cal type is identity re-scaling still required
    pct <- if (calibration.type == "identity") {
      scale.noise
    } else{
      calibration.type
    }

    # mean temperature in temperature units at each timepoint - use bioturbated signal
    mean.temperature <-  as.vector(ProxyConversion(proxy.value = out$proxy.bt,
                                                   calibration.type = pct, calibration = calibration,
                                                   slp.int.means = slp.int.means, slp.int.vcov = slp.int.vcov))

    sigma.meas.ind <- ProxyConversion(temperature = mean.temperature + sigma.meas.ind,
                                      calibration.type = pct, calibration = calibration,
                                      slp.int.means = slp.int.means, slp.int.vcov = slp.int.vcov) -
      ProxyConversion(temperature = mean.temperature,
                      calibration.type = pct, calibration = calibration,
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

  out$simulated.proxy = out$proxy.bt.sb.sampYM.b.n


  # Create smoothed climate signal -----
  if (is.na(plot.sig.res)) {
    out$clim.timepoints.ssr <- NA

  } else{
    out$clim.timepoints.ssr <- ChunkMatrix(timepoints, plot.sig.res, clim.signal)
  }

  # Add items to output list -----------
  out$timepoints = timepoints
  out$clim.signal.ann = rowSums(clim.signal[timepoints,  , drop = FALSE]) / ncol(clim.signal)
  
  #out$sed.acc.rate = sed.acc.rate
  

  # Organise output -------
  # simulated.proxy <-
  #   dplyr::tbl_df(out[c(
  #     "timepoints",
  #     "clim.signal.ann",
  #     "clim.timepoints.ssr",
  #     "proxy.bt",
  #     "proxy.bt.sb"#,
  #     #"sed.acc.rate",
  #     #"smoothing.width"
  #   )])
  #
  # simulated.proxy$proxy.bt.sb.sampY <- out$proxy.bt.sb.sampY[, 1, drop = TRUE]
  # simulated.proxy$proxy.bt.sb.sampYM <- out$proxy.bt.sb.sampYM[, 1, drop = TRUE]
  # simulated.proxy$proxy.bt.sb.inf.b <- out$proxy.bt.sb.inf.b[, 1, drop = TRUE]
  # simulated.proxy$proxy.bt.sb.inf.b.n <- out$proxy.bt.sb.inf.b.n[, 1, drop = TRUE]
  # simulated.proxy$proxy.bt.sb.sampYM.b <- out$proxy.bt.sb.sampYM.b[, 1, drop = TRUE]
  # simulated.proxy$proxy.bt.sb.sampYM.b.n <- out$proxy.bt.sb.sampYM.b.n[, 1, drop = TRUE]

  # if (all(is.finite(n.samples))) {
  #   simulated.proxy$simulated.proxy <- simulated.proxy$proxy.bt.sb.sampYM.b.n
  #   out$simulated.proxy <- out$proxy.bt.sb.sampYM.b.n
  # } else{
  #   simulated.proxy$simulated.proxy <- simulated.proxy$proxy.bt.sb.inf.b.n
  #   out$simulated.proxy <- out$proxy.bt.sb.inf.b.n
  # }


  # Add calibration uncertainty -------
  # If n.replicates > 1
  # First convert back to temperature units with fixed parameters
  # Then re-convert to proxy units with random parameters
  if (calibration.type != "identity"){

    if (n.replicates > 1){
      out$simulated.proxy.cal.err <-
        ProxyConversion(proxy.value = out$simulated.proxy,
                        calibration.type = calibration.type,
                        calibration = calibration,
                        point.or.sample = "point", n = 1,
                        slp.int.means = slp.int.means,
                        slp.int.vcov = slp.int.vcov)

      out$simulated.proxy.cal.err <-
        ProxyConversion(temperature = out$simulated.proxy.cal.err,
                        calibration.type = calibration.type,
                        calibration = calibration,
                        point.or.sample = "sample", n = n.replicates,
                        slp.int.means = slp.int.means, slp.int.vcov = slp.int.vcov)
    }else{
      out$simulated.proxy.cal.err <- out$simulated.proxy
    }

    # Do this in all cases, not just if n.replicates == 1
    out$reconstructed.climate <-
      ProxyConversion(proxy.value = out$simulated.proxy,
                      calibration.type = calibration.type,
                      calibration = calibration,
                      point.or.sample = "point", n = 1,
                      slp.int.means = slp.int.means, slp.int.vcov = slp.int.vcov)

  }else{
    out$simulated.proxy.cal.err <- out$simulated.proxy
    out$reconstructed.climate <- out$simulated.proxy
  }

  #simulated.proxy$simulated.proxy.cal.err <- out$simulated.proxy.cal.err[, 1, drop = TRUE]
  #simulated.proxy$reconstructed.climate <- out$reconstructed.climate[, 1, drop = TRUE]
  #print(out)
  #everything.2 <- MakePFMDataframe(out)
  #everything <- out
  everything <- as.data.frame(lapply(out, function(x) as.numeric(x))) 
  everything$replicate <- rep(1:n.replicates,
                              each = length(out$timepoints))
  
  everything <- gather(everything, stage, value, -timepoints, -replicate)  
  everything <- select(everything, timepoints, replicate, everything())
  
  smoothed.signal <- data.frame(
    timepoints = timepoints.smoothed,
    replicate = 1,
    value = clim.signal.smoothed
  )
  
  smoothed.signal$stage <- "clim.signal.smoothed"

  everything <- bind_rows(everything, smoothed.signal)
  
  
  slp.int.means <-
    if (is.null(slp.int.means) && calibration.type != "identity") {
      calibration <-
        if (calibration.type == "MgCa" & is.null(calibration)) {
          "Ten planktonic species_350-500"
        } else {
          calibration
        }

      cp <- data.frame(sedproxy::calibration.parameters)
      cfs.vcov <- cp[cp$calibration.type == calibration.type &
                       cp$calibration == calibration,]
      matrix(c(cfs.vcov$slope, cfs.vcov$intercept),
             ncol = 2,
             byrow = TRUE)
    } else{
      slp.int.means
    }

  calibration.pars <- list(calibration.type = calibration.type,
                           calibration = calibration,
                           slp.int.means = slp.int.means)

  #attr(simulated.proxy, "calibration.pars") <-  calibration.pars
  attr(everything, "calibration.pars") <-  calibration.pars

  out <- list(#simulated.proxy=simulated.proxy,
    IFA = IFA,
    smoothed.signal=smoothed.signal,
    everything = everything,
    calibration.pars = calibration.pars)

  class(out) <- "sedproxy.pfm"

  return(out)
}

GetIndices <- function(tp, sed.acc.rate, layer.width, bio.depth){

  # Get bioturbation weights --------

  bio.wdth <-

    first.tp <- max.min.windows[tp, "min"]
  last.tp <- max.min.windows[tp, "max"]
  bioturb.window <- first.tp:last.tp

  bioturb.weights <- BioturbationWeights(z = bioturb.window, focal.z = timepoints[tp],
                                         layer.width = layer.width[tp], sed.acc.rate = sed.acc.rate[tp],
                                         bio.depth = bio.depth)

  # Get bioturbation X seasonality weights matrix ---------
  habitat.weights <- habitat.weights[first.tp:last.tp - min.clim.signal.i+1, , drop = FALSE]
  habitat.weights <- habitat.weights / sum(habitat.weights)
  clim.sig.weights <- bioturb.weights * habitat.weights
  clim.sig.weights <- clim.sig.weights / sum(clim.sig.weights)

  # Check weights sum to 1, within tolerance
  weight.err <- abs(sum(clim.sig.weights) - 1)
  if ((weight.err < 1e-10) == FALSE) stop(paste0("weight.err = ", weight.err))

  # call sample once for all replicates together, then take means of
  # groups of n.samples
  # Get indices not values
  samp.indices <-  sample(length(clim.sig.window),
                          n.samples[tp] * n.replicates,
                          prob = clim.sig.weights,
                          replace = TRUE)

}
