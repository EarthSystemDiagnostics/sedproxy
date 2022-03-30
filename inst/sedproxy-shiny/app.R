library(ggplot2) 
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
#'   rate \emph{sed.acc.rate} and the bioturbation depth \emph{bio.depth} which defaults
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
#'   multiple replicate proxies are generated (\emph{n.replicates} > 1) each
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
#' @param scale.noise Scale noise to proxy units. Defaults to TRUE if
#' calibration.type is not "identity"
#' @param n.replicates Number of replicate proxy time-series to simulate from
#'   the climate signal
#' @param n.bd Number of multiples of the bioturbation width at which to truncate
#' the bioturbation filter
#' @param top.of.core The theoretical minimum age at the top of the core, ie.
#' the year the core was sampled, defaults to the start of clim.in
#' @inheritParams ProxyConversion
#' @return \code{ClimToProxyClim} returns an object of class "sedproxy.pfm", a list with three elements:
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
#' @importFrom dplyr rename
#' @importFrom rlang .data
#' @export
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
ClimToProxyClim <- function(clim.signal,
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
                            top.of.core = NULL,
                            n.bd = 3) {
  # Check inputs --------

  n.timepoints <- length(timepoints)

  if((length(n.samples) == 1 | length(n.samples)==n.timepoints)==FALSE)
    stop("n.sample must be either a single value, or a vector the same
         length as timepoints")

  if((length(sigma.meas) == 1 | length(sigma.meas)==n.timepoints)==FALSE)
    stop("sigma.meas must be either a single value, or a vector the same
         length as timepoints")

  if((length(sigma.ind) == 1 | length(sigma.ind)==n.timepoints)==FALSE)
    stop("sigma.ind must be either a single value, or a vector the same
         length as timepoints")

  if (all(is.finite(n.samples))==FALSE & all(is.infinite(n.samples))==FALSE)
    stop("n.samples cannot be a mix of finite and infinite")

  stopifnot(is.matrix(clim.signal))
  if (stats::is.ts(clim.signal)==FALSE)
    stop("Since version 0.3.1 of sedproxy, ClimToProxyClim requires clim.signal to be a ts object")
  stopifnot(length(sed.acc.rate) == n.timepoints |
              length(sed.acc.rate) == 1)

  if (any(is.function(habitat.weights),
          (is.vector(habitat.weights)&length(habitat.weights) == ncol(clim.signal)),
          (is.matrix(habitat.weights)&dim(habitat.weights) == dim(clim.signal))) == FALSE)
    stop("habitat.weights must be either a vector of weights with length = ncol(clim.signal),
         a matrix of weights of the same dimensions as the input climate signal, or a function.
         Function names should be given unquoted, e.g. dnorm, not \"dnorm\"")


  if (is.null(top.of.core)){
    top.of.core <- stats::time(clim.signal)[1]
  }else{
    if (top.of.core < stats::time(clim.signal)[1]) stop("top.of.core cannot be younger than the start of clim.signal")
  }


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
  max.clim.signal.i <- stats::end(clim.signal)[1]
  min.clim.signal.i <- stats::start(clim.signal)[1]

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

  max.min.windows <- cbind(max = ceiling(timepoints + n.bd * bio.depth.timesteps),
                           min = floor(timepoints - bio.depth.timesteps - layer.width.years / 2))

  max.ind <- max.min.windows[,"max"] >= max.clim.signal.i
  min.ind <- max.min.windows[,"min"] <  min.clim.signal.i



# Use Rapid or Slow version ----------------------

  if (length(bio.depth)==1 && length(sed.acc.rate)==1 &&
      length(layer.width)==1 && length(n.samples)==1 &&
      (is.function(habitat.weights) == FALSE) && is.vector(habitat.weights)){

    # Rapid ------
    message("Using Rapid version")

    # Find mixed layer points ------
    # keep points in the mixed layer as well as those below but still inside time-signal
    tpts.above.core.top <- timepoints < top.of.core
    valid.inds <- max.ind == FALSE & tpts.above.core.top == FALSE

    # identify mixed layer
    mixed.layer.inds <-  min.ind == TRUE & tpts.above.core.top == FALSE
    mixed.layer.inds <- mixed.layer.inds[valid.inds]

    if (any(max.ind))
      warning(paste0("One or more requested timepoints is too old. Bioturbation window(s) for timepoint(s) ",
                     paste(timepoints[max.ind], collapse = ", "),
                     " extend(s) beyond end of input climate signal. Returning pseudo-proxy for valid timepoints."))

    if (any(max.min.windows[,"min"] < min.clim.signal.i))
      warning(paste0("Timepoint(s) ",
                     paste(timepoints[mixed.layer.inds], collapse = ", "),
                     " are in the mixed layer"))

    if (any(tpts.above.core.top))
      warning(paste0("One or more requested timepoints is too recent. Timepoint(s) ",
                       paste(timepoints[tpts.above.core.top], collapse = ", "),
                       " are more recent than the top of the core."))

    # remove too old or young timepoints
    timepoints <- timepoints[valid.inds]
    n.timepoints <- length(timepoints)

    # adjusted timepoints for the mixed layer
    # in the mixed layer the bioturbation window is centred around the
    # bottom of the mixed layer
    timepoints.adj <- timepoints
    timepoints.adj[mixed.layer.inds] <- 1 + bio.depth.timesteps + layer.width.years / 2

    #max.min.windows <- max.min.windows[valid.inds, , drop = FALSE]

    # # reset mixed window for mixed layer points
    # max.min.windows[mixed.layer.inds, ] <-
    #   c((n.bd+1) * bio.depth.timesteps + layer.width.years / 2, 0)

    # Scale sigma.ind by n.samples and create combined error term
    sigma.ind.scl <- ifelse(is.finite(n.samples),
                            sigma.ind / sqrt(n.samples), 0)

    sigma.meas.ind <- sqrt(sigma.meas^2 + sigma.ind.scl^2)




    # Ensure seasonal productivities are weights
    habitat.weights <- habitat.weights / sum(habitat.weights)


    # Get relative bioturbation window ----------
    #browser()

    first.tp <- -bio.depth.timesteps - layer.width.years / 2
    last.tp <- n.bd * bio.depth.timesteps

    bioturb.window <- first.tp:last.tp

    # Get bioturbation weights --------
    bioturb.weights <- BioturbationWeights(z = bioturb.window, focal.z = 0,
                                           layer.width = layer.width,
                                           sed.acc.rate = sed.acc.rate,
                                           bio.depth = bio.depth)


    # Get bioturbation X no-seasonality weights matrix ---------
    biot.sig.weights <- bioturb.weights %o% rep(1, ncol(proxy.clim.signal))
    biot.sig.weights <- biot.sig.weights / sum(biot.sig.weights)

    # Get bioturbation X seasonality weights matrix ---------
    clim.sig.weights <- bioturb.weights %*% t(habitat.weights)
    clim.sig.weights <- clim.sig.weights / sum(clim.sig.weights)

    # Check weights sum to 1, within tolerance
    weight.err <- abs(sum(clim.sig.weights) - 1)
    if ((weight.err < 1e-10) == FALSE) stop(paste0("weight.err = ", weight.err))

    # Do sampling ------
    if (is.finite(n.samples)) {
      # call sample once for all replicates and timepoints together, then take means of
      # groups of n.samples
      # Get indices not values
      tot.n.samples <- n.samples * n.replicates * n.timepoints

      length.clim.sig.w <- length(clim.sig.weights)

      samp.indices <-  sample(length.clim.sig.w,
                              tot.n.samples,
                              prob = clim.sig.weights,
                              replace = TRUE)}

    # For each timepoint ------
    out <- sapply(1:n.timepoints, function(tp) {

      # Get portion of clim.signal corresponding to bioturbation window for this timepoint -------
      clim.sig.window <-
        proxy.clim.signal[(bioturb.window + timepoints.adj[tp] - min.clim.signal.i +
                            1), , drop = FALSE]

      # Calculate mean clim.signal -------

      # Just bioturbation
      proxy.bt <- sum(biot.sig.weights * clim.sig.window)

      # Bioturbation + seasonal bias
      proxy.bt.sb <- sum(clim.sig.weights * clim.sig.window)

      # Bioturbation + seasonal bias + aliasing
      if (is.infinite(n.samples)) {
        proxy.bt.sb.sampY <- rep(NA, n.replicates)
        proxy.bt.sb.sampYM <- rep(NA, n.replicates)
      } else if (is.finite(n.samples)) {
        # Use previously sampled indices
        # get indices for this timepoint
        samp.indices.tp <- samp.indices[((tp-1)*n.samples*n.replicates+1):(tp*n.samples*n.replicates)]

        # convert vector to matrix (cheap only attributes changed), then means
        # can be taken across columns to get per replicate means
        samp <- matrix(clim.sig.window[samp.indices.tp], nrow = n.samples)
        #proxy.bt.sb.sampYM <- apply(samp, 2, mean)
        proxy.bt.sb.sampYM <- colMeans(samp)

        #Get without seasonal aliasing (bioturbation aliasing only)
        clim.sig.window.ann <- colSums(t(clim.sig.window) * habitat.weights)
        col.indices <- (samp.indices.tp-1) %% nrow(clim.sig.window) + 1

        samp.bt <- matrix(clim.sig.window.ann[col.indices], ncol = n.samples)
        proxy.bt.sb.sampY <- rowMeans(samp.bt)
      }

      # Gather output ----------
      list(
        #smoothing.width = smoothing.width,
        proxy.bt = proxy.bt,
        proxy.bt.sb = proxy.bt.sb,
        proxy.bt.sb.sampY = proxy.bt.sb.sampY,
        proxy.bt.sb.sampYM = proxy.bt.sb.sampYM)
    })
  }else{

# Slow ----
    #browser()
    # Find mixed layer points ------
    # keep points in the mixed layer as well as those below but still inside time-signal
    tpts.above.core.top <- timepoints < top.of.core

    # identify mixed layer
    # find oldest timepoint in mixed layer

    oldest.in.mix <- which.max(timepoints[min.ind == TRUE])

    if (length(oldest.in.mix)!=0) {
      mixed.layer.inds <- timepoints <= timepoints[oldest.in.mix] & tpts.above.core.top == FALSE
      #mixed.layer.inds <- mixed.layer.inds[tpts.above.core.top == FALSE]
    } else {
      mixed.layer.inds <- rep(FALSE, n.timepoints)
    }


    valid.inds <- max.ind == FALSE & tpts.above.core.top == FALSE


    if (any(max.ind))
      warning(paste0("One or more requested timepoints is too old. Bioturbation window(s) for timepoint(s) ",
                     paste(timepoints[max.ind], collapse = ", "),
                     " extend(s) beyond end of input climate signal. Returning pseudo-proxy for valid timepoints."))

    if (any(mixed.layer.inds))
      warning(paste0("Timepoint(s) ",
                     paste(timepoints[mixed.layer.inds], collapse = ", "),
                     " are in the mixed layer"))

    if (any(tpts.above.core.top))
      warning(paste0("One or more requested timepoints is too recent. Timepoint(s) ",
                     paste(timepoints[tpts.above.core.top], collapse = ", "),
                     " are more recent than the top of the core."))


    timepoints <- timepoints[valid.inds]
    n.timepoints <- length(timepoints)
    mixed.layer.inds <- mixed.layer.inds[valid.inds]


    # Scale sigma.ind by n.samples and create combined error term
    sigma.ind.scl <- ifelse(is.finite(n.samples),
                            sigma.ind / sqrt(n.samples), 0)

    sigma.meas.ind <- sqrt(sigma.meas^2 + sigma.ind.scl^2)



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

    if (length(n.samples) > n.timepoints) n.samples <-    n.samples   [valid.inds]

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


    timepoints.adj <- timepoints

    max.min.windows <- max.min.windows[valid.inds, , drop = FALSE]


    # set mixed layer sed.acc.rate to the lowest
    if(any(mixed.layer.inds)){
      sed.acc.rate[mixed.layer.inds] <- min(sed.acc.rate[mixed.layer.inds])

    # reset mixed window for mixed layer points
    bio.depth.timesteps <- round(1000 * bio.depth / sed.acc.rate)
    layer.width.years <- ceiling(1000 * layer.width / sed.acc.rate)


    # adjusted timepoints for the mixed layer
    # in the mixed layer the bioturbation window is centred around the
    # bottom of the mixed layer

    timepoints.adj[mixed.layer.inds] <- 1 +
      bio.depth.timesteps[mixed.layer.inds] +
      layer.width.years[mixed.layer.inds] / 2


    max.min.windows[mixed.layer.inds, ] <-
    cbind(ceiling((n.bd+1) * bio.depth.timesteps[mixed.layer.inds] +
            layer.width.years[mixed.layer.inds] / 2), top.of.core)

    }

    # For each timepoint ------
    out <- sapply(1:n.timepoints, function(tp) {

      # Get bioturbation window ----------
      first.tp <- max.min.windows[tp, "min"]
      last.tp <- max.min.windows[tp, "max"]
      bioturb.window <- first.tp:last.tp

      # Get bioturbation weights --------
      bioturb.weights <- BioturbationWeights(z = bioturb.window, focal.z = timepoints.adj[tp],
                                             layer.width = layer.width[tp], sed.acc.rate = sed.acc.rate[tp],
                                             bio.depth = bio.depth)


      clim.sig.window <-  proxy.clim.signal[which(stats::time(clim.signal)%in%(first.tp:last.tp)), , drop = FALSE]



      # Get bioturbation X no-seasonality weights matrix ---------
      biot.sig.weights <- bioturb.weights %o% rep(1, ncol(proxy.clim.signal))
      biot.sig.weights <- biot.sig.weights / sum(biot.sig.weights)

      #browser()


      # Get bioturbation X seasonality weights matrix ---------
      habitat.weights <- habitat.weights[which(stats::time(clim.signal)%in%(first.tp:last.tp+1)), , drop = FALSE]
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

        # habitat weights rows need to sum to 1
        habitat.weights.r1 <- habitat.weights / rowSums(habitat.weights)

        clim.sig.window.ann <- rowSums(clim.sig.window * habitat.weights.r1)
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

  }

  # #return(out)
  RestructOut <- function(out, n.replicates){
    if (n.replicates == 1){
      tmp <- apply(out, 1, function(x) simplify2array(x))
      as.list(data.frame(apply(tmp, 2, as.vector)))
    }else{
      apply(out, 1, function(x) simplify2array(x))
    }
  }

  out <- RestructOut(out, n.replicates = n.replicates)
  # #out <- apply(out, 1, function(x) simplify2array(x))
  # out <- plyr::alply(out, 1, function(x) simplify2array(x), .dims = TRUE)
  #
  #  # remove extra attributes added by alply
  #  attr(out, "split_type") <- NULL
  #  attr(out, "split_labels") <- NULL

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

    out$proxy.bt.sb.sampYM.b <- sweep(out$proxy.bt.sb.sampYM, 2, bias, FUN = "+")
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
    out$proxy.bt.sb.sampYM.b <- sweep(out$proxy.bt.sb.sampYM, 2, bias, "*")
    out$proxy.bt.sb.sampYM.b.n <- out$proxy.bt.sb.sampYM.b * noise

    # set intermediate bias stages to NA if no bias modelled
    if (meas.bias == 0) {
      out$proxy.bt.sb.inf.b[,] <- NA
      out$proxy.bt.sb.sampYM.b[,] <- NA
    }
  }


  # Create smoothed climate signal -----
  if (is.na(plot.sig.res)) {
    out$clim.timepoints.ssr <- NA

  } else{
    out$clim.timepoints.ssr <- ChunkMatrix(timepoints, plot.sig.res, clim.signal)
  }

  # Add items to output list -----------
  out$timepoints = timepoints
  out$n.samples = n.samples
  out$clim.signal.ann = rowSums(
    # clim.signal[time(clim.signal) %in% timepoints, , drop = FALSE]
    clim.signal[match(timepoints, stats::time(clim.signal)), , drop = FALSE]
    ) / ncol(clim.signal)
  #out$sed.acc.rate = sed.acc.rate
  out$timepoints.smoothed = timepoints.smoothed
  out$clim.signal.smoothed = clim.signal.smoothed

  # Organise output -------
  simulated.proxy <-
    dplyr::as_tibble(out[c(
      "timepoints",
      "n.samples",
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


  smoothed.signal <- dplyr::as_tibble(out[c(
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

  simulated.proxy$simulated.proxy.cal.err <- out$simulated.proxy.cal.err[, 1, drop = TRUE]
  simulated.proxy$reconstructed.climate <- out$reconstructed.climate[, 1, drop = TRUE]

  everything <- MakePFMDataframe(out)

  slp.int.means <-
    if (is.null(slp.int.means) && calibration.type != "identity") {
      calibration <-
        if (calibration.type == "MgCa" & is.null(calibration)) {
          "Ten planktonic species_350-500"
        } else {
          calibration
        }

      cp <- data.frame(calibration.parameters)
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

  if (stats::is.ts(climate.matrix)) {
    rel.wind <- 1:width -round(width/2)
    strt <- stats::start(climate.matrix)[1]
    n.row <- nrow(climate.matrix)
    sapply(timepoints, function(tp){
      inds <- rel.wind + tp - strt + 1
      inds <- inds[inds > 0 & inds < n.row]
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
#' @importFrom rlang .data
#' @keywords internal
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
  df$n.samples <- PFM$n.samples
  df$replicate <- rep(1:ncol(PFM$proxy.bt.sb.inf.b), each = length(PFM$timepoints))
  df <- dplyr::as_tibble(df)
  #df <- tidyr::gather(df, stage, value, -timepoints, -n.samples, -replicate)
  df <- tidyr::pivot_longer(df, cols = tidyr::contains(c("proxy", "climate")),
                          names_to = "stage", values_to = "value")
  df <- dplyr::arrange(df, .data$replicate, .data$stage, .data$timepoints)
  df2 <- data.frame(
    replicate = 1,
    timepoints = PFM$timepoints,
    n.samples = PFM$n.samples,
    proxy.bt = PFM$proxy.bt,
    proxy.bt.sb = PFM$proxy.bt.sb,
    clim.signal.ann = PFM$clim.signal.ann,
    clim.timepoints.ssr = PFM$clim.timepoints.ssr,
    stringsAsFactors = FALSE)
  #df2 <- tidyr::gather(df2, stage, value, -timepoints, -n.samples, -replicate)
  df2 <- tidyr::pivot_longer(df2, cols = tidyr::contains(c("proxy", "clim")),
                             names_to = "stage", values_to = "value")
  df2 <- dplyr::arrange(df2, .data$replicate, .data$stage, .data$timepoints)
  df.smoothed <- data.frame(
    replicate = 1,
    timepoints = PFM$timepoints.smoothed,
    stage = "clim.signal.smoothed",
    value = PFM$clim.signal.smoothed,
    stringsAsFactors = FALSE)

  rtn <- dplyr::bind_rows(df, df2, df.smoothed)

  rtn <- droplevels(dplyr::filter(rtn, stats::complete.cases(.data$value)))
  rtn <- dplyr::left_join(rtn, dplyr::select(sedproxy::stages.key, stage, scale, .data$label), by = "stage")
  #rtn <- select(rtn, -plotting.colour, -plotting.alpha)

  return(rtn)
}

#' Plot forward modelled sedimentary proxies
#'
#' @param PFMs A dataframe of forward modelled proxies
#' @param stage.order Controls the order in which proxy stages are plotted,
#' either sequentially, "seq", or in order of variance, "var". Defaults to var.
#' @param plot.stages Proxy stages to be plotted, "default", "all", or a custom character vector
#' @param colr.palette Colours for the proxy stages
#' @param alpha.palette Alpha levels for the proxy stages
#' @param levl.labels Labels for the proxy stages
#' @param max.replicates Maximum number of replicates to plot at once
#'
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @return a ggplot object of class "gg" "ggplot"
#' @export PlotPFMs
#'
#' @examples
#' library(ggplot2)
#' set.seed(26052017)
#' clim.in <- ts(N41.t21k.climate[nrow(N41.t21k.climate):1,] - 273.15)
#'
#' PFM <- ClimToProxyClim(clim.signal = clim.in,
#'                        timepoints = round(N41.proxy$Published.age),
#'                        calibration.type = "identity",
#'                        habitat.weights = N41.G.ruber.seasonality,
#'                        sed.acc.rate = N41.proxy$Sed.acc.rate.cm.ka,
#'                        sigma.meas = 0.45,
#'                        sigma.ind = 0,
#'                        n.samples = Inf,
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
PlotPFMs <- function(PFMs,
                     stage.order = c("var", "seq"),
                     plot.stages = c("default"),
                     max.replicates = 5,
                     colr.palette = "default",
                     alpha.palette = "default",
                     levl.labels = "default"){

  PFMs.in <- PFMs
  if ("sedproxy.pfm" %in% class(PFMs.in)) PFMs <- PFMs.in$everything

  if(exists("replicate", where = PFMs)){
    rug.dat <- dplyr::filter(PFMs, stage %in% c("simulated.proxy", "observed.proxy"),
                             replicate == 1)
  }else{
    rug.dat <- dplyr::filter(PFMs, stage %in% c("simulated.proxy", "observed.proxy"))
    rug.dat$replicate <- 1
    PFMs$replicate <- 1
  }

  if(exists("Location", where = PFMs)==FALSE){
    PFMs$Location <- ""
  }

  if(exists("ID.no", where = PFMs)==FALSE){
    PFMs$ID.no <- ""
  }
  if(exists("Proxy", where = PFMs)==FALSE){
    PFMs$Proxy <- ""
  }

  # assign default asthetic mappings

  breaks <- sedproxy::stages.key$stage

  if (colr.palette[1] == "default")
    colr.palette  <-
    structure(sedproxy::stages.key$plotting.colour,
              .Names = sedproxy::stages.key$stage)

  if (alpha.palette[1] == "default") alpha.palette  <-
    structure(sedproxy::stages.key$plotting.alpha,
              .Names = sedproxy::stages.key$stage)

  if (levl.labels[1] == "default") levl.labels  <-
    structure(sedproxy::stages.key$label,
              .Names = sedproxy::stages.key$stage)

  cali.attr <- attr(PFMs, "calibration.pars")

  if (is.null(cali.attr)) {
    cali.attr <- list(calibration.type = "identity")
  }

  if (plot.stages[1] == "default") {
    if (cali.attr$calibration.type == "identity"){
      plotting.levels <- c(
        "clim.signal.monthly", "clim.signal.smoothed", "proxy.bt", "proxy.bt.sb",
        "proxy.bt.sb.sampYM",  "simulated.proxy", "observed.proxy"
      )}else{
        plotting.levels <- c(
          "clim.signal.monthly", "clim.signal.smoothed", "proxy.bt", "proxy.bt.sb",
          "proxy.bt.sb.sampYM",  "simulated.proxy",  "simulated.proxy.cal.err", "reconstructed.climate", "observed.proxy"
        )
      }
  } else if (plot.stages == "all") {
    plotting.levels <- sedproxy::stages.key$stage
    plotting.levels <- subset(plotting.levels, plotting.levels %in% c("clim.signal.ann", "clim.timepoints.ssr") == FALSE)
  } else{
    plotting.levels <- plot.stages
  }

  PFMs <- dplyr::filter(PFMs, stage %in% plotting.levels,
                        replicate <= max.replicates)


  # match scaling flag
  PFMs <- dplyr::left_join(PFMs, sedproxy::stages.key[, c("stage", "scale")])

  #set factor level ordering for stages
  stage.order <- match.arg(stage.order)
  switch(stage.order,
         seq = PFMs$stage <- factor(PFMs$stage, levels = plotting.levels, ordered = TRUE),
         var = {
           var.order <- tapply(PFMs$value, PFMs$stage, FUN = stats::var)
           var.order <- rank(var.order, ties.method = "first")
           var.order <- names(sort(var.order, decreasing = TRUE))
           PFMs$stage <- factor(PFMs$stage,
                                levels = var.order, ordered = TRUE)
         })


  p <- ggplot2::ggplot(data = PFMs, aes(x = .data$timepoints, y = .data$value,
                                        colour = stage, alpha = stage,
                                        linetype = as.factor(replicate))) +
    geom_line() +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), legend.position = "top") +

    guides(colour = guide_legend(
      ncol = 2,
      override.aes = list(alpha = 1))) +
    labs(x = expression("Timepoints"),
         y = expression("Proxy value")) +
    scale_linetype_manual(values = rep(1, 13*length(unique(PFMs$replicate))), guide = "none")+
    scale_alpha_manual(guide = "none")

  pal.df <- data.frame(
    colr.palette = colr.palette,
    colr.breaks = names(colr.palette),
    labels = levl.labels,
    alpha.palette = alpha.palette,
    alpha.breaks = names(alpha.palette)
  )

  pal.df <- dplyr::filter(pal.df, .data$colr.breaks %in% unique(PFMs$stage))

  if (is.null(colr.palette) == FALSE)
    p <- p + scale_colour_manual("", values = pal.df$colr.palette, breaks = pal.df$colr.breaks,
                                 labels = pal.df$labels)

  if (is.null(alpha.palette) == FALSE)
    p <- p + scale_alpha_manual("", values = pal.df$alpha.palette, breaks = pal.df$alpha.breaks,
                                labels = pal.df$labels)

  if (cali.attr$calibration.type != "identity"){
    p <- p + #facet_wrap(~scale, scales = "free_y") +
      facet_wrap( ~ scale, strip.position = "left", scales = "free_y") +
      theme(
        # remove the default y-axis title, "wt"
        axis.title.y = element_blank(),
        # replace the strip backgrounds with transparent
        strip.background = element_rect(fill = 'transparent', colour = 'transparent'),
        strip.placement = 'outside')
  }

  return(p)
}







#' Bioturbation weights
#' @description For a given focal depth (or time), this function returns the probability
#' that material collected from that depth was originally deposited at depth(s)
#' z. In other words, that the material would have been found at depth z if there
#' had been no bioturbation. It is the convolution of the depth solution from
#' Berger and Heath (1968) with a uniform distribution to account for the width
#' of the sediment layer from which samples
#' were picked/extracted. It is a probability density function.
#' @inheritParams ClimToProxyClim
#' @param z A vector of times or depths at which to evaluate the bioturbation weights
#' @param focal.z The depth (or time) for which source dates are wanted
#' @param scale Whether to scale depths by sediment accumulation rate to give
#' positions in terms of time. Defaults to time.
#' @return a numerical vector of weights.
#' @export
#' @references Berger, W. H., & Heath, G. R. (1968).
#' Vertical mixing in pelagic sediments.
#' Journal of Marine Research, 26(2), 134–143.
#' @examples
#' z <- 0:10000
#' w <- BioturbationWeights(z, focal.z = 4000, layer.width = 1, sed.acc.rate = 5, bio.depth = 10)
#' plot(z, w, "l")
BioturbationWeights <- function(z, focal.z, layer.width=1, sed.acc.rate, bio.depth, scale = c("time", "depth")){

  sed.acc.rate <- sed.acc.rate / 1000

  scale <- match.arg(scale)

  if (scale == "depth"){
    z <- z / sed.acc.rate
    focal.z <- focal.z / sed.acc.rate
  }


  lwy <- ceiling(layer.width / sed.acc.rate)
  #lwy[lwy < 1] <- 1 #lwy can be 0
  mdy <- ceiling(bio.depth / sed.acc.rate)

  if (lwy == 0 & mdy == 0) lwy <- 1

  C <- lwy/2
  lam <- 1/mdy

  z <- z - focal.z + mdy

  if (mdy <= 1){
    fz <- stats::dunif(z, -C, C)
  }else if (lwy == 0){
    fz <- stats::dexp(z, 1/mdy)
  }else{
    fz <- (z < -C) * 0 +
      (z >= -C & z <= C) * (lam*(1/lam-exp(-lam*C-lam*z)/lam))/(2*C)  +
      (z > C) * (lam*(exp(lam*C-lam*z)/lam-exp(-lam*C-lam*z)/lam))/(2*C)
  }
  if (sum(fz) == 0){fz}else{
    fz <- fz / sum(fz, na.rm = T)
  }

  return(fz)
}
# Objects
stages.key <-
  structure(
    list(
      stage = c(
        "timepoints",
        "clim.signal.ann",
        "clim.signal.smoothed",
        "clim.timepoints.ssr",
        "proxy.bt",
        "proxy.bt.sb",
        "proxy.bt.sb.inf.b",
        "proxy.bt.sb.inf.b.n",
        "proxy.bt.sb.sampY",
        "proxy.bt.sb.sampYM",
        "proxy.bt.sb.sampYM.b",
        "proxy.bt.sb.sampYM.b.n",
        "simulated.proxy",
        "simulated.proxy.cal.err",
        "reconstructed.climate",
        "observed.proxy"
      ),
      label = c(
        "Requested timepoints",
        "(1) Input climate",
        "(1) Input climate",
        "(1) Input climate",
        "(2) +Bioturbation",
        "(3) +Habitat bias",
        "(.) +Bias",
        "(5) +Ind. error",
        "(4) +Aliasing Y",
        "(4) +Aliasing YM",
        "(.) +Bias",
        "(5) +Ind. error",
        "(5) +Ind. error",
        "(6) +Calibration uncertainty",
        "(7) Reconstructed climate",
        "(*) Observed proxy"
      ),
      description = c(
        "Requested timepoints",
        "Input climate signal at requested timepoints at annual resolution",
        "Input climate signal at regular time intervals and resolution = plot.sig.res",
        "Input climate signal at requested timepoints, smoothed to resolution = plot.sig.res",
        "Climate signal after bioturbation",
        "Climate signal after bioturbation and habitat bias",
        "Climate signal after bioturbation, habitat bias, and bias",
        "Climate signal after bioturbation, habitat bias, and measurement error",
        "Climate signal after bioturbation, habitat bias, and aliasing of inter-annual variation",
        "Climate signal after bioturbation, habitat bias, and aliasing of inter-annual and intra-annual variation such as monthly temperatures or depth habitats",
        "Climate signal after bioturbation, habitat bias, and aliasing of inter-annual and intra-annual variation such as monthly temperatures or depth habitats, and bias",
        "Climate signal after bioturbation, habitat bias, aliasing, and measurement error",
        "Final simulated pseudo-proxy, this will be same as proxy.bt.sb.inf.b.n when n.samples = Inf, and proxy.bt.sb.sampYM.b.n when n.samples is finite",
        "Final simulated pseudo-proxy + uncertainty in true relationship to the climate variable (calibration error)",
        "Final pseudo-proxy calibrated to temperature",
        "True observed proxy (when supplied)"
      ),
      scale = c(
        "time",
        "climate",
        "climate",
        "climate",
        "proxy",
        "proxy",
        "proxy",
        "proxy",
        "proxy",
        "proxy",
        "proxy",
        "proxy",
        "proxy",
        "proxy",
        "climate",
        NA
      ),
      plot.order = c(1, 1, 2, 3,
                     4, 5, 10, 6, 7, 8, 9, 11, 12, 13, 14, 15),
      plotting.colour = c(
        "Black",
        "#018571",
        "#018571",
        "#018571",
        "Green",
        "Gold",
        "Pink",
        "#7570b3",
        "#d95f02",
        "#d95f02",
        "Pink",
        "#7570b3",
        "#7570b3",
        "Red",
        "Blue",
        "Red"
      ),
      plotting.alpha = c(1, 1, 1, 1, 1, 1, 1, 1, 0.5, 0.5,
                         0.5, 0.5, 1, 1, 0.5, 0.5)), .Names = c("stage", "label", "description",
                                                                                                                                                                                                            "scale", "plot.order", "plotting.colour", "plotting.alpha"), row.names = c(NA,
                                                                                                                                                                                                                                                                                   -16L), class = c("tbl_df", "tbl", "data.frame"))
# Functions
SimPowerlaw <- function(beta, N)
{
  N2 <- (3 ^ ceiling(log(N, base = 3)))
  df  <- 1 / N2
  f <- seq(from = df, to = 1 / 2, by = df)
  Filter <- sqrt(1 / (f ^ beta))
  Filter <- c(max(Filter), Filter, rev(Filter))
  x   <- rnorm(N2, 1)
  fx  <- fft(x)
  ffx <- fx * Filter
  result <- Re(fft(ffx, inverse = TRUE))[1:N]
  return(scale(result)[1:N])
}

# Define UI ----
ui <- fluidPage(
  titlePanel("sedproxy"),
  p(
    em("sedproxy"),
    "provides a forward model for sediment archived climate proxies. It is based
    on work described in Laepple and Huybers (2013). A manuscript is in review
    at Climate of the Past Discussions, which more fully describes the forward
    model and its applications", a("(Dolman and Laepple, in review).",
    href = "https://www.clim-past-discuss.net/cp-2018-13/", target = "_blank"),"
    Please contact Dr Andrew Dolman <andrew.dolman@awi.de>,
    or Dr Thomas Laepple <tlaepple@awi.de>, at the Alfred-Wegener-Institute,
    Helmholtz Centre for Polar and Marine Research,
    Germany, for more information.
    "
  ),
  p("This work was supported by German Federal Ministry of Education and Research
    (BMBF) as Research for Sustainability initiative", a("(FONA)", href = "https://www.fona.de/",  target = "_blank"),
    "through the", a("PalMod", href = "https://www.palmod.de/",  target = "_blank"), "project (FKZ: 01LP1509C).",
    br(), br(),
    a(img(src="PalMod_Logo_RGB.png", align = "top"),
      href = "https://www.palmod.de/",
      target = "_blank")), br(),
  p(
    "Reference: ",
    "Laepple, T., & Huybers, P. (2013): Reconciling discrepancies between Uk37
    and Mg/Ca reconstructions of Holocene marine temperature variability.
    Earth and Planetary Science Letters, 375: 418-429."
  ),
  sidebarPanel(tabsetPanel(
    tabPanel(
      "Model parameters",
      fluidRow(
        h4(
          "Update the parameter values below
          and then run the proxy forward model."
        ),
        column(12,
               actionButton("run.pfm", "Run forward model")),
        hr()
        ),
      fluidRow(
        h4("Setup input climate signal"),
        column(
          width = 12,
          sliderInput(
            "clim.signal.length",
            h5("Length of input climate signal [years]"),
            value = 25000,
            step = 1000,
            min = 5000,
            max = 100000
          )
        ),
        column(
          width = 6,
          numericInput(
            "clim.signal.beta",
            h5("Slope of the power spectrum of the input climate signal"),
            value = 1,
            step = 0.1,
            min = 0.1,
            max = 3
          )
        ),
        column(
          width = 6,
          numericInput(
            "seas.amp",
            h5("Amplitude of the seasonal cycle"),
            value = 5,
            step = 0.5,
            min = 0,
            max = 20
          )
        )
      ),
      fluidRow(
        h4("Control sampling"),
        column(width = 6,
               numericInput(
                 "seed",
                 h5("Set RNG seed"),
                 value = 1,
                 step = 1,
                 min = 1
               )),
        column(
          width = 6,
          numericInput(
            "n.replicates",
            h5("No. replicates"),
            value = 1,
            step = 1,
            min = 1,
            max = 100
          )
        )
      ),
      fluidRow(
        column(
          width = 6,
          numericInput(
            "t.res",
            h5("Core sampling resolution [years]"),
            value = 100,
            step = 100,
            min = 1,
            max = 10000
          )
        ),
        column(
          width = 6,
          numericInput(
            "n.samples",
            h5("No. samples per timepoint"),
            value = 30,
            step = 1,
            min = 1,
            max = 1000
          )
        )
      ),
      fluidRow(
        h4("Sedimentation parameters"),
        column(
          6,
          numericInput(
            "bio.depth",
            h5("Bioturbation depth [cm]"),
            value = 10,
            step = 1,
            min = 0,
            max = 30
          )
        ),
        column(
          6,
          numericInput(
            "sed.acc.rate",
            h5("Sediment accumulation rate [cm/ka]"),
            value = 50,
            step = 1,
            min = 0,
            max = 100
          )
        )
      ),
      fluidRow(h4("Proxy production weights (monthly)"),
               column(12,
                      fluidRow(
                        column(
                          12,
                          radioButtons(
                            "seas",
                            label = NULL,
                            choices = c("Uniform", "Custom"),
                            selected = "Uniform",
                            inline = TRUE
                          ),
                          conditionalPanel(
                            condition = "input.seas == 'Custom'",
                            textInput(
                              "mon.vec",
                              "Modify the 12 monthly weights",
                              "1,1,1,1,1,1,1,1,1,1,1,1"
                            ),
                            span(textOutput("habitat.weights.check"), style = "color:red")
                          )
                        )
                      ))),
      fluidRow(
        h4("Noise parameters"),
        column(
          6,
          numericInput(
            "sigma.meas",
            h5("Measurement noise"),
            value = 0.26,
            step = 0.01,
            min = 0,
            max = 1
          )
        ),
        column(
          6,
          numericInput(
            "sigma.ind",
            h5("Individual noise"),
            value = 2,
            step = 0.1,
            min = 0,
            max = 2
          )
        )),
        fluidRow(
        column(
          6,
          numericInput(
            "meas.bias",
            h5("Measurement bias"),
            value = 0,
            step = 0.1,
            min = 0,
            max = 2
          )
        )
      )

  ),
  tabPanel("Plot appearance",
           fluidRow(
             h4("Plot proxy stages:"),
             checkboxGroupInput(
               "stages",
               "Stages:",
               choices = list(
                 "Input climate" = "clim.signal.smoothed",
                 "Bioturbated climate" = "proxy.bt",
                 "Bioturbated + seasonally biased climate" =  "proxy.bt.sb",
                 "Sampled climate inc. aliasing" = "proxy.bt.sb.sampYM",
                 "Final pseudo-proxy" = "simulated.proxy"
               ),
               selected = list(
                 "clim.signal.smoothed",
                 "proxy.bt",
                 "proxy.bt.sb",
                 "proxy.bt.sb.sampYM",
                 "simulated.proxy"
               )
             )
           ))
    )),
  mainPanel(tabsetPanel(
    tabPanel("Plots",
             plotOutput("pfm.plot", height = "800px")),
    tabPanel("Numbers",
             dataTableOutput("pfm.str")),
    tabPanel("Placeholder", textOutput("habitat.weights"))
  ))
  )


# Define server logic ----
server <- function(input, output) {
  clim <- eventReactive(input$run.pfm, {
    set.seed(input$seed)
    ann <-
      SimPowerlaw(input$clim.signal.beta, input$clim.signal.length)
    mon <-
      cos(seq(pi, 3 * pi, length.out = 12)) * input$seas.amp / 2
    clim <- outer(ann, mon, "+")
    clim <- ts(clim, start = 1)
    return(clim)
  }, ignoreNULL = FALSE)
  timepoints <- eventReactive(input$run.pfm, {
    #res <- 100
    tp <- seq(1, input$clim.signal.length, by = input$t.res)
    # t.min <-
    #   ceiling(1000 * input$bio.depth / input$sed.acc.rate) + 1
    t.min <- min(tp)
    t.max <- input$clim.signal.length - 3 * t.min
    tp <- tp[tp > t.min & tp < t.max]
    return(tp)
  }, ignoreNULL = FALSE)
  seasprod <- eventReactive({
    input$mon.vec
    input$seas
  }, {
    if (input$seas == 'Custom')
    {
      v <- as.numeric(unlist(strsplit(input$mon.vec, ",")))
    } else{
      v <- rep(1, 12)
    }
    return(v)
  }, ignoreNULL = FALSE)
  output$habitat.weights.check <- renderText({
    if (length(seasprod()) != 12)
    {
      paste0("You entered ",
             length(seasprod()),
             " values; 12 are required.")
    }
  })
  pfm <- eventReactive(input$run.pfm, {
    pfm <- ClimToProxyClim(
      clim.signal = clim(),
      timepoints = timepoints(),
      plot.sig.res = 100,
      bio.depth = input$bio.depth,
      sed.acc.rate = input$sed.acc.rate,
      habitat.weights = seasprod(),
      n.samples = input$n.samples,
      n.replicates = input$n.replicates,
      sigma.meas = input$sigma.meas,
      sigma.ind = input$sigma.ind,
      meas.bias = input$meas.bias
    )
  }, ignoreNULL = FALSE)
  output$pfm.plot <- renderPlot({
    dat <- pfm()$everything
    #cal.pars <- attr(dat, "calibration.pars")
    #dat <- subset(dat, dat$stage %in% input$stages)
    #attr(dat, "calibration.pars") <- cal.pars
    if (nrow(dat) > 0) {
      PlotPFMs(dat) +
        ggplot2::labs(x = "Age [years]")
    }
  }, res = 72 * 2)
  output$pfm.str <- renderDataTable({
    round(pfm()$simulated.proxy, 5)
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)

