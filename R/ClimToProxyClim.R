#' Title
#'
#' @param clim.signal The "true" climate signal, i.e. model output or
#'   instrumental record. A years x 12 (months) matrix of temperatures.
#' @param timepoints The timepoints for which the proxy record is to be modelled
#' @param seas.prod The seasonal pattern of productivity for the organism(s)
#'   archived in the proxy. Either a vector of 12 values or a matrix of the same
#'   dimensions as clim.signal. Defaults to uniform seasonal distribution.
#' @param bio.depth Depth of the bioturbated layer in metres, defaults to 0.1 m
#' @param acc.rate Sediment accumulation rate in metres per year.
#'   Defaults to 5e-04 m per year (0.5 m / kyr). Either a single value, or vector
#'    of same length as "timepoints"
#' @param n.samples Number of e.g. foraminifera sampled per timepoint
#'
#' @return
#' @export
#'
#' @examples
ClimToProxyClim <- function(clim.signal,
                            timepoints,
                            seas.prod = rep(1, 12),
                            bio.depth = 0.1,
                            acc.rate = 5e-04,
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

  if (length(acc.rate) != n.timepoints) {
    acc.rate <- rep(acc.rate, n.timepoints)
  }

  # For each timepoint ------

  # # Prepare result objects
  #
  # proxy.sig.inf <- rep(NA, n.timepoints)
  # proxy.sig.samp <-
  #   matrix(rep(NA, n.timepoints * n.replicates), ncol = n.replicates)


  proxy.sig <- sapply(1:length(timepoints), function(x) {
    # Get bioturbation window ----------
    bio.depth.timesteps <- round(bio.depth / acc.rate[x])
    bioturb.window <- GetBioturbWindow(bio.depth.timesteps)

    # Get bioturbation weights --------
    bioturb.weights <-
      ImpulseResponse(-bioturb.window, bio.depth.timesteps, z0 = 0)
    bioturb.weights <- bioturb.weights / sum(bioturb.weights)

    # Check depth and time order match
    # plot(bioturb.weights, (bioturb.window), type = "l", ylim = rev(range(bioturb.window)))

    # get portion of clim.signal corresponding to bioturbation window -------

    sig.window.i.1 <- bioturb.window + timepoints[x]
    sig.window.i <- sig.window.i.1[sig.window.i.1 > 0]
    clim.sig.window <- clim.signal[sig.window.i,]

    # Get weights matrix ---------
    clim.sig.weights <- bioturb.weights %o% seas.prod
    clim.sig.weights <- clim.sig.weights[sig.window.i.1 > 0,]
    clim.sig.weights <- clim.sig.weights / sum(clim.sig.weights)

    # Check weights sum to 1, within tolerance
    stopifnot(abs(sum(clim.sig.weights) - 1) < 1e-10)


    # Calculate mean clim.signal -------
    proxy.sig.inf <- sum(clim.sig.weights * clim.sig.window)

    if (is.infinite(n.samples)) {
      proxy.sig.samp <- NA
    } else if (is.finite(n.samples)) {
      proxy.sig.samp <- replicate(n = n.replicates, expr = {
        mean(sample(
          clim.sig.window,
          n.samples,
          prob = clim.sig.weights,
          replace = TRUE
        ))
      })
    }
    list(
      timepoints = timepoints[x],
      proxy.sig.inf = proxy.sig.inf,
      proxy.sig.samp = proxy.sig.samp
    )
  }, simplify = FALSE)

  # return.obj <-
  #   list(
  #     timepoints = timepoints,
  #     proxy.sig.inf = proxy.sig.inf,
  #     proxy.sig.samp = proxy.sig.samp
  #   )

#list(timepoints=timepoints, proxy.sig.inf=proxy.sig.inf, proxy.sig.samp=proxy.sig.samp)
return(proxy.sig)
}


GetBioturbWindow <- function(bio.depth.timesteps) {
  # width of the bioturbation window (transfer function)
  # 4x biowidth contains most of the area of the impulse response
  bioturb.window <- seq(
    from = -1 * bio.depth.timesteps,
    to = 3 * bio.depth.timesteps,
    by = 1
  )
}
