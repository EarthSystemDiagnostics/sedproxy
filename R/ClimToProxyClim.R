#' Title
#'
#' @param clim.signal The "true" climate signal, i.e. model output or instrumental record. A years x 12 (months) matrix of temperatures.
#' @param timepoints The timepoints at which the proxy is to be sampled
#' @param seas.prod The seasonal pattern of productivity for the organism(s) archived in the proxy. Either a vector of 12 values or a matrix of the same dimensions as clim.signal
#' @param bio.depth Depth of the bioturbated layer in metres, defaults to 0.1 m
#' @param acc.rate rate of sediment accumulation (sedimentation rate) in metres
#'   per year. Defaults to 0.01 m (1 cm per year).
#' @param n_samples Number of e.g. foraminifera sampled per timepoint
#'
#' @return
#' @export
#'
#' @examples
ClimToProxyClim <- function(clim.signal,
                            timepoints,
                            seas.prod,
                            bio.depth = 0.1,
                            acc.rate = 0.012,
                            n_samples = Inf) {
  stopifnot(is.matrix(clim.signal))

  sig.years.i <- 1:nrow(clim.signal)

  # Get bioturbation window ----------
  time.step <- 1
  bio.depth.timesteps <- round(bio.depth / acc.rate)

  # width of the bioturbation window (transfer function)
  # 4x biowidth contains most of the area of the impulse response
  bioturb.window <- seq(
    from = -1 * bio.depth.timesteps,
    to = 3 * bio.depth.timesteps,
    by = 1
  )

  # Get bioturbation weights --------
  # use inverted bioturb.window, then reverse order of weights
  bioturb.weights <- ImpulseResponse(-bioturb.window, bio.depth.timesteps, z0 = 0) #impulse response in year unit
  bioturb.weights <- (bioturb.weights / sum(bioturb.weights))

  # Check depth and time order match
  #plot(bioturb.window, bioturb.weights)

  # Get weights matrix ---------
  clim.sig.weights <- bioturb.weights %o% seas.prod
  clim.sig.weights <- clim.sig.weights / sum(clim.sig.weights)

  # Check weights sum to 1, within tolerance
  stopifnot(abs(sum(clim.sig.weights) - 1) < 1e-04)

  clim.sig.windows <- lapply(timepoints, function(i) {
    x <- bioturb.window + i
    x <- x[x > 0]
    clim.signal[x, ]
    })

  proxy.sig.inf <- sapply(clim.sig.windows, function(x) sum(clim.sig.weights * x))

if (is.infinite(n_samples)) {
  proxy.sig.samp <- NA
  } else if (is.finite(n_samples)) {
    proxy.sig.samp <- sapply(clim.sig.windows, function(x)
      mean(sample(x, n_samples, prob = clim.sig.weights, replace = TRUE))
    )
    }

return.obj <- list(timepoints = timepoints, proxy.sig.inf=proxy.sig.inf, proxy.sig.samp=proxy.sig.samp)

return(return.obj)

}
