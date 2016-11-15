#' @title BioturbateTimeseries
#' @name BioturbateTimeseries
#' @param tss a ts object containing an equidistant timeseries, (e.g. model
#'   result or interpolated time series)
#' @param timepoints the timepoints at which to sample the time series
#' @param biowidth_timesteps width of the averaging window in no. of time points
#' @param n_samples Number of e.g. foraminifera sampled per timepoint
#' @description Resamples an equidistant sampled timeseries (e.g. model result)
#'   at the specified "timepoints" by simulating the foraminifera sampling
#' @return a list containing the resampled data values and times
#' @export BioturbateTimeseries
#'
#' @examples
#' tss <- ts(c(rep(0, 100), 1, rep(0, 99)), 0:199))
#' timepoints <- c(time(tss))
#' par(mfcol = c(1, 2))
#' plot(tss)
#' biowidth_timesteps = 10
#' result <- BioturbateTimeseries(tss, timepoints, biowidth_timesteps, n_samples = 100)
#' plot(result$time, result$data, type = "l")
#' par(mfcol = c(1, 1))
#'
#' tss<-ts(sin(seq(-pi, 9*pi, length.out = 1000)), start = 0, frequency = 1/100)
#' timepoints<-c(time(tss)[1:1000])
#' biowidth_timesteps=100
#' result<-BioturbateTimeseries(tss,timepoints,biowidth_timesteps, n_samples = 30)
#'
#' plot(tss)
#' lines(data~time, data = result, col = "Red")

BioturbateTimeseries <- function(tss,
                                 timepoints,
                                 biowidth_timesteps,
                                 n_samples = Inf)
{
  time_step <- 1 / frequency(tss)
  #width of the transfer function
  #(compromise 4x biowidth contains most of the area of the impulse response)
  z <- seq(
    from = -1 * biowidth_timesteps,
    to = 3 * biowidth_timesteps,
    by = 1
  )
  #Go in the index space; index.time contains the indices of the timepoints as index# in tss
  index.time <- (timepoints - time(tss)[1]) / time_step + 1
  #Only use points inside the timeseries
  index.time <-
    index.time[index.time > 0]
  N.tss <- length(tss) #length of the timeseries
  index.time <- index.time[index.time <= N.tss]

  #Index of the weight which corresponds to time 0
  index.time0 <-
    (biowidth_timesteps + 1)

  result <- list()
  result$data <- rep(NA, length(index.time))

  if (is.infinite(n_samples)) {
    if (length(z) >= length(tss)) {
      z <- (-1.3 * biowidth_timesteps):(1.3 * biowidth_timesteps)
      warning(
        "Bioturbation filter was approximated with a shorter filter as the timeseries is too low"
      )
    }
    # use inverted z, then reverse order of weights
    weights <-
      ImpulseResponse(z, biowidth_timesteps, z0 = index.time0-1) #impulse response in year unit
    weights <- (weights / sum(weights))

    bioturbated <- filter(tss, weights, circular = TRUE)
    result$data <- bioturbated[index.time]
    result$time <- c(time(tss))[index.time]
    } else if (is.finite(n_samples)) {
    for (i in 1:length(index.time))
    {
      #impulse response in year units; use inverted z
      # 1*biowidth into existing sediment, 3*biowidth into "future" sediment
      weights <-  ImpulseResponse(-z, biowidth_timesteps, z0 = 0)

      shift <- index.time[i] - index.time0
      #Shift of the weights that zero time overlaps with the sampling time

      #Put weights outside of the timeseries to zero
      weights.cut <- weights
      index.weight.shift <-
        (1:length(weights)) + shift #index to simplify the
      weights.cut[(index.weight.shift < 1)] <- 0
      weights.cut[(index.weight.shift) > N.tss] <- 0

      #create random numbers with the probabilities of the weights
      #rIndex <- rWeights(n_samples, weights.cut) + shift
      rIndex <-
        sample(
          index.weight.shift,
          size = n_samples,
          prob = weights.cut,
          replace = TRUE
        )
      #print(rIndex)
      result$data[i] <- mean(tss[rIndex])
    }
      result$time <- c(time(tss))[index.time]

  }

  return(result)
}
