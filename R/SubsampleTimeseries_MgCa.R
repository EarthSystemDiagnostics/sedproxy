#' SubsampleTimeseries_MgCa
#'
#' @param tss a pTs object containing an equidistant timeseries, (e.g. model
#'   result or interpolated time series), units = kyr
#' @param timepoints the timepoints at which to sample the time series
#' @param biowidth.time width of the averaging window in years
#' @param NForam Number of foraminifera sampled per timepoint
#' @description Resamples an equidistant sampled timeseries (e.g. model result)
#'   at the specified "timepoints" by simulating the foraminifera sampling
#' @return a list list(time, data)
#' @export SubsampleTimeseries_MgCa
#'
#' @examples
#' tss <- ts(c(rep(0, 100), 1, rep(0, 99)), (1:200) / 1000)
#' timepoints <- c(time(tss))
#' par(mfcol = c(1, 2))
#' plot(tss)
#' biowidth.time = 10
#' result <- SubsampleTimeseries_MgCa(tss, timepoints, biowidth.time, 100)
#' plot(result$time, result$data, type = "l")
#'
#' tss<-ts(sin(seq(-pi, 9*pi, length.out = 1000)), start = 0, frequency = 1/100)
#' timepoints<-c(time(tss)[1:1000])
#' biowidth.time=100
#' result<-SubsampleTimeseries_MgCa(tss,timepoints,biowidth.time, NForam = 30)
#'
#' plot(tss)
#' lines(data~time, data = result, col = "Red")

SubsampleTimeseries_MgCa <-
  function(tss,
           timepoints,
           biowidth.time,
           NForam = 30)

  {
    time_step <- 1/frequency(tss)
    #width of the transfer function
    #(compromise 4x biowidth contains most of the area of the impulse response)
    z <- seq(from = -1 * biowidth.time,
             to = 3 * biowidth.time,
             by = 1)
    #impulse response in year units; use inverted z
    weights <-
      ImpulseResponse(-z, biowidth.time, z0 = 0)

    #Go in the index space; index.time contains the indices of the timepoints as index# in tss
    index.time <- (timepoints - c(time(tss[1]))+1) / time_step
    #Only use points inside the timeseries
    index.time <-
      index.time[index.time > 0]
    N.tss <- length(tss) #length of the timeseries
    index.time <- index.time[index.time <= N.tss]

    tss <- c(tss) #remove timeseries attributes


     index.time0 <-
       (biowidth.time + 1) #Index of the weight which corresponds to time 0

    result <- list()
    result$data <- rep(NA, length(timepoints))
    for (i in 1:length(index.time))
    {
      shift <-
        index.time[i] - index.time0
      #Shift of the weights that zero time overlaps with the sampling time

      #Put weights outside of the timeseries to zero
      weights.cut <- weights
      index.weight.shift <-
        (1:length(weights)) + shift #index to simplify the
      weights.cut[(index.weight.shift < 1)] <- 0
      weights.cut[(index.weight.shift) > N.tss] <- 0

      #create random numbers with the probabilities of the weights
      rIndex <- rWeights(NForam, weights.cut) + shift
      #print(rIndex)
      result$data[i] <- mean(tss[rIndex])
    }
    result$time <- timepoints
    return(result)
  }
