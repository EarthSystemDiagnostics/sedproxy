#' Resample a seasonal cycle simulating the seasonal contribution
#' @description Calculate the (possibly stochastic) bias in the proxy record,
#' due to seasonality in the abundance of the sampled organism.
#' @inheritParams SubsampleTimeseries_MgCa
#' @param sst.clim ts object containing the monthly SST climatology (12 values)
#' @param proxy.clim ts object containing the monthly Foram population index (12 values)
#'
#' @return A list. result=list(time, data)
#'
#' data is a single value for each timepoint, indicating the offset in the proxy
#' temperature record due to the seasonal abundance pattern of the organism,
#' which preferentially samples those periods of the year when it is most
#' abundant. There is a stochastic component to this due to measuring a finite
#' number of individuals. If Nforam is set to \code{Inf}, then these biases will be
#' constant for a given sst.clim and proxy.clim.
#' @export
#'
#' @examples
#'
SeasonalityBias <-
  function(sst.clim, timepoints, proxy.clim, NForam = 30)
  {
    result <- list()
    sst.seasonal <- vector()
    sst.clim <- scale(sst.clim, scale = FALSE)
    # for (i in 1:length(timepoints))
    #   sst.seasonal[i] <- mean(sample(sst.clim, size = NForam, prob = proxy.clim, replace = TRUE))
    if (is.infinite(NForam)) {
      val <- sum(sst.clim * proxy.clim) / sum(proxy.clim)
      sst.seasonal <- rep(val, length(timepoints))
    }
    if (is.finite(NForam)) {
      sst.seasonal <-
        replicate(length(timepoints), mean(
          sample(
            sst.clim,
            size = NForam,
            prob = proxy.clim,
            replace = TRUE
          )
        ))
    }
    result$time <- timepoints
    result$data <- sst.seasonal
    return(result)
  }
