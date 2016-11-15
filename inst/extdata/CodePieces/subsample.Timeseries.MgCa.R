## tss<-pTs(c(rep(0,100),1,rep(0,99)),(1:200)/1000)
## timepoints<-c(time(tss)[30:170])
## biowidth.time=0.02*1000
## result<-subsample.Timeseries.MgCa(tss,timepoints,biowidth.time,100)
## result1<-subsample.Timeseries.MgCa.fast(tss,timepoints,biowidth.time,100)
## par(mfcol=c(1,2))
## plot(result$time,result$data)
## plot(result1$time,result1$data)
## lines(sampled$time,sampled$data,col="red")
## abline(v=100/1000)

## result<-Sample.Timeseries.Uk37(tss,timepoints,biowidth.time)

subsample.Timeseries.MgCa <-
  function(tss,
           timepoints,
           biowidth.time,
           NForam = 30)
    #Resample a equidistant sampled timeseries (e.g. model result) at the
    #"timepoints" by simulating the foraminifera sampling Input: tss: pTs object
    #containing the equidistant timeseries, units = kyr timepoints: vector with the
    #points in time biowidth: width of the averaging window in years NForam: Number
    #of foraminifera Returns: result=list(time,data)
  {
    z <-
      (-1 * biowidth.time):(3 * biowidth.time) #width of the transfer function (compromise 4x biowidth contains most of the area of the impulse response
    weights <-
      ImpulseResponse.BergerHeath(-z, biowidth.time, 0) #impulse response in year units; use inverted z

    #Go in the index space; index.time contains the indices of the timepoints as index# in tss
    index.time <- (timepoints - c(time(tss)[1])) * 1000 + 1
    index.time <-
      index.time[index.time > 0] #Only use points inside the timeseries
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
      result$data[i] <- mean(tss[rIndex])
    }
    result$time <- timepoints
    return(result)
  }
