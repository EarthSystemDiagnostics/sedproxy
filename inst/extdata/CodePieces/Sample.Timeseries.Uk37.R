Sample.Timeseries.Uk37<-function(tss,data.raw,bioWidth.depth=0.1)
{
# Sample an irregular timeseries from an annual timeseries mimicing the Uk37
# measurement process
# Inputs:
# tss: annual timeseries
# data.raw      : list(time,data,lat,lon) containing the raw timeseries
# bioWidth.depth: Bioturbation in m
# optionally: add noise with sd=noise on the measurement level
# Returns
#         list(time,data,lat,lon)
#
biowidth.time=round(bioWidth.depth/data.raw$meanAccum*1000) #Bioturbation width in time domain in yr unit
timepoints<-data.raw$time #timepoints at which we sample

z<-(-3*biowidth.time):(3*biowidth.time) #width of the transfer function (compromise 3x biowidth on each side contains most of the area

if (length(z) >= length(tss)) {
  z<-(-1.3*biowidth.time):(1.3*biowidth.time)
  warning("Bioturbation filter was approximated with a shorter filter as the timeseries is too low")
}
weights<-ImpulseResponse.BergerHeath(z,biowidth.time,0) #impulse response in year unit
weights<-weights/sum(weights)

    bioturbated<-filter(tss,weights,circular=TRUE)
    index.time<- (timepoints-c(time(tss[1])))*1000+1
    index.time<-index.time[index.time>0]
    result.data<-bioturbated[index.time]
    result.time<-c(time(tss))[index.time]

#Bring them together
    sampled<-list()
    index<-!is.na(result.data)
    sampled$data<-result.data[index]
    sampled$time<-result.time[index]
    sampled$lat<-data.raw$lat
    sampled$lon<-data.raw$lon
    sampled$name<-data.raw$name
    return(sampled)
}

