Sample.Timeseries.MgCa<-function(tss,data.raw,bioWidth.depth=0.1,return.both=FALSE)
{
# Sample an irregular timeseries from an annual timeseries; e.g. climate model output simulating the foram sample process
#
# Inputs:
# rss: annual timeseries
# data.raw      : list(time,data,lat,lon) containing the raw timeseries
# bioWidth.depth: Bioturbation in m
# return.both: If true; separatly return the interannual and the seasonal sampling
# Returns
#         list(time,data,lat,lon)
#
# Uses
#             subsampleTimeseries.MgCa
#              subsampleTimeseries.Seasonality


sampled.interannual<-subsample.Timeseries.MgCa(tss,data.raw$time,biowidth.time=bioWidth.depth/data.raw$meanAccum*1000,NForam=data.raw$other)

sampled.seasonality<-subsampleTimeseries.Seasonality(sst.clim=data.raw$sst.clim,timepoints=data.raw$time,proxy.clim=data.raw$foram.clim,NForam=data.raw$other)

#Bring them together
    sampled<-list()
    index<-!is.na(sampled.interannual$data)
    sampled$data<-sampled.interannual$data[index]+sampled.seasonality$data[index]
    sampled$time<-sampled.interannual$time[index]
    sampled$lat<-data.raw$lat
    sampled$lon<-data.raw$lon
    sampled$name<-data.raw$name
    sampled$ts.annual<-tss

if (return.both)
{
    sampled.woSeason<-sampled
    sampled.woSeason$data<-sampled.interannual$data[index]
    return(list(both=sampled,woSeason=sampled.woSeason))
} else return(sampled)
}
