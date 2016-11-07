
subsampleTimeseries.Seasonality<-function(sst.clim,timepoints,proxy.clim,NForam=30)
#Resample a seasonal cycle simulating the seasonal contribution to foraminifera sampling
#Input:
# sst.clim: pTs object containing the monthly SST climatology (12 values)
# timepoints: vector with the points in time
# proxy.clim: pTs object containing the monthly foram population (12 values)
#NForam: Number of foraminifera
#Returns: result=list(time,data)
{
    result<-list()
    sst.seasonal<-vector()
    sst.clim<-scale(sst.clim,scale=FALSE)
    for (i in 1:length(timepoints)) sst.seasonal[i]<-mean(c(sst.clim)[rSeasonality(NForam,c(proxy.clim))])

    result$time<-timepoints
    result$data<-sst.seasonal
    return(result)
}


