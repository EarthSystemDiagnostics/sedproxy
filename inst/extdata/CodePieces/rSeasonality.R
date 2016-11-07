
rSeasonality<-function(n=1000,seasonality)
### Draw random number between 1 and 12 with the distribution
### from seasonality
#Input
#n Number of random numbers to be drawn
#seasonality= vector of 12 numbers; defines the distribution
#based on rWeights
{

    return(rWeights(n,seasonality))
}
