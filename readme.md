# Sedproxy: Simulation of Sediment Archived Climate Proxy Records.

------------------------------


## Installation

**sedproxy** can be installed directly from bitbucket


```r
if (!require("devtools")) {
  install.packages("devtools")
}

devtools::install_bitbucket("ecus/sedproxy")
```

## Example data

`sedproxy` includes example data for a single sediment core and location: core number 41 (MD97-2141) in Shakun et al. (2012). The climate signal is taken from the [TraCE-21ka](http://www.cgd.ucar.edu/ccr/TraCE/) Simulation of Transient Climate Evolution over the last 21,000 years, using the grid cell closest to core MD97-2141. Seasonality of *G.ruber*, the Foraminifera for which test Mg/Ca ratios were measured, is taken from the model of Fraile et al (2008). Sediment accumulation rates were estimated from the depth and age data associated with core MD97-2141, with a minimum rate of 0.2 * the mean rate.


**The MD97-2141 core**


```r
library(tidyverse)
library(knitr)
library(sedproxy)
```



```r
knitr::kable(N41.proxy.details %>% tidyr::gather())
```



key               value                                            
----------------  -------------------------------------------------
Number            41                                               
ID.no             N41                                              
Core              MD97-2141                                        
Location          Sulu Sea                                         
Proxy             Mg/Ca                                            
Lat               8.78333                                          
Lon               121.2833                                         
Elevation         -3633.000000                                     
Reference         Rosenthal et al., 2003                           
Resolution        77.8894472361809                                 
Calibration.ref   Rosenthal and Lohman, 2002                       
Calibration       T = ln(MgCa/0.28)/0.095                          
Foram.sp          G. ruber                                         
Ref.14C           de Garidel-Thoron et al., 2001, Paleoceanography 
Notes             NA                                               
Geo.cluster       Sulu Sea                                         
Archive.type      Marine sediment                                  


**Modelled climate signal**



```r
knitr::kable(N41.t21k.climate[1:5,])
```

        1          2          3          4          5          6          7          8          9         10         11         12
---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------
 297.8568   297.3878   298.0602   299.2245   299.6761   300.2167   300.5483   299.9182   299.6400   299.6386   299.9301   299.3381
 297.9886   297.5341   297.8263   299.0083   299.7151   299.6743   300.1609   300.6685   299.7819   299.8497   299.7823   298.9699
 297.8427   297.7466   298.3599   299.1459   299.6056   300.0913   300.1500   300.1423   299.5427   299.5985   299.8105   298.9172
 297.6961   297.6659   298.4921   299.5137   300.0200   299.8971   300.4432   299.9803   299.6973   300.0648   299.7363   298.9907
 297.7691   297.3436   297.9452   299.1741   299.9940   299.8244   300.1444   300.3978   299.9496   300.1610   299.8150   298.9635

**Actual proxy record**



```r
kable(head(N41.proxy))
```



 Published.age   Published.temperature   Sed.acc.rate.m.yr
--------------  ----------------------  ------------------
      4334.286                   28.92           0.0003679
      4527.429                   29.20           0.0003675
      4575.714                   29.15           0.0003677
      4720.571                   28.55           0.0003677
      4913.714                   28.33           0.0003670
      4994.400                   29.44           0.0003667

## Function `ProxyToProxyClim`

`ProxyToProxyClim` is the main function in package `sedproxy`. It simulates a sediment archived proxy from an assumed true climate signal, the sediment accumulation rate, seasonality of the encoding organism/process, and the number of samples per timepoint.



```r
PFM <- ClimToProxyClim(clim.signal = N41.t21k.climate, 
                timepoints = N41.proxy$Published.age,
                seas.prod = N41.G.ruber.seasonality,
                sed.acc.rate = N41.proxy$Sed.acc.rate.m.yr,
                meas.noise = 0.46, n.samples = 30,
                n.replicates = 10)
```

```
## Warning in FUN(X[[i]], ...): Window extends below end of clim.signal
```



```r
names(PFM)
```

```
## [1] "simulated.proxy" "smoothed.signal" "everything"
```


```r
head(PFM$simulated.proxy)
```

```
## # A tibble: 6 × 14
##   timepoints clim.timepoints.100 clim.timepoints.50 proxy.bt proxy.bt.sb
##        <dbl>               <dbl>              <dbl>    <dbl>       <dbl>
## 1   4334.286            299.2912           299.2918 299.2828    299.0861
## 2   4527.429            299.3389           299.3384 299.3079    299.1140
## 3   4575.714            299.3670           299.3686 299.3143    299.1216
## 4   4720.571            299.3536           299.3513 299.3339    299.1414
## 5   4913.714            299.3090           299.3187 299.3384    299.1478
## 6   4994.400            299.2964           299.2696 299.3334    299.1442
## # ... with 9 more variables: sed.acc.rate <dbl>, smoothing.width <dbl>,
## #   proxy.bt.sb.sampY <dbl>, proxy.bt.sb.sampYM <dbl>,
## #   proxy.bt.sb.inf.b <dbl>, proxy.bt.sb.inf.b.n <dbl>,
## #   proxy.bt.sb.sampYM.b <dbl>, proxy.bt.sb.sampYM.b.n <dbl>,
## #   simulated.proxy <dbl>
```

**Simple plotting**


```r
plot.df <- PFM$simulated.proxy %>% 
  select(timepoints, clim.timepoints.50, proxy.bt, proxy.bt.sb,
         proxy.bt.sb.sampYM, proxy.bt.sb.sampYM.b.n) %>% 
  gather(Stage, Temperature, -timepoints) %>% 
  mutate(Age = timepoints,
         Temperature = Temperature - 273.15) 

plot.df %>% 
ggplot(aes(x = Age, y = Temperature, colour = Stage)) +
  geom_line()
```

![](readme_files/figure-html/default_plot-1.png)<!-- -->






**The 10 replicates of the final simulated proxy**


```r
head(PFM$everything$proxy.bt.sb.sampYM.b.n)
```

```
##          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]
## [1,] 298.9181 298.9139 299.2183 299.2184 299.5292 299.6990 298.9461
## [2,] 299.2619 298.3106 299.8199 299.2231 299.1808 298.4447 299.4475
## [3,] 297.8604 298.7895 299.4450 300.1607 299.1840 299.4941 298.4525
## [4,] 298.5160 298.9258 298.4471 299.0891 299.0630 299.7203 299.2809
## [5,] 299.7555 299.9854 299.9712 299.8376 298.6157 299.3750 299.9131
## [6,] 299.2639 298.3424 298.3271 299.1964 298.6711 298.7615 299.0948
##          [,8]     [,9]    [,10]
## [1,] 298.6312 298.4010 299.0853
## [2,] 299.4400 299.7222 299.2758
## [3,] 299.3774 298.2424 298.8785
## [4,] 299.0144 299.5697 299.3177
## [5,] 299.8295 298.6621 298.9932
## [6,] 300.0099 299.8169 299.0914
```



## Literature cited

Fraile, I., M. Schulz, S. Mulitza, and M. Kucera. “Predicting the Global Distribution of Planktonic Foraminifera Using a Dynamic Ecosystem Model.” Biogeosciences 5, no. 3 (2008): 891–911.


Shakun, Jeremy D., Peter U. Clark, Feng He, Shaun A. Marcott, Alan C. Mix, Zhengyu Liu, Bette Otto-Bliesner, Andreas Schmittner, and Edouard Bard. 2012. “Global Warming Preceded by Increasing Carbon Dioxide Concentrations during the Last Deglaciation.” Nature 484 (7392): 49–54. doi:10.1038/nature10915.


