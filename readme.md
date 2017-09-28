# Sedproxy: Simulation of Sediment Archived Climate Proxy Records.

------------------------------

## Introduction

`sedproxy` provides a forward model for sediment archived climate proxies. It is based on work described in Laepple and Huybers (2013). A manuscript is in preparation, Dolman and Laepple (in prep.), which will more fully describe the forward model and its applications. Please contact Dr Andrew Dolman <<andrew.dolman@awi.de>>, or Dr Thomas Laepple <<tlaepple@awi.de>>, at the Alfred-Wegener-Institute, Helmholtz Centre for Polar and Marine Research, Germany, for more information.

 
## Installation

**sedproxy** can be installed directly from bitbucket


```r
if (!require("devtools")) {
  install.packages("devtools")
}

devtools::install_bitbucket("ecus/sedproxy")
```

## Example data

`sedproxy` includes example data for a single sediment core and location: core number 41 in the Shakun et al. (2012) compilation (MD97-2141, Rosenthal et al. 2003). The climate signal is taken from the [TraCE-21ka](http://www.cgd.ucar.edu/ccr/TraCE/) Simulation of Transient Climate Evolution over the last 21,000 years, using the grid cell closest to core MD97-2141. Seasonality of *G.ruber*, the Foraminifera for which test Mg/Ca ratios were measured, is taken from the model of Fraile et al (2008). Sediment accumulation rates were estimated from the depth and age data associated with core MD97-2141, with a minimum rate of 0.2 * the mean rate.


**The MD97-2141 core**


```r
library(tidyverse)
library(knitr)
library(sedproxy)
```



```r
knitr::kable(N41.proxy.details %>% tidyr::gather(), format = "markdown")
```



|key             |value                                            |
|:---------------|:------------------------------------------------|
|Number          |41                                               |
|ID.no           |N41                                              |
|Core            |MD97-2141                                        |
|Location        |Sulu Sea                                         |
|Proxy           |Mg/Ca                                            |
|Lat             |8.78333                                          |
|Lon             |121.2833                                         |
|Elevation       |-3633.000000                                     |
|Reference       |Rosenthal et al., 2003                           |
|Resolution      |77.8894472361809                                 |
|Calibration.ref |Rosenthal and Lohman, 2002                       |
|Calibration     |T = ln(MgCa/0.28)/0.095                          |
|Foram.sp        |G. ruber                                         |
|Ref.14C         |de Garidel-Thoron et al., 2001, Paleoceanography |
|Notes           |NA                                               |
|Geo.cluster     |Sulu Sea                                         |
|Archive.type    |Marine sediment                                  |


**Modelled climate signal**

The first 5 rows:


```r
knitr::kable(N41.t21k.climate[1:5,], format = "markdown", digits = 2)
```



|      1|      2|      3|      4|      5|      6|      7|      8|      9|     10|     11|     12|
|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|
| 297.86| 297.39| 298.06| 299.22| 299.68| 300.22| 300.55| 299.92| 299.64| 299.64| 299.93| 299.34|
| 297.99| 297.53| 297.83| 299.01| 299.72| 299.67| 300.16| 300.67| 299.78| 299.85| 299.78| 298.97|
| 297.84| 297.75| 298.36| 299.15| 299.61| 300.09| 300.15| 300.14| 299.54| 299.60| 299.81| 298.92|
| 297.70| 297.67| 298.49| 299.51| 300.02| 299.90| 300.44| 299.98| 299.70| 300.06| 299.74| 298.99|
| 297.77| 297.34| 297.95| 299.17| 299.99| 299.82| 300.14| 300.40| 299.95| 300.16| 299.82| 298.96|

**Actual proxy record**

Core MD97-2141 (Rosenthal et al. 2003)


```r
kable(head(N41.proxy), format = "markdown")
```



| Published.age| Published.temperature| Sed.acc.rate.m.yr|
|-------------:|---------------------:|-----------------:|
|      4334.286|                 28.92|         0.0003679|
|      4527.429|                 29.20|         0.0003675|
|      4575.714|                 29.15|         0.0003677|
|      4720.571|                 28.55|         0.0003677|
|      4913.714|                 28.33|         0.0003670|
|      4994.400|                 29.44|         0.0003667|

*******

## Function `ClimToProxyClim`

`ClimToProxyClim` is the main function in package `sedproxy`. It simulates a sediment archived proxy from an assumed true climate signal, the sediment accumulation rate, seasonality of the encoding organism/process, and the number of samples per timepoint.



```r
set.seed(26052017)

PFM <- ClimToProxyClim(clim.signal = N41.t21k.climate[nrow(N41.t21k.climate):1,] - 273.15, 
                       timepoints = round(N41.proxy$Published.age),
                       proxy.calibration.type = "identity",
                seas.prod = N41.G.ruber.seasonality,
                sed.acc.rate = N41.proxy$Sed.acc.rate.m.yr,
                meas.noise = 0.46, n.samples = 30,
                n.replicates = 10)
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
## # A tibble: 6 x 17
##   timepoints clim.signal.ann clim.timepoints.1000 clim.timepoints.100
##        <dbl>           <dbl>                <dbl>               <dbl>
## 1       4334        28.02843             27.93745            27.91754
## 2       4527        27.64365             27.92579            27.90224
## 3       4576        28.27605             27.92421            27.89731
## 4       4721        28.18325             27.91341            27.94076
## 5       4914        27.91505             27.89494            27.90284
## 6       4994        28.02830             27.89069            27.88720
## # ... with 13 more variables: clim.timepoints.50 <dbl>,
## #   clim.timepoints.ssr <dbl>, proxy.bt <dbl>, proxy.bt.sb <dbl>,
## #   sed.acc.rate <dbl>, smoothing.width <dbl>, proxy.bt.sb.sampY <dbl>,
## #   proxy.bt.sb.sampYM <dbl>, proxy.bt.sb.inf.b <dbl>,
## #   proxy.bt.sb.inf.b.n <dbl>, proxy.bt.sb.sampYM.b <dbl>,
## #   proxy.bt.sb.sampYM.b.n <dbl>, simulated.proxy <dbl>
```

**Simple plotting**


```r
PFM$everything %>% 
  PlotPFMs()
```

```
## Scale for 'alpha' is already present. Adding another scale for 'alpha',
## which will replace the existing scale.
```

![](readme_files/figure-html/default_plot-1.png)<!-- -->


**The 10 replicates of the final simulated proxy**


```r
PFM$everything %>% 
  filter(stage == "simulated.proxy") %>% 
  spread(replicate, value) %>% 
  slice(1:5) %>% 
  kable(., digits = 2, format = "markdown")
```



| timepoints|stage           |     1|     2|     3|     4|     5|     6|     7|     8|     9|    10|
|----------:|:---------------|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|
|       4334|simulated.proxy | 28.59| 27.24| 27.85| 27.82| 28.75| 27.67| 28.24| 27.61| 28.45| 27.48|
|       4527|simulated.proxy | 27.39| 27.74| 27.18| 27.68| 28.97| 27.14| 26.70| 27.64| 28.06| 28.13|
|       4576|simulated.proxy | 28.52| 27.93| 27.88| 27.69| 27.36| 27.86| 27.95| 28.40| 27.17| 27.99|
|       4721|simulated.proxy | 27.72| 27.18| 28.37| 26.73| 27.36| 27.41| 28.22| 27.74| 27.84| 27.45|
|       4914|simulated.proxy | 27.44| 27.90| 27.40| 29.15| 27.55| 28.36| 28.29| 28.63| 28.08| 28.43|


### Proxy types

The initial input climate signal can be converted into "proxy units" if a `proxy.calibration.type` is specified. This simulates the Environment -> Sensor stage of the proxy system.



```r
set.seed(26052017)
PFM_2 <- ClimToProxyClim(clim.signal = N41.t21k.climate[nrow(N41.t21k.climate):1,] - 273.15, 
                         timepoints = round(N41.proxy$Published.age),
                         proxy.calibration.type = "MgCa",
                         seas.prod = N41.G.ruber.seasonality,
                         sed.acc.rate = N41.proxy$Sed.acc.rate.m.yr,
                         meas.noise = 0.46, n.samples = 30,
                         n.replicates = 1)
```

`PlotPFM` returns a ggplot2 object so it is easy to customize figures with new axes labels an such


```r
PFM_2$everything %>% 
  PlotPFMs(.) +
  labs(y = "Mg/Ca")
```

```
## Scale for 'alpha' is already present. Adding another scale for 'alpha',
## which will replace the existing scale.
```

![](readme_files/figure-html/MgCa_plot-1.png)<!-- -->









## Literature cited

Fraile, I., Schulz, M., Mulitza, S., & Kucera, M. (2008): Predicting the global distribution of planktonic foraminifera using a dynamic ecosystem model. Biogeosciences, 5: 891–911.

Laepple, T., & Huybers, P. (2013): Reconciling discrepancies between Uk37 and Mg/Ca reconstructions of Holocene marine temperature variability. Earth and Planetary Science Letters, 375: 418–429.

Rosenthal, Y., Oppo, D. W., & Linsley, B. K. (2003): The amplitude and phasing of climate change during the last deglaciation in the Sulu Sea, western equatorial Pacific. Geophys. Res. Lett., 30: 1428.

Shakun, J. D., Clark, P. U., He, F., Marcott, S. A., Mix, A. C., Liu, Z., Otto-Bliesner, B., Schmittner, A., & Bard, E. (2012): Global warming preceded by increasing carbon dioxide concentrations during the last deglaciation. Nature, 484: 49–54.


