library(sedproxy)
library(tidyverse)

PlotWeights <- function(weights){
  df <- as.data.frame(weights)
  df$time <- time(weights)
  df <- gather(df, habitat, value, -time)
  df$habitat <- factor(df$habitat, levels = colnames(weights), ordered = TRUE)

  df$value[df$value == 0] <- NA

  p <- ggplot(df, aes(x = habitat, y = time, fill = value))
  p <- p + geom_raster()
  p <- p + viridis::scale_fill_viridis(
    limits = c(0, max(df$value)),
    direction = -1)
  p <- p + theme_bw()
  return(p)
}


clim.in <- N41.t21k.climate[nrow(N41.t21k.climate):1,] - 273.15
clim.in <- ts(clim.in, start = -39)


wts.plafom <- matrix(rep(N41.G.ruber.seasonality, nrow(clim.in)),
                                     nrow = nrow(clim.in), byrow = TRUE)
colnames(wts.plafom) <- colnames(clim.in)

wts.norm <- dnorm(clim.in, mean = 27, sd = 5)


## With FAME 1.0

# install.packages("reticulate")
# requires a python installation

reticulate::source_python(system.file("extdata/py/forams_prod_l09.py", package = "sedproxy"))
#reticulate::source_python("py/fame-code-gmd-2017-251.py")

growth_rate_l09_array("ruber", c(280, 281))


# FAME <- function(m, taxon){
#   #n.row <- nrow(m)
#   n.col <- ncol(m)
#   c.names <- colnames(m)
#   v <- sapply(m+273.15, function(x) growth_rate_l09(taxon, x))
#   m2 <- matrix(v, ncol = n.col, byrow = FALSE)
#   #m2 <- m2 / sum(m2)
#   colnames(m2) <- c.names
#   return(m2)
# }
#
# system.time(
#   wts.fame <- FAME(clim.in[, ], "sacculifer")
#   )

wts.fame <- growth_rate_l09_array("ruber", clim.in + 273.15)
colnames(wts.fame) <- colnames(clim.in)

wts.fame.R <- growth_rate_l09_R("ruber", clim.in + 273.15)

table(wts.fame - wts.fame.R == 0)

PlotWeights(wts.plafom)
PlotWeights(wts.norm)
#PlotWeights((wts.fame))
PlotWeights((wts.fame.R))




# The input climate signal should be a time series object
# The Trace simulation runs to the year 1990 AD, therefore the start time for
# the input climate is -39 years BP

PLAFOM <- ClimToProxyClim(clim.signal = clim.in,
                       timepoints = round(N41.proxy$Published.age),
                       proxy.calibration.type = "identity",
                       proxy.prod.weights = wts.plafom,
                       sed.acc.rate = N41.proxy$Sed.acc.rate.cm.ka,
                       meas.noise = 0.46,
                       n.samples = 30,
                       n.replicates = 10)$everything %>%
  mutate(weights = "PLAFOM (modern static)")

FAME <- ClimToProxyClim(clim.signal = clim.in,
                       timepoints = round(N41.proxy$Published.age),
                       proxy.calibration.type = "identity",
                       proxy.prod.weights = wts.fame.R,
                       sed.acc.rate = N41.proxy$Sed.acc.rate.cm.ka,
                       meas.noise = 0.46,
                       n.samples = 30,
                       n.replicates = 10)$everything %>%
  mutate(weights = "FAME (dynamic)")


PFM <- bind_rows(PLAFOM = PLAFOM, FAME = FAME)


p.wts <- PFM %>% mutate(timepoints = timepoints/1000) %>%
  PlotPFMs(max.replicates = 1) +
  facet_wrap(~weights) +
  scale_x_continuous("Age [ka BP]")

ggsave("p.wts.png", p.wts, width = 6, height = 4, dpi = 300, units = "in")
