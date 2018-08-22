# Create example data
library("tidyverse")

#load(system.file("extdata/shak.t21k.list.rda", package = "ecusdata"))

N41.t21k.climate <- ecusdata::shak.t21k.list[["N41"]]

N41.proxy <- climproxyrecords::shakun.sed.acc %>%
  filter(Published.age < 21000,
         ID.no == "N41") %>%
  mutate(Sed.acc.rate.m.yr = ifelse(Sed.acc.rate.m.yr < 0.2 * mean(Sed.acc.rate.m.yr),
                                        0.2 * mean(Sed.acc.rate.m.yr),
                                        Sed.acc.rate.m.yr),
         Sed.acc.rate.cm.ka = round(Sed.acc.rate.m.yr * 1000*100, 2)) %>%
  select(Published.age, Published.temperature, Sed.acc.rate.cm.ka, Proxy.value)

N41.proxy.details <- climproxyrecords::shakun.metadata %>%
  filter(ID.no == "N41")

N41.G.ruber.seasonality <- ecusdata::shak.fraile.seasonality %>%
  ungroup() %>%
  filter(stats::complete.cases(V1),
         ID.no == "N41") %>%
  select(starts_with("V")) %>%
  gather() %>%
  .[["value"]]

devtools::use_data(N41.t21k.climate, N41.proxy, N41.proxy.details, N41.G.ruber.seasonality, overwrite = TRUE)


# Scussolini data
scussolini.tab1 <- ecusdata::scussolini.tab1
devtools::use_data(scussolini.tab1, overwrite = TRUE)
