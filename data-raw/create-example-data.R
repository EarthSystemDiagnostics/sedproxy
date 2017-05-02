# Create example data

library("tidyverse")

load(system.file("extdata/shak.t21k.list.rda", package = "ecusdata"))

N41.t21k.climate <- shak.t21k.list[["N41"]]

N41.proxy <- climproxyrecords::shakun.sed.acc %>%
  filter(Published.age < 21000,
         ID.no == "N41") %>%
  mutate(Sed.acc.rate.m.yr = ifelse(Sed.acc.rate.m.yr < 0.2 * mean(Sed.acc.rate.m.yr),
                                        0.2 * mean(Sed.acc.rate.m.yr),
                                        Sed.acc.rate.m.yr)) %>%
  select(Published.age, Published.temperature, Sed.acc.rate.m.yr)

N41.proxy.details <- climproxyrecords::shakun.metadata %>%
  filter(ID.no == "N41")

N41.G.ruber.seasonality <- ecusdata::shak.fraile.seasonality %>%
  ungroup() %>%
  filter(complete.cases(V1),
         ID.no == "N41") %>%
  select(starts_with("V")) %>%
  gather() %>%
  .[["value"]]

N41.sed.acc.rate <- climproxyrecords::shakun.sed.acc %>%
  filter(ID.no == "N41",
         Published.age < 21000) %>%
  select(Sed.acc.rate.m.yr) %>%
  mutate(Sed.acc.rate.m.yr = ifelse(Sed.acc.rate.m.yr < 0.2 * mean(Sed.acc.rate.m.yr),
                                    0.2 * mean(Sed.acc.rate.m.yr),
                                    Sed.acc.rate.m.yr)) %>%
  .[[1]]

devtools::use_data(N41.sed.acc.rate, N41.t21k.climate, N41.proxy, N41.proxy.details, N41.G.ruber.seasonality, overwrite = TRUE)
