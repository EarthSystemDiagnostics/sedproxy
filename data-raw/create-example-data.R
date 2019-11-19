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


# Rescaled Gisp2 from Löwemark, L., Konstantinou, K. I. and Steinke, S.: Bias in
# foraminiferal multispecies reconstructions of paleohydrographic conditions
# caused by foraminiferal abundance variations and bioturbational mixing: A
# model approach, Marine Geology, 256(1–4), 101–106,
# doi:10.1016/j.margeo.2008.10.005, 2008.

#Alley, R.B.. 2004. GISP2 Ice Core Temperature and Accumulation Data. IGBP
#PAGES/World Data Center for Paleoclimatology Data Contribution Series
##2004-013. NOAA/NGDC Paleoclimatology Program, Boulder CO, USA.

gisp2 <- climproxyrecords::alley.temperature %>% 
  mutate(Age = (Age * 1000)) %>% 
  filter(complete.cases(Temperature))

gisp2.ann.ts <- PaleoSpec::MakeEquidistant(gisp2$Age, gisp2$Temperature,
                                           time.target = seq(1, 50000, by = 1))

gisp2.ann <- data.frame(age.yr.bp = time(gisp2.ann.ts),
                        temperature = gisp2.ann.ts) %>%
  tbl_df() %>% 
  filter(complete.cases(temperature)) %>% 
  mutate(temperature.rescaled = -temperature / max(temperature),
         temperature.rescaled = temperature.rescaled - min(temperature.rescaled),
         temperature.rescaled = -2*temperature.rescaled / max(temperature.rescaled))

#plot(temperature.rescaled~age.yr.bp, data = gisp2.ann, type = "l")

usethis::use_data(gisp2.ann, overwrite = TRUE)
