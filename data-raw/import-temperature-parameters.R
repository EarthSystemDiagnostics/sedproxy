# import temperature preference parameters
library(tidyverse)
taxa.temperature.prefs <- readr::read_delim("data-raw/temperature-parameters.csv", delim = ";") %>%
  dplyr::mutate(Tmin = ifelse(is.na(Tmin), -Inf, Tmin),
                Tmax = ifelse(is.na(Tmax), Inf, Tmax),
                ID = paste(Taxon, Source, sep = " "))

usethis::use_data(taxa.temperature.prefs, overwrite = T)
