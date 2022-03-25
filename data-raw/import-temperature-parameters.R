# import temperature preference parameters
library(tidyverse)
taxa.temperature.prefs <- readr::read_delim("data-raw/temperature-parameters.csv", delim = ";") %>%
  dplyr::mutate(Tmin = ifelse(is.na(Tmin), -Inf, Tmin),
                Tmax = ifelse(is.na(Tmax), Inf, Tmax),
                ID = paste(Taxon, Source, sep = " "))

#usethis::use_data(taxa.temperature.prefs, overwrite = T)



l09_cnsts_dic <- list(sacculifer = c(0.3, 2155, 289.3, 304.4, 94385, 171209), 
                      bulloides = c(0.25, 6482, 280, 298.5, 295374, 159618),
                      pachy_d = c(0.13, 6584, 277.8, 293.7, 300239, 110583), 
                      siphonifera = c(0.17, 7105, 284.9, 301.8, 349991, 130852),
                      universa = c(0.23, 626, 289.8, 304, 30953, 249412),
                      pachy_s = c(0.37, 6584, 999, 279.8, 0, 59491),
                      dutertrei = c(0.1, 6876, 280.1, 298.8, 210746, 57855),
                      ruber = c(0.24, 1086, 292.3, 303.4, 53765, 165490))


l09_maxgrowth_dic <- list(sacculifer = 0.369316159643382, bulloides = 0.317097315352481, 
     siphonifera = 0.279187426094619, pachy_d = 0.109441954385795, 
     universa = 0.240849035558603, pachy_s = 0.0488642301133836, 
     dutertrei = 0.109516330927052, ruber = 0.260441142387609)




# Calibration parameters

## MgCa
### Create table combining estimates vcov matrices with publishd parameters from Anand et al. 2003

MgCa.foram.pars <- climproxycalibration::MgCa.foram.pars
tbl3 <- climproxycalibration::anand.table.3


MgCa.foram.pars.df <- plyr::ldply(MgCa.foram.pars, function(x) {
  data.frame(slope =  x$means[1],
             intercept = x$means[2],
             vcov = I(list(x$vcov)))}
  , .id = "Species") %>%
  mutate(Species = as.character(Species))


MgCa.foram.pars.df[1, "Species"] <- "Ten planktonic species_350-500"

MgCa.foram.pars.Anand <- left_join(MgCa.foram.pars.df, tbl3[c(1, 22:36),], by = "Species") %>%
  tbl_df() %>%
  rename(calibration = Species) %>%
  mutate(calibration.type = "MgCa") %>%
  mutate(slope = ifelse(is.na(A), slope, A),
         intercept = ifelse(is.na(B), intercept, log(B))) %>%
  select(calibration.type, calibration, slope, intercept, vcov, Notes)



Uk37.pars.df <- tibble(
  calibration.type = "Uk37",
  calibration = "Mueller global",
  slope = 0.033,
  intercept = 0.044,
  vcov = list(structure(
    structure(c(1.46053818728255e-07, -2.80815422746781e-06,
                -2.80815422746781e-06,  6.06536807458765e-05),
              .Dim = c(2L, 2L)),
    .Dimnames = list(c("slope", "intercept"),
                     c("slope", "intercept"))
  )))

calibration.parameters <- bind_rows(MgCa.foram.pars.Anand, Uk37.pars.df)


# Uk37.pars <- list(mueller.uk37 = structure(list(
#   means = structure(
#     # Set to Eq. 30, Table 1. Müller et al (1998)
#     c(0.033, 0.044),
#     .Names = c("slope", "intercept")
#   ),
#   vcov = structure(
#     structure(c(1.46053818728255e-07, -2.80815422746781e-06,
#                 -2.80815422746781e-06, 6.06536807458765e-05),
#               .Dim = c(2L, 2L)),
#     .Dimnames = list(c("slope", "intercept"),
#                      c("slope", "intercept"))
#   )
# ),
# .Names = c("means", "vcov")))


# CalibrationParameters <- list(MgCa = climproxycalibration::MgCa.foram.pars,
#                               Uk37 = Uk37.pars)


#devtools::use_data(Uk37.pars, overwrite == TRUE)
#devtools::use_data(MgCa.foram.pars, overwrite == TRUE)
#usethis::use_data(calibration.parameters, overwrite == TRUE)
#usethis::use_data(CalibrationParameters, overwrite == TRUE)


usethis::use_data(calibration.parameters, taxa.temperature.prefs, l09_cnsts_dic, l09_maxgrowth_dic,
                  internal = TRUE, overwrite = TRUE)

