library(dplyr)

## Create key and label lookup tables
stages.key <- dplyr::tibble(
  stage = c(
    "timepoints",
    "clim.signal.ann",
    "clim.signal.smoothed",
    "clim.timepoints.ssr",
    "proxy.bt",
    "proxy.bt.sb",
    "proxy.bt.sb.inf.b",
    "proxy.bt.sb.inf.b.n",
    "proxy.bt.sb.sampY",
    "proxy.bt.sb.sampYM",
    "proxy.bt.sb.sampYM.b",
    "proxy.bt.sb.sampYM.b.n",
    "simulated.proxy",
    "simulated.proxy.cal.err",
    "reconstructed.climate",
    "observed.proxy"
  ),
  label = c(
    "Requested timepoints",
    "(1) Input climate",
    "(1) Input climate",
    "(1) Input climate",
    "(2) +Bioturbation",
    "(3) +Habitat bias",
    "(.) +Bias",
    "(5) +Independent error",
    "(4) +Aliasing Y",
    "(4) +Aliasing YM",
    "(.) +Bias",
    "(5) +Independent error",
    "(5) +Independent error",
    "(6) +Calibration uncertainty",
    "(7) Reconstructed climate",
    "(*) Observed proxy"
  ),
  description = c(
    "Requested timepoints",
    "Input climate signal at requested timepoints at annual resolution",
    "Input climate signal at regular time intervals and resolution = plot.sig.res",
    "Input climate signal at requested timepoints, smoothed to resolution = plot.sig.res",
    "Climate signal after bioturbation",
    "Climate signal after bioturbation and habitat bias",
    "Climate signal after bioturbation, habitat bias, and bias",
    "Climate signal after bioturbation, habitat bias, and measurement error",
    "Climate signal after bioturbation, habitat bias, and aliasing of inter-annual variation",
    "Climate signal after bioturbation, habitat bias, and aliasing of inter-annual and intra-annual variation such as monthly temperatures or depth habitats",
    "Climate signal after bioturbation, habitat bias, and aliasing of inter-annual and intra-annual variation such as monthly temperatures or depth habitats, and bias",
    "Climate signal after bioturbation, habitat bias, aliasing, and measurement error",
    "Final simulated pseudo-proxy, this will be same as proxy.bt.sb.inf.b.n when n.samples = Inf, and proxy.bt.sb.sampYM.b.n when n.samples is finite",
    "Final simulated pseudo-proxy + uncertainty in true relationship to the climate variable (calibration error)",
    "Final pseudo-proxy calibrated to temperature",
    "True observed proxy (when supplied)"
  ),
  scale = c("time", rep("Climate units", 3), rep("Proxy units", 10), "Climate units", NA),
  plot.order = c(1, 1, 2, 3, 4, 5, 10, 6, 7, 8, 9, 11, 12, 13, 14, 15),
  plotting.colour = c("Black", "#018571", "#018571","#018571",
                      "Green", "Gold",
                      "Pink", "#7570b3",
                      "#d95f02", "#d95f02",
                      "Pink", "#7570b3", "#7570b3",
                      "Red",
                      "Blue",
                      "Red"),
  plotting.alpha = c(
    "timepoints" = 1,
    "clim.signal.ann" = 1,
    "clim.signal.smoothed" = 1,
    "clim.timepoints.ssr" = 1,
    "proxy.bt" = 1,
    "proxy.bt.sb" = 1,
    "proxy.bt.sb.inf.b" = 0.5,
    "proxy.bt.sb.inf.b.n" = 0.5,
    "proxy.bt.sb.sampY" = 0.5,
    "proxy.bt.sb.sampYM" = 1,
    "proxy.bt.sb.sampYM.b" = 0.5,
    "proxy.bt.sb.sampYM.b.n" = 0.5,
    "simulated.proxy" = 1,
    "simulated.proxy.cal.err" = 1,
    "reconstructed.climate" = 0.5,
    "observed.proxy" = 0.5
    #1, 1, 1, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 1, 1, 0.5, 0.5
    )
  #plotting.alpha = rep(1, 13)
)

#with(stages.key, factor(stage, ordered = TRUE, levels = stage[order(plot.order)]))
stage.labels <- stages.key$label
names(stage.labels) <- stages.key$stage

usethis::use_data(stages.key, stage.labels, overwrite = TRUE)

# import parameter description table
param.tab <- readr::read_delim("data-raw/parameter-descriptions.csv", delim = ";")
usethis::use_data(param.tab, overwrite = TRUE)


# create .rd format tables

# paste output into roxygen header of function definition
tb <- sedproxy::stages.key[, c("stage", "description")]
names(tb) <- c("Variable", "Description")
cat(sinew::tabular(tb))

# or like this
cat("#' \\describe{\n", paste0("#'    \\item{", tb$Variable, "}{", tb$Description, "}\n"),"#' }")



