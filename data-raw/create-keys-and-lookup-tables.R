## Create key and label lookup tables
stages.key <- dplyr::tibble(
  stage = c(
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
    "observed.proxy"
  ),
  label = c(
    "(1) Input climate",
    "(1) Input climate",
    "(1) Input climate",
    "(2) +Bioturbation",
    "(3) +Production bias",
    "(.) +Calibration bias",
    "(5) +Measurement error",
    "(4) +Aliasing Y",
    "(4) +Aliasing YM",
    "(.) +Calibration bias",
    "(5) +Measurement error",
    "(5) Simulated proxy",
    "(*) Observed proxy"
  ),
  description = c(
    "Input climate signal at requested timepoints at annual resolution",
    "Input climate signal at regular time intervals and resolution = smoothed.signal.res",
    "Input climate signal at requested timepoints, smoothed to resolution = smoothed.signal.res",
    "Climate signal after bioturbation",
    "Climate signal after bioturbation and production bias",
    "Climate signal after bioturbation, production bias, and calibration bias",
    "Climate signal after bioturbation, production bias, and measurement error",
    "Climate signal after bioturbation, production bias, and aliasing of inter-annual variation",
    "Climate signal after bioturbation, production bias, and aliasing of inter-annual and intra-annual variation such as monthly temperatures or depth habitats",
    "Climate signal after bioturbation, production bias, and aliasing of inter-annual and intra-annual variation such as monthly temperatures or depth habitats, and calibration bias",
    "Climate signal after bioturbation, production bias, aliasing, and measurement error",
    "Final simulated pseudo-proxy",
    "True observed proxy (when supplied)"
  ),
  plot.order = c(1, 2, 3, 4, 5, 10, 6, 7, 8, 9, 11, 12, 13),
  plotting.colour = c("#018571", "#018571","#018571",
                      "Green", "Gold",
                      "Pink", "#7570b3",
                      "#d95f02", "#d95f02",
                      "Pink", "#7570b3", "#7570b3",
                      "Red"),
  plotting.alpha = c(1, 1, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
  #plotting.alpha = rep(1, 13)
)

#with(stages.key, factor(stage, ordered = TRUE, levels = stage[order(plot.order)]))
stage.labels <- stages.key$label
names(stage.labels) <- stages.key$stage

devtools::use_data(stages.key, stage.labels, overwrite = TRUE)

