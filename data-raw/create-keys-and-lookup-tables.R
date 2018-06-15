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
    "(.) +Calibration bias",
    "(5) +Measurement error",
    "(4) +Aliasing Y",
    "(4) +Aliasing YM",
    "(.) +Calibration bias",
    "(5) +Measurement error",
    "(5) Simulated proxy",
    "(6) Reconstructed climate",
    "(*) Observed proxy"
  ),
  description = c(
    "Requested timepoints",
    "Input climate signal at requested timepoints at annual resolution",
    "Input climate signal at regular time intervals and resolution = smoothed.signal.res",
    "Input climate signal at requested timepoints, smoothed to resolution = smoothed.signal.res",
    "Climate signal after bioturbation",
    "Climate signal after bioturbation and habitat bias",
    "Climate signal after bioturbation, habitat bias, and calibration bias",
    "Climate signal after bioturbation, habitat bias, and measurement error",
    "Climate signal after bioturbation, habitat bias, and aliasing of inter-annual variation",
    "Climate signal after bioturbation, habitat bias, and aliasing of inter-annual and intra-annual variation such as monthly temperatures or depth habitats",
    "Climate signal after bioturbation, habitat bias, and aliasing of inter-annual and intra-annual variation such as monthly temperatures or depth habitats, and calibration bias",
    "Climate signal after bioturbation, habitat bias, aliasing, and measurement error",
    "Final simulated pseudo-proxy, this will be same as proxy.bt.sb.inf.b.n when n.samples = Inf, and proxy.bt.sb.sampYM.b.n when n.samples is finite",
    "True observed proxy (when supplied)",
    "Final pseudo-proxy calibrated to temperature"
  ),
  scale = c("time", rep("climate", 3), rep("proxy", 9), "climate", NA),
  plot.order = c(1, 1, 2, 3, 4, 5, 10, 6, 7, 8, 9, 11, 12, 13, 14),
  plotting.colour = c("Black", "#018571", "#018571","#018571",
                      "Green", "Gold",
                      "Pink", "#7570b3",
                      "#d95f02", "#d95f02",
                      "Pink", "#7570b3", "#7570b3",
                      "Darkblue",
                      "Red"),
  plotting.alpha = c(1, 1, 1, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
  #plotting.alpha = rep(1, 13)
)

#with(stages.key, factor(stage, ordered = TRUE, levels = stage[order(plot.order)]))
stage.labels <- stages.key$label
names(stage.labels) <- stages.key$stage

devtools::use_data(stages.key, stage.labels, overwrite = TRUE)

# import parameter description table
param.tab <- readr::read_delim("data-raw/parameter-descriptions.csv", delim = ";")
devtools::use_data(param.tab, overwrite = TRUE)


# create .rd format tables

# paste output into roxygen header of function definition
tb <- sedproxy::stages.key[, c("stage", "description")]
names(tb) <- c("Variable", "Description")
cat(sinew::tabular(tb))

# or like this
cat("#' \\describe{\n", paste0("#'    \\item{", tb$Variable, "}{", tb$Description, "}\n"),"#' }")



# calibration parameters

MgCa.foram.pars <- climproxycalibration::MgCa.foram.pars
devtools::use_data(MgCa.foram.pars, overwrite == TRUE)

UK37.pars <- list(mueller.uk37 = structure(list(
  means = structure(
    c(0.0328750614815548, 0.0686612340110185),
    .Names = c("intercept", "slope")
  ),
  vcov = structure(
    structure(c(1.46053818728255e-07, -2.80815422746781e-06,
                -2.80815422746781e-06, 6.06536807458765e-05),
              .Dim = c(2L, 2L)),
    .Dimnames = list(c("intercept", "slope"),
                     c("intercept", "slope"))
  )
),
.Names = c("means", "vcov")))

devtools::use_data(UK37.pars, overwrite == TRUE)

