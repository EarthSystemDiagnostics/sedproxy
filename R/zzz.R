# Create self contained Shiny app by pasting together
# the required functions and the ui and server files.
obj <- paste0(
  "library(ggplot2) \n",
  readr::read_file("R/ClimToProxyClim.R"),
  readr::read_file("R/PlotPFMs.R"),
  readr::read_file("R/BioturbationWeights.R"),
  readr::read_file("inst/ui-and-server.R"),
  sep = "\n",
  collapse = "\n"
)
readr::write_file(obj, file = "inst/sedproxy-shiny/app.R")

