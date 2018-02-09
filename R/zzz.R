# Create self contained Shiny app by pasting together
# the required functions and the ui and server files.
library(readr)
obj <- paste0(
  "library(ggplot2) \n",
  read_file("R/ClimToProxyClim.R"),
  read_file("R/PlotPFMs.R"),
  read_file("R/BioturbationWeights.R"),
  read_file("inst/ui-and-server.R"),
  sep = "\n",
  collapse = "\n"
)
write_file(obj, path = "inst/sedproxy-shiny/app.R")

