  # Create self contained Shiny app by pasting together
  # the required functions and the ui and server files.
  if (file.mtime("inst/extdata/ui-and-server.R") >
      file.mtime("inst/sedproxy-shiny/app.R") |
      file.exists("inst/sedproxy-shiny/app.R") == FALSE) {
    library(readr)
    obj <- paste0(
      "library(ggplot2) \n",
      read_file("R/ClimToProxyClim.R"),
      read_file("R/PlotPFMs.R"),
      read_file("R/ImpulseResponse.R"),
      read_file("inst/extdata/ui-and-server.R"),
      sep = "\n",
      collapse = "\n"
    )
    write_file(obj, path = "inst/sedproxy-shiny/app.R")
  }
