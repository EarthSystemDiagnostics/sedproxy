#' Run ShinySedproxy
#'
#' Launches the shiny app GUI interface to sedproxy
#'
#' @export
#' @importFrom readr read_file
#' @importFrom readr write_file
#' @examples
#' \dontrun{ShinySedproxy()}
ShinySedproxy <- function(){
  fl <- system.file("sedproxy-shiny/app.R", package = "sedproxy")
  shiny::runApp(fl)
}

