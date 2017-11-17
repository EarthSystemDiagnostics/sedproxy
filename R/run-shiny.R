#' Run ShinySedproxy
#'
#' Launches the shiny app GUI interface to sedproxy
#'
#' @export
#'
#' @examples
#' ShinySedproxy()
ShinySedproxy <- function(){
  fl <- system.file("ShinySedproxy/app.R", package = "sedproxy")
  shiny::runApp(fl)
}

