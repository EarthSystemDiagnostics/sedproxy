#' Run ShinySedproxy
#'
#' Launches the shiny app GUI interface to sedproxy
#'
#' @export
#'
#' @examples
#' ShinySedproxy()
ShinySedproxy <- function(){
  fl <- system.file("sedproxy-shiny/app.R", package = "sedproxy")
  shiny::runApp(fl)
}

