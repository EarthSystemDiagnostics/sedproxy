#' Run ShinySedproxy
#'
#' @export
#'
#' @examples
#' ShinySedproxy()
ShinySedproxy <- function(){
  fl <- system.file("ShinySedproxy/app.R", package = "sedproxy")
  shiny::runApp(fl)
}

