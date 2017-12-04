#' sedproxy: Simulation of Sediment Archived Climate Proxy Records
#'
#' The sedproxy package provides functions to simulate sediment archived proxies
#' @docType package
#' @name sedproxy
NULL

#' @title Description of proxy stages
#' @description A description of the proxy stages in the output of \code{ClimToProxyClim}
#' and default labels, colours and order for plotting
#' @format A data frame with 13 rows and 6 variables:
#' \tabular{ll}{
#'   \cr  \code{stage}  \tab proxy stage
#'   \cr  \code{label}  \tab label for proxy stage
#'   \cr  \code{description}  \tab description of proxy stage
#'   \cr  \code{plot.order}  \tab default plotting order of stages
#'   \cr  \code{plotting.colour}  \tab default colour for plotting
#'   \cr  \code{plotting.alpha}  \tab default alpha level for plotting
#'}
"stages.key"

#' @title Labels for proxy stages
#' @description Labels for proxy stages. For plotting.
#' @format A named character vector
"stage.labels"

