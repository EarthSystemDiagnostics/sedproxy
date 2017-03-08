#' Plot forward modelled sedimentary proxies
#'
#' @param PFMs A dataframe of forward modelled proxies
#' @param colr.palette Colours for the proxy stages
#' @param alpha.palette Alpha levels for the proxy stages
#' @param breaks Proxy stages for legend, in order
#' @param levl.labels Labels for the proxy stages, in order
#'
#' @return
#' @export
#'
#' @examples
PlotPFMs <- function(PFMs, colr.palette = NULL, alpha.palette = NULL, breaks, levl.labels){

  rug.dat <- filter(PFMs, Stage %in% c("proxy.sig.inf"),
                    replicate == 1)

  p <- ggplot(data = PFMs, aes(x = Age / 1000, y = Temperature,
                               colour = Stage, alpha = Stage)) +
    geom_rug(data = rug.dat, sides = "b", colour = "Darkgrey") +
    geom_line(aes(linetype = factor(replicate))) +
    facet_wrap(~ Location + ID.no, scales = "free_y",
               labeller = labeller(.multi_line = FALSE), ncol = 4) +
    theme(legend.position = "top") +
    guides(colour = guide_legend(label.position = "top",
                                 label.hjust = 1,
                                 nrow = 1,
                                 override.aes = list(alpha = 1))) +
    labs(title = unique(PFMs$Proxy),
         x=expression(Age~"[ka]"),
         y = expression(Temperature~"[°C]")) +
    scale_linetype(guide = FALSE) +
    scale_alpha_discrete(guide = FALSE)

  if (is.null(colr.palette) == FALSE)
    p <- p + scale_colour_manual("", values = colr.palette, breaks = breaks, labels = levl.labels)

  if (is.null(alpha.palette) == FALSE)
    p <- p + scale_alpha_manual("", values = levl.alphas, breaks = levl.names,
                                labels = levl.labels)

  return(p)
}



#' Plot all ggforce::facet_wrap_paginate pages
#'
#' @description A helper function to return all pages from
#' ggforce::facet_wrap_paginate. It automatically calculates
#' required number of pages.
#'
#' @param ggplot.obj A ggplot object
#' @param facets Faceting formula, see \code{\link{ggplot2::facet_wrap}}
#' @param nrow number of columns per page
#' @param ncol number of rows per page
#'
#' @return A list of ggplot objects
#' @export
#'
#' @examples
#' gg <- facet_wrap_paginate_auto(uk37.plots.2,
#'         facets = formula(~Location + Proxy + ID.no),
#'         nrow = 4, ncol = 3)
#' gg
facet_wrap_paginate_auto <- function(ggplot.obj, facets, nrow, ncol){
  n.pages <- ceiling(length(ggplot_build(ggplot.obj)$layout$panel_layout$PANEL) / (nrow * ncol))
  ggs <- lapply(1:n.pages, function(i) {
    ggplot.obj +
      ggforce::facet_wrap_paginate(facets,
                                   scales = "free_y",
                                   labeller = labeller(.multi_line = FALSE),
                                   ncol = ncol, nrow = nrow, page = i)}
  )
  return(ggs)
}