#' Plot forward modelled sedimentary proxies
#'
#' @param PFMs A dataframe of forward modelled proxies
#' @param breaks Proxy stages for legend, in order
#' @param colr.palette Colours for the proxy stages
#' @param alpha.palette Alpha levels for the proxy stages
#' @param levl.labels Labels for the proxy stages, in order
#'
#' @import ggplot2
#' @export PlotPFMs
#'
#' @examples
PlotPFMs <- function(PFMs,
                     breaks = c("clim.signal.ann", "clim.timepoints.1000", "clim.signal.smoothed", "clim.timepoints.50", "proxy.bt", "proxy.bt.sb",
                                "proxy.bt.sb.inf.b.n", "proxy.bt.sb.sampYM", "proxy.bt.sb.sampYM.b.n", "Actual proxy"),
                     #dfc27d
                     colr.palette = structure(c("#018571", "#018571","#018571", "#018571", "Green", "Gold", "#7570b3", "#d95f02", "#7570b3", "Red"),
                                               .Names = c("clim.signal.ann", "clim.timepoints.1000", "clim.signal.smoothed", "clim.timepoints.50", "proxy.bt", "proxy.bt.sb",
                                                          "proxy.bt.sb.inf.b.n", "proxy.bt.sb.sampYM",
                                                          "proxy.bt.sb.sampYM.b.n", "Actual proxy")),
                     alpha.palette = structure(c(1, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5),
                                                .Names = c("clim.signal.ann","clim.timepoints.1000", "clim.signal.smoothed", "clim.timepoints.50",
                                                           "proxy.bt", "proxy.bt.sb", "proxy.bt.sb.inf.b.n",
                                                           "proxy.bt.sb.sampYM", "proxy.bt.sb.sampYM.b.n", "Actual proxy")),
                     levl.labels = structure(c(" Modelled climate", " Modelled climate", " Modelled climate", " Modelled climate", "+Bioturbation",
                                                "+Seasonality", "+Measurement error",
                                                "+Finite sample", "+Measurement error", "Actual proxy"),
                                              .Names = c("clim.signal.ann", "clim.timepoints.1000", "clim.signal.smoothed", "clim.timepoints.50",
                                                         "proxy.bt", "proxy.bt.sb", "proxy.bt.sb.inf.b.n",
                                                         "proxy.bt.sb.sampYM", "proxy.bt.sb.sampYM.b.n", "Actual proxy"))
){

  if(exists("replicate", where = PFMs)){
    rug.dat <- dplyr::filter(PFMs, Stage %in% c("proxy.bt"),
                      replicate == 1)
  }else{
    rug.dat <- dplyr::filter(PFMs, Stage %in% c("proxy.bt"))
    rug.dat$replicate <- 1
    PFMs$replicate <- 1
  }

  if(exists("Location", where = PFMs)==FALSE){
    PFMs$Location <- ""
    }

  if(exists("ID.no", where = PFMs)==FALSE){
    PFMs$ID.no <- ""
  }

  p <- ggplot(data = PFMs, aes(x = Age / 1000, y = Temperature,
                               colour = Stage, alpha = Stage)) +
    geom_rug(data = rug.dat, sides = "b", colour = "Darkgrey") +
    geom_line(aes(linetype = factor(replicate))) +
    facet_wrap(~ Location + ID.no, scales = "free_y",
               labeller = labeller(.multi_line = FALSE), ncol = 4) +
    theme_bw() +
    theme(legend.position = "top", panel.grid.minor = element_blank()) +
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
    p <- p + scale_alpha_manual("", values = alpha.palette, breaks = breaks,
                                labels = levl.labels)

  return(p)
}
