#' Plot forward modelled sedimentary proxies
#'
#' @param PFMs A dataframe of forward modelled proxies
#' @param stage.order Controls the order in which proxy stages are plotted,
#' either sequentially, "seq", or in order of variance, "var". Defaults to var.
#' @param plot.stages Proxy stages to be plotted, "default", "all", or a custom character vector
#' @param colr.palette Colours for the proxy stages
#' @param alpha.palette Alpha levels for the proxy stages
#' @param levl.labels Labels for the proxy stages
#' @param max.replicates Maximum number of replicates to plot at once
#'
#' @import ggplot2
#' @importFrom dplyr filter
#' @export PlotPFMs
#'
#' @examples
#' library(ggplot2)
#' set.seed(26052017)
#' clim.in <- N41.t21k.climate[nrow(N41.t21k.climate):1,] - 273.15
#'
#' PFM <- ClimToProxyClim(clim.signal = clim.in,
#'                        timepoints = round(N41.proxy$Published.age),
#'                        proxy.calibration.type = "identity",
#'                        proxy.prod.weights = N41.G.ruber.seasonality,
#'                        sed.acc.rate = N41.proxy$Sed.acc.rate.cm.ka,
#'                        sigma.measurement = 0.45,
#'                        sigma.individual = 0,
#'                        n.samples = Inf,
#'                        smoothed.signal.res = 10, meas.bias = 1,
#'                        n.replicates = 10)
#'
#' PlotPFMs(PFM$everything, max.replicates = 1, stage.order = "seq") +
#'   facet_wrap(~stage)
#'
#' PlotPFMs(PFM$everything, max.replicates = 1, stage.order = "var")
#'
#' PlotPFMs(PFM$everything, stage.order = "var", plot.stages = "all")
#'
PlotPFMs <- function(PFMs,
                     stage.order = c("var", "seq"),
                     plot.stages = c("default"),
                     max.replicates = 5,
                     colr.palette = "default",
                     alpha.palette = "default",
                     levl.labels = "default"){

  if(exists("replicate", where = PFMs)){
    rug.dat <- dplyr::filter(PFMs, stage %in% c("simulated.proxy", "observed.proxy"),
                             replicate == 1)
  }else{
    rug.dat <- dplyr::filter(PFMs, stage %in% c("simulated.proxy", "observed.proxy"))
    rug.dat$replicate <- 1
    PFMs$replicate <- 1
  }

  if(exists("Location", where = PFMs)==FALSE){
    PFMs$Location <- ""
  }

  if(exists("ID.no", where = PFMs)==FALSE){
    PFMs$ID.no <- ""
  }
  if(exists("Proxy", where = PFMs)==FALSE){
    PFMs$Proxy <- ""
  }

  # assign default asthetic mappings

  breaks <- stages.key$stage

  if (colr.palette[1] == "default")
    colr.palette  <-
      structure(stages.key$plotting.colour,
                .Names = stages.key$stage)

  if (alpha.palette[1] == "default") alpha.palette  <-
      structure(stages.key$plotting.alpha,
                .Names = stages.key$stage)

  if (levl.labels[1] == "default") levl.labels  <-
      structure(stages.key$label,
                .Names = stages.key$stage)

  if (plot.stages[1] == "default") {
    plotting.levels <- c(
      "clim.signal.monthly", "clim.signal.smoothed", "proxy.bt", "proxy.bt.sb",
      "proxy.bt.sb.sampYM",  "simulated.proxy",  "simulated.proxy.cal.err", "reconstructed.climate", "observed.proxy"
      )
  } else if (plot.stages == "all") {
    plotting.levels <- stages.key$stage
    plotting.levels <- subset(plotting.levels, plotting.levels %in% c("clim.signal.ann", "clim.timepoints.ssr") == FALSE)
  } else{
    plotting.levels <- plot.stages
  }

  PFMs <- dplyr::filter(PFMs, stage %in% plotting.levels,
                        replicate <= max.replicates)


  # match scaling flag
  PFMs <- dplyr::left_join(PFMs, stages.key[, c("stage", "scale")])

  #set factor level ordering for stages
  stage.order <- match.arg(stage.order)
  switch(stage.order,
         seq = PFMs$stage <- factor(PFMs$stage, levels = plotting.levels, ordered = TRUE),
         var = {
           var.order <- tapply(PFMs$value, PFMs$stage, FUN = var)
           var.order <- rank(var.order, ties.method = "first")
           var.order <- names(sort(var.order, decreasing = TRUE))
           PFMs$stage <- factor(PFMs$stage,
                                levels = var.order, ordered = TRUE)
           })


  p <- ggplot2::ggplot(data = PFMs, aes(x = timepoints, y = value,
                               colour = stage, alpha = stage,
                               linetype = as.factor(replicate))) +
    #geom_rug(data = rug.dat, sides = "b", colour = "Darkgrey") +
    geom_line() +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), legend.position = "top") +
    # guides(colour = guide_legend(label.position = "top",
    #                              label.hjust = 1,
    #                              nrow = 1, byrow = TRUE,
    #                              override.aes = list(alpha = 1))) +
    guides(colour = guide_legend(#label.position = "top",
                                 #label.hjust = 1,
                                 ncol = 2,
                                 #byrow = TRUE,
                                 override.aes = list(alpha = 1))) +
    labs(x = expression("Timepoints"),
         y = expression("Proxy value")) +
    scale_linetype_manual(values = rep(1, 13*length(unique(PFMs$replicate))), guide = FALSE)+
    scale_alpha_manual(guide = FALSE)

  if (is.null(colr.palette) == FALSE)
    p <- p + scale_colour_manual("", values = colr.palette, breaks = names(colr.palette),
                                 labels = levl.labels)

  if (is.null(alpha.palette) == FALSE)
    p <- p + scale_alpha_manual("", values = alpha.palette, breaks = names(alpha.palette),
                                labels = levl.labels)

  p <- p + facet_wrap(~scale, scales = "free_y")

  return(p)
}






