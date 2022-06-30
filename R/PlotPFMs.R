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
#' @importFrom rlang .data
#' @return a ggplot object of class "gg" "ggplot"
#' @export PlotPFMs
#'
#' @examples
#' library(ggplot2)
#' set.seed(26052017)
#' clim.in <- ts(N41.t21k.climate[nrow(N41.t21k.climate):1,] - 273.15)
#'
#' PFM <- ClimToProxyClim(clim.signal = clim.in,
#'                        timepoints = round(N41.proxy$Published.age),
#'                        calibration.type = "identity",
#'                        habitat.weights = N41.G.ruber.seasonality,
#'                        sed.acc.rate = N41.proxy$Sed.acc.rate.cm.ka,
#'                        sigma.meas = 0.45,
#'                        sigma.ind = 0,
#'                        n.samples = Inf,
#'                        plot.sig.res = 10, meas.bias = 1,
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

  PFMs.in <- PFMs
  if ("sedproxy.pfm" %in% class(PFMs.in)) PFMs <- PFMs.in$everything

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

  breaks <- sedproxy::stages.key$stage

  if (colr.palette[1] == "default")
    colr.palette  <-
    structure(sedproxy::stages.key$plotting.colour,
              .Names = sedproxy::stages.key$stage)

  if (alpha.palette[1] == "default") alpha.palette  <-
    structure(sedproxy::stages.key$plotting.alpha,
              .Names = sedproxy::stages.key$stage)

  if (levl.labels[1] == "default") levl.labels  <-
    structure(sedproxy::stages.key$label,
              .Names = sedproxy::stages.key$stage)

  cali.attr <- attr(PFMs, "calibration.pars")

  if (is.null(cali.attr)) {
    cali.attr <- list(calibration.type = "identity")
  }

  if (plot.stages[1] == "default") {
    if (cali.attr$calibration.type == "identity"){
      plotting.levels <- c(
        "clim.signal.monthly", "clim.signal.smoothed", "proxy.bt", "proxy.bt.sb",
        "proxy.bt.sb.sampYM",  "simulated.proxy", "observed.proxy"
      )}else{
        plotting.levels <- c(
          "clim.signal.monthly", "clim.signal.smoothed", "proxy.bt", "proxy.bt.sb",
          "proxy.bt.sb.sampYM",  "simulated.proxy",  "simulated.proxy.cal.err", "reconstructed.climate", "observed.proxy"
        )
      }
  } else if (plot.stages == "all") {
    plotting.levels <- sedproxy::stages.key$stage
    plotting.levels <- subset(plotting.levels, plotting.levels %in% c("clim.signal.ann", "clim.timepoints.ssr") == FALSE)
  } else{
    plotting.levels <- plot.stages
  }

  PFMs <- dplyr::filter(PFMs, stage %in% plotting.levels,
                        replicate <= max.replicates)


  # match scaling flag
  PFMs <- dplyr::left_join(PFMs, sedproxy::stages.key[, c("stage", "scale")])

  #set factor level ordering for stages
  stage.order <- match.arg(stage.order)
  switch(stage.order,
         seq = PFMs$stage <- factor(PFMs$stage, levels = plotting.levels, ordered = TRUE),
         var = {
           var.order <- tapply(PFMs$value, PFMs$stage, FUN = stats::var)
           var.order <- rank(var.order, ties.method = "first")
           var.order <- names(sort(var.order, decreasing = TRUE))
           PFMs$stage <- factor(PFMs$stage,
                                levels = var.order, ordered = TRUE)
         })


  p <- ggplot2::ggplot(data = PFMs, aes(x = .data$timepoints, y = .data$value,
                                        colour = stage, alpha = stage,
                                        linetype = as.factor(replicate))) +
    geom_line() +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), legend.position = "top") +

    guides(colour = guide_legend(
      ncol = 2,
      override.aes = list(alpha = 1))) +
    labs(x = expression("Timepoints"),
         y = expression("Proxy value")) +
    scale_linetype_manual(values = rep(1, 13*length(unique(PFMs$replicate))), guide = "none")+
    scale_alpha_manual(guide = "none")

  pal.df <- data.frame(
    colr.palette = colr.palette,
    colr.breaks = names(colr.palette),
    labels = levl.labels,
    alpha.palette = alpha.palette,
    alpha.breaks = names(alpha.palette)
  )

  pal.df <- dplyr::filter(pal.df, .data$colr.breaks %in% unique(PFMs$stage))

  if (is.null(colr.palette) == FALSE)
    p <- p + scale_colour_manual("", values = pal.df$colr.palette, breaks = pal.df$colr.breaks,
                                 labels = pal.df$labels)

  if (is.null(alpha.palette) == FALSE)
    p <- p + scale_alpha_manual("", values = pal.df$alpha.palette, breaks = pal.df$alpha.breaks,
                                labels = pal.df$labels)

  if (cali.attr$calibration.type != "identity"){
    p <- p + #facet_wrap(~scale, scales = "free_y") +
      facet_wrap( ~ scale, strip.position = "left", scales = "free_y") +
      theme(
        # remove the default y-axis title, "wt"
        axis.title.y = element_blank(),
        # replace the strip backgrounds with transparent
        strip.background = element_rect(fill = 'transparent', colour = 'transparent'),
        strip.placement = 'outside')
  }

  return(p)
}







