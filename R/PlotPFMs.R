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
                     breaks = "default.breaks",

                     colr.palette = "default.colr.palette",
                     alpha.palette = "default.alpha.palette",
                     levl.labels = "default.levl.labels"){

  # set factor level ordering for Stage
  PFMs$Stage <- factor(PFMs$Stage, levels = rev(c(
    "Observed proxy",
    "proxy.bt.sb.sampYM.b.n",
    "proxy.bt.sb.inf.b.n",
    "proxy.bt.sb.sampYM.b",
    "proxy.bt.sb.inf.b",
    "proxy.bt.sb.sampYM",
    "proxy.bt.sb",
    "proxy.bt",
    "clim.signal.ann",
    "clim.timepoints.1000",
    "clim.timepoints.100",
    "clim.timepoints.50",
    "clim.timepoints.ssr",
    "clim.signal.smoothed"
  )), ordered = TRUE)


  if(exists("replicate", where = PFMs)){
    rug.dat <- dplyr::filter(PFMs, Stage %in% c("proxy.bt.sb.inf.b.n", "proxy.bt.sb.sampYM.b.n", "Observed proxy"),
                             replicate == 1)
  }else{
    rug.dat <- dplyr::filter(PFMs, Stage %in% c("proxy.bt.sb.inf.b.n", "proxy.bt.sb.sampYM.b.n", "Observed proxy"))
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
  if (breaks == "default.breaks"){
    breaks <-
      c(
        "clim.signal.ann",
        "clim.timepoints.1000",
        "clim.timepoints.100",
        "clim.timepoints.50",
        "clim.timepoints.ssr",
        "clim.signal.smoothed",
        "proxy.bt",
        "proxy.bt.sb",
        "proxy.bt.sb.inf.b",
        "proxy.bt.sb.sampYM.b",
        "proxy.bt.sb.sampYM",
        "proxy.bt.sb.inf.b.n",
        "proxy.bt.sb.sampYM.b.n",
        "simulated.proxy",
        "Observed proxy"
      )}
  
  if (colr.palette == "default.colr.palette") colr.palette  <-
      structure(c("#018571", "#018571","#018571", "#018571","#018571", "#018571",
                  "Green", "Gold",
                  "White", "White",
                  "#d95f02",
                  "#7570b3", "#7570b3", "#7570b3",
                  "Red"),
                .Names = breaks)
  
  if (alpha.palette == "default.alpha.palette") alpha.palette  <-
      structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5),
                .Names = breaks)
  
  if (levl.labels == "default.levl.labels") levl.labels  <-
      structure(c(rep(" Input climate", 6),
                  "+Bioturbation",
                  "+Seasonality",
                  "+Finite sample",
                  rep("+Bias", 2),
                  rep("+Measurement error", 3),
                  "Observed proxy"),
                .Names = breaks)
  
    p <- ggplot(data = PFMs, aes(x = Age / 1000, y = Temperature,
                               colour = Stage, alpha = Stage)) +
    geom_rug(data = rug.dat, sides = "b", colour = "Darkgrey") +
    geom_line(aes(group = paste0(Stage, replicate))) +
    #geom_line() +
    facet_wrap(~ Location + ID.no, scales = "free_y",
               labeller = labeller(.multi_line = FALSE), ncol = 4) +
    theme_bw() +
    theme(legend.position = "top", panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(label.position = "top",
                                 label.hjust = 1,
                                 nrow = 1,
                                 override.aes = list(alpha = 1))) +
    labs(title = unique(PFMs$Proxy),
         x = expression(Age~"[ka]"),
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


# PlotPFMs_2 <- function(PFMs,
#                        breaks = "default.breaks",
#                        
#                        colr.palette = "default.colr.palette",
#                        alpha.palette = "default.alpha.palette",
#                        levl.labels = "default.levl.labels"){
#   
#   
#   
#   if(exists("replicate", where = PFMs)){
#     rug.dat <- dplyr::filter(PFMs, stage %in% c("simulated.proxy", "Observed proxy"),
#                              replicate == 1)
#   }else{
#     rug.dat <- dplyr::filter(PFMs, stage %in% c("simulated.proxy", "Observed proxy"))
#     rug.dat$replicate <- 1
#     PFMs$replicate <- 1
#   }
#   
#   if(exists("Location", where = PFMs)==FALSE){
#     PFMs$Location <- ""
#   }
#   
#   if(exists("ID.no", where = PFMs)==FALSE){
#     PFMs$ID.no <- ""
#   }
#   if(exists("Proxy", where = PFMs)==FALSE){
#     PFMs$Proxy <- ""
#   }
#   
#   # assign default asthetic mappings
#   if (breaks == "default.breaks"){
#     breaks <-
#       c(
#         "clim.signal.ann",
#         "clim.timepoints.1000",
#         "clim.timepoints.100",
#         "clim.timepoints.50",
#         "clim.timepoints.ssr",
#         "clim.signal.smoothed",
#         "proxy.bt",
#         "proxy.bt.sb",
#         "proxy.bt.sb.inf.b",
#         "proxy.bt.sb.sampYM.b",
#         "proxy.bt.sb.sampYM",
#         "proxy.bt.sb.inf.b.n",
#         "proxy.bt.sb.sampYM.b.n",
#         "simulated.proxy",
#         "Observed proxy"
#       )}
#   
#   if (colr.palette == "default.colr.palette") colr.palette  <-
#       structure(c("#018571", "#018571","#018571", "#018571","#018571", "#018571",
#                   "Green", "Gold",
#                   "White", "White",
#                   "#d95f02",
#                   "#7570b3", "#7570b3", "#7570b3",
#                   "Red"),
#                 .Names = breaks)
#   
#   if (alpha.palette == "default.alpha.palette") alpha.palette  <-
#       structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5),
#                 .Names = breaks)
#   
#   if (levl.labels == "default.levl.labels") levl.labels  <-
#       structure(c(rep(" Input climate", 6),
#                   "+Bioturbation",
#                   "+Seasonality",
#                   "+Finite sample",
#                   rep("+Bias", 2),
#                   rep("+Measurement error", 3),
#                   "Observed proxy"),
#                 .Names = breaks)
#   
#   
#   # Filter stages 
#   #if ("proxy.bt" %in% PFMs$stage == TRUE){
#   plotting.levels <- c("clim.signal.smoothed", "proxy.bt",
#                        "proxy.bt.sb", "proxy.bt.sb.sampYM",
#                        "simulated.proxy", "proxy.bt.sb.sampYM.b.n")
#   # }else if ()
#   
#   
#   PFMs <- subset(PFMs, PFMs$stage %in% plotting.levels)
#   
#   #set factor level ordering for stage
#   PFMs$stage <- factor(PFMs$stage, levels = rev(plotting.levels), ordered = TRUE)
#   
#   
#   
#   p <- ggplot(data = PFMs, aes(x = timepoints, y = value,
#                                colour = stage, alpha = stage)) +
#     # geom_rug(data = rug.dat, sides = "b", colour = "Darkgrey") +
#     geom_line(aes(group = paste0(stage, replicate))) +
#     theme_bw() +
#     theme(legend.position = "top", panel.grid.minor = element_blank()) +
#     guides(colour = guide_legend(label.position = "top",
#                                  label.hjust = 1,
#                                  nrow = 1,
#                                  override.aes = list(alpha = 1))) +
#     labs(x = expression("timepoints"),
#          y = expression("Proxy value")) +
#     
#     scale_linetype(guide = FALSE) +
#     scale_alpha_discrete(guide = FALSE)
#   
#   if (is.null(colr.palette) == FALSE)
#     p <- p + scale_colour_manual("", values = colr.palette, breaks = breaks, labels = levl.labels)
#   
#   if (is.null(alpha.palette) == FALSE)
#     p <- p + scale_alpha_manual("", values = alpha.palette, breaks = breaks,
#                                 labels = levl.labels)
#   
#   return(p)
# }
# 
# 
# PFM$everything %>% 
#   filter(stage %in% c("clim.signal.smoothed", "proxy.bt", "proxy.bt.sb", "proxy.bt.sb.sampYM",
#                       "simulated.proxy")) %>%
#   data.frame(.) %>% 
#   PlotPFMs_2()



