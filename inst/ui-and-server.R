# Objects
stages.key <- structure(list(stage = c("timepoints", "clim.signal.ann", "clim.signal.smoothed",
                         "clim.timepoints.ssr", "proxy.bt", "proxy.bt.sb", "proxy.bt.sb.inf.b",
                         "proxy.bt.sb.inf.b.n", "proxy.bt.sb.sampY", "proxy.bt.sb.sampYM",
                         "proxy.bt.sb.sampYM.b", "proxy.bt.sb.sampYM.b.n", "simulated.proxy",
                         "observed.proxy"), label = c("Requested timepoints", "(1) Input climate",
                                                      "(1) Input climate", "(1) Input climate", "(2) +Bioturbation",
                                                      "(3) +Production bias", "(.) +Calibration bias", "(5) +Measurement error",
                                                      "(4) +Aliasing Y", "(4) +Aliasing YM", "(.) +Calibration bias",
                                                      "(5) +Measurement error", "(5) Simulated proxy", "(*) Observed proxy"
                         ), description = c("Requested timepoints", "Input climate signal at requested timepoints at annual resolution",
                                            "Input climate signal at regular time intervals and resolution = smoothed.signal.res",
                                            "Input climate signal at requested timepoints, smoothed to resolution = smoothed.signal.res",
                                            "Climate signal after bioturbation", "Climate signal after bioturbation and production bias",
                                            "Climate signal after bioturbation, production bias, and calibration bias",
                                            "Climate signal after bioturbation, production bias, and measurement error",
                                            "Climate signal after bioturbation, production bias, and aliasing of inter-annual variation",
                                            "Climate signal after bioturbation, production bias, and aliasing of inter-annual and intra-annual variation such as monthly temperatures or depth habitats",
                                            "Climate signal after bioturbation, production bias, and aliasing of inter-annual and intra-annual variation such as monthly temperatures or depth habitats, and calibration bias",
                                            "Climate signal after bioturbation, production bias, aliasing, and measurement error",
                                            "Final simulated pseudo-proxy, this will be same as proxy.bt.sb.inf.b.n when n.samples = Inf, and proxy.bt.sb.sampYM.b.n when n.samples is finite",
                                            "True observed proxy (when supplied)"), plot.order = c(1, 1,
                                                                                                   2, 3, 4, 5, 10, 6, 7, 8, 9, 11, 12, 13), plotting.colour = c("Black",
                                                                                                                                                                "#018571", "#018571", "#018571", "Green", "Gold", "Pink", "#7570b3",
                                                                                                                                                                "#d95f02", "#d95f02", "Pink", "#7570b3", "#7570b3", "Red"), plotting.alpha = c(1,
                                                                                                                                                                                                                                               1, 1, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)), .Names = c("stage",
                                                                                                                                                                                                                                                                                                               "label", "description", "plot.order", "plotting.colour", "plotting.alpha"
                                                                                                                                                                                                                                               ), class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA,
                                                                                                                                                                                                                                                                                                          -14L))

# Functions
SimPowerlaw <- function(beta, N)
{
  N2 <- (3 ^ ceiling(log(N, base = 3)))
  df  <- 1 / N2
  f <- seq(from = df, to = 1 / 2, by = df)
  Filter <- sqrt(1 / (f ^ beta))
  Filter <- c(max(Filter), Filter, rev(Filter))
  x   <- rnorm(N2, 1)
  fx  <- fft(x)
  ffx <- fx * Filter
  result <- Re(fft(ffx, inverse = TRUE))[1:N]
  return(scale(result)[1:N])
}

# Define UI ----
ui <- fluidPage(
  titlePanel("sedproxy"),
  p(
    em("sedproxy"),
    "is a forward model for sediment archived climate proxies.
    It is based on work described in Laepple and Huybers (2013).
    A manuscript is in preparation, Dolman and Laepple (in prep.),
    which will more fully describe the forward model and its applications.
    Please contact Dr Andrew Dolman <andrew.dolman@awi.de>,
    or Dr Thomas Laepple <tlaepple@awi.de>, at the Alfred-Wegener-Institute,
    Helmholtz Centre for Polar and Marine Research,
    Germany, for more information.
    "
  ),
  p("This work is supported by the BMBF funded PalMod project."),
  p(
    "Reference: ",
    "Laepple, T., & Huybers, P. (2013): Reconciling discrepancies between Uk37
    and Mg/Ca reconstructions of Holocene marine temperature variability.
    Earth and Planetary Science Letters, 375: 418-429."
  ),
  sidebarPanel(tabsetPanel(
    tabPanel(
      "Model parameters",
      fluidRow(
        h4(
          "Update the parameter values below
          and then run the proxy forward model."
        ),
        column(12,
               actionButton("run.pfm", "Run forward model")),
        hr()
        ),
      fluidRow(
        h4("Setup input climate signal"),
        column(
          width = 12,
          sliderInput(
            "clim.signal.length",
            h5("Length of input climate signal [years]"),
            value = 10000,
            step = 1000,
            min = 5000,
            max = 100000
          )
        ),
        column(
          width = 6,
          numericInput(
            "clim.signal.beta",
            h5("Slope of the power spectrum of the input climate signal"),
            value = 1,
            step = 0.1,
            min = 0.1,
            max = 3
          )
        ),
        column(
          width = 6,
          numericInput(
            "seas.amp",
            h5("Amplitude of the seasonal cycle"),
            value = 5,
            step = 0.5,
            min = 0,
            max = 20
          )
        )
      ),
      fluidRow(
        h4("Control sampling"),
        column(width = 6,
               numericInput(
                 "seed",
                 h5("Set RNG seed"),
                 value = 1,
                 step = 1,
                 min = 1
               )),
        column(
          width = 6,
          numericInput(
            "n.replicates",
            h5("No. replicates"),
            value = 1,
            step = 1,
            min = 1,
            max = 100
          )
        )
      ),
      fluidRow(
        column(
          width = 6,
          numericInput(
            "t.res",
            h5("Core sampling resolution [years]"),
            value = 100,
            step = 100,
            min = 1,
            max = 10000
          )
        ),
        column(
          width = 6,
          numericInput(
            "n.samples",
            h5("No. samples per timepoint"),
            value = 30,
            step = 1,
            min = 1,
            max = 1000
          )
        )
      ),
      fluidRow(
        h4("Sedimentation parameters"),
        column(
          6,
          numericInput(
            "bio.depth",
            h5("Bioturbation depth [cm]"),
            value = 10,
            step = 1,
            min = 0,
            max = 30
          )
        ),
        column(
          6,
          numericInput(
            "sed.acc.rate",
            h5("Sediment accumulation rate [cm/ka]"),
            value = 50,
            step = 1,
            min = 0,
            max = 100
          )
        )
      ),
      fluidRow(h4("Proxy production weights (monthly)"),
               column(12,
                      fluidRow(
                        column(
                          12,
                          radioButtons(
                            "seas",
                            label = NULL,
                            choices = c("Uniform", "Custom"),
                            selected = "Uniform",
                            inline = TRUE
                          ),
                          conditionalPanel(
                            condition = "input.seas == 'Custom'",
                            textInput(
                              "mon.vec",
                              "Modify the 12 monthly weights",
                              "1,1,1,1,1,1,1,1,1,1,1,1"
                            ),
                            span(textOutput("proxy.prod.weights.check"), style = "color:red")
                          )
                        )
                      ))),
      fluidRow(
        h4("Noise parameters"),
        column(
          6,
          numericInput(
            "meas.noise",
            h5("Measurement noise"),
            value = 0.46,
            step = 0.01,
            min = 0,
            max = 1
          )
        ),
        column(
          6,
          numericInput(
            "meas.bias",
            h5("Measurement bias"),
            value = 0,
            step = 0.1,
            min = 0,
            max = 2
          )
        )
      )

  ),
  tabPanel("Plot appearance",
           fluidRow(
             h4("Plot proxy stages:"),
             checkboxGroupInput(
               "stages",
               "Stages:",
               choices = list(
                 "Input climate" = "clim.signal.smoothed",
                 "Bioturbated climate" = "proxy.bt",
                 "Bioturbated + seasonally biased climate" =  "proxy.bt.sb",
                 "Sampled climate inc. aliasing" = "proxy.bt.sb.sampYM",
                 "Final pseudo-proxy" = "simulated.proxy"
               ),
               selected = list(
                 "clim.signal.smoothed",
                 "proxy.bt",
                 "proxy.bt.sb",
                 "proxy.bt.sb.sampYM",
                 "simulated.proxy"
               )
             )
           ))
    )),
  mainPanel(tabsetPanel(
    tabPanel("Plots",
             plotOutput("pfm.plot", height = "800px")),
    tabPanel("Numbers",
             dataTableOutput("pfm.str")),
    tabPanel("Placeholder", textOutput("proxy.prod.weights"))
  ))
  )


# Define server logic ----
server <- function(input, output) {
  clim <- eventReactive(input$run.pfm, {
    set.seed(input$seed)
    ann <-
      SimPowerlaw(input$clim.signal.beta, input$clim.signal.length)
    mon <-
      cos(seq(pi, 3 * pi, length.out = 12)) * input$seas.amp / 2
    clim <- outer(ann, mon, "+")
    clim <- ts(clim, start = 1)
    return(clim)
  }, ignoreNULL = FALSE)
  timepoints <- eventReactive(input$run.pfm, {
    #res <- 100
    tp <- seq(1, input$clim.signal.length, by = input$t.res)
    t.min <-
      ceiling(1000 * input$bio.depth / input$sed.acc.rate) + 1
    t.max <- input$clim.signal.length - 3 * t.min
    tp <- tp[tp > t.min & tp < t.max]
    return(tp)
  }, ignoreNULL = FALSE)
  seasprod <- eventReactive({
    input$mon.vec
    input$seas
  }, {
    if (input$seas == 'Custom')
    {
      v <- as.numeric(unlist(strsplit(input$mon.vec, ",")))
    } else{
      v <- rep(1, 12)
    }
    return(v)
  }, ignoreNULL = FALSE)
  output$proxy.prod.weights.check <- renderText({
    if (length(seasprod()) != 12)
    {
      paste0("You entered ",
             length(seasprod()),
             " values; 12 are required.")
    }
  })
  pfm <- eventReactive(input$run.pfm, {
    pfm <- ClimToProxyClim(
      clim.signal = clim(),
      timepoints = timepoints(),
      smoothed.signal.res = 100,
      bio.depth = input$bio.depth,
      sed.acc.rate = input$sed.acc.rate,
      proxy.prod.weights = seasprod(),
      n.samples = input$n.samples,
      n.replicates = input$n.replicates,
      meas.noise = input$meas.noise,
      meas.bias = input$meas.bias
    )
  }, ignoreNULL = FALSE)
  output$pfm.plot <- renderPlot({
    dat <- pfm()$everything
    dat <- subset(dat, dat$stage %in% input$stages)
    if (nrow(dat) > 0) {
      PlotPFMs(dat) +
        ggplot2::labs(x = "Age [years]")
    }
  }, res = 72 * 2)
  output$pfm.str <- renderDataTable({
    round(pfm()$simulated.proxy, 5)
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)
