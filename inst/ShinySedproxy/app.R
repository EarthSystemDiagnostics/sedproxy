# Packages ---
library(sedproxy)

# Functions ----
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
    "Laepple, T., & Huybers, P. (2013): Reconciling discrepancies between Uk37
    and Mg/Ca reconstructions of Holocene marine temperature variability.
    Earth and Planetary Science Letters, 375: 418–429."
  ),
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        p(
          "Update the parameter values below
          and then run the proxy forward model."
        ),
        actionButton("run.pfm", "Run forward model"),
        hr()
        ),
      fluidRow(
        h4("Setup input climate"),
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
            h5("Bioturbation depth [m]"),
            value = 0.1,
            step = 0.01,
            min = 0,
            max = 1
          )
        ),
        column(
          6,
          numericInput(
            "sed.acc.rate",
            h5("Sediment accumulation rate [m/ka]"),
            value = 5e-04,
            step = 0.01 / 100,
            min = 0,
            max = 1
          )
        )
      ),
      fluidRow(h4("Proxy production weights (monthly)"),
               column(
                 12,
                 fluidRow(
                   radioButtons(
                     "seas",
                     label = NULL,
                     choices = c("Uniform", "Custom"),
                     selected = "Uniform",
                     inline = TRUE
                   ),
                   conditionalPanel(
                     condition = "input.seas == 'Custom'",
                     column(
                       4,
                       sliderInput("Jan", "Jan", 0, 1, 0.5, 0.1),
                       sliderInput("Feb", "Feb", 0, 1, 0.5, 0.1),
                       sliderInput("Mar", "Mar", 0, 1, 0.5, 0.1),
                       sliderInput("Apr", "Apr", 0, 1, 0.5, 0.1)
                     ),
                     column(
                       4,
                       sliderInput("May", "May", 0, 1, 0.5, 0.1),
                       sliderInput("Jun", "Jun", 0, 1, 0.5, 0.1),
                       sliderInput("Jul", "Jul", 0, 1, 0.5, 0.1),
                       sliderInput("Aug", "Aug", 0, 1, 0.5, 0.1)
                     ),
                     column(
                       4,
                       sliderInput("Sep", "Sep", 0, 1, 0.5, 0.1),
                       sliderInput("Oct", "Oct", 0, 1, 0.5, 0.1),
                       sliderInput("Nov", "Nov", 0, 1, 0.5, 0.1),
                       sliderInput("Dec", "Dec", 0, 1, 0.5, 0.1)
                     )
                   )
                 )
               )),
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
    mainPanel(tabsetPanel(
      tabPanel("Plots",
               plotOutput("pfm.plot", height = "800px")),
      tabPanel("Numbers",
               dataTableOutput("pfm.str")),
      tabPanel("Placeholder", textOutput("seas.prod"))
    ))
  )
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
    return(clim)
  }, ignoreNULL = FALSE)
  timepoints <- eventReactive(input$run.pfm, {
    #res <- 100
    tp <- seq(1, input$clim.signal.length, by = input$t.res)
    t.min <- ceiling(input$bio.depth / input$sed.acc.rate) + 1
    t.max <- input$clim.signal.length - 3 * t.min
    tp <- tp[tp > t.min & tp < t.max]
    return(tp)
  }, ignoreNULL = FALSE)
  seasprod <- eventReactive({
    input$run.pfm
    # input$Jan
    # input$Feb
    # input$Mar
    # input$Apr
    # input$May
    # input$Jun
    # input$Jul
    # input$Aug
    # input$Sep
    # input$Oct
    # input$Nov
    # input$Dec
  }, {
    c(
      input$Jan,
      input$Feb,
      input$Mar,
      input$Apr,
      input$May,
      input$Jun,
      input$Jul,
      input$Aug,
      input$Sep,
      input$Oct,
      input$Nov,
      input$Dec
    )
  }, ignoreNULL = FALSE)
  pfm <- eventReactive(input$run.pfm, {
    pfm <- ClimToProxyClim(
      clim.signal = clim(),
      timepoints = timepoints(),
      smoothed.signal.res = 100,
      bio.depth = input$bio.depth,
      sed.acc.rate = input$sed.acc.rate,
      seas.prod = seasprod(),
      n.samples = input$n.samples,
      n.replicates = input$n.replicates,
      meas.noise = input$meas.noise,
      meas.bias = input$meas.bias
    )
  }, ignoreNULL = FALSE)
  output$pfm.plot <- renderPlot({
    PlotPFMs(pfm()$everything) +
      ggplot2::labs(x = "Age [years]")
  }, res = 72 * 2)
  output$pfm.str <- renderDataTable({
    round(pfm()$simulated.proxy, 5)
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)
