# Define UI ----
ui <- fluidPage(titlePanel("sedproxy"),
                sidebarLayout(
                  sidebarPanel(
                    fluidRow(
                      title = "Start",
                      h4("Help text"),
                      helpText(
                        "Note: help text isn't a true widget,",
                        "but it provides an easy way to add text to",
                        "accompany other widgets."
                      ),
                      actionButton("goButton", "Run forward model"),
                      br(),
                      br()
                    ),
                    fluidRow(
                      h3("Setup input climate"),
                      sliderInput(
                        "clim.signal.length",
                        h4("Length of input climate signal [years]"),
                        value = 1000,
                        step = 1000,
                        min = 1000,
                        max = 10000
                      ),
                      column(width = 6,
                      numericInput(
                        "clim.signal.beta",
                        h4("Slope of the power spectrum of the input climate signal"),
                        value = 1,
                        step = 0.1,
                        min = 0.1,
                        max = 3
                      )),
                      column(width = 6,
                      numericInput(
                        "seas.amp",
                        h4("Amplitude of the seasonal cycle"),
                        value = 1,
                        step = 0.1,
                        min = 0.5,
                        max = 3
                      )
                    )),
                    fluidRow(
                      h3("Sedimentation parameters"),
                      column(6, 
                             numericInput(
                        "bio.depth",
                        h4("Bioturbation depth [m]"),
                        value = 0.1,
                        step = 0.01,
                        min = 0,
                        max = 1
                      )),
                      column(6,
                      numericInput(
                        "sed.acc.rate",
                        h4("Sediment accumulation rate [m/ka]"),
                        value = 5e-04,
                        step = 0.01 / 100,
                        min = 0,
                        max = 1
                      ))
                    ),
                    fluidRow(
                      h3("Noise parameters"),
                      numericInput(
                        "seed",
                        h4("Set RNG seed"),
                        value = 1,
                        step = 1,
                        min = 1
                      ),
                      column(6, 
                      numericInput(
                        "meas.noise",
                        h4("Measurement noise"),
                        value = 0.46,
                        step = 0.01,
                        min = 0,
                        max = 1
                      )),
                      column(6,
                      numericInput(
                        "meas.bias",
                        h4("Measurement bias"),
                        value = 0,
                        step = 0.1,
                        min = 0,
                        max = 2
                      ))
                    ),
                    fluidRow(
                      title = "Samples",
                      numericInput(
                        "n.samples",
                        h4("Number of samples per timepoint"),
                        value = 30,
                        step = 1,
                        min = 1,
                        max = Inf
                      ),
                      numericInput(
                        "n.replicates",
                        h4("Number of replicates"),
                        value = 1,
                        step = 1,
                        min = 1,
                        max = 100
                      ),
                      numericInput(
                        "smoothed.signal.res",
                        h4("Smoothed signal resolution for plotting [years]"),
                        value = 100,
                        step = 10,
                        min = 1,
                        max = 10000
                      )
                    )
                  ),
                  mainPanel(width = 6,
                            tabsetPanel(
                              tabPanel(
                                "Plots",
                                plotOutput("clim.plot", width = "100%", height = "400px")
                              ),
                              tabPanel("Placeholder"),
                              tabPanel("Placeholder")
                            ))
                ))



# Define server logic ----
server <- function(input, output) {
  input.clim <- reactive({
    set.seed(input$seed)
    ann <- PaleoSpec::SimPowerlaw(input$clim.signal.beta, input$clim.signal.length)
    mon <- cos(seq(pi, 3*pi, length.out = 12)) * input$seas.amp/2
    clim <- outer(ann, mon, "+")
    return(clim)
  })
  output$clim.plot <- renderPlot({
    plot(input.clim()[,1], type = "l")
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)