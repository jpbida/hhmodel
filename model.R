# Install and load the solver package
# install.packages("deSolve")
# install.packages("shiny")
library(deSolve)
library(shiny)

# Define the user interfact for HH Model Constants ----
# We have sliderbars for constants and a plot for the
# four differential equations
ui <- fluidPage(
  # App title ----
  titlePanel("Hodgkin-Huxley Constants"),
  sidebarLayout(
    # Sidebars to collect constant values
    sidebarPanel(
      sliderInput("Cm", "Cm:",
        min = 0, max = 5,
        value = 1, step = 0.1
      ),
      sliderInput("gNa", "gNa:",
        min = 60, max = 200,
        value = 120, step = 1
      ),
      sliderInput("gK", "gK:",
        min = 0, max = 80,
        value = 36, step = 0.1
      ),
      sliderInput("gL", "gL:",
        min = 0, max = 5,
        value = 0.3, step = 0.01
      ),
      sliderInput("ENa", "ENa:",
        min = 0, max = 100,
        value = 50, step = 0.1
      ),
      sliderInput("EK", "EK:",
        min = -100, max = 0,
        value = -77, step = 0.1
      ),
      sliderInput("EL", "EL:",
        min = -100, max = 0,
        value = -54.387, step = 0.01
      ),
    ),
    mainPanel(
      plotOutput(outputId = "HHPlot")
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  ## Every time the constants are updated run the ODE
  ## solver and plot the results
  output$HHPlot <- renderPlot({
    # 1. Define Model Parameters
    params <- c(
      Cm    = input$Cm, # uF/cm2
      gNa   = input$gNa, # mS/cm2
      gK    = input$gK, # mS/cm2
      gL    = input$gL, # mS/cm2
      ENa   = input$ENa, # mV
      EK    = input$EK, # mV
      EL    = input$EL # mV
    )

    # 2. Rate Functions for Gating Variables
    alpha_m <- function(V) 0.1 * (V + 40) / (1 - exp(-(V + 40) / 10))
    beta_m <- function(V) 4.0 * exp(-(V + 65) / 18)
    alpha_h <- function(V) 0.07 * exp(-(V + 65) / 20)
    beta_h <- function(V) 1 / (1 + exp(-(V + 35) / 10))
    alpha_n <- function(V) 0.01 * (V + 55) / (1 - exp(-(V + 55) / 10))
    beta_n <- function(V) 0.125 * exp(-(V + 65) / 80)

    # 3. The ODE System
    hh_model <- function(t, state, parameters) {
      with(as.list(c(state, parameters)), {
        # TODO: Make a user interface for this feature
        # External stimulus: 10 uA from t=10 to t=60
        I_ext <- if (t >= 10 && t <= 60) 10 else 0

        # Calculate currents
        INa <- gNa * m^3 * h * (V - ENa)
        IK <- gK * n^4 * (V - EK)
        IL <- gL * (V - EL)

        # Differential equations
        dV <- (I_ext - INa - IK - IL) / Cm
        dm <- alpha_m(V) * (1 - m) - beta_m(V) * m
        dh <- alpha_h(V) * (1 - h) - beta_h(V) * h
        dn <- alpha_n(V) * (1 - n) - beta_n(V) * n

        return(list(c(dV, dm, dh, dn)))
      })
    }

    # 4. Initial Conditions (Resting State)
    V_start <- -65
    initial_state <- c(
      V = V_start,
      m = alpha_m(V_start) / (alpha_m(V_start) + beta_m(V_start)),
      h = alpha_h(V_start) / (alpha_h(V_start) + beta_h(V_start)),
      n = alpha_n(V_start) / (alpha_n(V_start) + beta_n(V_start))
    )

    # 5. Solve the Model
    times <- seq(0, 100, by = 0.01)
    out <- ode(y = initial_state, times = times, func = hh_model, parms = params)

    # 6. Plotting
    par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
    plot(out[, "time"], out[, "V"],
      type = "l", lwd = 2, col = "black",
      main = "Hodgkin-Huxley Simulation ", ylab = "Membrane Potential (mV)", xlab = ""
    )
    grid()

    plot(out[, "time"], out[, "m"],
      type = "l", col = "blue", ylim = c(0, 1),
      ylab = "Gating Variables", xlab = "Time (ms)"
    )
    lines(out[, "time"], out[, "h"], col = "red")
    lines(out[, "time"], out[, "n"], col = "green")
    legend("topright", legend = c("m", "h", "n"), col = c("blue", "red", "green"), lty = 1, bty = "n")
    grid()
  })
}

# Create Shiny app ----
shinyApp(ui, server)
