#LIBRARIES----------------
library(shiny)

#UI---------------
ui <- fluidPage(# Application title
  titlePanel("2AG Fluorescence Simulation"),
  
  # Sidebar with inputs
  sidebarLayout(
    sidebarPanel(width=2,
                 selectInput("set_id_selection", "Experiment data", choices = c()),
                 
                 # General Input
                 h4("Pre-ATP Simulation Settings"),
                 sliderInput("signal_fraction", "Signaling pool saturation (%) :", min = 0, max = 100, value = 100, step = 1),
                 
                 # Subsection header for Pre-ATP simulation
                 h4("Pre-ATP Simulation Settings"),
                 numericInput("run_simulation_preATP", "Simulation duration (minutes):", 3*60),
                 # numericInput("simulation_dt_preATP", "Simulation timestep (minutes):", 1),
                 checkboxInput("preATP_in_percentage", "Display results in Percentage:", T),
                 
                 # Subsection header for Post-ATP simulation
                 h4("Post-ATP Simulation Settings"),
                 numericInput("run_simulation_postATP", "Simulation duration (minutes):", 30),
                 # numericInput("simulation_dt_postATP", "Simulation timestep (minutes):", 0.1),
                 checkboxInput("show_C_HEK", "Show [2AG] in HEK cells:", T),
                 
                 # reset button
                 actionButton("resetButton", "Reset", class = "btn-primary")
    ),
    
    # Main panel for displaying plots
    mainPanel(
      plotOutput("plot_combined"),
      # Description for the plots
      p("A) Pre-ATP value; simulating main and signaling pool replenishment in N2A cells"),
      p("B) Post-ATP value; simulating signal output when N2A was ATP-stimulated at a given signaling pool saturation. Control: at 100%"),
      
      h4("Simulation parameter multiplier (Default = 1)"),
      uiOutput("multipliersUI")
    )
  ))