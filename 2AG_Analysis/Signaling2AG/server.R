#INPUTS------------
mainwd <- "path_to_folder/" # change to path to shiny folder!
sourceFunction <- "sourceFunctions.R"

#LIBRARIES-------------
library(shiny)
library(cowplot)
source(paste0(mainwd, sourceFunction))

#SERVER-----------------
server <- function(input, output, session) {
  # Input Operations #
  # Load simulation components
  updateSelectInput(session, "set_id_selection", choices = setID_selection)
  
  # Generate UI for parameter multipliers dynamically
  output$multipliersUI <- renderUI({
    numParams <- nrow(pars) # Example number of parameters
    nCols <- 3 # Adjust the number of columns as needed
    nRows <- ceiling(numParams / nCols)
    
    # Create a matrix to hold UI elements for grid layout, fill by column
    inputsMatrix <- matrix(list(), nrow = nRows, ncol = nCols, byrow = FALSE)
    
    for (i in 1:numParams) {
      rowIndex <- ((i - 1) %% nRows) + 1
      colIndex <- ((i - 1) %/% nRows) + 1
      inputsMatrix[rowIndex, colIndex] <- list(
        numericInput(inputId = paste0("multiplier", i),
                     label = pars$names[i],
                     value = 1)
      )
    }
    
    # Generate UI with fluid rows and columns
    numericUI <- lapply(1:nrow(inputsMatrix), function(row) {
      fluidRow(
        lapply(1:nCols, function(col) {
          input <- inputsMatrix[row, col]
          if (!is.null(input[[1]])) {
            column(width = 12 / nCols, input[[1]])
          }
        })
      )
    })
    
    do.call(tagList, numericUI)
  })
  
  # Combine plots based on inputs
  output$plot_combined <- renderPlot({
    # data selection
    data_set <- subset(simPackage[['data']], setID==input$set_id_selection) %>%
      mutate(time = time/60)
    
    # Collect multipliers from inputs
    multipliers <- sapply(1:nrow(pars), function(i) input[[paste0("multiplier", i)]])

    # Pre-ATP
    preATP_simulation <- runSimulation_preATP(multipliers, input$signal_fraction/100, 
                                              simulation_tmax = input$run_simulation_preATP, 
                                              # simulation_dt = input$simulation_dt_preATP, 
                                              as_fraction = input$preATP_in_percentage)
    
    # Post-ATP
    postATP_simulation <- runSimulation_postATP(multipliers, input$signal_fraction/100, 
                                                simulation_tmax = input$run_simulation_postATP, 
                                                # simulation_dt = input$simulation_dt_postATP, 
                                                data_set = data_set, show_concentration=input$show_C_HEK)
    
    # combined plot
    plot_combined <- plot_grid(preATP_simulation, postATP_simulation, ncol=2, scale=0.95, labels=c("A)", "B)"))

    return(plot_combined)
  })
  
  observeEvent(input$resetButton, {
    sapply(c(1:nrow(pars)), function(i){
      updateNumericInput(session, inputId = paste0("multiplier", i),
                         value = 1)
    })
    updateSliderInput(session, inputId = "signal_fraction", value = 100)
  })
}
