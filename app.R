#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(vroom)
library(tidyverse)
library(tibble)

datamatrix <- read.table("02_data_integration/datamatrix_tmm.txt",
                         row.names = 1,
                         sep = "\t")
log2(datamatrix + 1)

ui <- fluidPage(
  theme = shinytheme("yeti"),
  title = "Zebrafish Heart Transcriptome Through Development",
  
  
  selectInput("dataset", label = "Dataset", choices = ls("package:datasets")),
  verbatimTextOutput("summary"),
  tableOutput("table")
)
server <- function(input, output, session) {
  # Create a reactive expression
  dataset <- reactive({
    get(input$dataset, "package:datasets")
  })
  
  output$summary <- renderPrint({
    # Use a reactive expression by calling it like a function
    summary(dataset())
  })
  
  output$table <- renderTable({
    dataset()
  })
}
shinyApp(ui, server)