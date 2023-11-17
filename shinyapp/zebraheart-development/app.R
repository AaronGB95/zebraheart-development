#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)

NJData <- data.frame(
  "Year" = c(2000, 2000, 2001, 2001, 2002, 2002, 2003, 2003),
  "WUCode" = c("IR", "PS", "IR", "PS","IR", "PS","IR", "PS"),
  "Annual.Value" = c(12, 14, 19, 7, 11, 13, 20, 17)
)

ui <- fluidPage(
  titlePanel("Subsetting Dataset"),
  sidebarLayout(
    sidebarPanel(
      selectInput("codeInput1", label = "Choose Year", choices = unique(NJData$Year)),
      selectInput("codeInput2", label = "Choose WU Type", choices = unique(NJData$WUCode))
    ),
    mainPanel(
      plotOutput("plot")
    )
  )
)

server <- function(input, output) {
  
  dataset <- reactive({
    subset(
      NJData,
      (Year == input$codeInput1 & WUCode == input$codeInput2)
    )
  })
  
  output$plot <- renderPlot({
    ggplot(
      dataset(),
      aes_string(
        x = "Year",
        y = "Annual.Value",
        fill = "Year"
      )
    ) + geom_boxplot()
  })
  
}

shinyApp(ui = ui, server = server)