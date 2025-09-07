library(shiny)
library(tidyverse)
library(reshape2)
library(DT)

# Source the functions
source("functions.R")

# UI
ui <- fluidPage(
  titlePanel("scRNA Marker Gene Comparison"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Upload Dataset 1 (CSV)",
                accept = c(".csv")),
      fileInput("file2", "Upload Dataset 2 (CSV)",
                accept = c(".csv")),
      br(),
      actionButton("analyze", "Analyze Overlap", 
                   class = "btn-primary"),
      br(), br(),
      downloadButton("downloadData", "Download Results")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Heatmap", 
                 plotOutput("heatmap", height = "600px")),
        tabPanel("Data", 
                 DTOutput("dataTable")),
        tabPanel("Help", 
                 includeMarkdown("help.md"))
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  values <- reactiveValues(
    data1 = NULL,
    data2 = NULL,
    overlap = NULL
  )
  
  observeEvent(input$file1, {
    values$data1 <- read_csv(input$file1$datapath)
  })
  
  observeEvent(input$file2, {
    values$data2 <- read_csv(input$file2$datapath)
  })
  
  observeEvent(input$analyze, {
    req(values$data1, values$data2)
    
    values$overlap <- computeMarkersOverlap(
      values$data1,
      values$data2,
      denom = "Sample_1"
    )
  })
  
  output$heatmap <- renderPlot({
    req(values$overlap)
    heatmapMarkersOverlap(values$overlap)
  })
  
  output$dataTable <- renderDT({
    req(values$overlap)
    values$overlap
  })
  
  output$downloadData <- downloadHandler(
    filename = "marker_overlap.csv",
    content = function(file) {
      write_csv(values$overlap, file)
    }
  )
}

shinyApp(ui = ui, server = server)