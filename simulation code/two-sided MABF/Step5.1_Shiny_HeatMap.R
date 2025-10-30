# This is the Shiny app for Step4.0_Data_Visualization_HeatMap
# This Shiny will visualize TPR and FPR of MABF methods using heatmap
# Directly using helper functions result in slightly different counts than using this Shiny app. I suggest using the helper functions to generate/save heatmaps.
library(shiny)
library(ggplot2)
library(dplyr)
library(scales)
library(stringr)
# Load helper functions (FPR_heatmap, TPR_heatmap, etc.)
source("./helper functions_Shiny.R")

# Set the path to the directory containing RDS files
folder_path <- "./MABF4ROC/4TSFS" 
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

# ========== UI ==========
ui <- fluidPage(
  titlePanel("MABF Simulation Results Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("metric", "Choose Metric", choices = c(
        "False Positive Rate" = "FPR",
        "True Positive Rate"  = "TPR",
        "False Success Rate II" = "FSR2",
        "True Success Rate I"   = "TSR1"
      )),
      selectInput("method", "MABF Method", choices = c("FEMABF", "BFbMA", "iBF", "iBF2", "EUBF")),
      selectInput("effect_size", "True Effect Size", choices = c("0", "0.2", "0.5")),
      uiOutput("pcutoff_ui"),  # dynamic pcutoff
      radioButtons("output_format", "Show As", choices = c("Percentage", "Absolute Count")),
      checkboxInput("collapsed", "Collapse Across Datasets", value = FALSE),
      numericInput("cutoff", "Bayes Factor Cutoff", value = 1, min = 0, step = 0.1)
    ),
    
    mainPanel(
      plotOutput("heatmapPlot", height = "800px")
    )
  )
)

# ========== Server ==========
server <- function(input, output, session) {
  
  # Dynamically show/hide pcutoff input
  output$pcutoff_ui <- renderUI({
    if (input$metric %in% c("FSR2", "TSR1")) {
      selectInput("pcutoff", "P-value Cutoff", choices = c(0.05, 0.01), selected = 0.05)
    }
  })
  
  output$heatmapPlot <- renderPlot({
    method      <- input$method
    format      <- input$output_format
    collapsed   <- input$collapsed
    cutoff      <- input$cutoff
    effect_size <- input$effect_size
    metric      <- input$metric
    label_type  <- ifelse(format == "Percentage", "percent", "count")
    
    if (metric == "FPR") {
      if (collapsed) {
        return(FPR_heatmap_collapsed(method_name = method, cutoff = cutoff, label_type = label_type, effect_size = effect_size))
      } else {
        return(FPR_heatmap(method_name = method, cutoff = cutoff, label_type = label_type, effect_size = effect_size))
      }
    }
    
    if (metric == "TPR") {
      if (collapsed) {
        return(TPR_heatmap_collapsed(method_name = method, cutoff = cutoff, label_type = label_type, effect_size = effect_size))
      } else {
        return(TPR_heatmap(method_name = method, cutoff = cutoff, label_type = label_type, effect_size = effect_size))
      }
    }
    
    if (metric == "FSR2") {
      pcutoff <- input$pcutoff
      if (collapsed) {
        return(FSR2_heatmap_collapsed(method_name = method, cutoff = cutoff, pcutoff = pcutoff, label_type = label_type, effect_size = effect_size))
      } else {
        return(FSR2_heatmap(method_name = method, cutoff = cutoff, pcutoff = pcutoff, label_type = label_type, effect_size = effect_size))
      }
    }
    
    if (metric == "TSR1") {
      pcutoff <- input$pcutoff
      if (collapsed) {
        return(TSR1_heatmap_collapsed(method_name = method, cutoff = cutoff, pcutoff = pcutoff, label_type = label_type, effect_size = effect_size))
      } else {
        return(TSR1_heatmap(method_name = method, cutoff = cutoff, pcutoff = pcutoff, label_type = label_type, effect_size = effect_size))
      }
    }
  })
}

shinyApp(ui, server)