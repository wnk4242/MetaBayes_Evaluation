
library(shiny)
library(readr)
library(dplyr)
library(ggplot2)


# Set the path to the directory containing RDS files
folder_path <- "./MABFanalyses/matrix-wise/rates4Viz/"
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

# Combine rate outcomes of different methods in order to compare them using AUC
# Combine datasets of BFMA methods
rates_0.2null <- rbind(BFbMA_TFPS4Viz_0.2null, EUBF_TFPS4Viz_0.2null, FEMABF_TFPS4Viz_0.2null, iBF_TFPS4Viz_0.2null)
# Change threshold_p_r to threshold
FEMA_TFPS4Viz_0.2null <- FEMA_TFPS4Viz_0.2null %>% 
  rename(threshold = threshold_p_r)
# Change threshold_BF to threshold
rates_0.2null <- rates_0.2null  %>% 
  rename(threshold = threshold_BF)
# Combine rows of MABF methods and FEMA
rates_0.2null <- rbind(FEMA_TFPS4Viz_0.2null, rates_0.2null)

# Convert character to factor
rates_0.2null <- rates_0.2null %>%
  mutate(true_es = 0.2) %>% 
  mutate(across(c(method, true_es, orig.n, `QRP level`, `PB level`, `rep number`,rep.n), as.factor)) %>% 
  relocate(true_es, .after = `true effect`) %>% 
  relocate(group, .before = `method`) %>% 
  rename(scenario = group)
# Remove and add columns
df <- rates_0.2null %>% 
  select(-`true effect`,-true_es, -threshold_ES) %>% 
  mutate(bias_level = `PB level`) %>% 
  select(-`QRP level`, -`PB level`) %>% 
  relocate(bias_level, .before = `rep number`) %>% 
  mutate(TPR = TP/(TP+FN),
         FPR = FP/(FP+TN),
         
         TSR1 = TS1/Successes1, #True success rate; Of all the claimed successful replications (i.e., original study and replications are consistent), what proportion of them are true successes (i.e., both original and replications reflects the true underlying effect)
         FSR1 = FS1/Successes1, #False success rate; Of all the claimed successful replications (i.e., original study and replications are consistent),  what proportion of them are false successes (i.e., neither original and replications reflects the true underlying effect)
         TFR1 = TF1/Failures1, #True failure rate; Of all the claimed failed replications (i.e., original and replications are not consistent), what proportion of them are true failures (i.e., only the original study reflects the true underlying effect)
         FFR1 = FF1/Failures1, #False failure rate; Of all the claimed failed replications (i.e., original and replications are not consistent), what proportion of them are false failures (i.e., only the replications reflect the true underlying effect)
         TSR2 = TS2/Successes2,
         FSR2 = FS2/Successes2,
         TFR2 = TF2/Failures2,
         FFR2 = FF2/Failures2,
         
         SAR1 = TS1/(TS1 + TF1), #Successful Affirmation Rate; When original p-value is correct, what's the probability that MABF is also correct
         FAR1 = TF1/(TS1 + TF1), #Failed Affirmation Rate; When original p-value is correct, what's the probability that MABF is not correct
         FCR1 = FS1/(FS1 + FF1), #Failed Correction Rate; When original p-value is not correct, what's the probability that MABF is not correct, either
         SCR1 = FF1/(FS1 + FF1), #Successful Correction Rate; When original p-value is not correct, what's the probability that MABF is correct
         SAR2 = TS2/(TS2 + TF2), 
         FAR2 = TF2/(TS2 + TF2), 
         FCR2 = FS2/(FS2 + FF2), 
         SCR2 = FF2/(FS2 + FF2)) %>% 
  select(-c(9:26)) %>% 
  filter(threshold != "Inf")

df_selected <- df
remove(df)

# Load data
#df <- read_csv("df_selected.csv")

# Pre-compute dropdown options
methods <- unique(df_selected$method)
orig_n_values <- unique(df_selected$orig.n)
bias_levels <- unique(df_selected$bias_level)
rep_nums <- unique(df_selected$`rep number`)
rep_n_values <- unique(df_selected$rep.n)
threshold_p_levels <- unique(df_selected$threshold_p)

# Define rate columns available for plotting
rate_choices <- c("FSR2","TSR1", "TPR", "FPR", 
                  "FSR1", "TFR1", "FFR1",
                  "TSR2",  "TFR2", "FFR2",
                  "SAR1", "FAR1", "SCR1", "FCR1",
                  "SAR2", "FAR2", "SCR2", "FCR2")

# UI
ui <- fluidPage(
  titlePanel("Interactive Rate Plot"),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,  # sidebar takes 4 columns
      selectInput("method", "Select method(s):",
                  choices = methods,
                  selected = methods[1],
                  multiple = TRUE),
      selectInput("orig_n", "Original study sample size:", choices = orig_n_values, selected = orig_n_values[1]),
      selectInput("bias", "Bias level:", choices = bias_levels, selected = bias_levels[1]),
      selectInput("rep_num", "Number of replications:", choices = rep_nums, selected = rep_nums[1]),
      selectInput("rep_n", "Replication sample size:", choices = rep_n_values, selected = rep_n_values[1]),
      selectInput("threshold_p", "Threshold p level:", choices = threshold_p_levels, selected = threshold_p_levels[1]),
      sliderInput("threshold", "Select threshold upper limit:",
                  min = min(df_selected$threshold, na.rm = TRUE),
                  max = max(df_selected$threshold, na.rm = TRUE),
                  value = min(df_selected$threshold, na.rm = TRUE),
                  step = 0.01,
                  animate = TRUE),
      selectInput("x_rate", "Select x-axis rate:", choices = rate_choices, selected = "FPR"),
      selectInput("y_rate", "Select y-axis rate:", choices = rate_choices, selected = "TPR"),
      actionButton("swap_axes", "Swap X and Y axes")
    ),
    mainPanel(
      width = 8,
      plotOutput("ratePlot", height = "500px", brush = brushOpts(id = "ratePlot_brush")),
      h4("Zoomed-in view"),
      plotOutput("zoomPlot", height = "500px")
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive filtered data with caching
  filtered_data <- reactive({
    req(input$threshold)
    df_selected %>%
      filter(
        method %in% input$method,
        orig.n == input$orig_n,
        bias_level == input$bias,
        `rep number` == input$rep_num,
        rep.n == input$rep_n,
        threshold_p == input$threshold_p,
        threshold <= input$threshold
      )
  }) %>% bindCache(
    input$method, input$orig_n, input$bias,
    input$rep_num, input$rep_n, input$threshold_p,
    input$threshold
  )
  
  # Swap x and y axes when button is clicked
  observeEvent(input$swap_axes, {
    current_x <- isolate(input$x_rate)
    current_y <- isolate(input$y_rate)
    
    updateSelectInput(session, inputId = "x_rate", selected = current_y)
    updateSelectInput(session, inputId = "y_rate", selected = current_x)
  })
  
  # Main plot
  output$ratePlot <- renderPlot({
    data <- filtered_data()
    
    ggplot(data, aes_string(x = input$x_rate, y = input$y_rate, color = "method")) +
      geom_step(linewidth = 0.5) +
      geom_point(size = 0.1) +
      labs(x = input$x_rate, y = input$y_rate, color = "Method") +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      theme_minimal()
  })
  
  # Zoomed-in plot based on brush
  output$zoomPlot <- renderPlot({
    data <- filtered_data()
    brush <- input$ratePlot_brush
    
    if (!is.null(brush)) {
      xmin <- brush$xmin
      xmax <- brush$xmax
      ymin <- brush$ymin
      ymax <- brush$ymax
      
      zoomed_data <- data %>%
        filter(
          !!sym(input$x_rate) >= xmin, !!sym(input$x_rate) <= xmax,
          !!sym(input$y_rate) >= ymin, !!sym(input$y_rate) <= ymax
        )
      
      # Manually set color palette based on all selected methods
      selected_methods <- input$method
      palette <- scales::hue_pal()(length(methods))  # consistent color palette
      names(palette) <- methods  # assign method names to colors
      
      ggplot(zoomed_data, aes_string(x = input$x_rate, y = input$y_rate, color = "method")) +
        geom_step(linewidth = 0.5) +
        geom_point(size = 0.1) +
        labs(x = paste("Zoomed", input$x_rate),
             y = paste("Zoomed", input$y_rate),
             color = "Method") +
        scale_color_manual(values = palette) +  # force consistent color mapping
        scale_x_continuous(limits = c(xmin, xmax), expand = c(0, 0)) +
        scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0)) +
        theme_minimal()
    }
  })
  
}



# Run the app
shinyApp(ui = ui, server = server)



