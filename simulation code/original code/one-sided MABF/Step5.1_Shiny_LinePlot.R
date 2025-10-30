library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(DT)

# # Set the path to the directory containing RDS files
# folder_path <- "./MABFanalyses/row-wise"
# # List all RDS files in the directory
# rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# # Read in all RDS files into the workspace
# for (file_path in rds_files) {
#   # Extract the base name without the extension
#   file_name <- tools::file_path_sans_ext(basename(file_path))
#   # Create a variable with the name of the file and assign the data from readRDS
#   assign(file_name, readRDS(file_path), envir = .GlobalEnv)
# }
# 
# 
# 
# # --- dataset groups ---
# datasets_0.2 <- c(
#   "MABF_rates_0.2_rowise_BFcutoff1",
#   "MABF_rates_0.2_rowise_BFcutoff3",
#   "MABF_rates_0.2_rowise_BFcutoff10",
#   "MABF_rates_0.2_rowise_BFcutoff30"
# )
# 
# datasets_0.5 <- c(
#   "MABF_rates_0.5_rowise_BFcutoff1",
#   "MABF_rates_0.5_rowise_BFcutoff3",
#   "MABF_rates_0.5_rowise_BFcutoff10",
#   "MABF_rates_0.5_rowise_BFcutoff30"
# )
# 
# datasets_null <- c(
#   "MABF_rates_null_rowise_BFcutoff1",
#   "MABF_rates_null_rowise_BFcutoff3",
#   "MABF_rates_null_rowise_BFcutoff10",
#   "MABF_rates_null_rowise_BFcutoff30"
# )
# 
# # --- reshape function ---
# reshape_dataset <- function(df, name) {
#   if (grepl("null", name)) {
#     df <- df %>% rename(ADR = ADR_null_anecdotal)
#     dv_cols <- c("ADR","TNR","FPR","TSR2","FSR2","TFR2","FFR2")
#   } else {
#     df <- df %>% rename(ADR = ADR_true_anecdotal)
#     dv_cols <- c("ADR","TPR","FNR","TSR1","FSR1","TFR1","FFR1")
#   }
#   
#   df %>%
#     pivot_longer(
#       cols = all_of(dv_cols),
#       names_to = "Outcome",
#       values_to = "Rate"
#     ) %>%
#     mutate(
#       Outcome = recode(Outcome,
#                        "TSR1"="TSR","TSR2"="TSR",
#                        "FSR1"="FSR","FSR2"="FSR",
#                        "TFR1"="TFR","TFR2"="TFR",
#                        "FFR1"="FFR","FFR2"="FFR"),
#       dataset = name
#     )
# }
# 
# # --- process dataset groups ---
# process_group <- function(dataset_names) {
#   bind_rows(lapply(dataset_names, function(nm) {
#     reshape_dataset(get(nm), nm)
#   })) %>%
#     mutate(
#       rep.number = factor(rep.number, levels = c(2, 5, 10)),
#       rep.n      = factor(rep.n, levels = c(40, 100, 400)),
#       orig.n     = factor(orig.n, levels = c(20, 50, 200)),
#       bias.level = factor(bias.level, levels = c("low", "medium", "high")),
#       orig.alpha = factor(orig.alpha, levels = c(0.01, 0.05)),
#       BFcutoff   = str_extract(dataset, "(?<=BFcutoff)\\d+"),
#       BFcutoff   = factor(BFcutoff, levels = c("1","3","10","30"))
#     ) %>%
#     select( -QRP.level, -PB.level, -SAR1, -FAR1, -FCR1, -SCR1,-TF1, -FS1, -TP, -FN, -TS1, -FF1, -TF1, -FS1, -starts_with("AD_"), -starts_with("ADR_"),-starts_with("ADOdds_"))
# }
# 
# process_group2 <- function(dataset_names) {
#   bind_rows(lapply(dataset_names, function(nm) {
#     reshape_dataset(get(nm), nm)
#   })) %>%
#     mutate(
#       rep.number = factor(rep.number, levels = c(2, 5, 10)),
#       rep.n      = factor(rep.n, levels = c(40, 100, 400)),
#       orig.n     = factor(orig.n, levels = c(20, 50, 200)),
#       bias.level = factor(bias.level, levels = c("low", "medium", "high")),
#       orig.alpha = factor(orig.alpha, levels = c(0.01, 0.05)),
#       BFcutoff   = str_extract(dataset, "(?<=BFcutoff)\\d+"),
#       BFcutoff   = factor(BFcutoff, levels = c("1","3","10","30"))
#     ) %>%
#     select( -QRP.level, -PB.level, -SAR2, -FAR2, -FCR2, -SCR2,-TF2, -FS2, -TN, -FP, -TS2, -FF2, -TF2, -FS2, -starts_with("AD_"), -starts_with("ADR_"),-starts_with("ADOdds_"))
# }
# 
# all_data_0.2  <- process_group(datasets_0.2)
# all_data_0.5  <- process_group(datasets_0.5)
# all_data_null <- process_group2(datasets_null)

load("./MABFanalyses/row-wise/LinePlotData.RData")

# --- helper for summary table (and plot data) ---
summarise_plot_data <- function(data, outcome_var, x_var, cutoff_filter, method_filter) {
  df <- data %>% filter(Outcome == outcome_var,
                        BFcutoff %in% cutoff_filter)
  if (length(method_filter) > 0) {
    df <- df %>% filter(method %in% method_filter)
  }
  
  df %>%
    group_by(method, !!sym(x_var), Outcome, BFcutoff, true.effect) %>%
    summarise(
      mean_rate   = round(mean(Rate, na.rm = TRUE), 2),
      sd_rate     = round(sd(Rate, na.rm = TRUE), 2),
      prop_missing = round(mean(is.na(Rate)), 4),  
      .groups     = "drop"
    )%>%
    mutate(mean_rate = round(mean_rate, 2))
}

# --- plotting function ---
plot_outcome_avg <- function(df, outcome_var, x_var) {
  cutoff_vals <- unique(df$BFcutoff)
  eff_vals    <- unique(df$true.effect)
  
  # APA-style labels
  x_labels <- c(
    "rep.number" = "Number of Replications",
    "rep.n"      = "Replication Sample Size",
    "orig.n"     = "Original Sample Size",
    "bias.level" = "Bias Level",
    "orig.alpha" = "Original α Level"
  )
  
  ggplot(df, aes_string(x = x_var, y = "mean_rate", color = "method", group = "method")) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
    labs(
      title = paste0(outcome_var, " across ", x_labels[x_var],
                     " (BF cutoff = ", cutoff_vals,
                     ", true effect = ", eff_vals, ")"),
      x = x_labels[x_var],
      y = paste0("Mean ", outcome_var),
      color = "Method"
    ) +
    theme_minimal(base_size = 14)
}

# --- Shiny app ---
ui <- fluidPage(
  titlePanel("MABF Simulation Results"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset_group", "Choose Dataset:",
                  choices = c("Null"="null", "Effect size 0.2"="0.2", "Effect size 0.5"="0.5")),
      selectInput("outcome", "Choose Outcome:",
                  choices = c("ADR","TPR","FNR","TNR","FPR","TSR","FSR","TFR","FFR")),
      selectInput("xvar", "X-axis Variable:",
                  choices = c("rep.number","rep.n","orig.n","bias.level","orig.alpha")),
      selectInput("cutoff", "Bayes Factor Cutoff:", choices = c("1","3","10","30")),
      checkboxGroupInput("methods", "Methods:",
                         choices = c("FEMABF","EUBF","BFbMA","iBF"),
                         selected = c("FEMABF","EUBF","BFbMA","iBF")),
      checkboxInput("show_boxplot", "Show Boxplots", value = TRUE)
    ),
    mainPanel(
      plotOutput("plot"),
      h3("Exact Values"),
      DTOutput("table")
    )
  )
)

server <- function(input, output) {
  
  # raw data (needed for boxplots)
  reactive_raw <- reactive({
    data <- switch(input$dataset_group,
                   "0.2"  = all_data_0.2,
                   "0.5"  = all_data_0.5,
                   "null" = all_data_null)
    
    data %>%
      filter(Outcome == input$outcome,
             BFcutoff %in% input$cutoff,
             method %in% input$methods)
  })
  
  # summary data (needed for means table + line plot)
  reactive_data <- reactive({
    df <- reactive_raw()
    df %>%
      group_by(method, !!sym(input$xvar), Outcome, BFcutoff, true.effect) %>%
      summarise(
        mean_rate   = round(mean(Rate, na.rm = TRUE), 2),
        sd_rate     = round(sd(Rate, na.rm = TRUE), 2),
        prop_missing = round(mean(is.na(Rate)), 4),  
        .groups     = "drop"
      )
  })
  
  output$plot <- renderPlot({
    df_raw  <- reactive_raw()
    df_mean <- reactive_data()
    
    # APA-style labels
    x_labels <- c(
      "rep.number" = "Number of Replications",
      "rep.n"      = "Replication Sample Size",
      "orig.n"     = "Original Sample Size",
      "bias.level" = "Bias Level",
      "orig.alpha" = "Original α Level"
    )
    
    p <- ggplot() 
    
    # Add boxplots if checkbox is ticked
    if (input$show_boxplot) {
      p <- p + geom_boxplot(data = df_raw,
                            aes_string(x = input$xvar, y = "Rate", fill = "method",
                                       group = paste0("interaction(method,", input$xvar, ")")),
                            alpha = 0.2, width = 0.3, position = position_dodge(width = 0.6),
                            outlier.shape = NA, color = NA)
    }
    
    # Add mean line + points
    p <- p +
      geom_line(data = df_mean,
                aes_string(x = input$xvar, y = "mean_rate", color = "method", group = "method"),
                size = 1, position = position_dodge(width = 0.6)) +
      geom_point(data = df_mean,
                 aes_string(x = input$xvar, y = "mean_rate", color = "method"),
                 size = 3, position = position_dodge(width = 0.6)) +
      scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
      labs(
        title = paste0(input$outcome, " across ", x_labels[input$xvar],
                       " (BF cutoff = ", unique(df_mean$BFcutoff),
                       ", true effect = ", unique(df_mean$true.effect), ")"),
        x = x_labels[input$xvar],
        y = paste0("Rate / Mean ", input$outcome),
        color = "Method",
        fill = "Method"
      ) +
      theme_minimal(base_size = 14)
    
    p
  })
  
  
  output$table <- renderDT({
    df <- reactive_data()
    datatable(df, options = list(pageLength = 25), rownames = FALSE)
  })
}


shinyApp(ui, server)


