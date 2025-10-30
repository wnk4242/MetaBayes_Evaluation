#This Shiny app visualizes performance metrics (such as TPR, FSR2, etc.) under varying conditions in the shape of curves
##############
######Version 1: Only row and column facets
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
rates_0.2null <- rbind(BFbMA_TFPS4Viz_0.2null,iBF_TFPS4Viz_0.2null)
# Change threshold_p_r to threshold
#FEMA_TSFS4RepCurve_0.2null <- FEMA_TSFS4RepCurve_0.2null %>% 
#  rename(threshold = threshold_p_r)
# Change threshold_BF to threshold
rates_0.2null <- rates_0.2null  %>% 
  rename(threshold = threshold_BF)
# Combine rows of MABF methods and FEMA
#rates_0.2null <- rbind(FEMA_TSFS4RepCurve_0.2null, rates_0.2null)
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
#remove(df)
# === Preprocessing ===
df_selected <- df_selected %>%
  rename(rep_num = `rep number`) %>%
  mutate(
    orig.n = factor(orig.n, levels = c(20, 50 ,200)),
    rep.n = factor(rep.n, levels = c(40, 100, 400)),
    bias_level = factor(bias_level, levels = c("low", "medium", "high")),
    method = factor(method),
    rep_num = factor(rep_num, levels = c(2, 5, 10)),
    threshold_p = factor(threshold_p, levels = c(0.01, 0.05, 0.1)),
    thresh_bin = ntile(threshold, 50)
  )

rate_choices <- c("TPR", "FPR", 
                  "TSR1", "FSR1", "TFR1", "FFR1",
                  "TSR2", "FSR2", "TFR2", "FFR2",
                  "SAR1", "FAR1", "SCR1", "FCR1",
                  "SAR2", "FAR2", "SCR2", "FCR2")
facet_vars <- c("None", "orig.n", "bias_level", "rep_num", "rep.n", "threshold_p")

# === UI ===
ui <- fluidPage(
  titlePanel(textOutput("plotTitle")),
  sidebarLayout(
    sidebarPanel(
      selectInput("method", "Select method(s):",
                  choices = unique(df_selected$method),
                  selected = unique(df_selected$method)[1],
                  multiple = TRUE),
      selectInput("x_rate", "X-axis rate:", choices = rate_choices, selected = "FCR2"),
      selectInput("y_rate", "Y-axis rate:", choices = rate_choices, selected = "SAR1"),
      actionButton("swap_xy", "Swap X and Y axes"),
      br(), br(), br(),
      selectInput("facet_row", "Facet Rows:", choices = facet_vars, selected = "rep.n"),
      selectInput("facet_col", "Facet Columns:", choices = facet_vars, selected = "bias_level"),
      actionButton("swap_facets", "Swap Row/Column Facets")
    ),
    mainPanel(
      width = 8,
      plotOutput("ratePlot", height = "600px")
    )
  )
)

# === Server ===
swap_xy_state <- reactiveVal(FALSE)
swap_facets_state <- reactiveVal(FALSE)
server <- function(input, output, session) {
  
  observeEvent(input$swap_xy, {
    swap_xy_state(!swap_xy_state())
  })
  observeEvent(input$swap_facets, {
    swap_facets_state(!swap_facets_state())
  })
  
  filtered_data <- reactive({
    group_vars <- c("method", "thresh_bin")
    row_var <- if (swap_facets_state()) input$facet_col else input$facet_row
    col_var <- if (swap_facets_state()) input$facet_row else input$facet_col
    
    if (row_var != "None") group_vars <- c(group_vars, row_var)
    if (col_var != "None" && col_var != row_var) group_vars <- c(group_vars, col_var)
    
    df_selected %>%
      filter(method %in% input$method) %>%
      group_by(across(all_of(group_vars))) %>%
      summarize(
        TPR = mean(TPR, na.rm = TRUE),
        FPR = mean(FPR, na.rm = TRUE),
        
        TSR1 = mean(TSR1, na.rm = TRUE),
        FSR1 = mean(FSR1, na.rm = TRUE),
        TFR1 = mean(TFR1, na.rm = TRUE),
        FFR1 = mean(FFR1, na.rm = TRUE),
        TSR2 = mean(TSR2, na.rm = TRUE),
        FSR2 = mean(FSR2, na.rm = TRUE),
        TFR2 = mean(TFR2, na.rm = TRUE),
        FFR2 = mean(FFR2, na.rm = TRUE),
        
        SAR1 = mean(SAR1, na.rm = TRUE),
        FAR1 = mean(FAR1, na.rm = TRUE),
        FCR1 = mean(FCR1, na.rm = TRUE),
        SCR1 = mean(SCR1, na.rm = TRUE),
        SAR2 = mean(SAR2, na.rm = TRUE),
        FAR2 = mean(FAR2, na.rm = TRUE),
        FCR2 = mean(FCR2, na.rm = TRUE),
        SCR2 = mean(SCR2, na.rm = TRUE),
        
        .groups = "drop"
      )
  })
  
  output$plotTitle <- renderText({
    x_label <- if (swap_xy_state()) input$y_rate else input$x_rate
    y_label <- if (swap_xy_state()) input$x_rate else input$y_rate
    paste(y_label, "vs", x_label)
  })
  
  output$ratePlot <- renderPlot({
    data <- filtered_data()
    
    xvar <- if (swap_xy_state()) input$y_rate else input$x_rate
    yvar <- if (swap_xy_state()) input$x_rate else input$y_rate
    
    row_facet <- if (swap_facets_state()) input$facet_col else input$facet_row
    col_facet <- if (swap_facets_state()) input$facet_row else input$facet_col
    
    rename_label <- function(var) {
      named <- c(
        "rep_num" = "Number of Replications",
        "rep.n" = "Replication Sample Size",
        "bias_level" = "Bias Level",
        "orig.n" = "Original Sample Size",
        "threshold_p" = "Original Alpha Level"
      )
      if (var %in% names(named)) named[[var]] else var
    }
    
    facet_title <- ""
    if (row_facet != "None" && col_facet != "None") {
      facet_title <- paste("Method Performance Across",
                           rename_label(row_facet), "and",
                           rename_label(col_facet))
    } else if (row_facet != "None") {
      facet_title <- paste("Method Performance Across",
                           rename_label(row_facet))
    } else if (col_facet != "None") {
      facet_title <- paste("Method Performance Across",
                           rename_label(col_facet))
    }
    
    p <- ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
      geom_path(linewidth = 1) +
      geom_point(size = 0.3, alpha = 0.8) +
      labs(
        x = paste0(xvar, " (Lower is Better)"),
        y = paste0(yvar, " (Higher is Better)"),
        title = facet_title,
        color = "Method"
      ) +
      theme_minimal()
    
    if (row_facet != "None" || col_facet != "None") {
      facet_formula <- as.formula(paste(
        ifelse(row_facet == "None", ".", row_facet),
        "~",
        ifelse(col_facet == "None", ".", col_facet)
      ))
      p <- p + facet_grid(facet_formula)
    }
    
    p
  })
}

# Run app
shinyApp(ui = ui, server = server)



############################
#Version 2: Add a third facet, but repeated lables and legends
library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)

# Load and prepare data
folder_path <- "./MABFanalyses/matrix-wise/rates4Viz/"
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
for (file_path in rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

rates_0.2null <- rbind(BFbMA_TFPS4Viz_0.2null, iBF_TFPS4Viz_0.2null)
rates_0.2null <- rates_0.2null %>% rename(threshold = threshold_BF)

rates_0.2null <- rates_0.2null %>%
  mutate(true_es = 0.2) %>% 
  mutate(across(c(method, true_es, orig.n, `QRP level`, `PB level`, `rep number`, rep.n), as.factor)) %>% 
  relocate(true_es, .after = `true effect`) %>% 
  relocate(group, .before = `method`) %>% 
  rename(scenario = group)

df <- rates_0.2null %>% 
  select(-`true effect`, -true_es, -threshold_ES) %>% 
  mutate(bias_level = `PB level`) %>% 
  select(-`QRP level`, -`PB level`) %>% 
  relocate(bias_level, .before = `rep number`) %>% 
  mutate(
    TPR = TP / (TP + FN),
    FPR = FP / (FP + TN),
    TSR1 = TS1 / Successes1,
    FSR1 = FS1 / Successes1,
    TFR1 = TF1 / Failures1,
    FFR1 = FF1 / Failures1,
    TSR2 = TS2 / Successes2,
    FSR2 = FS2 / Successes2,
    TFR2 = TF2 / Failures2,
    FFR2 = FF2 / Failures2,
    SAR1 = TS1 / (TS1 + TF1),
    FAR1 = TF1 / (TS1 + TF1),
    FCR1 = FS1 / (FS1 + FF1),
    SCR1 = FF1 / (FS1 + FF1),
    SAR2 = TS2 / (TS2 + TF2),
    FAR2 = TF2 / (TS2 + TF2),
    FCR2 = FS2 / (FS2 + FF2),
    SCR2 = FF2 / (FS2 + FF2)
  ) %>% 
  select(-c(9:26)) %>% 
  filter(threshold != "Inf")

df_selected <- df %>%
  rename(rep_num = `rep number`) %>%
  mutate(
    orig.n = factor(orig.n, levels = c(20, 50 ,200)),
    rep.n = factor(rep.n, levels = c(40, 100, 400),
                   labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400")),
    bias_level = factor(bias_level, levels = c("low", "medium", "high")),
    method = factor(method),
    rep_num = factor(rep_num, levels = c(2, 5, 10),
                     labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10")),
    threshold_p = factor(threshold_p, levels = c(0.01, 0.05, 0.1)),
    thresh_bin = ntile(threshold, 50)
  )

rate_choices <- c("TPR", "FPR", 
                  "TSR1", "FSR1", "TFR1", "FFR1",
                  "TSR2", "FSR2", "TFR2", "FFR2",
                  "SAR1", "FAR1", "SCR1", "FCR1",
                  "SAR2", "FAR2", "SCR2", "FCR2")
facet_vars <- c("None", "orig.n", "bias_level", "rep_num", "rep.n", "threshold_p")

# === UI ===
ui <- fluidPage(
  titlePanel(textOutput("plotTitle")),
  sidebarLayout(
    sidebarPanel(
      selectInput("method", "Select method(s):",
                  choices = unique(df_selected$method),
                  selected = unique(df_selected$method)[1],
                  multiple = TRUE),
      selectInput("x_rate", "X-axis rate:", choices = rate_choices, selected = "FCR2"),
      selectInput("y_rate", "Y-axis rate:", choices = rate_choices, selected = "SAR1"),
      actionButton("swap_xy", "Swap X and Y axes"),
      br(), br(), br(),
      selectInput("facet_row", "Facet Rows:", choices = facet_vars, selected = "rep.n"),
      selectInput("facet_col", "Facet Columns:", choices = facet_vars, selected = "rep_num"),
      selectInput("facet_3rd", "Third Facet (Right):", choices = facet_vars, selected = "bias_level"),
      actionButton("swap_facets", "Swap Row/Column Facets")
    ),
    mainPanel(
      width = 8,
      plotOutput("ratePlot", height = "800px")
    )
  )
)

# === Server ===
swap_xy_state <- reactiveVal(FALSE)
swap_facets_state <- reactiveVal(FALSE)

server <- function(input, output, session) {
  observeEvent(input$swap_xy, {
    swap_xy_state(!swap_xy_state())
  })
  
  observeEvent(input$swap_facets, {
    swap_facets_state(!swap_facets_state())
  })
  
  filtered_data <- reactive({
    group_vars <- c("method", "thresh_bin")
    row_var <- if (swap_facets_state()) input$facet_col else input$facet_row
    col_var <- if (swap_facets_state()) input$facet_row else input$facet_col
    third_var <- input$facet_3rd
    
    if (row_var != "None") group_vars <- c(group_vars, row_var)
    if (col_var != "None" && col_var != row_var) group_vars <- c(group_vars, col_var)
    if (third_var != "None" && !(third_var %in% group_vars)) group_vars <- c(group_vars, third_var)
    
    df_selected %>%
      filter(method %in% input$method) %>%
      group_by(across(all_of(group_vars))) %>%
      summarize(across(all_of(rate_choices), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  })
  
  output$plotTitle <- renderText({
    x_label <- if (swap_xy_state()) input$y_rate else input$x_rate
    y_label <- if (swap_xy_state()) input$x_rate else input$y_rate
    paste(y_label, "vs", x_label)
  })
  
  output$ratePlot <- renderPlot({
    data <- filtered_data()
    
    xvar <- if (swap_xy_state()) input$y_rate else input$x_rate
    yvar <- if (swap_xy_state()) input$x_rate else input$y_rate
    
    row_facet <- if (swap_facets_state()) input$facet_col else input$facet_row
    col_facet <- if (swap_facets_state()) input$facet_row else input$facet_col
    third_facet <- input$facet_3rd
    
    rename_label <- function(var) {
      named <- c(
        "rep_num" = "Number of Replications",
        "rep.n" = "Replication Sample Size",
        "bias_level" = "Bias Level",
        "orig.n" = "Original Sample Size",
        "threshold_p" = "Original Alpha Level"
      )
      if (var %in% names(named)) named[[var]] else var
    }
    
    facet_formula <- as.formula(paste(
      ifelse(row_facet == "None", ".", row_facet),
      "~",
      ifelse(col_facet == "None", ".", col_facet)
    ))
    
    if (third_facet == "None" || third_facet %in% c(row_facet, col_facet)) {
      ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
        geom_path(linewidth = 1) +
        geom_point(size = 0.3, alpha = 0.8) +
        labs(
          x = paste0(rename_label(xvar), " (Lower is Better)"),
          y = paste0(rename_label(yvar), " (Higher is Better)"),
          title = paste("Performance across", rename_label(row_facet), "and", rename_label(col_facet)),
          color = "Method"
        ) +
        facet_grid(facet_formula, labeller = labeller(rep.n = label_parsed, rep_num = label_parsed)) +
        theme_minimal() +
        theme(legend.position = "bottom")
    } else {
      third_levels <- unique(data[[third_facet]])
      
      plots <- lapply(third_levels, function(level_value) {
        df_subset <- data[data[[third_facet]] == level_value, , drop = FALSE]
        if (nrow(df_subset) == 0) return(NULL)
        
        ggplot(df_subset, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
          geom_path(linewidth = 1) +
          geom_point(size = 0.3, alpha = 0.8) +
          labs(
            x = paste0(rename_label(xvar), " (Lower is Better)"),
            y = paste0(rename_label(yvar), " (Higher is Better)"),
            title = paste(rename_label(third_facet), "=", level_value),
            color = "Method"
          ) +
          facet_grid(facet_formula, labeller = labeller(rep.n = label_parsed, rep_num = label_parsed)) +
          theme_minimal()
      })
      
      plots <- Filter(Negate(is.null), plots)
      
      wrap_plots(plots, ncol = length(plots)) +
        plot_annotation(
          theme = theme(
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12)
          )
        ) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
    }
  })
}

# Run app
shinyApp(ui = ui, server = server)



#######################


#Version 3: Basically what I want but the stacked facets are fixed to N_rep and n_rep
library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rlang)

# Load and prepare data
folder_path <- "./MABFanalyses/matrix-wise/rates4Viz/"
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
for (file_path in rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

rates_0.2null <- rbind(BFbMA_TFPS4Viz_0.2null, iBF_TFPS4Viz_0.2null)
rates_0.2null <- rates_0.2null %>% rename(threshold = threshold_BF)

rates_0.2null <- rates_0.2null %>%
  mutate(true_es = 0.2) %>% 
  mutate(across(c(method, true_es, orig.n, `QRP level`, `PB level`, `rep number`, rep.n), as.factor)) %>% 
  relocate(true_es, .after = `true effect`) %>% 
  relocate(group, .before = `method`) %>% 
  rename(scenario = group)

df <- rates_0.2null %>% 
  select(-`true effect`, -true_es, -threshold_ES) %>% 
  mutate(bias_level = `PB level`) %>% 
  select(-`QRP level`, -`PB level`) %>% 
  relocate(bias_level, .before = `rep number`) %>% 
  mutate(
    TPR = TP / (TP + FN),
    FPR = FP / (FP + TN),
    TSR1 = TS1 / Successes1,
    FSR1 = FS1 / Successes1,
    TFR1 = TF1 / Failures1,
    FFR1 = FF1 / Failures1,
    TSR2 = TS2 / Successes2,
    FSR2 = FS2 / Successes2,
    TFR2 = TF2 / Failures2,
    FFR2 = FF2 / Failures2,
    SAR1 = TS1 / (TS1 + TF1),
    FAR1 = TF1 / (TS1 + TF1),
    FCR1 = FS1 / (FS1 + FF1),
    SCR1 = FF1 / (FS1 + FF1),
    SAR2 = TS2 / (TS2 + TF2),
    FAR2 = TF2 / (TS2 + TF2),
    FCR2 = FS2 / (FS2 + FF2),
    SCR2 = FF2 / (FS2 + FF2)
  ) %>% 
  select(-c(9:26)) %>% 
  filter(threshold != "Inf")

df_selected <- df %>%
  rename(rep_num = `rep number`) %>%
  mutate(
    orig.n = factor(orig.n, levels = c(20, 50 ,200)),
    rep.n = factor(rep.n, levels = c(40, 100, 400)),
    rep_num = factor(rep_num, levels = c(2, 5, 10)),
    bias_level = factor(bias_level, levels = c("low", "medium", "high")),
    method = factor(method),
    threshold_p = factor(threshold_p, levels = c(0.01, 0.05, 0.1)),
    thresh_bin = ntile(threshold, 50),
    
    # ✅ stacked facet label using atop()
    facet_col_combined = forcats::fct_relevel(
      paste0("atop(n[rep]==", rep.n, ", N[rep]==", rep_num, ")"),
      paste0("atop(n[rep]==", c(40, 100, 400), ", N[rep]==", 2, ")"),
      paste0("atop(n[rep]==", c(40, 100, 400), ", N[rep]==", 5, ")"),
      paste0("atop(n[rep]==", c(40, 100, 400), ", N[rep]==", 10, ")")
    )
  )

rate_choices <- c("TPR", "FPR", 
                  "TSR1", "FSR1", "TFR1", "FFR1",
                  "TSR2", "FSR2", "TFR2", "FFR2",
                  "SAR1", "FAR1", "SCR1", "FCR1",
                  "SAR2", "FAR2", "SCR2", "FCR2")

facet_vars <- c("None", "orig.n", "bias_level", "threshold_p")

# === UI ===
ui <- fluidPage(
  titlePanel(textOutput("plotTitle")),
  sidebarLayout(
    sidebarPanel(
      selectInput("method", "Select method(s):",
                  choices = unique(df_selected$method),
                  selected = unique(df_selected$method)[1],
                  multiple = TRUE),
      selectInput("x_rate", "X-axis rate:", choices = rate_choices, selected = "FCR2"),
      selectInput("y_rate", "Y-axis rate:", choices = rate_choices, selected = "SAR1"),
      actionButton("swap_xy", "Swap X and Y axes"),
      br(), br(), br(),
      selectInput("facet_row", "Facet Rows:", choices = facet_vars, selected = "bias_level"),
      selectInput("facet_3rd", "Third Facet (Vertical Panels):", choices = facet_vars, selected = "None"),
      actionButton("swap_facets", "Swap Row/Column Facets")
    ),
    mainPanel(
      width = 8,
      plotOutput("ratePlot", height = "800px")
    )
  )
)

# === Server ===
swap_xy_state <- reactiveVal(FALSE)
swap_facets_state <- reactiveVal(FALSE)

server <- function(input, output, session) {
  observeEvent(input$swap_xy, {
    swap_xy_state(!swap_xy_state())
  })
  
  observeEvent(input$swap_facets, {
    swap_facets_state(!swap_facets_state())
  })
  
  filtered_data <- reactive({
    group_vars <- c("method", "thresh_bin", "facet_col_combined")
    row_var <- input$facet_row
    third_var <- input$facet_3rd
    
    if (row_var != "None") group_vars <- c(group_vars, row_var)
    if (third_var != "None" && !(third_var %in% group_vars)) group_vars <- c(group_vars, third_var)
    
    df_selected %>%
      filter(method %in% input$method) %>%
      group_by(across(all_of(group_vars))) %>%
      summarize(across(all_of(rate_choices), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  })
  
  output$plotTitle <- renderText({
    x_label <- if (swap_xy_state()) input$y_rate else input$x_rate
    y_label <- if (swap_xy_state()) input$x_rate else input$y_rate
    paste(y_label, "vs", x_label)
  })
  
  output$ratePlot <- renderPlot({
    data <- filtered_data()
    
    xvar <- if (swap_xy_state()) input$y_rate else input$x_rate
    yvar <- if (swap_xy_state()) input$x_rate else input$y_rate
    
    row_facet <- input$facet_row
    third_facet <- input$facet_3rd
    
    if (third_facet == "None" || third_facet == row_facet) {
      ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
        geom_path(linewidth = 1) +
        geom_point(size = 0.3, alpha = 0.8) +
        labs(
          x = paste0(xvar, " (Lower is Better)"),
          y = paste0(yvar, " (Higher is Better)"),
          color = "Method"
        ) +
        facet_grid(
          rows = vars(!!sym(row_facet)),
          cols = vars(facet_col_combined),
          labeller = label_parsed,
          switch = "both"
        ) +
        theme_minimal() +
        theme(
          strip.placement = "outside",
          legend.position = "bottom"
        )
    } else {
      third_levels <- unique(data[[third_facet]])
      
      plots <- lapply(third_levels, function(level_value) {
        df_subset <- data[data[[third_facet]] == level_value, , drop = FALSE]
        if (nrow(df_subset) == 0) return(NULL)
        
        ggplot(df_subset, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
          geom_path(linewidth = 1) +
          geom_point(size = 0.3, alpha = 0.8) +
          labs(
            x = paste0(xvar, " (Lower is Better)"),
            y = paste0(yvar, " (Higher is Better)"),
            title = paste(third_facet, "=", level_value),
            color = "Method"
          ) +
          facet_grid(
            rows = vars(!!sym(row_facet)),
            cols = vars(facet_col_combined),
            labeller = label_parsed,
            switch = "both"
          ) +
          theme_minimal() +
          theme(strip.placement = "outside")
      })
      
      plots <- Filter(Negate(is.null), plots)
      
      wrap_plots(plots, ncol = 1) +  # vertical stacking of third facet
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
    }
  })
}

# Run app
shinyApp(ui = ui, server = server)





##################
#Version 4: we can select stacked facets but this version runs very slow
library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rlang)
library(forcats)

# === Load and Prepare Data (same as before) ===
folder_path <- "./MABFanalyses/matrix-wise/rates4Viz/"
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
for (file_path in rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

rates_0.2null <- rbind(BFbMA_TFPS4Viz_0.2null, iBF_TFPS4Viz_0.2null)
rates_0.2null <- rates_0.2null %>% rename(threshold = threshold_BF)

rates_0.2null <- rates_0.2null %>%
  mutate(true_es = 0.2) %>% 
  mutate(across(c(method, true_es, orig.n, `QRP level`, `PB level`, `rep number`, rep.n), as.factor)) %>% 
  relocate(true_es, .after = `true effect`) %>% 
  relocate(group, .before = `method`) %>% 
  rename(scenario = group)

df <- rates_0.2null %>% 
  select(-`true effect`, -true_es, -threshold_ES) %>% 
  mutate(bias_level = `PB level`) %>% 
  select(-`QRP level`, -`PB level`) %>% 
  relocate(bias_level, .before = `rep number`) %>% 
  mutate(
    TPR = TP / (TP + FN),
    FPR = FP / (FP + TN),
    TSR1 = TS1 / Successes1,
    FSR1 = FS1 / Successes1,
    TFR1 = TF1 / Failures1,
    FFR1 = FF1 / Failures1,
    TSR2 = TS2 / Successes2,
    FSR2 = FS2 / Successes2,
    TFR2 = TF2 / Failures2,
    FFR2 = FF2 / Failures2,
    SAR1 = TS1 / (TS1 + TF1),
    FAR1 = TF1 / (TS1 + TF1),
    FCR1 = FS1 / (FS1 + FF1),
    SCR1 = FF1 / (FS1 + FF1),
    SAR2 = TS2 / (TS2 + TF2),
    FAR2 = TF2 / (TS2 + TF2),
    FCR2 = FS2 / (FS2 + FF2),
    SCR2 = FF2 / (FS2 + FF2)
  ) %>% 
  select(-c(9:26)) %>% 
  filter(threshold != "Inf")

df_selected <- df %>%
  rename(rep_num = `rep number`) %>%
  mutate(
    orig.n = factor(orig.n, levels = c(20, 50 ,200)),
    rep.n = factor(rep.n, levels = c(40, 100, 400)),
    rep_num = factor(rep_num, levels = c(2, 5, 10)),
    bias_level = factor(bias_level, levels = c("low", "medium", "high")),
    method = factor(method),
    threshold_p = factor(threshold_p, levels = c(0.01, 0.05, 0.1)),
    thresh_bin = ntile(threshold, 50)
  )

rate_choices <- c("TPR", "FPR", 
                  "TSR1", "FSR1", "TFR1", "FFR1",
                  "TSR2", "FSR2", "TFR2", "FFR2",
                  "SAR1", "FAR1", "SCR1", "FCR1",
                  "SAR2", "FAR2", "SCR2", "FCR2")

facet_vars <- c("None", "orig.n", "bias_level", "rep.n", "rep_num")
facet_col_vars <- c("rep.n", "rep_num", "orig.n", "bias_level")

# === UI ===
ui <- fluidPage(
  titlePanel(textOutput("plotTitle")),
  sidebarLayout(
    sidebarPanel(
      selectInput("method", "Select method(s):",
                  choices = unique(df_selected$method),
                  selected = unique(df_selected$method)[1],
                  multiple = TRUE),
      selectInput("x_rate", "X-axis rate:", choices = rate_choices, selected = "FSR2"),
      selectInput("y_rate", "Y-axis rate:", choices = rate_choices, selected = "TSR1"),
      actionButton("swap_xy", "Swap X and Y axes"),
      br(), br(),
      selectInput("facet_row", "Facet Rows:", choices = facet_vars, selected = "bias_level"),
      selectInput("facet_3rd", "Third Facet (Vertical Panels):", choices = facet_vars, selected = "None"),
      selectInput("facet_col_var1", "Facet Column Variable 1 (Bottom Stack):",
                  choices = facet_col_vars, selected = "rep.n"),
      selectInput("facet_col_var2", "Facet Column Variable 2 (Above):",
                  choices = facet_col_vars, selected = "rep_num"),
      actionButton("swap_facets", "Swap Row/Column Facets")
    ),
    mainPanel(
      width = 8,
      plotOutput("ratePlot", height = "800px")
    )
  )
)

# === Server ===
swap_xy_state <- reactiveVal(FALSE)
swap_facets_state <- reactiveVal(FALSE)

server <- function(input, output, session) {
  observeEvent(input$swap_xy, {
    swap_xy_state(!swap_xy_state())
  })
  
  observeEvent(input$swap_facets, {
    swap_facets_state(!swap_facets_state())
  })
  
  pretty_var_labels <- list(
    "rep.n" = "n[rep]",
    "rep_num" = "N[rep]",
    "orig.n" = "n[orig]",
    "bias_level" = "Bias"
  )
  
  filtered_data <- reactive({
    row_var <- input$facet_row
    third_var <- input$facet_3rd
    col_vars <- c(input$facet_col_var1, input$facet_col_var2)
    
    group_vars <- c("method", "thresh_bin")
    if (!is.null(col_vars)) group_vars <- c(group_vars, col_vars)
    if (row_var != "None") group_vars <- c(group_vars, row_var)
    if (third_var != "None" && !(third_var %in% group_vars)) group_vars <- c(group_vars, third_var)
    
    df_filtered <- df_selected %>%
      filter(method %in% input$method)
    
    if (!is.null(col_vars) && length(col_vars) > 0) {
      df_filtered <- df_filtered %>%
        mutate(facet_col_combined = purrr::pmap_chr(across(all_of(col_vars)), function(...) {
          vals <- list(...)
          labels <- mapply(function(v, val) {
            pretty <- pretty_var_labels[[v]]
            if (is.null(pretty)) pretty <- v
            paste0(pretty, "==", val)
          }, v = col_vars, val = vals, SIMPLIFY = TRUE)
          paste0("atop(", paste(labels, collapse = ", "), ")")
        })) %>%
        mutate(facet_col_combined = forcats::fct_reorder(
          facet_col_combined,
          as.numeric(as.character(rep_num)) * 1000 + as.numeric(as.character(rep.n))
        ))
    } else {
      df_filtered$facet_col_combined <- "atop()"
    }
    
    df_filtered %>%
      group_by(across(all_of(c(group_vars, "facet_col_combined")))) %>%
      summarize(across(all_of(rate_choices), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  })
  
  output$plotTitle <- renderText({
    x_label <- if (swap_xy_state()) input$y_rate else input$x_rate
    y_label <- if (swap_xy_state()) input$x_rate else input$y_rate
    paste(y_label, "vs", x_label)
  })
  
  output$ratePlot <- renderPlot({
    data <- filtered_data()
    
    xvar <- if (swap_xy_state()) input$y_rate else input$x_rate
    yvar <- if (swap_xy_state()) input$x_rate else input$y_rate
    
    row_facet <- input$facet_row
    third_facet <- input$facet_3rd
    
    if (third_facet == "None" || third_facet == row_facet) {
      ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
        geom_path(linewidth = 1) +
        geom_point(size = 0.3, alpha = 0.8) +
        labs(
          x = paste0(xvar, " (Lower is Better)"),
          y = paste0(yvar, " (Higher is Better)"),
          color = "Method"
        ) +
        facet_grid(
          rows = vars(!!sym(row_facet)),
          cols = vars(facet_col_combined),
          labeller = label_parsed,
          switch = "both"
        ) +
        theme_minimal() +
        theme(
          strip.placement = "outside",
          legend.position = "bottom"
        )
    } else {
      third_levels <- unique(data[[third_facet]])
      plots <- lapply(third_levels, function(level_value) {
        df_subset <- data[data[[third_facet]] == level_value, , drop = FALSE]
        if (nrow(df_subset) == 0) return(NULL)
        
        ggplot(df_subset, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
          geom_path(linewidth = 1) +
          geom_point(size = 0.3, alpha = 0.8) +
          labs(
            x = paste0(xvar, " (Lower is Better)"),
            y = paste0(yvar, " (Higher is Better)"),
            title = paste(third_facet, "=", level_value),
            color = "Method"
          ) +
          facet_grid(
            rows = vars(!!sym(row_facet)),
            cols = vars(facet_col_combined),
            labeller = label_parsed,
            switch = "both"
          ) +
          theme_minimal() +
          theme(strip.placement = "outside")
      })
      
      plots <- Filter(Negate(is.null), plots)
      wrap_plots(plots, ncol = 1) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
    }
  })
}

# Run the app
shinyApp(ui = ui, server = server)



#######
#version 5: this will precalculate all combinations of variables, which occupy large memory.
library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rlang)
library(forcats)

# === Load and Prepare Data ===
folder_path <- "./MABFanalyses/matrix-wise/rates4Viz/"
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
for (file_path in rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

# Combine datasets of BFMA methods
rates_0.2null <- rbind(BFbMA_TFPS4Viz_0.2null, EUBF_TFPS4Viz_0.2null, FEMABF_TFPS4Viz_0.2null, iBF_TFPS4Viz_0.2null, iBF2_TFPS4Viz_0.2null)
# Change threshold_p_r to threshold
FEMA_TFPS4Viz_0.2null <- FEMA_TFPS4Viz_0.2null %>% 
  rename(threshold = threshold_p_r)
# Change threshold_BF to threshold
rates_0.2null <- rates_0.2null  %>% 
  rename(threshold = threshold_BF)
# Combine rows of MABF methods and FEMA
rates_0.2null <- rbind(FEMA_TFPS4Viz_0.2null, rates_0.2null)

rates_0.2null <- rates_0.2null %>%
  mutate(true_es = 0.2) %>% 
  mutate(across(c(method, true_es, orig.n, `QRP level`, `PB level`, `rep number`, rep.n), as.factor)) %>% 
  relocate(true_es, .after = `true effect`) %>% 
  relocate(group, .before = `method`) %>% 
  rename(scenario = group)

df <- rates_0.2null %>% 
  select(-`true effect`, -true_es, -threshold_ES) %>% 
  mutate(bias_level = `PB level`) %>% 
  select(-`QRP level`, -`PB level`) %>% 
  relocate(bias_level, .before = `rep number`) %>% 
  mutate(
    TPR = TP / (TP + FN),
    FPR = FP / (FP + TN),
    TSR1 = TS1 / Successes1,
    FSR1 = FS1 / Successes1,
    TFR1 = TF1 / Failures1,
    FFR1 = FF1 / Failures1,
    TSR2 = TS2 / Successes2,
    FSR2 = FS2 / Successes2,
    TFR2 = TF2 / Failures2,
    FFR2 = FF2 / Failures2,
    SAR1 = TS1 / (TS1 + TF1),
    FAR1 = TF1 / (TS1 + TF1),
    FCR1 = FS1 / (FS1 + FF1),
    SCR1 = FF1 / (FS1 + FF1),
    SAR2 = TS2 / (TS2 + TF2),
    FAR2 = TF2 / (TS2 + TF2),
    FCR2 = FS2 / (FS2 + FF2),
    SCR2 = FF2 / (FS2 + FF2)
  ) %>% 
  select(-c(9:26)) %>% 
  filter(threshold != "Inf")

# Preprocess base dataframe
df_base <- df %>%
  rename(rep_num = `rep number`) %>%
  mutate(
    orig.n = factor(orig.n, levels = c(20, 50 ,200)),
    rep.n = factor(rep.n, levels = c(40, 100, 400)),
    rep_num = factor(rep_num, levels = c(2, 5, 10)),
    bias_level = factor(bias_level, levels = c("low", "medium", "high")),
    method = factor(method),
    threshold_p = factor(threshold_p, levels = c(0.01, 0.05, 0.1)),
    thresh_bin = ntile(threshold, 50)
  )

# === Precompute stacked facet labels ===
facet_col_vars <- c("rep.n", "rep_num", "orig.n", "bias_level", "threshold_p")
facet_col_labels <- list(
  rep.n = "n[rep]",
  rep_num = "N[rep]",
  orig.n = "n[orig]",
  bias_level = "Bias",
  threshold_p = "p[threshold]"
)

facet_col_pairs <- purrr::set_names(
  unlist(lapply(facet_col_vars, function(v1) lapply(facet_col_vars, function(v2) c(v1, v2))), recursive = FALSE),
  unlist(lapply(facet_col_vars, function(v1) paste0(v1, "_", facet_col_vars)))
)

df_selected_long <- purrr::map_dfr(names(facet_col_pairs), function(key) {
  vars <- facet_col_pairs[[key]]
  var1 <- vars[1]
  var2 <- vars[2]
  
  df_tmp <- df_base %>%
    mutate(
      facet_key = key,
      facet_col_combined = paste0("atop(", facet_col_labels[[var1]], "==", .data[[var1]], ", ", facet_col_labels[[var2]], "==", .data[[var2]], ")")
    ) %>%
    arrange(.data[[var1]], .data[[var2]]) %>%
    mutate(facet_col_combined = factor(facet_col_combined, levels = unique(facet_col_combined)))
})

# Add UI and server for Shiny
rate_choices <- c("TPR", "FPR", 
                  "TSR1", "FSR1", "TFR1", "FFR1",
                  "TSR2", "FSR2", "TFR2", "FFR2",
                  "SAR1", "FAR1", "SCR1", "FCR1",
                  "SAR2", "FAR2", "SCR2", "FCR2")

facet_vars <- c("None", facet_col_vars)

ui <- fluidPage(
  titlePanel(textOutput("plotTitle")),
  sidebarLayout(
    sidebarPanel(
      selectInput("method", "Select method(s):",
                  choices = unique(df_selected_long$method),
                  selected = unique(df_selected_long$method)[1],
                  multiple = TRUE),
      selectInput("x_rate", "X-axis rate:", choices = rate_choices, selected = "FSR2"),
      selectInput("y_rate", "Y-axis rate:", choices = rate_choices, selected = "TSR1"),
      actionButton("swap_xy", "Swap X and Y axes"),
      br(), br(),
      selectInput("facet_row", "Facet Rows:", choices = facet_vars, selected = "bias_level"),
      selectInput("facet_3rd", "Third Facet (Vertical Panels):", choices = facet_vars, selected = "None"),
      selectInput("facet_col_var1", "Facet Column Variable 1 (Bottom):",
                  choices = facet_col_vars, selected = "rep.n"),
      selectInput("facet_col_var2", "Facet Column Variable 2 (Above):",
                  choices = facet_col_vars, selected = "rep_num")
    ),
    mainPanel(
      width = 8,
      plotOutput("ratePlot", height = "800px")
    )
  )
)

swap_xy_state <- reactiveVal(FALSE)

server <- function(input, output, session) {
  observeEvent(input$swap_xy, {
    swap_xy_state(!swap_xy_state())
  })
  
  filtered_data <- reactive({
    key <- paste(input$facet_col_var1, input$facet_col_var2, sep = "_")
    row_var <- input$facet_row
    third_var <- input$facet_3rd
    
    df_filtered <- df_selected_long %>%
      filter(facet_key == key, method %in% input$method) %>%
      mutate(facet_row_combined = case_when(
        row_var == "rep.n" ~ paste0("n[rep]==", rep.n),
        row_var == "rep_num" ~ paste0("N[rep]==", rep_num),
        row_var == "orig.n" ~ paste0("n[orig]==", orig.n),
        row_var == "bias_level" ~ paste0("Bias==\"", bias_level, "\""),
        row_var == "threshold_p" ~ paste0("p[threshold]==", threshold_p),
        TRUE ~ NA_character_
      )) %>%
      mutate(facet_row_combined = factor(facet_row_combined, levels = unique(facet_row_combined)))
    
    group_vars <- c("method", "thresh_bin", "facet_col_combined", "facet_row_combined")
    if (third_var != "None" && !(third_var %in% group_vars)) group_vars <- c(group_vars, third_var)
    
    df_filtered %>%
      group_by(across(all_of(group_vars))) %>%
      summarize(across(all_of(rate_choices), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  })
  
  output$plotTitle <- renderText({
    x_label <- if (swap_xy_state()) input$y_rate else input$x_rate
    y_label <- if (swap_xy_state()) input$x_rate else input$y_rate
    paste(y_label, "vs", x_label)
  })
  
  output$ratePlot <- renderPlot({
    data <- filtered_data()
    xvar <- if (swap_xy_state()) input$y_rate else input$x_rate
    yvar <- if (swap_xy_state()) input$x_rate else input$y_rate
    third_facet <- input$facet_3rd
    
    if (third_facet == "None" || third_facet == input$facet_row) {
      ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
        geom_path(linewidth = 1) +
        geom_point(size = 0.3, alpha = 0.8) +
        labs(
          x = paste0(xvar, " (Lower is Better)"),
          y = paste0(yvar, " (Higher is Better)"),
          color = "Method"
        ) +
        facet_grid(
          rows = vars(facet_row_combined),
          cols = vars(facet_col_combined),
          labeller = label_parsed,
          switch = "both"
        ) +
        theme_minimal() +
        theme(
          strip.placement = "outside",
          legend.position = "bottom"
        )
    } else {
      third_levels <- unique(data[[third_facet]])
      plots <- lapply(third_levels, function(level_value) {
        df_subset <- data[data[[third_facet]] == level_value, , drop = FALSE]
        if (nrow(df_subset) == 0) return(NULL)
        
        ggplot(df_subset, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
          geom_path(linewidth = 1) +
          geom_point(size = 0.3, alpha = 0.8) +
          labs(
            x = paste0(xvar, " (Lower is Better)"),
            y = paste0(yvar, " (Higher is Better)"),
            title = paste(third_facet, "=", level_value),
            color = "Method"
          ) +
          facet_grid(
            rows = vars(facet_row_combined),
            cols = vars(facet_col_combined),
            labeller = label_parsed,
            switch = "both"
          ) +
          theme_minimal() +
          theme(strip.placement = "outside")
      })
      
      plots <- Filter(Negate(is.null), plots)
      wrap_plots(plots, ncol = 1) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
    }
  })
}

shinyApp(ui = ui, server = server)
###################
# Version 6: This is a good version, but lacking some user friendly features
library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rlang)
library(forcats)

# === Load and Prepare Data ===
folder_path <- "./MABFanalyses/matrix-wise/rates4Viz/"
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
for (file_path in rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

# Combine datasets of BFMA methods
rates_0.2null <- rbind(BFbMA_TFPS4Viz_0.2null, EUBF_TFPS4Viz_0.2null, FEMABF_TFPS4Viz_0.2null, iBF_TFPS4Viz_0.2null, iBF2_TFPS4Viz_0.2null)
# Change threshold_p_r to threshold
FEMA_TFPS4Viz_0.2null <- FEMA_TFPS4Viz_0.2null %>% 
  rename(threshold = threshold_p_r)
# Change threshold_BF to threshold
rates_0.2null <- rates_0.2null  %>% 
  rename(threshold = threshold_BF)
# Combine rows of MABF methods and FEMA
rates_0.2null <- rbind(FEMA_TFPS4Viz_0.2null, rates_0.2null)

rates_0.2null <- rates_0.2null %>%
  mutate(true_es = 0.2) %>% 
  mutate(across(c(method, true_es, orig.n, `QRP level`, `PB level`, `rep number`, rep.n), as.factor)) %>% 
  relocate(true_es, .after = `true effect`) %>% 
  relocate(group, .before = `method`) %>% 
  rename(scenario = group)

df <- rates_0.2null %>% 
  select(-`true effect`, -true_es, -threshold_ES) %>% 
  mutate(bias_level = `PB level`) %>% 
  select(-`QRP level`, -`PB level`) %>% 
  relocate(bias_level, .before = `rep number`) %>% 
  mutate(
    TPR = TP / (TP + FN),
    FPR = FP / (FP + TN),
    TSR1 = TS1 / Successes1,
    FSR1 = FS1 / Successes1,
    TFR1 = TF1 / Failures1,
    FFR1 = FF1 / Failures1,
    TSR2 = TS2 / Successes2,
    FSR2 = FS2 / Successes2,
    TFR2 = TF2 / Failures2,
    FFR2 = FF2 / Failures2,
    SAR1 = TS1 / (TS1 + TF1),
    FAR1 = TF1 / (TS1 + TF1),
    FCR1 = FS1 / (FS1 + FF1),
    SCR1 = FF1 / (FS1 + FF1),
    SAR2 = TS2 / (TS2 + TF2),
    FAR2 = TF2 / (TS2 + TF2),
    FCR2 = FS2 / (FS2 + FF2),
    SCR2 = FF2 / (FS2 + FF2)
  ) %>% 
  select(-c(9:26)) %>% 
  filter(threshold != "Inf")

# Preprocess base dataframe
df_base <- df %>%
  rename(rep_num = `rep number`) %>%
  mutate(
    orig.n = factor(orig.n, levels = c(20, 50 ,200)),
    rep.n = factor(rep.n, levels = c(40, 100, 400)),
    rep_num = factor(rep_num, levels = c(2, 5, 10)),
    bias_level = factor(bias_level, levels = c("low", "medium", "high")),
    method = factor(method),
    threshold_p = factor(threshold_p, levels = c(0.01, 0.05, 0.1)),
    thresh_bin = ntile(threshold, 50)
  )

# Label mapping
facet_col_vars <- c("rep.n", "rep_num", "orig.n", "bias_level", "threshold_p")
facet_col_labels <- list(
  rep.n = "n[rep]",
  rep_num = "N[rep]",
  orig.n = "n[orig]",
  bias_level = "Bias",
  threshold_p = "p[threshold]"
)

# === UI ===
rate_choices <- c("TPR", "FPR", 
                  "TSR1", "FSR1", "TFR1", "FFR1",
                  "TSR2", "FSR2", "TFR2", "FFR2",
                  "SAR1", "FAR1", "SCR1", "FCR1",
                  "SAR2", "FAR2", "SCR2", "FCR2")

facet_vars <- c("None", facet_col_vars)

ui <- fluidPage(
  titlePanel(textOutput("plotTitle")),
  sidebarLayout(
    sidebarPanel(
      selectInput("method", "Select method(s):",
                  choices = unique(df_base$method),
                  selected = unique(df_base$method)[1],
                  multiple = TRUE),
      selectInput("x_rate", "X-axis rate:", choices = rate_choices, selected = "FSR2"),
      selectInput("y_rate", "Y-axis rate:", choices = rate_choices, selected = "TSR1"),
      actionButton("swap_xy", "Swap X and Y axes"),
      br(), br(),
      selectInput("facet_row", "Facet Rows:", choices = facet_vars, selected = "bias_level"),
      selectInput("facet_3rd", "Third Facet (Vertical Panels):", choices = facet_vars, selected = "None"),
      selectInput("facet_col_var1", "Facet Column Variable 1 (Bottom):",
                  choices = facet_col_vars, selected = "rep.n"),
      selectInput("facet_col_var2", "Facet Column Variable 2 (Above):",
                  choices = facet_col_vars, selected = "rep_num")
    ),
    mainPanel(
      width = 8,
      plotOutput("ratePlot", height = "800px")
    )
  )
)

# === Server ===
swap_xy_state <- reactiveVal(FALSE)

server <- function(input, output, session) {
  observeEvent(input$swap_xy, {
    swap_xy_state(!swap_xy_state())
  })
  
  filtered_data <- reactive({
    var1 <- input$facet_col_var1
    var2 <- input$facet_col_var2
    row_var <- input$facet_row
    third_var <- input$facet_3rd
    
    df_tmp <- df_base %>%
      filter(method %in% input$method) %>%
      mutate(
        facet_col_combined = paste0("atop(", facet_col_labels[[var1]], "==", .data[[var1]], ", ",
                                    facet_col_labels[[var2]], "==", .data[[var2]], ")"),
        facet_row_combined = case_when(
          row_var == "rep.n" ~ paste0("n[rep]==", rep.n),
          row_var == "rep_num" ~ paste0("N[rep]==", rep_num),
          row_var == "orig.n" ~ paste0("n[orig]==", orig.n),
          row_var == "bias_level" ~ paste0("Bias==\"", bias_level, "\""),
          row_var == "threshold_p" ~ paste0("p[threshold]==", threshold_p),
          TRUE ~ NA_character_
        )
      ) %>%
      mutate(
        facet_col_combined = forcats::fct_inorder(facet_col_combined),
        facet_row_combined = forcats::fct_inorder(facet_row_combined)
      )
    
    group_vars <- c("method", "thresh_bin", "facet_col_combined", "facet_row_combined")
    if (third_var != "None") group_vars <- c(group_vars, third_var)
    
    df_tmp %>%
      group_by(across(all_of(group_vars))) %>%
      summarize(across(all_of(rate_choices), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  })
  
  output$plotTitle <- renderText({
    x_label <- if (swap_xy_state()) input$y_rate else input$x_rate
    y_label <- if (swap_xy_state()) input$x_rate else input$y_rate
    paste(y_label, "vs", x_label)
  })
  
  output$ratePlot <- renderPlot({
    data <- filtered_data()
    xvar <- if (swap_xy_state()) input$y_rate else input$x_rate
    yvar <- if (swap_xy_state()) input$x_rate else input$y_rate
    third_facet <- input$facet_3rd
    
    if (third_facet == "None" || third_facet == input$facet_row) {
      ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
        geom_path(linewidth = 1) +
        geom_point(size = 0.3, alpha = 0.8) +
        labs(
          x = paste0(xvar, " (Lower is Better)"),
          y = paste0(yvar, " (Higher is Better)"),
          color = "Method"
        ) +
        facet_grid(
          rows = vars(facet_row_combined),
          cols = vars(facet_col_combined),
          labeller = label_parsed,
          switch = "both"
        ) +
        theme_minimal() +
        theme(
          strip.placement = "outside",
          legend.position = "bottom"
        )
    } else {
      third_levels <- unique(data[[third_facet]])
      plots <- lapply(third_levels, function(level_value) {
        df_subset <- data[data[[third_facet]] == level_value, , drop = FALSE]
        if (nrow(df_subset) == 0) return(NULL)
        
        ggplot(df_subset, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
          geom_path(linewidth = 1) +
          geom_point(size = 0.3, alpha = 0.8) +
          labs(
            x = paste0(xvar, " (Lower is Better)"),
            y = paste0(yvar, " (Higher is Better)"),
            title = paste(third_facet, "=", level_value),
            color = "Method"
          ) +
          facet_grid(
            rows = vars(facet_row_combined),
            cols = vars(facet_col_combined),
            labeller = label_parsed,
            switch = "both"
          ) +
          theme_minimal() +
          theme(strip.placement = "outside")
      })
      
      plots <- Filter(Negate(is.null), plots)
      wrap_plots(plots, ncol = 1) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
    }
  })
}

shinyApp(ui = ui, server = server)




#######
# This version works well and include a third facet
library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rlang)
library(forcats)

# === Load and Prepare Data ===
folder_path <- "./MABFanalyses/matrix-wise/rates4Viz/"
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
for (file_path in rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

# Combine datasets of BFMA methods
rates_0.2null <- rbind(BFbMA_TFPS4Viz_0.2null, EUBF_TFPS4Viz_0.2null, FEMABF_TFPS4Viz_0.2null, iBF_TFPS4Viz_0.2null, iBF2_TFPS4Viz_0.2null)
# Change threshold_p_r to threshold
FEMA_TFPS4Viz_0.2null <- FEMA_TFPS4Viz_0.2null %>% 
  rename(threshold = threshold_p_r)
# Change threshold_BF to threshold
rates_0.2null <- rates_0.2null  %>% 
  rename(threshold = threshold_BF)
# Combine rows of MABF methods and FEMA
rates_0.2null <- rbind(FEMA_TFPS4Viz_0.2null, rates_0.2null)

rates_0.2null <- rates_0.2null %>%
  mutate(true_es = 0.2) %>% 
  mutate(across(c(method, true_es, orig.n, `QRP level`, `PB level`, `rep number`, rep.n), as.factor)) %>% 
  relocate(true_es, .after = `true effect`) %>% 
  relocate(group, .before = `method`) %>% 
  rename(scenario = group)

df <- rates_0.2null %>% 
  select(-`true effect`, -true_es, -threshold_ES) %>% 
  mutate(bias_level = `PB level`) %>% 
  select(-`QRP level`, -`PB level`) %>% 
  relocate(bias_level, .before = `rep number`) %>% 
  mutate(
    TPR = TP / (TP + FN),
    FPR = FP / (FP + TN),
    TSR1 = TS1 / Successes1,
    FSR1 = FS1 / Successes1,
    TFR1 = TF1 / Failures1,
    FFR1 = FF1 / Failures1,
    TSR2 = TS2 / Successes2,
    FSR2 = FS2 / Successes2,
    TFR2 = TF2 / Failures2,
    FFR2 = FF2 / Failures2,
    SAR1 = TS1 / (TS1 + TF1),
    FAR1 = TF1 / (TS1 + TF1),
    FCR1 = FS1 / (FS1 + FF1),
    SCR1 = FF1 / (FS1 + FF1),
    SAR2 = TS2 / (TS2 + TF2),
    FAR2 = TF2 / (TS2 + TF2),
    FCR2 = FS2 / (FS2 + FF2),
    SCR2 = FF2 / (FS2 + FF2)
  ) %>% 
  select(-c(9:26)) %>% 
  filter(threshold != "Inf")

# Preprocess base dataframe
df_base <- df %>%
  rename(rep_num = `rep number`) %>%
  mutate(
    orig.n = factor(orig.n, levels = c(20, 50 ,200)),
    rep.n = factor(rep.n, levels = c(40, 100, 400)),
    rep_num = factor(rep_num, levels = c(2, 5, 10)),
    bias_level = factor(bias_level, levels = c("low", "medium", "high")),
    method = factor(method),
    threshold_p = factor(threshold_p, levels = c(0.01, 0.05, 0.1)),
    thresh_bin = ntile(threshold, 50)
  )

# Label mapping
facet_col_vars <- c("None", "rep.n", "rep_num", "orig.n", "bias_level", "threshold_p")
facet_col_labels <- list(
  rep.n = "n[rep]",
  rep_num = "N[rep]",
  orig.n = "n[orig]",
  bias_level = "Bias",
  threshold_p = "p[threshold]"
)

# === UI ===
rate_choices <- c("TPR", "FPR", 
                  "TSR1", "FSR1", "TFR1", "FFR1",
                  "TSR2", "FSR2", "TFR2", "FFR2",
                  "SAR1", "FAR1", "SCR1", "FCR1",
                  "SAR2", "FAR2", "SCR2", "FCR2")

facet_vars <- c("None", facet_col_vars[-1])

ui <- fluidPage(
  titlePanel(textOutput("plotTitle")),
  sidebarLayout(
    sidebarPanel(
      selectInput("method", "Select method(s):",
                  choices = unique(df_base$method),
                  selected = unique(df_base$method)[1],
                  multiple = TRUE),
      selectInput("x_rate", "X-axis rate:", choices = rate_choices, selected = "FSR2"),
      selectInput("y_rate", "Y-axis rate:", choices = rate_choices, selected = "TSR1"),
      actionButton("swap_xy", "Swap X and Y axes"),
      br(), br(),
      selectInput("facet_row", "Facet Rows:", choices = facet_vars, selected = "bias_level"),
      selectInput("facet_3rd", "Third Facet (Vertical Panels):", choices = facet_vars, selected = "None"),
      selectInput("facet_col_var1", "Facet Column Variable 1 (Bottom):",
                  choices = facet_col_vars, selected = "rep.n"),
      selectInput("facet_col_var2", "Facet Column Variable 2 (Above):",
                  choices = facet_col_vars, selected = "rep_num"),
      actionButton("make_plot", "Make the graph"),
      downloadButton("downloadPlot", "Export Graph as PNG")
    ),
    mainPanel(
      width = 8,
      plotOutput("ratePlot", height = "800px")
    )
  )
)

# === Server ===
server <- function(input, output, session) {
  observeEvent(input$swap_xy, {
    current_x <- isolate(input$x_rate)
    current_y <- isolate(input$y_rate)
    
    updateSelectInput(session, "x_rate", selected = current_y)
    updateSelectInput(session, "y_rate", selected = current_x)
  })
  
  filtered_data <- eventReactive(input$make_plot, {
    var1 <- input$facet_col_var1
    var2 <- input$facet_col_var2
    row_var <- input$facet_row
    third_var <- input$facet_3rd
    
    df_tmp <- df_base %>%
      filter(method %in% input$method)
    
    if (var1 != "None" & var2 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0("atop(", facet_col_labels[[var1]], "==", .data[[var1]], ", ",
                                    facet_col_labels[[var2]], "==", .data[[var2]], ")")
      )
    } else if (var1 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0(facet_col_labels[[var1]], "==", .data[[var1]])
      )
    } else if (var2 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0(facet_col_labels[[var2]], "==", .data[[var2]])
      )
    } else {
      df_tmp <- df_tmp %>% mutate(facet_col_combined = "''")
    }
    
    df_tmp <- df_tmp %>% mutate(
      facet_row_combined = case_when(
        row_var == "rep.n" ~ paste0("n[rep]==", rep.n),
        row_var == "rep_num" ~ paste0("N[rep]==", rep_num),
        row_var == "orig.n" ~ paste0("n[orig]==", orig.n),
        row_var == "bias_level" ~ paste0("Bias==\"", bias_level, "\""),
        row_var == "threshold_p" ~ paste0("p[threshold]==", threshold_p),
        TRUE ~ NA_character_
      ),
      facet_col_combined = fct_inorder(facet_col_combined),
      facet_row_combined = fct_inorder(facet_row_combined)
    )
    
    group_vars <- c("method", "thresh_bin", "facet_col_combined", "facet_row_combined")
    if (third_var != "None") group_vars <- c(group_vars, third_var)
    
    df_tmp %>%
      group_by(across(all_of(group_vars))) %>%
      summarize(across(all_of(rate_choices), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  })
  
  output$plotTitle <- renderText({
    paste0("Y-axis: ", input$y_rate, "   |   X-axis: ", input$x_rate)
  })
  
  plot_event <- eventReactive(input$make_plot, {
    data <- filtered_data()
    req(nrow(data) > 0)
    
    xvar <- input$x_rate
    yvar <- input$y_rate
    third_facet <- input$facet_3rd
    
    if (third_facet == "None" || third_facet == input$facet_row) {
      ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
        geom_path(linewidth = 1) +
        geom_point(size = 0.3, alpha = 0.8) +
        labs(
          x = paste0(xvar, " (Lower is Better)"),
          y = paste0(yvar, " (Higher is Better)"),
          color = "Method"
        ) +
        facet_grid(
          rows = vars(facet_row_combined),
          cols = vars(facet_col_combined),
          labeller = label_parsed,
          switch = "both"
        ) +
        theme_minimal() +
        theme(
          strip.placement = "outside",
          legend.position = "bottom"
        )
    } else {
      third_levels <- unique(data[[third_facet]])
      plots <- lapply(third_levels, function(level_value) {
        df_subset <- data[data[[third_facet]] == level_value, , drop = FALSE]
        if (nrow(df_subset) == 0) return(NULL)
        
        ggplot(df_subset, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
          geom_path(linewidth = 1) +
          geom_point(size = 0.3, alpha = 0.8) +
          labs(
            x = paste0(xvar, " (Lower is Better)"),
            y = paste0(yvar, " (Higher is Better)"),
            title = paste(third_facet, "=", level_value),
            color = "Method"
          ) +
          facet_grid(
            rows = vars(facet_row_combined),
            cols = vars(facet_col_combined),
            labeller = label_parsed,
            switch = "both"
          ) +
          theme_minimal() +
          theme(strip.placement = "outside")
      })
      
      plots <- Filter(Negate(is.null), plots)
      wrap_plots(plots, ncol = 1) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
    }
  })
  
  output$ratePlot <- renderPlot({
    plot_event()
  })
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("replication_curve_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(file, plot = plot_event(), device = "png", width = 12, height = 8, dpi = 300, bg = "white")
    }
  )
}

shinyApp(ui = ui, server = server)


##########
# This version doesn't contain the third facet
library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rlang)
library(forcats)

# === Load and Prepare Data ===
folder_path <- "./MABFanalyses/matrix-wise/rates4Viz/"
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
for (file_path in rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

# Combine datasets of BFMA methods
rates_0.2null <- rbind(BFbMA_TFPS4Viz_0.2null, EUBF_TFPS4Viz_0.2null, FEMABF_TFPS4Viz_0.2null, iBF_TFPS4Viz_0.2null, iBF2_TFPS4Viz_0.2null)
# Change threshold_p_r to threshold
FEMA_TFPS4Viz_0.2null <- FEMA_TFPS4Viz_0.2null %>% 
  rename(threshold = threshold_p_r)
# Change threshold_BF to threshold
rates_0.2null <- rates_0.2null  %>% 
  rename(threshold = threshold_BF)
# Combine rows of MABF methods and FEMA
rates_0.2null <- rbind(FEMA_TFPS4Viz_0.2null, rates_0.2null)

rates_0.2null <- rates_0.2null %>%
  mutate(true_es = 0.2) %>% 
  mutate(across(c(method, true_es, orig.n, `QRP level`, `PB level`, `rep number`, rep.n), as.factor)) %>% 
  relocate(true_es, .after = `true effect`) %>% 
  relocate(group, .before = `method`) %>% 
  rename(scenario = group)

df <- rates_0.2null %>% 
  select(-`true effect`, -true_es, -threshold_ES) %>% 
  mutate(bias_level = `PB level`) %>% 
  select(-`QRP level`, -`PB level`) %>% 
  relocate(bias_level, .before = `rep number`) %>% 
  mutate(
    TPR = TP / (TP + FN),
    FPR = FP / (FP + TN),
    TSR1 = TS1 / Successes1,
    FSR1 = FS1 / Successes1,
    TFR1 = TF1 / Failures1,
    FFR1 = FF1 / Failures1,
    TSR2 = TS2 / Successes2,
    FSR2 = FS2 / Successes2,
    TFR2 = TF2 / Failures2,
    FFR2 = FF2 / Failures2,
    SAR1 = TS1 / (TS1 + TF1),
    FAR1 = TF1 / (TS1 + TF1),
    FCR1 = FS1 / (FS1 + FF1),
    SCR1 = FF1 / (FS1 + FF1),
    SAR2 = TS2 / (TS2 + TF2),
    FAR2 = TF2 / (TS2 + TF2),
    FCR2 = FS2 / (FS2 + FF2),
    SCR2 = FF2 / (FS2 + FF2)
  ) %>% 
  select(-c(9:26)) %>% 
  filter(threshold != "Inf")

# Preprocess base dataframe
df_base <- df %>%
  rename(rep_num = `rep number`) %>%
  mutate(
    orig.n = factor(orig.n, levels = c(20, 50 ,200)),
    rep.n = factor(rep.n, levels = c(40, 100, 400)),
    rep_num = factor(rep_num, levels = c(2, 5, 10)),
    bias_level = factor(bias_level, levels = c("low", "medium", "high")),
    method = factor(method),
    threshold_p = factor(threshold_p, levels = c(0.01, 0.05, 0.1)),
    thresh_bin = ntile(threshold, 50)
  )

# Label mapping
facet_col_vars <- c("None", "rep.n", "rep_num", "orig.n", "bias_level", "threshold_p")
facet_col_labels <- list(
  rep.n = "n[rep]",
  rep_num = "N[rep]",
  orig.n = "n[orig]",
  bias_level = "Bias",
  threshold_p = "p[threshold]"
)

# === UI ===
rate_choices <- c("TPR", "FPR", 
                  "TSR1", "FSR1", "TFR1", "FFR1",
                  "TSR2", "FSR2", "TFR2", "FFR2",
                  "SAR1", "FAR1", "SCR1", "FCR1",
                  "SAR2", "FAR2", "SCR2", "FCR2")

facet_vars <- c("None", facet_col_vars[-1])

ui <- fluidPage(
  titlePanel(textOutput("plotTitle")),
  sidebarLayout(
    sidebarPanel(
      selectInput("method", "Select method(s):",
                  choices = unique(df_base$method),
                  selected = unique(df_base$method)[1],
                  multiple = TRUE),
      selectInput("x_rate", "X-axis rate:", choices = rate_choices, selected = "FSR2"),
      selectInput("y_rate", "Y-axis rate:", choices = rate_choices, selected = "TSR1"),
      actionButton("swap_xy", "Swap X and Y axes"),
      br(), br(),
      selectInput("facet_row", "Facet Rows:", choices = facet_vars, selected = "bias_level"),
      selectInput("facet_col_var1", "Facet Column Variable 1 (Bottom):",
                  choices = facet_col_vars, selected = "rep.n"),
      selectInput("facet_col_var2", "Facet Column Variable 2 (Above):",
                  choices = facet_col_vars, selected = "rep_num"),
      actionButton("make_plot", "Make the graph"),
      downloadButton("downloadPlot", "Export Graph as PNG")
    ),
    mainPanel(
      width = 8,
      plotOutput("ratePlot", height = "800px")
    )
  )
)

# === Server ===
server <- function(input, output, session) {
  observeEvent(input$swap_xy, {
    current_x <- isolate(input$x_rate)
    current_y <- isolate(input$y_rate)
    
    updateSelectInput(session, "x_rate", selected = current_y)
    updateSelectInput(session, "y_rate", selected = current_x)
  })
  
  filtered_data <- eventReactive(input$make_plot, {
    var1 <- input$facet_col_var1
    var2 <- input$facet_col_var2
    row_var <- input$facet_row
    
    df_tmp <- df_base %>%
      filter(method %in% input$method)
    
    if (var1 != "None" & var2 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0("atop(", facet_col_labels[[var1]], "==", .data[[var1]], ", ",
                                    facet_col_labels[[var2]], "==", .data[[var2]], ")")
      )
    } else if (var1 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0(facet_col_labels[[var1]], "==", .data[[var1]])
      )
    } else if (var2 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0(facet_col_labels[[var2]], "==", .data[[var2]])
      )
    } else {
      df_tmp <- df_tmp %>% mutate(facet_col_combined = "''")
    }
    
    df_tmp <- df_tmp %>% mutate(
      facet_row_combined = case_when(
        row_var == "rep.n" ~ paste0("n[rep]==", rep.n),
        row_var == "rep_num" ~ paste0("N[rep]==", rep_num),
        row_var == "orig.n" ~ paste0("n[orig]==", orig.n),
        row_var == "bias_level" ~ paste0("Bias==\"", bias_level, "\""),
        row_var == "threshold_p" ~ paste0("p[threshold]==", threshold_p),
        TRUE ~ NA_character_
      ),
      facet_col_combined = fct_inorder(facet_col_combined),
      facet_row_combined = fct_inorder(facet_row_combined)
    )
    
    group_vars <- c("method", "thresh_bin", "facet_col_combined", "facet_row_combined")
    
    df_tmp %>%
      group_by(across(all_of(group_vars))) %>%
      summarize(across(all_of(rate_choices), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  })
  
  output$plotTitle <- renderText({
    paste0("Y-axis: ", input$y_rate, "   |   X-axis: ", input$x_rate)
  })
  
  plot_event <- eventReactive(input$make_plot, {
    data <- filtered_data()
    req(nrow(data) > 0)
    
    xvar <- input$x_rate
    yvar <- input$y_rate
    
    ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
      geom_path(linewidth = 1) +
      geom_point(size = 0.3, alpha = 0.8) +
      labs(
        x = paste0(xvar, " (Lower is Better)"),
        y = paste0(yvar, " (Higher is Better)"),
        color = "Method"
      ) +
      facet_grid(
        rows = vars(facet_row_combined),
        cols = vars(facet_col_combined),
        labeller = label_parsed,
        switch = "both"
      ) +
      theme_minimal() +
      theme(
        strip.placement = "outside",
        legend.position = "bottom"
      )
  })
  
  output$ratePlot <- renderPlot({
    plot_event()
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("replication_curve_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(file, plot = plot_event(), device = "png", width = 12, height = 8, dpi = 300, bg = "white")
    }
  )
}

shinyApp(ui = ui, server = server)

##############
#This version include descriptive drop-down menu options and a description of the study
library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rlang)
library(forcats)

# === Load and Prepare Data ===
folder_path <- "./MABFanalyses/matrix-wise/rates4Viz/"
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
for (file_path in rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

rates_0.2null <- rbind(BFbMA_TFPS4Viz_0.2null, EUBF_TFPS4Viz_0.2null, FEMABF_TFPS4Viz_0.2null, iBF_TFPS4Viz_0.2null, iBF2_TFPS4Viz_0.2null)
FEMA_TFPS4Viz_0.2null <- FEMA_TFPS4Viz_0.2null %>% rename(threshold = threshold_p_r)
rates_0.2null <- rates_0.2null %>% rename(threshold = threshold_BF)
rates_0.2null <- rbind(FEMA_TFPS4Viz_0.2null, rates_0.2null)

rates_0.2null <- rates_0.2null %>%
  mutate(true_es = 0.2) %>%
  mutate(across(c(method, true_es, orig.n, `QRP level`, `PB level`, `rep number`, rep.n), as.factor)) %>%
  relocate(true_es, .after = `true effect`) %>%
  relocate(group, .before = `method`) %>%
  rename(scenario = group)

df <- rates_0.2null %>%
  select(-`true effect`, -true_es, -threshold_ES) %>%
  mutate(bias_level = `PB level`) %>%
  select(-`QRP level`, -`PB level`) %>%
  relocate(bias_level, .before = `rep number`) %>%
  mutate(
    TPR = TP / (TP + FN),
    FPR = FP / (FP + TN),
    TSR1 = TS1 / Successes1,
    FSR2 = FS2 / Successes2,
    SAR1 = TS1 / (TS1 + TF1),
    FCR2 = FS2 / (FS2 + FF2)
  ) %>%
  select(scenario, method, bias_level, rep.n, rep_num = `rep number`, orig.n, threshold, threshold_p, TPR, FPR, TSR1, FSR2, SAR1, FCR2) %>%
  filter(threshold != "Inf")

df_base <- df %>%
  mutate(
    orig.n = factor(orig.n, levels = c(20, 50 ,200)),
    rep.n = factor(rep.n, levels = c(40, 100, 400)),
    rep_num = factor(rep_num, levels = c(2, 5, 10)),
    bias_level = factor(bias_level, levels = c("low", "medium", "high")),
    method = factor(method),
    threshold_p = factor(threshold_p, levels = c(0.01, 0.05, 0.1)),
    thresh_bin = ntile(threshold, 50)
  )

facet_col_labels <- list(
  rep.n = "n[rep]",
  rep_num = "N[rep]",
  orig.n = "n[orig]",
  bias_level = "Bias",
  threshold_p = "p[threshold]"
)

rate_choices <- c("TPR", "FPR", "TSR1", "FSR2", "SAR1", "FCR2")

facet_vars <- setNames(
  c("None", "rep.n", "rep_num", "orig.n", "bias_level", "threshold_p"),
  c("None",
    "Replication Study Sample Size",
    "Number of Replications",
    "Original Study Sample Size",
    "Level of Bias",
    "Original Study Alpha Level")
)

# === UI ===
ui <- fluidPage(
  titlePanel(textOutput("plotTitle")),
  sidebarLayout(
    sidebarPanel(
      selectInput("method", "Select method(s):",
                  choices = unique(df_base$method),
                  selected = unique(df_base$method)[1],
                  multiple = TRUE),
      
      selectInput("x_rate", "X-axis rate:", choices = rate_choices, selected = "FSR2"),
      tags$small(textOutput("xrate_desc"), style = "color: gray; font-style: italic; margin-bottom: 10px;"),
      
      selectInput("y_rate", "Y-axis rate:", choices = rate_choices, selected = "TSR1"),
      tags$small(textOutput("yrate_desc"), style = "color: gray; font-style: italic; margin-bottom: 10px;"),
      
      actionButton("swap_xy", "Swap X and Y axes"),
      br(), br(),
      
      selectInput("facet_row", "Facet Rows:", choices = facet_vars, selected = "bias_level"),
      selectInput("facet_3rd", "Third Facet (Vertical Panels):", choices = facet_vars, selected = "None"),
      selectInput("facet_col_var1", "Facet Column Variable 1 (Bottom):",
                  choices = facet_vars, selected = "rep.n"),
      selectInput("facet_col_var2", "Facet Column Variable 2 (Above):",
                  choices = facet_vars, selected = "rep_num"),
      
      actionButton("make_plot", "Make the graph"),
      downloadButton("downloadPlot", "Export Graph as PNG"),
      
      tags$div(
        style = "margin-top: 30px; padding: 15px; background-color: #f7f7f7; border: 1px solid #ddd; border-radius: 5px;",
        HTML("<strong>About this study:</strong><br/>
        This app is part of a dissertation project that simulates a large-scale dataset to compare six Bayesian statistical methods (FEMA, FEMABF, BFbMA, EUBF, iBF, iBF2) for combining replication studies.<br/><br/>
        The goal is to evaluate how well these methods determine whether original findings are credible.<br/><br/>
        The app visualizes performance metrics (such as TPR, FSR2, etc.) under varying conditions:
        <ul style='margin-bottom: 0;'>
        <li>Original and replication sample sizes</li>
        <li>Number of replications</li>
        <li>Different levels of publication bias</li>
        <li>Threshold values used for decision making</li>
        </ul>
        It helps explore how often methods succeed or fail to support true or false hypotheses across these conditions.")
      )
    ),
    mainPanel(
      width = 8,
      plotOutput("ratePlot", height = "800px")
    )
  )
)

# === Server ===
server <- function(input, output, session) {
  rate_descriptions <- c(
    TPR = "True Positive Rate: proportion of true effects correctly identified.",
    FPR = "False Positive Rate: proportion of null effects wrongly identified as true.",
    TSR1 = "True Success Rate: percent of replications that confirm the original when the effect exists.",
    FSR2 = "False Success Rate: percent of replications that confirm the original when the effect does not exist.",
    SAR1 = "Successful Affirmation Rate: replication confirms the direction of the original when it's correct.",
    FCR2 = "Failed Correction Rate: replication fails to correct the direction of a wrong original finding."
  )
  
  output$xrate_desc <- renderText({
    rate_descriptions[[input$x_rate]]
  })
  
  output$yrate_desc <- renderText({
    rate_descriptions[[input$y_rate]]
  })
  
  observeEvent(input$swap_xy, {
    current_x <- isolate(input$x_rate)
    current_y <- isolate(input$y_rate)
    
    updateSelectInput(session, "x_rate", selected = current_y)
    updateSelectInput(session, "y_rate", selected = current_x)
  })
  
  filtered_data <- eventReactive(input$make_plot, {
    var1 <- input$facet_col_var1
    var2 <- input$facet_col_var2
    row_var <- input$facet_row
    third_var <- input$facet_3rd
    
    df_tmp <- df_base %>%
      filter(method %in% input$method)
    
    if (var1 != "None" & var2 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0("atop(", facet_col_labels[[var1]], "==", .data[[var1]], ", ",
                                    facet_col_labels[[var2]], "==", .data[[var2]], ")")
      )
    } else if (var1 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0(facet_col_labels[[var1]], "==", .data[[var1]])
      )
    } else if (var2 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0(facet_col_labels[[var2]], "==", .data[[var2]])
      )
    } else {
      df_tmp <- df_tmp %>% mutate(facet_col_combined = "''")
    }
    
    df_tmp <- df_tmp %>% mutate(
      facet_row_combined = case_when(
        row_var == "rep.n" ~ paste0("n[rep]==", rep.n),
        row_var == "rep_num" ~ paste0("N[rep]==", rep_num),
        row_var == "orig.n" ~ paste0("n[orig]==", orig.n),
        row_var == "bias_level" ~ paste0("Bias==\"", bias_level, "\""),
        row_var == "threshold_p" ~ paste0("p[threshold]==", threshold_p),
        TRUE ~ NA_character_
      ),
      facet_col_combined = fct_inorder(facet_col_combined),
      facet_row_combined = fct_inorder(facet_row_combined)
    )
    
    group_vars <- c("method", "thresh_bin", "facet_col_combined", "facet_row_combined")
    if (third_var != "None") group_vars <- c(group_vars, third_var)
    
    df_tmp %>%
      group_by(across(all_of(group_vars))) %>%
      summarize(across(all_of(rate_choices), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  })
  
  output$plotTitle <- renderText({
    paste0("Y-axis: ", input$y_rate, "   |   X-axis: ", input$x_rate)
  })
  
  plot_event <- eventReactive(input$make_plot, {
    data <- filtered_data()
    req(nrow(data) > 0)
    
    xvar <- input$x_rate
    yvar <- input$y_rate
    third_facet <- input$facet_3rd
    
    if (third_facet == "None" || third_facet == input$facet_row) {
      ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
        geom_path(linewidth = 1) +
        geom_point(size = 0.3, alpha = 0.8) +
        labs(
          x = paste0(xvar, " (Lower is Better)"),
          y = paste0(yvar, " (Higher is Better)"),
          color = "Method"
        ) +
        facet_grid(
          rows = vars(facet_row_combined),
          cols = vars(facet_col_combined),
          labeller = label_parsed,
          switch = "both"
        ) +
        theme_minimal() +
        theme(
          strip.placement = "outside",
          legend.position = "bottom"
        )
    } else {
      third_levels <- unique(data[[third_facet]])
      plots <- lapply(third_levels, function(level_value) {
        df_subset <- data[data[[third_facet]] == level_value, , drop = FALSE]
        if (nrow(df_subset) == 0) return(NULL)
        
        ggplot(df_subset, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
          geom_path(linewidth = 1) +
          geom_point(size = 0.3, alpha = 0.8) +
          labs(
            x = paste0(xvar, " (Lower is Better)"),
            y = paste0(yvar, " (Higher is Better)"),
            title = paste(third_facet, "=", level_value),
            color = "Method"
          ) +
          facet_grid(
            rows = vars(facet_row_combined),
            cols = vars(facet_col_combined),
            labeller = label_parsed,
            switch = "both"
          ) +
          theme_minimal() +
          theme(strip.placement = "outside")
      })
      
      plots <- Filter(Negate(is.null), plots)
      wrap_plots(plots, ncol = 1) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
    }
  })
  
  output$ratePlot <- renderPlot({
    plot_event()
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("replication_curve_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(file, plot = plot_event(), device = "png", width = 12, height = 8, dpi = 300, bg = "white")
    }
  )
}

shinyApp(ui = ui, server = server)

#######
#This version is the published Shiny version with two tabs and descriptive performance metrics
library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rlang)
library(forcats)

# === Load and Prepare Data ===
folder_path <- "./MABFanalyses/matrix-wise/rates4Viz/"
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
for (file_path in rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

rates_0.2null <- rbind(BFbMA_TFPS4Viz_0.2null, EUBF_TFPS4Viz_0.2null, FEMABF_TFPS4Viz_0.2null, iBF_TFPS4Viz_0.2null)
FEMA_TFPS4Viz_0.2null <- FEMA_TFPS4Viz_0.2null %>% rename(threshold = threshold_p_r)
rates_0.2null <- rates_0.2null %>% rename(threshold = threshold_BF)
rates_0.2null <- rbind(FEMA_TFPS4Viz_0.2null, rates_0.2null)

rates_0.2null <- rates_0.2null %>%
  mutate(true_es = 0.2) %>%
  mutate(across(c(method, true_es, orig.n, `QRP level`, `PB level`, `rep number`, rep.n), as.factor)) %>%
  relocate(true_es, .after = `true effect`) %>%
  relocate(group, .before = `method`) %>%
  rename(scenario = group)

df <- rates_0.2null %>%
  select(-`true effect`, -true_es, -threshold_ES) %>%
  mutate(bias_level = `PB level`) %>%
  select(-`QRP level`, -`PB level`) %>%
  relocate(bias_level, .before = `rep number`) %>%
  mutate(
    TPR = TP / (TP + FN),
    FPR = FP / (FP + TN),
    TSR1 = TS1 / Successes1,
    
    FSR2 = FS2 / Successes2,
    
    SAR1 = TS1 / (TS1 + TF1),
    
    FCR2 = FS2 / (FS2 + FF2)
    
  ) %>%
  select(scenario, method, bias_level, rep.n, rep_num = `rep number`, orig.n, threshold, threshold_p, TPR, FPR, TSR1, FSR2, SAR1, FCR2) %>%
  filter(threshold != "Inf")

df_base <- df %>%
  mutate(
    orig.n = factor(orig.n, levels = c(20, 50 ,200)),
    rep.n = factor(rep.n, levels = c(40, 100, 400)),
    rep_num = factor(rep_num, levels = c(2, 5, 10)),
    bias_level = factor(bias_level, levels = c("low", "medium", "high")),
    method = factor(method),
    threshold_p = factor(threshold_p, levels = c(0.01, 0.05, 0.1)),
    thresh_bin = ntile(threshold, 50)
  )

# Label mapping
facet_col_vars <- c("None", "rep.n", "rep_num", "orig.n", "bias_level", "threshold_p")
facet_col_labels <- list(
  rep.n = "n[rep]",
  rep_num = "N[rep]",
  orig.n = "n[orig]",
  bias_level = "Bias",
  threshold_p = "p[threshold]"
)

# === UI ===
rate_choices <- c("TPR", "FPR", 
                  "TSR1", "FSR2", 
                  "SAR1", "FCR2")

facet_vars <- c("None", facet_col_vars[-1])

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .custom-sidebar {
        background-image: url('https://www.transparenttextures.com/patterns/cubes.png');
        background-color: #ffffff;
        padding: 20px;
        border-radius: 5px;
        height: 100%;
        border: 1px solid #ddd;
      }
    "))
  ),
  
  titlePanel("Replication Curve Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      class = "custom-sidebar",
      
      selectInput("method", "Select method(s):",
                  choices = unique(df_base$method),
                  selected = unique(df_base$method)[1],
                  multiple = TRUE),
      
      selectInput("x_rate", "X-axis rate:",
                  choices = c(
                    "True Positive Rate" = "TPR",
                    "False Positive Rate" = "FPR",
                    "True Success Rate (when effect exists)" = "TSR1",
                    
                    "False Success Rate (when effect is spurious)" = "FSR2",
                    
                    "Successful Affirmation Rate (when effect exists)" = "SAR1",
                    
                    "Failed Correction Rate (when effect is spurious)" = "FCR2"
                  ),
                  selected = "FSR2"),
      
      selectInput("y_rate", "Y-axis rate:",
                  choices = c(
                    "True Positive Rate" = "TPR",
                    "False Positive Rate" = "FPR",
                    "True Success Rate (when effect exists)" = "TSR1",
                    
                    "False Success Rate (when effect is spurious)" = "FSR2",
                    
                    "Successful Affirmation Rate (when effect exists)" = "SAR1",
                    
                    "Failed Correction Rate (when effect is spurious)" = "FCR2"
                  ),
                  selected = "TSR1"),
      
      actionButton("swap_xy", "Swap X and Y axes"),
      br(), br(),
      
      selectInput("facet_row", "Covariate 1:",
                  choices = c(
                    "None" = "None",
                    "Replication Sample Size" = "rep.n",
                    "Number of Replications" = "rep_num",
                    "Original Sample Size" = "orig.n",
                    "Bias Level (p-hacking and publication bias)" = "bias_level",
                    "Original Alpha Level" = "threshold_p"
                  ),
                  selected = "bias_level"),
      
      selectInput("facet_3rd", "Covariate 2:",
                  choices = c(
                    "None" = "None",
                    "Replication Sample Size" = "rep.n",
                    "Number of Replications" = "rep_num",
                    "Original Sample Size" = "orig.n",
                    "Bias Level (p-hacking and publication bias)" = "bias_level",
                    "Original Alpha Level" = "threshold_p"
                  ),
                  selected = "None"),
      
      selectInput("facet_col_var1", "Explanatory Variable 1:",
                  choices = c(
                    "None" = "None",
                    "Replication Sample Size" = "rep.n",
                    "Number of Replications" = "rep_num",
                    "Original Sample Size" = "orig.n",
                    "Bias Level (p-hacking and publication bias)" = "bias_level",
                    "Original Alpha Level" = "threshold_p"
                  ),
                  selected = "rep.n"),
      
      selectInput("facet_col_var2", "Explanatory Variable 2:",
                  choices = c(
                    "None" = "None",
                    "Replication Sample Size" = "rep.n",
                    "Number of Replications" = "rep_num",
                    "Original Sample Size" = "orig.n",
                    "Bias Level (p-hacking and publication bias)" = "bias_level",
                    "Original Alpha Level" = "threshold_p"
                  ),
                  selected = "rep_num"),
      
      actionButton("make_plot", "Make the graph"),
      downloadButton("downloadPlot", "Export Graph as PNG")
    ),
    
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Study Overview",
                 tags$div(style = "padding: 30px; max-width: 900px; text-align: left;",
                          HTML("
              <p>Recent concerns about the replication crisis in psychology have underscored the need for better frameworks to assess whether original findings can be reliably replicated. This study aims to compare the performance of several meta-analytic Bayes factors (FEMABF, BFbMA, EUBF, and iBF) both relative to one another and against the traditional fixed-effect meta-analysis in evaluating replication success.</p>
              <p>To do so, we use a large-scale simulation to generate original studies across varied scenarios by manipulating true effect size, research environment (levels of p-hacking and publication bias), and original sample size. Replication studies are then simulated by varying replication sample size and number of replications, yielding 243 distinct scenarios and over 60 million data points.</p> 
              <p>This Shiny app visualizes MABF model performance using ROC-alike curves (called replication curve) to compare how well each method distinguishes between true and null effects across these conditions.</p>
            ")
                 )
        ),
        
        tabPanel("Visualization",
                 plotOutput("ratePlot", height = "800px")
        )
      )
    )
  )
)




# === Server ===
server <- function(input, output, session) {
  observeEvent(input$swap_xy, {
    current_x <- isolate(input$x_rate)
    current_y <- isolate(input$y_rate)
    
    updateSelectInput(session, "x_rate", selected = current_y)
    updateSelectInput(session, "y_rate", selected = current_x)
  })
  
  filtered_data <- eventReactive(input$make_plot, {
    var1 <- input$facet_col_var1
    var2 <- input$facet_col_var2
    row_var <- input$facet_row
    third_var <- input$facet_3rd
    
    df_tmp <- df_base %>%
      filter(method %in% input$method)
    
    if (var1 != "None" & var2 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0("atop(", facet_col_labels[[var1]], "==", .data[[var1]], ", ",
                                    facet_col_labels[[var2]], "==", .data[[var2]], ")")
      )
    } else if (var1 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0(facet_col_labels[[var1]], "==", .data[[var1]])
      )
    } else if (var2 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0(facet_col_labels[[var2]], "==", .data[[var2]])
      )
    } else {
      df_tmp <- df_tmp %>% mutate(facet_col_combined = "''")
    }
    
    df_tmp <- df_tmp %>% mutate(
      facet_row_combined = case_when(
        row_var == "rep.n" ~ paste0("n[rep]==", rep.n),
        row_var == "rep_num" ~ paste0("N[rep]==", rep_num),
        row_var == "orig.n" ~ paste0("n[orig]==", orig.n),
        row_var == "bias_level" ~ paste0("Bias==\"", bias_level, "\""),
        row_var == "threshold_p" ~ paste0("p[threshold]==", threshold_p),
        TRUE ~ NA_character_
      ),
      facet_col_combined = fct_inorder(facet_col_combined),
      facet_row_combined = fct_inorder(facet_row_combined)
    )
    
    group_vars <- c("method", "thresh_bin", "facet_col_combined", "facet_row_combined")
    if (third_var != "None") group_vars <- c(group_vars, third_var)
    
    df_tmp %>%
      group_by(across(all_of(group_vars))) %>%
      summarize(across(all_of(rate_choices), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  })
  
  output$plotTitle <- renderText({
    "Replication Curve Analysis"
  })
  
  plot_event <- eventReactive(input$make_plot, {
    data <- filtered_data()
    req(nrow(data) > 0)
    
    xvar <- input$x_rate
    yvar <- input$y_rate
    third_facet <- input$facet_3rd
    
    if (third_facet == "None" || third_facet == input$facet_row) {
      ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
        geom_path(linewidth = 1) +
        geom_point(size = 0.3, alpha = 0.8) +
        labs(
          x = paste0(xvar, " (Lower is Better)"),
          y = paste0(yvar, " (Higher is Better)"),
          color = "Method"
        ) +
        facet_grid(
          rows = vars(facet_row_combined),
          cols = vars(facet_col_combined),
          labeller = label_parsed,
          switch = "both"
        ) +
        theme_minimal() +
        theme(
          strip.placement = "outside",
          legend.position = "bottom"
        )
    } else {
      third_levels <- unique(data[[third_facet]])
      plots <- lapply(third_levels, function(level_value) {
        df_subset <- data[data[[third_facet]] == level_value, , drop = FALSE]
        if (nrow(df_subset) == 0) return(NULL)
        
        ggplot(df_subset, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
          geom_path(linewidth = 1) +
          geom_point(size = 0.3, alpha = 0.8) +
          labs(
            x = paste0(xvar, " (Lower is Better)"),
            y = paste0(yvar, " (Higher is Better)"),
            title = paste(third_facet, "=", level_value),
            color = "Method"
          ) +
          facet_grid(
            rows = vars(facet_row_combined),
            cols = vars(facet_col_combined),
            labeller = label_parsed,
            switch = "both"
          ) +
          theme_minimal() +
          theme(strip.placement = "outside")
      })
      
      plots <- Filter(Negate(is.null), plots)
      wrap_plots(plots, ncol = 1) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
    }
  })
  
  output$ratePlot <- renderPlot({
    plot_event()
  })
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("replication_curve_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(file, plot = plot_event(), device = "png", width = 12, height = 8, dpi = 300, bg = "white")
    }
  )
}

shinyApp(ui = ui, server = server)

#######
#######This is the dual dataset version (you can select between true effect = 0.2 or 0.5)
library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rlang)
library(forcats)

# === Load and Prepare Data Function ===
prepare_dataset <- function(rates_data, FEMA_data, true_es_value) {
  FEMA_data <- FEMA_data %>% rename(threshold = threshold_p_r)
  rates_data <- rates_data %>% rename(threshold = threshold_BF)
  rates <- rbind(FEMA_data, rates_data)
  
  rates <- rates %>%
    mutate(true_es = true_es_value) %>%
    mutate(across(c(method, true_es, orig.n, `QRP level`, `PB level`, `rep number`, rep.n), as.factor)) %>%
    relocate(true_es, .after = `true effect`) %>%
    relocate(group, .before = `method`) %>%
    rename(scenario = group)
  
  df <- rates %>%
    select(-`true effect`, -true_es, -threshold_ES) %>%
    mutate(bias_level = `PB level`) %>%
    select(-`QRP level`, -`PB level`) %>%
    relocate(bias_level, .before = `rep number`) %>%
    mutate(
      TPR = TP / (TP + FN),
      FPR = FP / (FP + TN),
      TSR1 = TS1 / Successes1,
      FSR2 = FS2 / Successes2,
      SAR1 = TS1 / (TS1 + TF1),
      FCR2 = FS2 / (FS2 + FF2)
    ) %>%
    select(scenario, method, bias_level, rep.n, rep_num = `rep number`, orig.n, threshold, threshold_p, 
           TPR, FPR, TSR1, FSR2, SAR1, FCR2) %>%
    filter(threshold != "Inf")
  
  df %>%
    mutate(
      orig.n = factor(orig.n, levels = c(20, 50 ,200)),
      rep.n = factor(rep.n, levels = c(40, 100, 400)),
      rep_num = factor(rep_num, levels = c(2, 5, 10)),
      bias_level = factor(bias_level, levels = c("low", "medium", "high")),
      method = factor(method),
      threshold_p = factor(threshold_p, levels = c(0.01, 0.05, 0.1)),
      thresh_bin = ntile(threshold, 50)
    )
}

# === Load Data ===
folder_path <- "./MABFanalyses/matrix-wise/rates4Viz/"
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
for (file_path in rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

df_0.2null <- prepare_dataset(
  rbind(BFbMA_TFPS4Viz_0.2null, EUBF_TFPS4Viz_0.2null, FEMABF_TFPS4Viz_0.2null, iBF_TFPS4Viz_0.2null),
  FEMA_TFPS4Viz_0.2null,
  0.2
)

df_0.5null <- prepare_dataset(
  rbind(BFbMA_TFPS4Viz_0.5null, EUBF_TFPS4Viz_0.5null, FEMABF_TFPS4Viz_0.5null, iBF_TFPS4Viz_0.5null),
  FEMA_TFPS4Viz_0.5null,
  0.5
)

# Label mapping
facet_col_vars <- c("None", "rep.n", "rep_num", "orig.n", "bias_level", "threshold_p")
facet_col_labels <- list(
  rep.n = "n[rep]",
  rep_num = "N[rep]",
  orig.n = "n[orig]",
  bias_level = "Bias",
  threshold_p = "p[threshold]"
)

# === UI ===
rate_choices <- c("TPR", "FPR", "TSR1", "FSR2", "SAR1", "FCR2")

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .custom-sidebar {
        background-image: url('https://www.transparenttextures.com/patterns/cubes.png');
        background-color: #ffffff;
        padding: 20px;
        border-radius: 5px;
        height: 100%;
        border: 1px solid #ddd;
      }
    "))
  ),
  
  titlePanel("Replication Curve Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      class = "custom-sidebar",
      
      # --- Updated dropdown menu ---
      selectInput("dataset", "Select true effect:", 
                  choices = c("0.2" = "0.2null", "0.5" = "0.5null"), 
                  selected = "0.2null"),
      
      selectInput("method", "Select method(s):", choices = NULL, multiple = TRUE),
      
      selectInput("x_rate", "X-axis rate:",
                  choices = c(
                    "True Positive Rate" = "TPR",
                    "False Positive Rate" = "FPR",
                    "True Success Rate (when effect exists)" = "TSR1",
                    "False Success Rate (when effect is spurious)" = "FSR2",
                    "Successful Affirmation Rate (when effect exists)" = "SAR1",
                    "Failed Correction Rate (when effect is spurious)" = "FCR2"
                  ),
                  selected = "FSR2"),
      
      selectInput("y_rate", "Y-axis rate:",
                  choices = c(
                    "True Positive Rate" = "TPR",
                    "False Positive Rate" = "FPR",
                    "True Success Rate (when effect exists)" = "TSR1",
                    "False Success Rate (when effect is spurious)" = "FSR2",
                    "Successful Affirmation Rate (when effect exists)" = "SAR1",
                    "Failed Correction Rate (when effect is spurious)" = "FCR2"
                  ),
                  selected = "TSR1"),
      
      actionButton("swap_xy", "Swap X and Y axes"),
      br(), br(),
      
      selectInput("facet_row", "Covariate 1:",
                  choices = c("None" = "None",
                              "Replication Sample Size" = "rep.n",
                              "Number of Replications" = "rep_num",
                              "Original Sample Size" = "orig.n",
                              "Bias Level (p-hacking and publication bias)" = "bias_level",
                              "Original Alpha Level" = "threshold_p"),
                  selected = "bias_level"),
      
      selectInput("facet_3rd", "Covariate 2:",
                  choices = c("None" = "None",
                              "Replication Sample Size" = "rep.n",
                              "Number of Replications" = "rep_num",
                              "Original Sample Size" = "orig.n",
                              "Bias Level (p-hacking and publication bias)" = "bias_level",
                              "Original Alpha Level" = "threshold_p"),
                  selected = "None"),
      
      selectInput("facet_col_var1", "Explanatory Variable 1:",
                  choices = c("None" = "None",
                              "Replication Sample Size" = "rep.n",
                              "Number of Replications" = "rep_num",
                              "Original Sample Size" = "orig.n",
                              "Bias Level (p-hacking and publication bias)" = "bias_level",
                              "Original Alpha Level" = "threshold_p"),
                  selected = "rep.n"),
      
      selectInput("facet_col_var2", "Explanatory Variable 2:",
                  choices = c("None" = "None",
                              "Replication Sample Size" = "rep.n",
                              "Number of Replications" = "rep_num",
                              "Original Sample Size" = "orig.n",
                              "Bias Level (p-hacking and publication bias)" = "bias_level",
                              "Original Alpha Level" = "threshold_p"),
                  selected = "rep_num"),
      
      actionButton("make_plot", "Make the graph"),
      downloadButton("downloadPlot", "Export Graph as PNG")
    ),
    
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Study Overview",
                 tags$div(style = "padding: 30px; max-width: 900px; text-align: left;",
                          HTML("
              <p>Recent concerns about the replication crisis in psychology have underscored the need for better frameworks to assess whether original findings can be reliably replicated. This study aims to compare the performance of several meta-analytic Bayes factors (FEMABF, BFbMA, EUBF, and iBF) both relative to one another and against the traditional fixed-effect meta-analysis in evaluating replication success.</p>
              <p>To do so, we use a large-scale simulation to generate original studies across varied scenarios by manipulating true effect size, research environment (levels of p-hacking and publication bias), and original sample size. Replication studies are then simulated by varying replication sample size and number of replications, yielding 243 distinct scenarios and over 60 million data points.</p> 
              <p>This Shiny app visualizes MABF model performance using ROC-alike curves (called replication curve) to compare how well each method distinguishes between true and null effects across these conditions.</p>
            "))
        ),
        
        tabPanel("Visualization",
                 plotOutput("ratePlot", height = "800px")
        )
      )
    )
  )
)


# === Server ===
server <- function(input, output, session) {
  
  # update method choices based on dataset
  observeEvent(input$dataset, {
    df_current <- if (input$dataset == "0.2null") df_0.2null else df_0.5null
    updateSelectInput(session, "method",
                      choices = unique(df_current$method),
                      selected = unique(df_current$method)[1])
  })
  
  observeEvent(input$swap_xy, {
    current_x <- isolate(input$x_rate)
    current_y <- isolate(input$y_rate)
    updateSelectInput(session, "x_rate", selected = current_y)
    updateSelectInput(session, "y_rate", selected = current_x)
  })
  
  filtered_data <- eventReactive(input$make_plot, {
    df_base <- if (input$dataset == "0.2null") df_0.2null else df_0.5null
    
    var1 <- input$facet_col_var1
    var2 <- input$facet_col_var2
    row_var <- input$facet_row
    third_var <- input$facet_3rd
    
    df_tmp <- df_base %>% filter(method %in% input$method)
    
    if (var1 != "None" & var2 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0("atop(", facet_col_labels[[var1]], "==", .data[[var1]], ", ",
                                    facet_col_labels[[var2]], "==", .data[[var2]], ")")
      )
    } else if (var1 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0(facet_col_labels[[var1]], "==", .data[[var1]])
      )
    } else if (var2 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0(facet_col_labels[[var2]], "==", .data[[var2]])
      )
    } else {
      df_tmp <- df_tmp %>% mutate(facet_col_combined = "''")
    }
    
    df_tmp <- df_tmp %>% mutate(
      facet_row_combined = case_when(
        row_var == "rep.n" ~ paste0("n[rep]==", rep.n),
        row_var == "rep_num" ~ paste0("N[rep]==", rep_num),
        row_var == "orig.n" ~ paste0("n[orig]==", orig.n),
        row_var == "bias_level" ~ paste0("Bias==\"", bias_level, "\""),
        row_var == "threshold_p" ~ paste0("p[threshold]==", threshold_p),
        TRUE ~ NA_character_
      ),
      facet_col_combined = fct_inorder(facet_col_combined),
      facet_row_combined = fct_inorder(facet_row_combined)
    )
    
    group_vars <- c("method", "thresh_bin", "facet_col_combined", "facet_row_combined")
    if (third_var != "None") group_vars <- c(group_vars, third_var)
    
    df_tmp %>%
      group_by(across(all_of(group_vars))) %>%
      summarize(across(all_of(rate_choices), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  })
  
  plot_event <- eventReactive(input$make_plot, {
    data <- filtered_data()
    req(nrow(data) > 0)
    
    xvar <- input$x_rate
    yvar <- input$y_rate
    third_facet <- input$facet_3rd
    
    if (third_facet == "None" || third_facet == input$facet_row) {
      ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
        geom_path(linewidth = 1) +
        geom_point(size = 0.3, alpha = 0.8) +
        labs(x = paste0(xvar, " (Lower is Better)"),
             y = paste0(yvar, " (Higher is Better)"),
             color = "Method") +
        facet_grid(rows = vars(facet_row_combined),
                   cols = vars(facet_col_combined),
                   labeller = label_parsed,
                   switch = "both") +
        theme_minimal() +
        theme(strip.placement = "outside", legend.position = "bottom")
    } else {
      third_levels <- unique(data[[third_facet]])
      plots <- lapply(third_levels, function(level_value) {
        df_subset <- data[data[[third_facet]] == level_value, , drop = FALSE]
        if (nrow(df_subset) == 0) return(NULL)
        
        ggplot(df_subset, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
          geom_path(linewidth = 1) +
          geom_point(size = 0.3, alpha = 0.8) +
          labs(x = paste0(xvar, " (Lower is Better)"),
               y = paste0(yvar, " (Higher is Better)"),
               title = paste(third_facet, "=", level_value),
               color = "Method") +
          facet_grid(rows = vars(facet_row_combined),
                     cols = vars(facet_col_combined),
                     labeller = label_parsed,
                     switch = "both") +
          theme_minimal() +
          theme(strip.placement = "outside")
      })
      
      plots <- Filter(Negate(is.null), plots)
      wrap_plots(plots, ncol = 1) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
    }
  })
  
  output$ratePlot <- renderPlot({ plot_event() })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste0("replication_curve_", Sys.Date(), ".png") },
    content = function(file) {
      ggsave(file, plot = plot_event(), device = "png", width = 12, height = 8, dpi = 300, bg = "white")
    }
  )
}

shinyApp(ui = ui, server = server)



#######Use this version for each method individually. This version superimpose different levels of Nrep in one graph for a method
library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rlang)
library(forcats)

# === Load and Prepare Data Function ===
prepare_dataset <- function(rates_data, FEMA_data, true_es_value) {
  FEMA_data <- FEMA_data %>% rename(threshold = threshold_p_r)
  rates_data <- rates_data %>% rename(threshold = threshold_BF)
  rates <- rbind(FEMA_data, rates_data)
  
  rates <- rates %>%
    mutate(true_es = true_es_value) %>%
    mutate(across(c(method, true_es, orig.n, `QRP level`, `PB level`, `rep number`, rep.n), as.factor)) %>%
    relocate(true_es, .after = `true effect`) %>%
    relocate(group, .before = `method`) %>%
    rename(scenario = group)
  
  df <- rates %>%
    select(-`true effect`, -true_es, -threshold_ES) %>%
    mutate(bias_level = `PB level`) %>%
    select(-`QRP level`, -`PB level`) %>%
    relocate(bias_level, .before = `rep number`) %>%
    mutate(
      TPR = TP / (TP + FN),
      FPR = FP / (FP + TN),
      TSR1 = TS1 / Successes1,
      FSR2 = FS2 / Successes2,
      FFR1 = FF1 / Failures1,
      TFR2 = TF2 / Failures2
      
    ) %>%
    select(scenario, method, bias_level, rep.n, rep_num = `rep number`, orig.n, threshold, threshold_p, 
           TPR, FPR, TSR1, FSR2, FFR1, TFR2) %>%
    filter(threshold != "Inf")
  
  df %>%
    mutate(
      orig.n = factor(orig.n, levels = c(20, 50 ,200)),
      rep.n = factor(rep.n, levels = c(40, 100, 400)),
      rep_num = factor(rep_num, levels = c(2, 5, 10)),
      bias_level = factor(bias_level, levels = c("low", "medium", "high")),
      method = factor(method),
      threshold_p = factor(threshold_p, levels = c(0.1, 0.05, 0.01)),
      thresh_bin = ntile(threshold, 5)
    )
}

# === Load Data ===
folder_path <- "./MABFanalyses/matrix-wise/rates4Viz/"
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
for (file_path in rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

df_0.2null <- prepare_dataset(
  rbind(BFbMA_TFPS4Viz_0.2null, EUBF_TFPS4Viz_0.2null, FEMABF_TFPS4Viz_0.2null, iBF_TFPS4Viz_0.2null),
  FEMA_TFPS4Viz_0.2null,
  0.2
)

df_0.5null <- prepare_dataset(
  rbind(BFbMA_TFPS4Viz_0.5null, EUBF_TFPS4Viz_0.5null, FEMABF_TFPS4Viz_0.5null, iBF_TFPS4Viz_0.5null),
  FEMA_TFPS4Viz_0.5null,
  0.5
)

# Label mapping
facet_col_vars <- c("None", "rep.n", "rep_num", "orig.n", "bias_level", "threshold_p")
facet_col_labels <- list(
  rep.n = "n[rep]",
  rep_num = "N[rep]",
  orig.n = "n[orig]",
  bias_level = "Bias",
  threshold_p = "p[threshold]"
)

# === UI ===
rate_choices <- c("TPR", "FPR", "TSR1", "FSR2", "FFR1", "TFR2")

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .custom-sidebar {
        background-image: url('https://www.transparenttextures.com/patterns/cubes.png');
        background-color: #ffffff;
        padding: 20px;
        border-radius: 5px;
        height: 100%;
        border: 1px solid #ddd;
      }
    "))
  ),
  
  titlePanel("Replication Curve Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      class = "custom-sidebar",
      
      # --- Updated dropdown menu ---
      selectInput("dataset", "Select true effect:", 
                  choices = c("0.2" = "0.2null", "0.5" = "0.5null"), 
                  selected = "0.2null"),
      
      selectInput("method", "Select method(s):", choices = NULL, multiple = TRUE),
      
      selectInput("x_rate", "X-axis rate:",
                  choices = c(
                    "True Positive Rate" = "TPR",
                    "False Positive Rate" = "FPR",
                    "True Success Rate (when effect exists)" = "TSR1",
                    "False Success Rate (when effect is spurious)" = "FSR2",
                    "False Failure Rate (When effect exists)" = "FFR1",
                    "True Failure Rate (When effect is spurisous)" = "TFR2"
                  ),
                  selected = "FSR2"),
      
      selectInput("y_rate", "Y-axis rate:",
                  choices = c(
                    "True Positive Rate" = "TPR",
                    "False Positive Rate" = "FPR",
                    "True Success Rate (when effect exists)" = "TSR1",
                    "False Success Rate (when effect is spurious)" = "FSR2",
                    "False Failure Rate (When effect exists)" = "FFR1",
                    "True Failure Rate (When effect is spurisous)" = "TFR2"
                  ),
                  selected = "TSR1"),
      
      actionButton("swap_xy", "Swap X and Y axes"),
      br(), br(),
      
      selectInput("facet_row", "Covariate 1:",
                  choices = c("None" = "None",
                              "Replication Sample Size" = "rep.n",
                              "Number of Replications" = "rep_num",
                              "Original Sample Size" = "orig.n",
                              "Bias Level (p-hacking and publication bias)" = "bias_level",
                              "Original Alpha Level" = "threshold_p"),
                  selected = "bias_level"),
      
      selectInput("facet_3rd", "Covariate 2:",
                  choices = c("None" = "None",
                              "Replication Sample Size" = "rep.n",
                              "Number of Replications" = "rep_num",
                              "Original Sample Size" = "orig.n",
                              "Bias Level (p-hacking and publication bias)" = "bias_level",
                              "Original Alpha Level" = "threshold_p"),
                  selected = "None"),
      
      selectInput("facet_col_var1", "Explanatory Variable 1:",
                  choices = c("None" = "None",
                              "Replication Sample Size" = "rep.n",
                              "Number of Replications" = "rep_num",
                              "Original Sample Size" = "orig.n",
                              "Bias Level (p-hacking and publication bias)" = "bias_level",
                              "Original Alpha Level" = "threshold_p"),
                  selected = "rep.n"),
      
      selectInput("facet_col_var2", "Explanatory Variable 2:",
                  choices = c("None" = "None",
                              "Replication Sample Size" = "rep.n",
                              "Number of Replications" = "rep_num",
                              "Original Sample Size" = "orig.n",
                              "Bias Level (p-hacking and publication bias)" = "bias_level",
                              "Original Alpha Level" = "threshold_p"),
                  selected = "rep_num"),
      
      actionButton("make_plot", "Make the graph"),
      downloadButton("downloadPlot", "Export Graph as PNG")
    ),
    
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Study Overview",
                 tags$div(style = "padding: 30px; max-width: 900px; text-align: left;",
                          HTML("
              <p>Recent concerns about the replication crisis in psychology have underscored the need for better frameworks to assess whether original findings can be reliably replicated. This study aims to compare the performance of several meta-analytic Bayes factors (FEMABF, BFbMA, EUBF, and iBF) both relative to one another and against the traditional fixed-effect meta-analysis in evaluating replication success.</p>
              <p>To do so, we use a large-scale simulation to generate original studies across varied scenarios by manipulating true effect size, research environment (levels of p-hacking and publication bias), and original sample size. Replication studies are then simulated by varying replication sample size and number of replications, yielding 243 distinct scenarios and over 60 million data points.</p> 
              <p>This Shiny app visualizes MABF model performance using ROC-alike curves (called replication curve) to compare how well each method distinguishes between true and null effects across these conditions.</p>
            "))
        ),
        
        tabPanel("Visualization",
                 plotOutput("ratePlot", height = "800px")
        )
      )
    )
  )
)


# === Server ===
server <- function(input, output, session) {
  
  # update method choices based on dataset
  observeEvent(input$dataset, {
    df_current <- if (input$dataset == "0.2null") df_0.2null else df_0.5null
    updateSelectInput(session, "method",
                      choices = unique(df_current$method),
                      selected = unique(df_current$method)[1])
  })
  
  observeEvent(input$swap_xy, {
    current_x <- isolate(input$x_rate)
    current_y <- isolate(input$y_rate)
    updateSelectInput(session, "x_rate", selected = current_y)
    updateSelectInput(session, "y_rate", selected = current_x)
  })
  
  filtered_data <- eventReactive(input$make_plot, {
    df_base <- if (input$dataset == "0.2null") df_0.2null else df_0.5null
    
    var1 <- input$facet_col_var1
    var2 <- input$facet_col_var2
    row_var <- input$facet_row
    third_var <- input$facet_3rd
    
    df_tmp <- df_base %>% filter(method %in% input$method)
    
    # --- build facet labels, but exclude var2 if both selected ---
    if (var1 != "None" & var2 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0(facet_col_labels[[var1]], "==", .data[[var1]])
      )
    } else if (var1 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0(facet_col_labels[[var1]], "==", .data[[var1]])
      )
    } else if (var2 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0(facet_col_labels[[var2]], "==", .data[[var2]])
      )
    } else {
      df_tmp <- df_tmp %>% mutate(facet_col_combined = "''")
    }
    
    df_tmp <- df_tmp %>% mutate(
      facet_row_combined = case_when(
        row_var == "rep.n" ~ paste0("n[rep]==", rep.n),
        row_var == "rep_num" ~ paste0("N[rep]==", rep_num),
        row_var == "orig.n" ~ paste0("n[orig]==", orig.n),
        row_var == "bias_level" ~ paste0("Bias==\"", bias_level, "\""),
        row_var == "threshold_p" ~ paste0("p[threshold]==", threshold_p),
        TRUE ~ NA_character_
      ),
      facet_col_combined = fct_inorder(facet_col_combined),
      facet_row_combined = fct_inorder(facet_row_combined)
    )
    
    # group vars keep rep_num, rep.n, etc.
    group_vars <- c("method", "thresh_bin", "facet_col_combined", "facet_row_combined")
    if (var1 != "None") group_vars <- c(group_vars, var1)
    if (var2 != "None") group_vars <- c(group_vars, var2)
    if (row_var != "None") group_vars <- c(group_vars, row_var)
    if (third_var != "None") group_vars <- c(group_vars, third_var)
    
    df_tmp %>%
      group_by(across(all_of(group_vars))) %>%
      summarize(across(all_of(rate_choices), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  })
  
  plot_event <- eventReactive(input$make_plot, {
    data <- filtered_data()
    req(nrow(data) > 0)
    
    xvar <- input$x_rate
    yvar <- input$y_rate
    third_facet <- input$facet_3rd
    
    facet_var <- input$facet_col_var1
    group_var <- input$facet_col_var2
    if (facet_var == "None" & group_var != "None") {
      facet_var <- group_var
      group_var <- "None"
    }
    
    p <- ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]]))
    
    if (group_var != "None") {
      # use group_var as color/shape/linetype
      p <- p +
        geom_path(aes(color = .data[[group_var]], 
                      linetype = .data[[group_var]],
                      group = interaction(method, .data[[group_var]])),
                  linewidth = 0.8, alpha = 0.6) +
        geom_point(aes(color = .data[[group_var]], 
                       shape = .data[[group_var]]), 
                   size = 2, stroke = 0.8) +
        labs(
          color    = ifelse(group_var == "rep_num", expression(N[rep]), group_var),
          shape    = ifelse(group_var == "rep_num", expression(N[rep]), group_var),
          linetype = ifelse(group_var == "rep_num", expression(N[rep]), group_var)
        )
    } else {
      # default coloring by method
      p <- p +
        geom_path(aes(color = method, group = method),
                  linewidth = 0.8, alpha = 0.6) +
        geom_point(aes(color = method, shape = method), 
                   size = 2, stroke = 0.8) +
        labs(color = "Method", shape = "Method")
    }
    
    # Faceting logic
    if (facet_var != "None" || input$facet_row != "None") {
      p <- p +
        facet_grid(
          rows = if (input$facet_row != "None") vars(facet_row_combined) else NULL,
          cols = if (facet_var != "None") vars(facet_col_combined) else NULL,
          labeller = label_parsed,
          switch = "both"
        )
    }
    
    
    p + theme_minimal() +
      theme(strip.placement = "outside", legend.position = "bottom")
  })
  
  output$ratePlot <- renderPlot({ plot_event() })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste0("replication_curve_", Sys.Date(), ".png") },
    content = function(file) {
      ggsave(file, plot = plot_event(), device = "png", width = 12, height = 8, dpi = 300, bg = "white")
    }
  )
}



shinyApp(ui = ui, server = server)


#####This version uses the van Ravenzwaaij and Ioannidis (2019) ROC curve style (no 'steps' like in the ROC curves)
library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rlang)
library(forcats)

# === Load and Prepare Data Function ===
prepare_dataset <- function(rates_data, FEMA_data, true_es_value) {
  FEMA_data <- FEMA_data %>% rename(threshold = threshold_p_r)
  rates_data <- rates_data %>% rename(threshold = threshold_BF)
  rates <- rbind(FEMA_data, rates_data)
  
  rates <- rates %>%
    mutate(true_es = true_es_value) %>%
    mutate(across(c(method, true_es, orig.n, `QRP level`, `PB level`, `rep number`, rep.n), as.factor)) %>%
    relocate(true_es, .after = `true effect`) %>%
    relocate(group, .before = `method`) %>%
    rename(scenario = group)
  
  df <- rates %>%
    select(-`true effect`, -true_es, -threshold_ES) %>%
    mutate(bias_level = `PB level`) %>%
    select(-`QRP level`, -`PB level`) %>%
    relocate(bias_level, .before = `rep number`) %>%
    mutate(
      TPR = TP / (TP + FN),
      FPR = FP / (FP + TN),
      TSR1 = TS1/Successes1,
      FSR1 = FS1/Successes1,
      TFR1 = TF1/Failures1,
      FFR1 = FF1/Failures1,
      
      TSR2 = TS2/Successes2,
      FSR2 = FS2/Successes2,
      TFR2 = TF2/Failures2,
      FFR2 = FF2/Failures2,
      
      SAR1 = TS1/(TS1 + TF1), #Successful Affirmation Rate
      FAR1 = TF1/(TS1 + TF1), #Failed Affirmation Rate
      FCR1 = FS1/(FS1 + FF1), #Failed Correction Rate
      SCR1 = FF1/(FS1 + FF1), #Successful Correction Rate
      
      SAR2 = TS2/(TS2 + TF2), #Successful Affirmation Rate
      FAR2 = TF2/(TS2 + TF2), #Failed Affirmation Rate
      FCR2 = FS2/(FS2 + FF2), #Failed Correction Rate
      SCR2 = FF2/(FS2 + FF2) #Successful Correction Rate
      
    ) %>%
    select(scenario, method, bias_level, rep.n, rep_num = `rep number`, orig.n, threshold, threshold_p, 
           TPR, FPR, TSR1, FSR1, TFR1, FFR1, TSR2, FSR2, TFR2, FFR2,
           SAR1, FAR1, FCR1, SCR1, SAR2, FAR2, FCR2, SCR2) %>%
    filter(threshold != "Inf")
  
  df %>%
    mutate(
      orig.n = factor(orig.n, levels = c(20, 50 ,200)),
      rep.n = factor(rep.n, levels = c(40, 100, 400)),
      rep_num = factor(rep_num, levels = c(2, 5, 10)),
      bias_level = factor(bias_level, levels = c("low", "medium", "high")),
      method = factor(method),
      threshold_p = factor(threshold_p, levels = c(0.1, 0.05, 0.01)),
      thresh_bin = ntile(threshold, 5)
    )
}

# === Load Data ===
folder_path <- "./MABFanalyses/matrix-wise/rates4Viz/"
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
for (file_path in rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

df_0.2null <- prepare_dataset(
  rbind(BFbMA_TFPS4Viz_0.2null, EUBF_TFPS4Viz_0.2null, FEMABF_TFPS4Viz_0.2null, iBF_TFPS4Viz_0.2null),
  FEMA_TFPS4Viz_0.2null,
  0.2
)

df_0.5null <- prepare_dataset(
  rbind(BFbMA_TFPS4Viz_0.5null, EUBF_TFPS4Viz_0.5null, FEMABF_TFPS4Viz_0.5null, iBF_TFPS4Viz_0.5null),
  FEMA_TFPS4Viz_0.5null,
  0.5
)

# Label mapping
facet_col_vars <- c("None", "rep.n", "rep_num", "orig.n", "bias_level", "threshold_p")
facet_col_labels <- list(
  rep.n = "n[rep]",
  rep_num = "N[rep]",
  orig.n = "n[orig]",
  bias_level = "Bias",
  threshold_p = "p[threshold]"
)


# === UI ===
rate_choices <- c("TPR", "FPR", "TSR1", "FSR1", "TFR1", "FFR1", "TSR2", "FSR2", "TFR2", "FFR2",
                  "SAR1", "FAR1", "FCR1", "SCR1", "SAR2", "FAR2", "FCR2", "SCR2")

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .custom-sidebar {
        background-image: url('https://www.transparenttextures.com/patterns/cubes.png');
        background-color: #ffffff;
        padding: 20px;
        border-radius: 5px;
        height: 100%;
        border: 1px solid #ddd;
      }
    "))
  ),
  
  titlePanel("Replication Curve Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      class = "custom-sidebar",
      
      selectInput("dataset", "Select true effect:", 
                  choices = c("0.2" = "0.2null", "0.5" = "0.5null"), 
                  selected = "0.2null"),
      
      selectInput("method", "Select method(s):", choices = NULL, multiple = TRUE),
      
      selectInput("x_rate", "X-axis rate:",
                  choices = c("True Positive Rate" = "TPR",
                              "False Positive Rate" = "FPR",
                              "True Success Rate (when effect exists)" = "TSR1",
                              "False Success Rate (when effect exists)" = "FSR1",
                              "True Failure Rate (when effect exists)" = "TFR1", 
                              "False Failure Rate (when effect exists)" = "FFR1", 
                              "True Success Rate (When effect is spurisous)" = "TSR2", 
                              "False Success Rate (When effect is spurisous)" = "FSR2", 
                              "True Failure Rate (When effect is spurisous)" = "TFR2", 
                              "False Failure Rate (When effect is spurisous)" = "FFR2",
                              "Successful Affirmation Rate (when effect exists)" = "SAR1", 
                              "Failed Affirmation Rate (when effect exists)" = "FAR1", 
                              "Failed Correction Rate (when effect exists)" = "FCR1", 
                              "Successful Correction Rate (when effect exists)" = "SCR1", 
                              "Successful Affirmation Rate (When effect is spurisous)" = "SAR2", 
                              "Failed Affirmation Rate (When effect is spurisous)" = "FAR2", 
                              "Failed Correction Rate (When effect is spurisous)" = "FCR2", 
                              "Successful Correction Rate (When effect is spurisous)" = "SCR2"),
                  selected = "FSR2"),
      
      selectInput("y_rate", "Y-axis rate:",
                  choices = c("True Positive Rate" = "TPR",
                              "False Positive Rate" = "FPR",
                              "True Success Rate (when effect exists)" = "TSR1",
                              "False Success Rate (when effect exists)" = "FSR1",
                              "True Failure Rate (when effect exists)" = "TFR1", 
                              "False Failure Rate (when effect exists)" = "FFR1", 
                              "True Success Rate (When effect is spurisous)" = "TSR2", 
                              "False Success Rate (When effect is spurisous)" = "FSR2", 
                              "True Failure Rate (When effect is spurisous)" = "TFR2", 
                              "False Failure Rate (When effect is spurisous)" = "FFR2",
                              "Successful Affirmation Rate (when effect exists)" = "SAR1", 
                              "Failed Affirmation Rate (when effect exists)" = "FAR1", 
                              "Failed Correction Rate (when effect exists)" = "FCR1", 
                              "Successful Correction Rate (when effect exists)" = "SCR1", 
                              "Successful Affirmation Rate (When effect is spurisous)" = "SAR2", 
                              "Failed Affirmation Rate (When effect is spurisous)" = "FAR2", 
                              "Failed Correction Rate (When effect is spurisous)" = "FCR2", 
                              "Successful Correction Rate (When effect is spurisous)" = "SCR2"),
                  selected = "TSR1"),
      
      actionButton("swap_xy", "Swap X and Y axes"),
      br(), br(),
      
      selectInput("facet_row", "Covariate 1:",
                  choices = c("None" = "None",
                              "Replication Sample Size" = "rep.n",
                              "Number of Replications" = "rep_num",
                              "Original Sample Size" = "orig.n",
                              "Bias Level (p-hacking and publication bias)" = "bias_level",
                              "Original Alpha Level" = "threshold_p"),
                  selected = "bias_level"),
      
      selectInput("facet_3rd", "Covariate 2:",
                  choices = c("None" = "None",
                              "Replication Sample Size" = "rep.n",
                              "Number of Replications" = "rep_num",
                              "Original Sample Size" = "orig.n",
                              "Bias Level (p-hacking and publication bias)" = "bias_level",
                              "Original Alpha Level" = "threshold_p"),
                  selected = "None"),
      
      selectInput("facet_col_var1", "Explanatory Variable 1:",
                  choices = c("None" = "None",
                              "Replication Sample Size" = "rep.n",
                              "Number of Replications" = "rep_num",
                              "Original Sample Size" = "orig.n",
                              "Bias Level (p-hacking and publication bias)" = "bias_level",
                              "Original Alpha Level" = "threshold_p"),
                  selected = "rep.n"),
      
      selectInput("facet_col_var2", "Explanatory Variable 2:",
                  choices = c("None" = "None",
                              "Replication Sample Size" = "rep.n",
                              "Number of Replications" = "rep_num",
                              "Original Sample Size" = "orig.n",
                              "Bias Level (p-hacking and publication bias)" = "bias_level",
                              "Original Alpha Level" = "threshold_p"),
                  selected = "rep_num"),
      checkboxInput("fix_axes", "Fix axes to 0–1 for ROC line", value = FALSE),
      
      actionButton("make_plot", "Make the graph"),
      downloadButton("downloadPlot", "Export Graph as PNG")
    ),
    
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Study Overview",
                 tags$div(style = "padding: 30px; max-width: 900px; text-align: left;",
                          HTML("
              <p>Recent concerns about the replication crisis in psychology have underscored the need for better frameworks to assess whether original findings can be reliably replicated. This study aims to compare the performance of several meta-analytic Bayes factors (FEMABF, BFbMA, EUBF, and iBF) both relative to one another and against the traditional fixed-effect meta-analysis in evaluating replication success.</p>
              <p>To do so, we use a large-scale simulation to generate original studies across varied scenarios by manipulating true effect size, research environment (levels of p-hacking and publication bias), and original sample size. Replication studies are then simulated by varying replication sample size and number of replications, yielding 243 distinct scenarios and over 60 million data points.</p> 
              <p>This Shiny app visualizes MABF model performance using ROC-alike curves (called replication curve) to compare how well each method distinguishes between true and null effects across these conditions.</p>
            "))
        ),
        
        tabPanel("Visualization",
                 plotOutput("ratePlot", height = "800px", click = "plot_click"),
                 verbatimTextOutput("click_info")   # nearest symbol info
        )
        
      )
    )
  )
)

server <- function(input, output, session) {
  
  observeEvent(input$dataset, {
    df_current <- if (input$dataset == "0.2null") df_0.2null else df_0.5null
    updateSelectInput(session, "method",
                      choices = unique(df_current$method),
                      selected = unique(df_current$method)[4])
  })
  
  observeEvent(input$swap_xy, {
    current_x <- isolate(input$x_rate)
    current_y <- isolate(input$y_rate)
    updateSelectInput(session, "x_rate", selected = current_y)
    updateSelectInput(session, "y_rate", selected = current_x)
  })
  
  # --- Build the ROC curve data for all thresholds in [3, 30] ---
  filtered_data <- eventReactive(input$make_plot, {
    df_base <- if (input$dataset == "0.2null") df_0.2null else df_0.5null
    
    var1      <- input$facet_col_var1
    var2      <- input$facet_col_var2
    row_var   <- input$facet_row
    third_var <- input$facet_3rd
    
    df_tmp <- df_base %>% 
      filter(method %in% input$method)
    
    # facet labels
    if (var1 != "None" & var2 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0(facet_col_labels[[var1]], "==", .data[[var1]])
      )
    } else if (var1 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0(facet_col_labels[[var1]], "==", .data[[var1]])
      )
    } else if (var2 != "None") {
      df_tmp <- df_tmp %>% mutate(
        facet_col_combined = paste0(facet_col_labels[[var2]], "==", .data[[var2]])
      )
    } else {
      df_tmp <- df_tmp %>% mutate(facet_col_combined = "''")
    }
    
    df_tmp <- df_tmp %>%
      mutate(
        facet_row_combined = dplyr::case_when(
          row_var == "rep.n"       ~ paste0("n[rep]==", rep.n),
          row_var == "rep_num"     ~ paste0("N[rep]==", rep_num),
          row_var == "orig.n"      ~ paste0("n[orig]==", orig.n),
          row_var == "bias_level"  ~ paste0("Bias==\"", bias_level, "\""),
          row_var == "threshold_p" ~ paste0("p[threshold]==", threshold_p),
          TRUE ~ "''"   # empty string parsed, no NA facet
        ),
        facet_col_combined = forcats::fct_inorder(facet_col_combined),
        facet_row_combined = forcats::fct_inorder(facet_row_combined),
        threshold_num      = as.numeric(threshold)
      )
    
    # average across everything not in facets or the color grouping
    group_cols_curve <- c("method", "facet_col_combined", "facet_row_combined")
    if (var2 != "None") group_cols_curve <- c(group_cols_curve, var2)
    
    df_curve <- df_tmp %>%
      filter(threshold_num >= 0.01, threshold_num <= 30) %>%
      filter(!(threshold_p %in% c("0.1"))) %>%
      group_by(across(all_of(group_cols_curve)), threshold_num) %>%
      summarize(across(all_of(rate_choices), ~mean(.x, na.rm = TRUE)), .groups = "drop")
    
    df_curve
  })
  
  rate_labels <- c(
    "TPR"  = "TPR",
    "FPR"  = "FPR",
    "TSR1" = "TSR[1]",
    "FSR1" = "FSR[1]",
    "TFR1" = "TFR[1]",
    "FFR1" = "FFR[1]",
    "TSR2" = "TSR[2]",
    "FSR2" = "FSR[2]",
    "TFR2" = "TFR[2]",
    "FFR2" = "FFR[2]",
    "SAR1" = "SAR[1]",
    "FAR1" = "FAR[1]",
    "FCR1" = "FCR[1]",
    "SCR1" = "SCR[1]",
    "SAR2" = "SAR[2]",
    "FAR2" = "FAR[2]",
    "FCR2" = "FCR[2]",
    "SCR2" = "SCR[2]"
  )
  
  
  plot_event <- eventReactive(input$make_plot, {
    data <- filtered_data()
    req(nrow(data) > 0)
    
    xvar <- input$x_rate
    yvar <- input$y_rate
    
    facet_var <- input$facet_col_var1
    group_var <- input$facet_col_var2
    if (facet_var == "None" & group_var != "None") {
      facet_var <- group_var
      group_var <- "None"
    }
    
    group_cols_curve <- c("method", "facet_col_combined", "facet_row_combined")
    if (group_var != "None") group_cols_curve <- c(group_cols_curve, group_var)
    
    cutoff_vals <- c(1, 3, 10, 30)
    df_markers <- data %>%
      group_by(across(all_of(group_cols_curve))) %>%
      dplyr::group_modify(~{
        d <- .x
        rows <- lapply(cutoff_vals, function(cv) {
          d[which.min(abs(d$threshold_num - cv)), , drop = FALSE] %>%
            dplyr::mutate(cutoff = cv)
        })
        dplyr::bind_rows(rows)
      }) %>%
      ungroup()
    
    cutoff_shapes <- c("1" = 16, "3" = 17, "10" = 15, "30" = 7)
    
    # dynamic mappings with tidy eval
    if (group_var != "None") {
      col_mapping <- aes(
        x = !!sym(xvar), y = !!sym(yvar),
        color = !!sym(group_var),
        group = interaction(method, !!sym(group_var))
      )
      marker_mapping <- aes(
        x = !!sym(xvar), y = !!sym(yvar),
        color = !!sym(group_var),
        shape = factor(cutoff, levels = cutoff_vals)
      )
      color_lab <- if (group_var == "rep_num") expression(N[rep]) else group_var
    } else {
      col_mapping <- aes(
        x = !!sym(xvar), y = !!sym(yvar),
        color = method,
        group = method
      )
      marker_mapping <- aes(
        x = !!sym(xvar), y = !!sym(yvar),
        color = method,
        shape = factor(cutoff, levels = cutoff_vals)
      )
      color_lab <- "Method"
    }
    
    p <- ggplot() +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_line(data = df_markers, mapping = col_mapping, linewidth = 0.5) +
      geom_point(data = df_markers, mapping = marker_mapping, size = 2.8, stroke = 1) +
      scale_shape_manual(values = cutoff_shapes) +
      labs(
        color = color_lab,
        shape = "BF cutoff",
        x = parse(text = rate_labels[[xvar]]),
        y = parse(text = rate_labels[[yvar]])
      )
    
    
    
    # Faceting logic
    if (facet_var != "None" || input$facet_row != "None") {
      p <- p +
        facet_grid(
          rows = if (input$facet_row != "None") vars(facet_row_combined) else NULL,
          cols = if (facet_var != "None") vars(facet_col_combined) else NULL,
          labeller = label_parsed,
          switch = "both"
        )
    }
    
    
    p + theme_minimal(base_size = 12) +
      theme(strip.placement = "outside", legend.position = "bottom") 
    # Add fixed axes if checkbox selected
    if (input$fix_axes) {
      p <- p + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
    }
    
    p
  })
  # Show nearest symbol (marker) info
  # Show nearest symbol (marker) info
  output$click_info <- renderPrint({
    req(input$plot_click)
    req(filtered_data())  # use your current filtered data
    
    # regenerate marker dataset (same logic you used in plot_event)
    data <- filtered_data()
    cutoff_vals <- c(1, 3, 10, 30)
    
    facet_var <- input$facet_col_var1
    group_var <- input$facet_col_var2
    group_cols_curve <- c("method", "facet_col_combined", "facet_row_combined")
    if (group_var != "None") group_cols_curve <- c(group_cols_curve, group_var)
    
    df_markers <- data %>%
      group_by(across(all_of(group_cols_curve))) %>%
      dplyr::group_modify(~{
        d <- .x
        rows <- lapply(cutoff_vals, function(cv) {
          d[which.min(abs(d$threshold_num - cv)), , drop = FALSE] %>%
            dplyr::mutate(cutoff = cv)
        })
        dplyr::bind_rows(rows)
      }) %>%
      ungroup()
    
    # find nearest symbol to the click
    click <- input$plot_click
    df_markers$dist <- (df_markers[[input$x_rate]] - click$x)^2 + 
      (df_markers[[input$y_rate]] - click$y)^2
    nearest <- df_markers[which.min(df_markers$dist), ]
    
    # display info
    cat("Nearest symbol:\n",
        paste0(input$x_rate, " = ", round(nearest[[input$x_rate]], 4)), "\n",
        paste0(input$y_rate, " = ", round(nearest[[input$y_rate]], 4)))
  })
  
  
  
  
  
  output$ratePlot <- renderPlot({ plot_event() })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste0("replication_curve_", Sys.Date(), ".png") },
    content = function(file) {
      ggsave(file, plot = plot_event(), device = "png",
             width = 12, height = 8, dpi = 300, bg = "white")
    }
  )
}

shinyApp(ui = ui, server = server)

