#####
# library(shiny)
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(stringr)
# library(DT)
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
#####
# load("./MABFanalyses/row-wise/LinePlotData.RData")
# 
# # --- helper for summary table (and plot data) ---
# summarise_plot_data <- function(data, outcome_var, x_var, cutoff_filter, method_filter) {
#   df <- data %>% filter(Outcome == outcome_var,
#                         BFcutoff %in% cutoff_filter)
#   if (length(method_filter) > 0) {
#     df <- df %>% filter(method %in% method_filter)
#   }
#   
#   df %>%
#     group_by(method, !!sym(x_var), Outcome, BFcutoff, true.effect) %>%
#     summarise(
#       mean_rate   = round(mean(Rate, na.rm = TRUE), 2),
#       sd_rate     = round(sd(Rate, na.rm = TRUE), 2),
#       prop_missing = round(mean(is.na(Rate)), 4),  
#       .groups     = "drop"
#     )%>%
#     mutate(mean_rate = round(mean_rate, 2))
# }
# 
# # --- plotting function ---
# plot_outcome_avg <- function(df, outcome_var, x_var) {
#   cutoff_vals <- unique(df$BFcutoff)
#   eff_vals    <- unique(df$true.effect)
#   
#   # APA-style labels
#   x_labels <- c(
#     "rep.number" = "Number of Replications",
#     "rep.n"      = "Replication Sample Size",
#     "orig.n"     = "Original Sample Size",
#     "bias.level" = "Bias Level",
#     "orig.alpha" = "Original α Level"
#   )
#   
#   ggplot(df, aes_string(x = x_var, y = "mean_rate", color = "method", group = "method")) +
#     geom_line(size = 1) +
#     geom_point(size = 3) +
#     scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
#     labs(
#       title = paste0(outcome_var, " across ", x_labels[x_var],
#                      " (BF cutoff = ", cutoff_vals,
#                      ", true effect = ", eff_vals, ")"),
#       x = x_labels[x_var],
#       y = paste0("Mean ", outcome_var),
#       color = "Method"
#     ) +
#     theme_minimal(base_size = 14)
# }
# 
# # --- Shiny app ---
# ui <- fluidPage(
#   titlePanel("MABF Simulation Results"),
#   sidebarLayout(
#     sidebarPanel(
#       selectInput("dataset_group", "Choose Dataset:",
#                   choices = c("Null"="null", "Effect size 0.2"="0.2", "Effect size 0.5"="0.5")),
#       selectInput("outcome", "Choose Outcome:",
#                   choices = c("ADR","TPR","FNR","TNR","FPR","TSR","FSR","TFR","FFR")),
#       selectInput("xvar", "X-axis Variable:",
#                   choices = c("rep.number","rep.n","orig.n","bias.level","orig.alpha")),
#       selectInput("cutoff", "Bayes Factor Cutoff:", choices = c("1","3","10","30")),
#       checkboxGroupInput("methods", "Methods:",
#                          choices = c("FEMABF","EUBF","BFbMA","iBF"),
#                          selected = c("FEMABF","EUBF","BFbMA","iBF")),
#       checkboxInput("show_boxplot", "Show Boxplots", value = TRUE)
#     ),
#     mainPanel(
#       plotOutput("plot"),
#       h3("Exact Values"),
#       DTOutput("table")
#     )
#   )
# )
# 
# server <- function(input, output) {
#   
#   # raw data (needed for boxplots)
#   reactive_raw <- reactive({
#     data <- switch(input$dataset_group,
#                    "0.2"  = all_data_0.2,
#                    "0.5"  = all_data_0.5,
#                    "null" = all_data_null)
#     
#     data %>%
#       filter(Outcome == input$outcome,
#              BFcutoff %in% input$cutoff,
#              method %in% input$methods)
#   })
#   
#   # summary data (needed for means table + line plot)
#   reactive_data <- reactive({
#     df <- reactive_raw()
#     df %>%
#       group_by(method, !!sym(input$xvar), Outcome, BFcutoff, true.effect) %>%
#       summarise(
#         mean_rate   = round(mean(Rate, na.rm = TRUE), 2),
#         sd_rate     = round(sd(Rate, na.rm = TRUE), 2),
#         prop_missing = round(mean(is.na(Rate)), 4),  
#         .groups     = "drop"
#       )
#   })
#   
#   output$plot <- renderPlot({
#     df_raw  <- reactive_raw()
#     df_mean <- reactive_data()
#     
#     # APA-style labels
#     x_labels <- c(
#       "rep.number" = "Number of Replications",
#       "rep.n"      = "Replication Sample Size",
#       "orig.n"     = "Original Sample Size",
#       "bias.level" = "Bias Level",
#       "orig.alpha" = "Original α Level"
#     )
#     
#     p <- ggplot() 
#     
#     # Add boxplots if checkbox is ticked
#     if (input$show_boxplot) {
#       p <- p + geom_boxplot(data = df_raw,
#                             aes_string(x = input$xvar, y = "Rate", fill = "method",
#                                        group = paste0("interaction(method,", input$xvar, ")")),
#                             alpha = 0.2, width = 0.3, position = position_dodge(width = 0.6),
#                             outlier.shape = NA, color = NA)
#     }
#     
#     # Add mean line + points
#     p <- p +
#       geom_line(data = df_mean,
#                 aes_string(x = input$xvar, y = "mean_rate", color = "method", group = "method"),
#                 size = 1, position = position_dodge(width = 0.6)) +
#       geom_point(data = df_mean,
#                  aes_string(x = input$xvar, y = "mean_rate", color = "method"),
#                  size = 3, position = position_dodge(width = 0.6)) +
#       scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
#       labs(
#         title = paste0(input$outcome, " across ", x_labels[input$xvar],
#                        " (BF cutoff = ", unique(df_mean$BFcutoff),
#                        ", true effect = ", unique(df_mean$true.effect), ")"),
#         x = x_labels[input$xvar],
#         y = paste0("Rate / Mean ", input$outcome),
#         color = "Method",
#         fill = "Method"
#       ) +
#       theme_minimal(base_size = 14)
#     
#     p
#   })
#   
#   
#   output$table <- renderDT({
#     df <- reactive_data()
#     datatable(df, options = list(pageLength = 25), rownames = FALSE)
#   })
# }
# 
# 
# shinyApp(ui, server)
# 


#####
######
#Updated version: this version show theta = 0, 0.2, 0.5 separately
library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(DT)
library(grid)

load("./MABFanalyses/row-wise/LinePlotData.RData")

ui <- fluidPage(
  titlePanel("MABF Simulation Results"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "dataset_group",
        "Choose Dataset:",
        choices = c("Null" = "null", "Effect size 0.2" = "0.2", "Effect size 0.5" = "0.5")
      ),
      selectInput(
        "outcome",
        "Choose Outcome:",
        choices = c("ADR","TPR","FNR","TNR","FPR","TSR","FSR","TFR","FFR")
      ),
      selectInput(
        "cutoff",
        "Bayes Factor Cutoff:",
        choices = c("1","3","10","30")
      ),
      checkboxGroupInput(
        "methods",
        "Methods:",
        choices = c("FEMABF","EUBF","BFbMA","iBF"),
        selected = c("FEMABF","EUBF","BFbMA","iBF")
      ),
      selectInput(
        "interval_type",
        "Display interval:",
        choices = c(
          "Boxplot" = "boxplot",
          "MCSE error bar" = "mcse",
          "95% CI based on MCSE" = "ci95"
        ),
        selected = "mcse"
      ),
      hr(),
      h4("Panel layout"),
      selectInput(
        "facet_row",
        "Panel rows:",
        choices = c("rep.number", "rep.n", "orig.n", "bias.level", "orig.alpha"),
        selected = "rep.number"
      ),
      selectInput(
        "facet_col",
        "Panel columns:",
        choices = c("rep.number", "rep.n", "orig.n", "bias.level", "orig.alpha"),
        selected = "rep.n"
      )
    ),
    mainPanel(
      plotOutput("plot", height = "900px"),
      h3("Exact Values"),
      DTOutput("table")
    )
  )
)

server <- function(input, output, session) {
  
  design_vars <- c("rep.number", "rep.n", "orig.n", "bias.level", "orig.alpha")
  
  observeEvent(input$facet_row, {
    current_col <- isolate(input$facet_col)
    new_choices <- design_vars[design_vars != input$facet_row]
    new_selected <- if (!is.null(current_col) && current_col %in% new_choices) current_col else new_choices[1]
    
    updateSelectInput(
      session,
      "facet_col",
      choices = setNames(new_choices, new_choices),
      selected = new_selected
    )
  }, ignoreInit = TRUE)
  
  observeEvent(input$facet_col, {
    current_row <- isolate(input$facet_row)
    new_choices <- design_vars[design_vars != input$facet_col]
    new_selected <- if (!is.null(current_row) && current_row %in% new_choices) current_row else new_choices[1]
    
    updateSelectInput(
      session,
      "facet_row",
      choices = setNames(new_choices, new_choices),
      selected = new_selected
    )
  }, ignoreInit = TRUE)
  
  reactive_raw <- reactive({
    data <- switch(
      input$dataset_group,
      "0.2"  = all_data_0.2,
      "0.5"  = all_data_0.5,
      "null" = all_data_null
    )
    
    data %>%
      filter(
        Outcome == input$outcome,
        BFcutoff %in% input$cutoff,
        method %in% input$methods
      ) %>%
      mutate(
        rep.number = factor(rep.number, levels = c(2, 5, 10)),
        rep.n      = factor(rep.n, levels = c(40, 100, 400)),
        orig.n     = factor(orig.n, levels = c(20, 50, 200)),
        bias.level = factor(bias.level, levels = c("low", "medium", "high")),
        orig.alpha = factor(orig.alpha, levels = c(0.01, 0.05)),
        method     = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF"))
      )
  })
  
  reactive_data <- reactive({
    df <- reactive_raw()
    
    facet_row_levels <- levels(df[[input$facet_row]])
    facet_col_levels <- levels(df[[input$facet_col]])
    
    df %>%
      group_by(
        method,
        .data[[input$facet_row]],
        .data[[input$facet_col]],
        Outcome,
        BFcutoff,
        true.effect
      ) %>%
      summarise(
        mean_rate    = mean(Rate, na.rm = TRUE),
        sd_rate      = sd(Rate, na.rm = TRUE),
        mcse_rate    = sd(Rate, na.rm = TRUE) / sqrt(500),
        ci_lower     = pmax(0, mean_rate - 1.96 * mcse_rate),
        ci_upper     = pmin(1, mean_rate + 1.96 * mcse_rate),
        prop_missing = mean(is.na(Rate)),
        .groups      = "drop"
      ) %>%
      rename(
        facet_row_value = !!sym(input$facet_row),
        facet_col_value = !!sym(input$facet_col)
      ) %>%
      mutate(
        facet_row_value = factor(facet_row_value, levels = facet_row_levels),
        facet_col_value = factor(facet_col_value, levels = facet_col_levels),
        method = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF")),
        mean_label = sprintf("%.2f", mean_rate)
      )
  })
  
  output$plot <- renderPlot({
    df_raw  <- reactive_raw()
    df_mean <- reactive_data()
    
    req(nrow(df_mean) > 0)
    
    y_label <- switch(
      input$interval_type,
      "boxplot" = paste0(input$outcome, " Proportion"),
      "mcse"    = paste0("Mean ", input$outcome, " ± MCSE"),
      "ci95"    = paste0("Mean ", input$outcome, " with 95% CI")
    )
    
    df_raw_plot <- df_raw %>%
      mutate(
        facet_row_value = factor(.data[[input$facet_row]], levels = levels(df_raw[[input$facet_row]])),
        facet_col_value = factor(.data[[input$facet_col]], levels = levels(df_raw[[input$facet_col]]))
      )
    
    row_labeller <- function(x) {
      if (input$facet_row == "rep.number") {
        paste0("N[rep] == ", x)
      } else if (input$facet_row == "rep.n") {
        paste0("n[rep] == ", x)
      } else if (input$facet_row == "orig.n") {
        paste0("n[orig] == ", x)
      } else if (input$facet_row == "orig.alpha") {
        paste0("alpha[orig] == ", x)
      } else {
        as.character(x)
      }
    }
    
    col_labeller <- function(x) {
      if (input$facet_col == "rep.number") {
        paste0("N[rep] == ", x)
      } else if (input$facet_col == "rep.n") {
        paste0("n[rep] == ", x)
      } else if (input$facet_col == "orig.n") {
        paste0("n[orig] == ", x)
      } else if (input$facet_col == "orig.alpha") {
        paste0("alpha[orig] == ", x)
      } else {
        as.character(x)
      }
    }
    
    p <- ggplot()
    
    if (input$interval_type == "boxplot") {
      p <- p +
        geom_boxplot(
          data = df_raw_plot,
          aes(x = method, y = Rate, fill = method),
          alpha = 0.25,
          width = 0.5,
          outlier.shape = NA,
          color = NA
        ) +
        geom_text(
          data = df_mean,
          aes(
            x = method,
            y = pmin(1, mean_rate + 0.03),
            label = mean_label,
            color = method
          ),
          size = 3.5,
          vjust = 0,
          show.legend = FALSE
        )
    }
    
    if (input$interval_type == "mcse") {
      p <- p +
        geom_errorbar(
          data = df_mean,
          aes(
            x = method,
            ymin = pmax(0, mean_rate - mcse_rate),
            ymax = pmin(1, mean_rate + mcse_rate),
            color = method
          ),
          width = 0.18,
          size = 0.8
        ) +
        geom_point(
          data = df_mean,
          aes(x = method, y = mean_rate, color = method),
          shape = 95,
          size = 8
        ) +
        geom_text(
          data = df_mean,
          aes(
            x = method,
            y = pmin(1, mean_rate + mcse_rate + 0.03),
            label = mean_label,
            color = method
          ),
          size = 3.5,
          vjust = 0,
          show.legend = FALSE
        )
    }
    
    if (input$interval_type == "ci95") {
      p <- p +
        geom_errorbar(
          data = df_mean,
          aes(
            x = method,
            ymin = ci_lower,
            ymax = ci_upper,
            color = method
          ),
          width = 0.18,
          size = 0.8
        ) +
        geom_point(
          data = df_mean,
          aes(x = method, y = mean_rate, color = method),
          shape = 95,
          size = 8
        ) +
        geom_text(
          data = df_mean,
          aes(
            x = method,
            y = pmin(1, ci_upper + 0.03),
            label = mean_label,
            color = method
          ),
          size = 3.5,
          vjust = 0,
          show.legend = FALSE
        )
    }
    
    p <- p +
      facet_grid(
        rows = vars(facet_row_value),
        cols = vars(facet_col_value),
        labeller = labeller(
          facet_row_value = row_labeller,
          facet_col_value = col_labeller,
          .default = label_parsed
        )
      ) +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.1)
      ) +
      labs(
        title = paste0(
          input$outcome, " by Method",
          " (BF cutoff = ", unique(df_mean$BFcutoff),
          ", true effect = ", unique(df_mean$true.effect), ")"
        ),
        x = "Method",
        y = y_label,
        color = "Method",
        fill = "Method"
      ) +
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 4, linetype = 0)),
        fill  = guide_legend(override.aes = list(shape = 15, size = 5))
      ) +
      theme_minimal(base_size = 14) +
      theme(
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.8),
        panel.spacing = unit(1, "lines"),
        panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey94", linewidth = 0.2),
        strip.background = element_rect(fill = "grey90", color = "grey40", linewidth = 0.6)
      )
    
    p
  })
  
  output$table <- renderDT({
    df <- reactive_data() %>%
      mutate(
        mean_rate    = round(mean_rate, 2),
        sd_rate      = round(sd_rate, 3),
        mcse_rate    = round(mcse_rate, 4),
        ci_lower     = round(ci_lower, 4),
        ci_upper     = round(ci_upper, 4),
        prop_missing = round(prop_missing, 4)
      )
    
    datatable(df, options = list(pageLength = 25), rownames = FALSE)
  })
}

shinyApp(ui, server)





######
#This version show theta = 0, 0.2, 0.5 in the same 3 x 3 grid 
library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(DT)
library(grid)

load("./MABFanalyses/row-wise/LinePlotData.RData")

ui <- fluidPage(
  titlePanel("MABF Simulation Results"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "outcome",
        "Choose Outcome:",
        choices = c("ADR","TPR","FNR","TNR","FPR","TSR","FSR","TFR","FFR")
      ),
      selectInput(
        "cutoff",
        "Bayes Factor Cutoff:",
        choices = c("1","3","10","30")
      ),
      checkboxGroupInput(
        "methods",
        "Methods:",
        choices = c("FEMABF","EUBF","BFbMA","iBF"),
        selected = c("FEMABF","EUBF","BFbMA","iBF")
      ),
      checkboxGroupInput(
        "datasets",
        "Datasets:",
        choices = c("Null" = "null", "Effect size 0.2" = "0.2", "Effect size 0.5" = "0.5"),
        selected = c("null", "0.2", "0.5")
      ),
      selectInput(
        "interval_type",
        "Display interval:",
        choices = c(
          "MCSE error bar" = "mcse",
          "95% CI based on MCSE" = "ci95",
          "Boxplot" = "boxplot"
        ),
        selected = "mcse"
      ),
      checkboxInput(
        "show_interval",
        "Show interval",
        value = TRUE
      ),
      hr(),
      h4("Panel layout"),
      selectInput(
        "facet_row",
        "Panel rows:",
        choices = c("rep.number", "rep.n", "orig.n", "bias.level", "orig.alpha"),
        selected = "rep.number"
      ),
      selectInput(
        "facet_col",
        "Panel columns:",
        choices = c("rep.number", "rep.n", "orig.n", "bias.level", "orig.alpha"),
        selected = "rep.n"
      )
    ),
    mainPanel(
      plotOutput("plot", height = "900px", width = "600px"),
      h3("Exact Values"),
      DTOutput("table")
    )
  )
)

server <- function(input, output, session) {
  
  design_vars <- c("rep.number", "rep.n", "orig.n", "bias.level", "orig.alpha")
  
  observeEvent(input$facet_row, {
    current_col <- isolate(input$facet_col)
    new_choices <- design_vars[design_vars != input$facet_row]
    new_selected <- if (!is.null(current_col) && current_col %in% new_choices) current_col else new_choices[1]
    
    updateSelectInput(
      session,
      "facet_col",
      choices = setNames(new_choices, new_choices),
      selected = new_selected
    )
  }, ignoreInit = TRUE)
  
  observeEvent(input$facet_col, {
    current_row <- isolate(input$facet_row)
    new_choices <- design_vars[design_vars != input$facet_col]
    new_selected <- if (!is.null(current_row) && current_row %in% new_choices) current_row else new_choices[1]
    
    updateSelectInput(
      session,
      "facet_row",
      choices = setNames(new_choices, new_choices),
      selected = new_selected
    )
  }, ignoreInit = TRUE)
  
  reactive_raw <- reactive({
    req(length(input$datasets) > 0)
    
    data_list <- list()
    
    if ("null" %in% input$datasets) {
      data_list[["null"]] <- all_data_null %>% mutate(dataset_group = "null")
    }
    if ("0.2" %in% input$datasets) {
      data_list[["0.2"]] <- all_data_0.2 %>% mutate(dataset_group = "0.2")
    }
    if ("0.5" %in% input$datasets) {
      data_list[["0.5"]] <- all_data_0.5 %>% mutate(dataset_group = "0.5")
    }
    
    bind_rows(data_list) %>%
      filter(
        Outcome == input$outcome,
        BFcutoff %in% input$cutoff,
        method %in% input$methods
      ) %>%
      mutate(
        rep.number = factor(rep.number, levels = c(2, 5, 10)),
        rep.n      = factor(rep.n, levels = c(40, 100, 400)),
        orig.n     = factor(orig.n, levels = c(20, 50, 200)),
        bias.level = factor(bias.level, levels = c("low", "medium", "high")),
        orig.alpha = factor(orig.alpha, levels = c(0.01, 0.05)),
        method     = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF")),
        dataset_group = factor(dataset_group,
                               levels = c("null", "0.2", "0.5"),
                               labels = c("Null", "0.2", "0.5"))
      )
  })
  
  reactive_data <- reactive({
    df <- reactive_raw()
    
    facet_row_levels <- levels(df[[input$facet_row]])
    facet_col_levels <- levels(df[[input$facet_col]])
    
    df %>%
      group_by(
        method,
        dataset_group,
        .data[[input$facet_row]],
        .data[[input$facet_col]],
        Outcome,
        BFcutoff
      ) %>%
      summarise(
        mean_rate    = mean(Rate, na.rm = TRUE),
        sd_rate      = sd(Rate, na.rm = TRUE),
        mcse_rate    = sd(Rate, na.rm = TRUE) / sqrt(500),
        ci_lower     = pmax(0, mean_rate - 1.96 * mcse_rate),
        ci_upper     = pmin(1, mean_rate + 1.96 * mcse_rate),
        prop_missing = mean(is.na(Rate)),
        .groups      = "drop"
      ) %>%
      rename(
        facet_row_value = !!sym(input$facet_row),
        facet_col_value = !!sym(input$facet_col)
      ) %>%
      mutate(
        facet_row_value = factor(facet_row_value, levels = facet_row_levels),
        facet_col_value = factor(facet_col_value, levels = facet_col_levels),
        method = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF")),
        dataset_group = factor(dataset_group, levels = c("Null", "0.2", "0.5")),
        mean_label = sprintf("%.2f", mean_rate)
      )
  })
  
  output$plot <- renderPlot({
    df_raw  <- reactive_raw()
    df_mean <- reactive_data()
    
    req(nrow(df_mean) > 0)
    
    y_label <- switch(
      input$interval_type,
      "boxplot" = paste0(input$outcome, " Proportion"),
      "mcse"    = paste0("Mean ", input$outcome, " ± MCSE"),
      "ci95"    = paste0("Mean ", input$outcome, " with 95% CI")
    )
    
    df_raw_plot <- df_raw %>%
      mutate(
        facet_row_value = factor(.data[[input$facet_row]], levels = levels(df_raw[[input$facet_row]])),
        facet_col_value = factor(.data[[input$facet_col]], levels = levels(df_raw[[input$facet_col]]))
      )
    
    row_labeller <- function(x) {
      if (input$facet_row == "rep.number") {
        paste0("N[rep] == ", x)
      } else if (input$facet_row == "rep.n") {
        paste0("n[rep] == ", x)
      } else if (input$facet_row == "orig.n") {
        paste0("n[orig] == ", x)
      } else if (input$facet_row == "orig.alpha") {
        paste0("alpha[orig] == ", x)
      } else {
        as.character(x)
      }
    }
    
    col_labeller <- function(x) {
      if (input$facet_col == "rep.number") {
        paste0("N[rep] == ", x)
      } else if (input$facet_col == "rep.n") {
        paste0("n[rep] == ", x)
      } else if (input$facet_col == "orig.n") {
        paste0("n[orig] == ", x)
      } else if (input$facet_col == "orig.alpha") {
        paste0("alpha[orig] == ", x)
      } else {
        as.character(x)
      }
    }
    
    dodge_width <- 0.55
    
    p <- ggplot()
    
    if (isTRUE(input$show_interval) && input$interval_type == "boxplot") {
      p <- p +
        geom_boxplot(
          data = df_raw_plot,
          aes(x = method, y = Rate, fill = method),
          alpha = 0.15,
          width = 0.45,
          outlier.shape = NA,
          color = NA
        )
    }
    
    if (isTRUE(input$show_interval) && input$interval_type == "mcse") {
      p <- p +
        geom_errorbar(
          data = df_mean,
          aes(
            x = method,
            ymin = pmax(0, mean_rate - mcse_rate),
            ymax = pmin(1, mean_rate + mcse_rate),
            color = method,
            group = interaction(method, dataset_group)
          ),
          width = 0.28,
          linewidth = 1.1,
          position = position_dodge(width = dodge_width)
        )
    }
    
    if (isTRUE(input$show_interval) && input$interval_type == "ci95") {
      p <- p +
        geom_errorbar(
          data = df_mean,
          aes(
            x = method,
            ymin = ci_lower,
            ymax = ci_upper,
            color = method,
            group = interaction(method, dataset_group)
          ),
          width = 0.28,
          linewidth = 1.1,
          position = position_dodge(width = dodge_width)
        )
    }
    
    p <- p +
      geom_point(
        data = df_mean,
        aes(
          x = method,
          y = mean_rate,
          color = method,
          shape = dataset_group,
          group = interaction(method, dataset_group)
        ),
        size = 3.2,
        stroke = 0.9,
        position = position_dodge(width = dodge_width)
      ) +
      facet_grid(
        rows = vars(facet_row_value),
        cols = vars(facet_col_value),
        labeller = labeller(
          facet_row_value = row_labeller,
          facet_col_value = col_labeller,
          .default = label_parsed
        )
      ) +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.1)
      ) +
      scale_color_manual(
        values = c(
          "BFbMA" = "blue2",
          "EUBF" = "#7CAE00",
          "FEMABF" = "firebrick2",
          "iBF" = "burlywood4"
        ),
        name = "Method"
      ) +
      scale_shape_manual(
        values = c("Null" = 3, "0.2" = 4, "0.5" = 8),
        name = "Dataset"
      ) +
      labs(
        title = paste0(
          input$outcome, " by Method and Dataset",
          " (BF cutoff = ", unique(df_mean$BFcutoff), ")"
        ),
        x = "Method",
        y = y_label,
        color = "Method",
        fill = "Method"
      ) +
      guides(
        color = guide_legend(order = 1, override.aes = list(shape = 16, size = 4, linetype = 0)),
        shape = guide_legend(order = 2),
        fill  = guide_legend(order = 3, override.aes = list(shape = 15, size = 5))
      ) +
      theme_minimal(base_size = 14) +
      theme(
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.8),
        panel.spacing = unit(1, "lines"),
        panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey94", linewidth = 0.2),
        strip.background = element_rect(fill = "grey90", color = "grey40", linewidth = 0.6)
      )
    
    p
  })
  
  output$table <- renderDT({
    df <- reactive_data() %>%
      mutate(
        mean_rate    = round(mean_rate, 2),
        sd_rate      = round(sd_rate, 3),
        mcse_rate    = round(mcse_rate, 4),
        ci_lower     = round(ci_lower, 4),
        ci_upper     = round(ci_upper, 4),
        prop_missing = round(prop_missing, 4)
      )
    
    datatable(df, options = list(pageLength = 25), rownames = FALSE)
  })
}

shinyApp(ui, server)