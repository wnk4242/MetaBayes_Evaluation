
#####This version uses the van Ravenzwaaij and Ioannidis (2019) ROC curve style (no 'steps' like in the ROC curves)
library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rlang)
library(forcats)

# === Load and prepare data function ===
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

# === Load data ===
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

