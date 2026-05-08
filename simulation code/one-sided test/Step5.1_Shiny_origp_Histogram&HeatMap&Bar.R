# This Shiny app visualize FPR and TPR of original study p values using histogram, heatmap, and bar plot.
# The code is based on Step4.0_Data_Visualization_origp_BarPlotLineplot (A code) and Step4.0_Data_Visualization_origp_GridHistogram&HeatMap (B code).
# The bar plot show slight difference in rate values than A code, due to the different method of calculation.
library(shiny)
library(ggplot2)
library(dplyr)
library(purrr)
library(stringr)


# Load and prepare data globally for speed
folder_path <- "./OGDG/10000"
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
for (file_path in rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}
df_lists <- df_lists_10000
remove(df_lists_10000)

# Helper to extract data per condition
extract_condition_df <- function(effect_size) {
  effect_conditions <- list(
    "0" = c("0_0.02_20_none_low", "0_0.02_50_none_low", "0_0.02_200_none_low",
            "0_0.02_20_medium_medium", "0_0.02_50_medium_medium", "0_0.02_200_medium_medium",
            "0_0.02_20_high_high", "0_0.02_50_high_high", "0_0.02_200_high_high"),
    "0.2" = c("0.2_0.05_20_none_low", "0.2_0.05_50_none_low", "0.2_0.05_200_none_low",
              "0.2_0.05_20_medium_medium", "0.2_0.05_50_medium_medium", "0.2_0.05_200_medium_medium",
              "0.2_0.05_20_high_high", "0.2_0.05_50_high_high", "0.2_0.05_200_high_high"),
    "0.5" = c("0.5_0.125_20_none_low", "0.5_0.125_50_none_low", "0.5_0.125_200_none_low",
              "0.5_0.125_20_medium_medium", "0.5_0.125_50_medium_medium", "0.5_0.125_200_medium_medium",
              "0.5_0.125_20_high_high", "0.5_0.125_50_high_high", "0.5_0.125_200_high_high")
  )
  
  conditions_to_plot <- effect_conditions[[as.character(effect_size)]]
  
  map_dfr(conditions_to_plot, function(cond_name) {
    parts <- str_split(cond_name, "_")[[1]]
    effect <- as.numeric(parts[1])
    orig_n <- as.integer(parts[3])
    p_hack <- parts[4]
    pub_bias <- parts[5]
    
    bias_level <- case_when(
      p_hack == "none" & pub_bias == "low" ~ "Low",
      p_hack == "medium" & pub_bias == "medium" ~ "Medium",
      p_hack == "high" & pub_bias == "high" ~ "High",
      TRUE ~ "other"
    )
    
    df_lists[[cond_name]] %>%
      mutate(
        p_value = p,
        orig_n = factor(orig_n, levels = c(20, 50, 200)),
        bias_level = factor(bias_level, levels = c("Low", "Medium", "High")),
        is_significant = p < 0.05,
        metric_label = case_when(
          effect == 0 & p < 0.05 ~ "False Positive",
          effect == 0 & p >= 0.05 ~ "True Negative",
          effect > 0 & p < 0.05 ~ "True Positive",
          effect > 0 & p >= 0.05 ~ "False Negative"
        ),
        effect_size = effect
      )
  })
}

# UI
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
  
  titlePanel("Impact of p-Hacking and Publication Bias on Study Results"),
  
  sidebarLayout(
    sidebarPanel(
      class = "custom-sidebar",
      selectInput("plot_type", "Choose Plot Type", choices = c("Histogram", "Heatmap", "Bar Plot")),
      selectInput("metric", "Choose Metric", choices = c(
        "False Positive Rate" = "FPR",
        "True Positive Rate" = "TPR",
        "True Negative Rate" = "TNR",
        "False Negative Rate" = "FNR"
      )),
      selectInput("effect_size", "Effect Size", choices = c("0", "0.2", "0.5"))
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Study Overview",
                 tags$div(style = "padding: 30px; max-width: 900px;",
                          HTML("
              <p>Recent concerns about the replication crisis in psychology have underscored the need for better frameworks to assess whether original findings can be reliably replicated. This study aims to compare the performance of several meta-analytic Bayes factors (FEMABF, BFbMA, EUBF, and iBF) both relative to one another and against the traditional fixed-effect meta-analysis in evaluating replication success.</p>
              <p>To do so, we use a large-scale simulation to generate original studies across varied scenarios by manipulating true effect size, research environment (levels of p-hacking and publication bias), and original sample size. Replication studies are then simulated by varying replication sample size and number of replications, yielding 243 distinct scenarios and over 60 million data points.</p>
              <p>This Shiny app visualizes the impact of p-hacking and publicatin bias on original study results across various scenarios. Specifically, the underlying effect sizes—expressed as standardized mean differences—are set to 0, 0.2, and 0.5. Original sample sizes for both the treatment and control groups are set to 20, 50, and 200 participants. In addition, the simulated research environment is shaped by three levels of p‐hacking and publication bias (i.e., none, medium, and high). </p>
            ")
                 )
        ),
        
        tabPanel("Visualization",
                 plotOutput("ratePlot", height = "1000px")
        )
      )
    )
  )
)

# Server
server <- function(input, output) {
  plot_data <- reactive({
    extract_condition_df(as.numeric(input$effect_size))
  })
  
  output$ratePlot <- renderPlot({
    df <- plot_data()
    
    if (input$plot_type == "Histogram") {
      
      # 1. Define which label to count as "positive" based on user-selected metric
      metric_map <- list(
        FPR = "False Positive",
        TPR = "True Positive",
        FNR = "False Negative",
        TNR = "True Negative"
      )
      metric_label <- metric_map[[input$metric]]
      
      # 2. Set color scheme for the histogram bars
      fill_colors <- c(
        "False Positive" = "firebrick",
        "True Negative" = "lightblue",
        "True Positive" = "green",
        "False Negative" = "lightblue"
      )
      
      # 3. Keep only rows that match valid color categories
      df_filtered <- df %>% filter(metric_label %in% names(fill_colors))
      
      # 4. Create the summary dataframe with rate annotations per facet
      annotation_df <- df_filtered %>%
        group_by(bias_level, orig_n) %>%
        summarise(
          rate = mean(metric_label == !!metric_label),
          .groups = "drop"
        ) %>%
        mutate(
          label = paste0(input$metric, " = ", round(rate, 3)),
          #x = -0.1,
          #y = 3000
          x = median(df_filtered$d, na.rm = TRUE),  # middle
          y = max(table(cut(df_filtered$d, breaks = 200))) * 0.9  # ~90% height
        )
      
      # 5. Plot with histogram + annotation
      ggplot(df_filtered, aes(x = d, fill = metric_label)) +
        geom_histogram(bins = 50, position = "stack", alpha = 0.4, color = "gray30") +
        geom_text(
          data = annotation_df,
          aes(x = x, y = y, label = label),
          inherit.aes = FALSE,
          size = 5, fontface = "italic"
        ) +
        scale_fill_manual(values = fill_colors) +
        facet_grid(bias_level ~ orig_n, labeller = label_both) +
        labs(
          title = paste("Histogram of", input$metric, "by Effect Size =", input$effect_size),
          x = "Observed Effect Size (d)", y = "Count",
          fill = "Result Type"
        ) +
        coord_cartesian(ylim = c(0, 3000)) +
        theme_minimal() +
        theme(strip.text = element_text(size = 12, face = "bold"),
              legend.position = "bottom")
      
    } else if (input$plot_type == "Heatmap") {
      summary_df <- df %>%
        group_by(bias_level, orig_n) %>%
        summarise(
          rate = case_when(
            input$metric == "FPR" ~ mean(metric_label == "False Positive"),
            input$metric == "TPR" ~ mean(metric_label == "True Positive"),
            input$metric == "FNR" ~ mean(metric_label == "False Negative"),
            input$metric == "TNR" ~ mean(metric_label == "True Negative")
          ),
          .groups = "drop"
        )
      
      fill_label <- case_when(
        input$metric == "FPR" ~ "False Positive Rate",
        input$metric == "TPR" ~ "True Positive Rate",
        input$metric == "FNR" ~ "False Negative Rate",
        input$metric == "TNR" ~ "True Negative Rate"
      )
      
      ggplot(summary_df, aes(x = orig_n, y = bias_level, fill = rate)) +
        geom_tile(color = "white") +
        geom_text(aes(label = round(rate, 3)), size = 5) +
        scale_fill_gradient(low = "#fef6d6", high = "#fb972c", name = fill_label) +
        labs(
          title = paste(fill_label, "by Bias Level and Sample Size"),
          x = "Original Sample Size", y = "Bias Level"
        ) +
        theme_minimal(base_size = 14) +
        theme(panel.grid = element_blank())
    } else if (input$plot_type == "Bar Plot") {
      if (input$metric == "FPR") {
        df %>%
          filter(effect_size == 0) %>%
          group_by(orig_n, bias_level) %>%
          summarise(FPR = mean(metric_label == "False Positive"), .groups = "drop") %>%
          ggplot(aes(x = orig_n, y = FPR, fill = bias_level)) +
          geom_col(position = position_dodge(width = 0.9)) +
          geom_text(
            aes(label = scales::percent(FPR, accuracy = 1)),
            position = position_dodge(width = 0.9),
            vjust = -0.3,
            size = 5
          ) +
          labs(
            #title = "False Positive Rate When Effect Size = 0",
            x = "Original Sample Size",
            y = "False Positive Rate",
            fill = "Bias level"
          ) +
          scale_y_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, by = 0.1),
            labels = scales::percent_format(accuracy = 1)
          ) +
          theme_minimal(base_size = 14)
      } else if (input$metric == "TPR") {
        df %>%
          filter(effect_size > 0) %>%
          group_by(effect_size, orig_n, bias_level) %>%
          summarise(TPR = mean(metric_label == "True Positive"), .groups = "drop") %>%
          mutate(effect_size = factor(effect_size)) %>%
          ggplot(aes(x = orig_n, y = TPR, fill = bias_level)) +
          geom_col(position = position_dodge(width = 0.9)) +
          geom_text(
            aes(label = scales::percent(TPR, accuracy = 1)),
            position = position_dodge(width = 0.9),
            vjust = -0.3,
            size = 5
          ) +
          
          labs(
            #title = "True Positive Rate When Effect Size > 0",
            x = "Original Sample Size",
            y = "True Positive Rate",
            fill = "Bias level"
          ) +
          scale_y_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, by = 0.1),
            labels = scales::percent_format(accuracy = 1)
          ) +
          theme_minimal(base_size = 14)
        
      }
    }
    
  })
}

shinyApp(ui, server)
