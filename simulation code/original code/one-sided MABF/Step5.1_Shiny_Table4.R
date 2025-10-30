# Stacked bar plot + Summary Table + APA Table (LaTeX, dynamic metrics)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(knitr)
library(kableExtra)
library(scales)
library(htmltools)
library(purrr)
library(tibble)

# === Data Preparation Function ===
load_and_prepare_data <- function(cutoff_folder = "c1") {
  load_group <- function(folder_path) {
    rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
    lapply(rds_files, function(file_path) {
      data <- readRDS(file_path)
      name <- tools::file_path_sans_ext(basename(file_path))
      assign(name, data, envir = .GlobalEnv)
    })
  }
  
  get_data_combined <- function(pcutoff, true_es) {
    folder <- paste0("./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=",
                     pcutoff, ", EScutoff_o=0/", cutoff_folder, "/")
    load_group(folder)
    method_list <- list(
      BFbMA = get(paste0("rates_BFbMA_", true_es, "null_", cutoff_folder)),
      EUBF = get(paste0("rates_EUBF_", true_es, "null_", cutoff_folder)),
      FEMABF = get(paste0("rates_FEMABF_", true_es, "null_", cutoff_folder)),
      iBF = get(paste0("rates_iBF_", true_es, "null_", cutoff_folder))
    )
    bind_rows(lapply(names(method_list), function(method) {
      method_list[[method]] %>%
        mutate(method = method) %>%
        relocate(method, .before = everything())
    })) %>%
      mutate(p_cutoff = as.numeric(pcutoff),
             true_es = as.numeric(true_es)) %>%
      relocate(true_es, p_cutoff)
  }
  
  bind_rows(
    get_data_combined("0.01", "0.2"),
    get_data_combined("0.01", "0.5"),
    get_data_combined("0.05", "0.2"),
    get_data_combined("0.05", "0.5")
  )
}

# === UI ===
ui <- fluidPage(
  titlePanel("MABF Stacked Bar Plot Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("bf_cutoff", "Bayes Factor Cutoff",
                  choices = c("3" = "c1", "10" = "c2"), selected = "c1"),
      selectInput("dataset_choice", "Data Scenario",
                  choices = c("Null" = "null", "0.2", "0.5"), selected = "null"),
      br(),
      downloadButton("download_plot", "Download Plot (PNG)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Stacked Plot", plotOutput("piePlot", height = "800px")),
        tabPanel("Summary Table", DT::dataTableOutput("summaryTable")),
        tabPanel("APA Table",
                 htmlOutput("apa_html"),
                 br(),
                 downloadButton("download_apa_tex", "Download APA Table (LaTeX)"),
                 br(),
                 verbatimTextOutput("apa_latex"))
      )
    )
  )
)

# === Server ===
server <- function(input, output) {
  rates_all <- reactive({
    load_and_prepare_data(input$bf_cutoff)
  })
  
  filtered_data <- reactive({
    use_rate <- if (input$dataset_choice == "null") "FP" else "TP"
    true_es_num <- if (input$dataset_choice == "null") 0.2 else as.numeric(input$dataset_choice)
    
    metric_cols <- if (use_rate == "TP") {
      c("ADTE", "TP", "FN")
    } else {
      c("ADNE", "TN", "FP")
    }
    
    rates_all() %>%
      filter(true_es == true_es_num) %>%
      mutate(
        rep.number = factor(rep.number, levels = c(2, 5, 10),
                            labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10")),
        rep.n = factor(rep.n, levels = c(40, 100, 400),
                       labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400")),
        method = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF"))
      ) %>%
      select(method, rep.n, rep.number, all_of(metric_cols)) %>%
      pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "value") %>%
      mutate(
        metric_label = dplyr::recode(metric,
                                     ADTE = "AD", TP = "TP", FN = "FN",
                                     ADNE = "AD", TN = "TN", FP = "FP"),
        metric_label = factor(metric_label, levels = c("AD", "TP", "FN", "TN", "FP"))
      ) %>%
      group_by(method, rep.n, rep.number, metric_label) %>%
      summarise(value = sum(value), .groups = "drop") %>%
      group_by(method, rep.n, rep.number) %>%
      mutate(proportion = value / sum(value))
  })
  
 
  # --- Plot object for reuse ---
  plot_obj <- reactive({
    plot_data <- filtered_data()
    ggplot(plot_data, aes(x = method, y = proportion, fill = metric_label)) +
      geom_bar(stat = "identity", position = "stack") +
      facet_grid(rep.n ~ rep.number,
                 labeller = labeller(rep.number = label_parsed,
                                     rep.n = label_parsed)) +
      scale_fill_manual(values = c("AD" = "red", "FP" = "green", "TN" = "blue",
                                   "TP" = "blue", "FN" = "green"),
                        name = "Metric") +
      labs(x = "MABF Method", y = "Proportion") +
      theme_gray(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom")
  })
  
  output$piePlot <- renderPlot({
    plot_obj()
  })
  # --- Download plot ---
  output$download_plot <- downloadHandler(
    filename = function() {
      cutoff_lbl <- ifelse(input$bf_cutoff == "c1", "3", "10")
      scen_lbl   <- input$dataset_choice
      paste0("StackedPlot_1sided_BF", cutoff_lbl, "_Scenario", scen_lbl, ".png")
    },
    content = function(file) {
      ggsave(file, plot = plot_obj(), device = "png", width = 10, height = 8, dpi = 300)
    }
  )
  
  
  # --- Summary Table ---
  output$summaryTable <- DT::renderDataTable({
    df <- filtered_data() %>%
      select(method, rep.number, rep.n, metric_label, value, proportion) %>%
      pivot_wider(names_from = metric_label,
                  values_from = c(value, proportion),
                  names_glue = "{metric_label}_{.value}") %>%
      replace(is.na(.), 0) %>%
      mutate(rep.number = gsub("N\\[rep\\] == ", "", as.character(rep.number)),
             rep.n      = gsub("n\\[rep\\] == ", "", as.character(rep.n))) %>%
      mutate(across(ends_with("_proportion"),
                    ~ formatC(as.numeric(.), format = "f", digits = 4)))
    
    colnames(df) <- dplyr::recode(colnames(df),
                                  "AD_value" = "AD (n)", "TP_value" = "TP (n)",
                                  "FN_value" = "FN (n)", "TN_value" = "TN (n)", "FP_value" = "FP (n)",
                                  "AD_proportion" = "AD (%)", "TP_proportion" = "TP (%)",
                                  "FN_proportion" = "FN (%)", "TN_proportion" = "TN (%)", "FP_proportion" = "FP (%)"
    )
    df
  }, options = list(pageLength = 20, autoWidth = TRUE))
  
  # --- APA Table ---
  apa_tbl <- reactive({
    methods <- c("BFbMA", "EUBF", "FEMABF", "iBF")
    
    # choose metrics based on scenario
    metrics <- if (input$dataset_choice == "null") {
      c("AD","TN","FP")
    } else {
      c("AD","TP","FN")
    }
    
    expected_cols <- unlist(lapply(methods, function(m) paste0(m, "_", metrics)))
    
    df <- filtered_data() %>%
      filter(metric_label %in% metrics) %>%
      group_by(method, rep.number, rep.n, metric_label) %>%
      summarise(rate = mean(proportion), .groups = "drop") %>%
      mutate(rep.number = gsub("N\\[rep\\] == ", "", as.character(rep.number)),
             rep.n      = gsub("n\\[rep\\] == ", "", as.character(rep.n)),
             rep.number = factor(rep.number, levels = c("2","5","10"), ordered = TRUE),
             rep.n      = factor(rep.n, levels = c("40","100","400"), ordered = TRUE)) %>%
      arrange(rep.number, rep.n, method, metric_label) %>%
      unite(".col", method, metric_label, sep = "_") %>%
      pivot_wider(names_from = ".col", values_from = rate) %>%
      { 
        mis <- setdiff(expected_cols, names(.))
        if (length(mis) > 0) {
          for (cn in mis) .[[cn]] <- 0
        }
        .
      } %>%
      select(rep.number, rep.n, all_of(expected_cols)) %>%
      mutate(rep.number = as.character(rep.number),
             rep.n      = as.character(rep.n))
    df
  })
  
  apa_kable <- reactive({
    df <- apa_tbl()
    methods <- c("BFbMA", "EUBF", "FEMABF", "iBF")
    metrics <- if (input$dataset_choice == "null") c("AD","TN","FP") else c("AD","TP","FN")
    expected_cols <- unlist(lapply(methods, function(m) paste0(m, "_", metrics)))
    
    
    df_fmt <- df %>%
      mutate(across(
        all_of(expected_cols),
        ~ formatC(as.numeric(.) * 100, format = "f", digits = 2)
      ))
    
    cutoff_label <- ifelse(input$bf_cutoff == "c1", 3, 10)
    cols <- names(df_fmt)
    hdr_list <- as.list(rep("", length(cols)))
    names(hdr_list) <- cols
    hdr_list[[1]] <- paste0("\\textbf{BF cutoff = ", cutoff_label, "}")
    hdr_tbl <- tibble::as_tibble_row(hdr_list)
    
    df_nested <- df_fmt %>%
      mutate(rep.number = ifelse(duplicated(rep.number), "", rep.number))
    
    bind_rows(hdr_tbl, df_nested)
  })
  
  output$apa_html <- renderUI({
    methods <- c("BFbMA", "EUBF", "FEMABF", "iBF")
    metrics <- if (input$dataset_choice == "null") c("AD","TN","FP") else c("AD","TP","FN")
    header_methods <- c(" " = 2, setNames(rep(length(metrics), length(methods)), methods))
    header_rates   <- c(" " = 2, rep(paste0(metrics, " (%)"), times = length(methods)))
    
    tab <- apa_kable()
    kable(tab, format = "html",
          col.names = c("$N_{rep}$", "$n_{rep}$", rep("", ncol(tab) - 2)),
          escape = FALSE) %>%
      kable_styling(full_width = FALSE, bootstrap_options = c("condensed", "hover")) %>%
      add_header_above(header_rates, bold = FALSE) %>%
      add_header_above(header_methods, bold = FALSE, line = FALSE) %>%
      as.character() %>% HTML()
  })
  
  output$apa_latex <- renderText({
    methods <- c("BFbMA", "EUBF", "FEMABF", "iBF")
    metrics <- if (input$dataset_choice == "null") c("AD","TN","FP") else c("AD","TP","FN")
    header_methods <- c(" " = 2, setNames(rep(length(metrics), length(methods)), methods))
    header_rates   <- c("$N_{\\\\text{rep}}$", "$n_{\\\\text{rep}}$",
                        rep(paste0(metrics, " (\\\\%)"), times = length(methods)))
    
    tab <- apa_kable()
    
    # Dynamic caption
    cutoff_lbl <- ifelse(input$bf_cutoff == "c1", "3", "10")
    caption_text <- if (input$dataset_choice == "null") {
      paste0("Proportions of anecdotal evidence, false positives, and true negatives across MABF methods (one-sided test)",
             "when the Bayes factor cutoff for anecdotal evidence cases is ", cutoff_lbl, ".")
    } else {
      paste0("Proportions of anecdotal evidence, true positives, and false negatives across MABF methods (one-sided test)",
             "when the Bayes factor cutoff for anecdotal evidence cases is ", cutoff_lbl, ".")
    }
    
    as.character(
      kable(tab, format = "latex", booktabs = TRUE,
            col.names = rep("", ncol(tab)),
            escape = FALSE, linesep = "",
            align = c("l","r", rep("r", ncol(tab) - 2)),
            caption = caption_text) %>%
        kable_styling(latex_options = "hold_position", full_width = FALSE) %>%
        add_header_above(header_rates,   bold = FALSE, line = FALSE, escape = FALSE) %>%
        add_header_above(header_methods, bold = FALSE, line = TRUE,  escape = FALSE)
    )
  })
  
  output$download_apa_tex <- downloadHandler(
    filename = function() paste0("APA_table_BF", ifelse(input$bf_cutoff=="c1", "3", "10"), ".tex"),
    content = function(file) {
      writeLines(output$apa_latex(), file)
    }
  )
}


# === Run ===
shinyApp(ui = ui, server = server)
