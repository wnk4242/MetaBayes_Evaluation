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
  
  # --- Plot ---
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
      paste0("StackedPlot_2sided_BF", cutoff_lbl, "_Scenario", scen_lbl, ".png")
    },
    content = function(file) {
      ggsave(file, plot = plot_obj(), device = "png", width = 10, height = 8, dpi = 300)
    }
  )
  
  
  # --- Summary table ---
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
    
    # ✅ Remove % signs, just numbers with 2 decimals
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
    
    # ✅ Dynamic caption
    cutoff_lbl <- ifelse(input$bf_cutoff == "c1", "3", "10")
    caption_text <- if (input$dataset_choice == "null") {
      paste0("Proportions of anecdotal evidence, false positives, and true negatives across MABF methods (two-sided test)",
             "when the Bayes factor cutoff for anecdotal evidence cases is ", cutoff_lbl, ".")
    } else {
      paste0("Proportions of anecdotal evidence, true positives, and false negatives across MABF methods (two-sided test)",
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

#####This version combine all theta=0, 0.2, 0.5 in one graph
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
library(shiny)

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
      br(),
      downloadButton("download_plot", "Download Plot (PNG)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Stacked Plot", plotOutput("piePlot", height = "700px"))
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
    df_all <- rates_all() %>%
      mutate(
        rep.number_num = rep.number,
        rep.n_num      = rep.n,
        rep.number = factor(rep.number,
                            levels = c(2, 5, 10),
                            labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10")),
        rep.n = factor(rep.n,
                       levels = c(40, 100, 400),
                       labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400")),
        method = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF"))
      )
    
    # theta = 0
    df_null <- df_all %>%
      filter(true_es == 0.2) %>%   # keep your original null convention
      select(method, rep.n, rep.number, rep.n_num, rep.number_num, ADNE, TN, FP) %>%
      pivot_longer(cols = c(ADNE, TN, FP), names_to = "metric", values_to = "value") %>%
      mutate(
        theta = factor("theta == 0",
                       levels = c("theta == 0", "theta == 0.2", "theta == 0.5")),
        metric_label = recode(metric,
                              ADNE = "AD",
                              TN   = "TN",
                              FP   = "FP")
      )
    
    # theta = 0.2
    df_02 <- df_all %>%
      filter(true_es == 0.2) %>%
      select(method, rep.n, rep.number, rep.n_num, rep.number_num, ADTE, TP, FN) %>%
      pivot_longer(cols = c(ADTE, TP, FN), names_to = "metric", values_to = "value") %>%
      mutate(
        theta = factor("theta == 0.2",
                       levels = c("theta == 0", "theta == 0.2", "theta == 0.5")),
        metric_label = recode(metric,
                              ADTE = "AD",
                              TP   = "TP",
                              FN   = "FN")
      )
    
    # theta = 0.5
    df_05 <- df_all %>%
      filter(true_es == 0.5) %>%
      select(method, rep.n, rep.number, rep.n_num, rep.number_num, ADTE, TP, FN) %>%
      pivot_longer(cols = c(ADTE, TP, FN), names_to = "metric", values_to = "value") %>%
      mutate(
        theta = factor("theta == 0.5",
                       levels = c("theta == 0", "theta == 0.2", "theta == 0.5")),
        metric_label = recode(metric,
                              ADTE = "AD",
                              TP   = "TP",
                              FN   = "FN")
      )
    
    bind_rows(df_null, df_02, df_05) %>%
      mutate(
        metric_label = factor(metric_label, levels = c("AD", "TN", "FP", "TP", "FN")),
        design = factor(
          paste0("atop(", as.character(rep.n), ", ", as.character(rep.number), ")"),
          levels = c(
            "atop(n[rep] == 40, N[rep] == 2)",
            "atop(n[rep] == 40, N[rep] == 5)",
            "atop(n[rep] == 40, N[rep] == 10)",
            "atop(n[rep] == 100, N[rep] == 2)",
            "atop(n[rep] == 100, N[rep] == 5)",
            "atop(n[rep] == 100, N[rep] == 10)",
            "atop(n[rep] == 400, N[rep] == 2)",
            "atop(n[rep] == 400, N[rep] == 5)",
            "atop(n[rep] == 400, N[rep] == 10)"
          )
        )
      ) %>%
      group_by(theta, design, method, metric_label) %>%
      summarise(value = sum(value), .groups = "drop") %>%
      group_by(theta, design, method) %>%
      mutate(proportion = value / sum(value)) %>%
      ungroup()
  })
  
  # --- Plot object for reuse ---
  plot_obj <- reactive({
    plot_data <- filtered_data()
    
    ggplot(plot_data, aes(x = method, y = proportion, fill = metric_label)) +
      geom_bar(stat = "identity", position = "stack") +
      facet_grid(
        theta ~ design,
        labeller = labeller(
          theta = label_parsed,
          design = label_parsed
        )
      ) +
      scale_fill_manual(
        values = c(
          "AD" = "grey70",
          "FP" = "darkslategray4",
          "TN" = "cadetblue3",
          "TP" = "wheat3",
          "FN" = "lightsalmon4"
        ),
        name = "Classification",
        labels = c(
          "AD" = "Anecdotal Evidence",
          "FP" = "False Positive",
          "TN" = "True Negative",
          "TP" = "True Positive",
          "FN" = "False Negative"
        )
      ) +
      labs(x = "MABF Method", y = "Proportion") +
      theme_gray(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        strip.text = element_text(size = 11),
        panel.border = element_rect(color = "grey90", fill = NA, linewidth = 0.8),
        strip.background = element_rect(fill = "grey90", color = "grey70", linewidth = 0.6),
        plot.background = element_rect(color = "grey80", fill = "white", linewidth = 0.8)
      )
  })
  
  output$piePlot <- renderPlot({
    plot_obj()
  })
  
  # --- Download plot ---
  output$download_plot <- downloadHandler(
    filename = function() {
      cutoff_lbl <- ifelse(input$bf_cutoff == "c1", "3", "10")
      paste0("StackedPlot_CombinedTheta_1sided_BF", cutoff_lbl, ".png")
    },
    content = function(file) {
      ggsave(file, plot = plot_obj(), device = "png", width = 12, height = 11, dpi = 300)
    }
  )
}

# === Run ===
shinyApp(ui = ui, server = server)


#####This version collapse over rep.n and rep.number
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
library(shiny)

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
    folder <- paste0(
      "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=",
      pcutoff, ", EScutoff_o=0/", cutoff_folder, "/"
    )
    
    load_group(folder)
    
    method_list <- list(
      BFbMA  = get(paste0("rates_BFbMA_", true_es, "null_", cutoff_folder)),
      EUBF   = get(paste0("rates_EUBF_", true_es, "null_", cutoff_folder)),
      FEMABF = get(paste0("rates_FEMABF_", true_es, "null_", cutoff_folder)),
      iBF    = get(paste0("rates_iBF_", true_es, "null_", cutoff_folder))
    )
    
    bind_rows(lapply(names(method_list), function(method) {
      method_list[[method]] %>%
        mutate(method = method) %>%
        relocate(method, .before = everything())
    })) %>%
      mutate(
        p_cutoff = as.numeric(pcutoff),
        true_es  = as.numeric(true_es)
      ) %>%
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
      selectInput(
        "bf_cutoff", "Bayes Factor Cutoff",
        choices = c("3" = "c1", "10" = "c2"),
        selected = "c1"
      ),
      br(),
      downloadButton("download_plot", "Download Plot (PNG)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Stacked Plot", plotOutput("piePlot", height = "650px"))
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
    df_all <- rates_all() %>%
      mutate(
        method = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF"))
      )
    
    # theta = 0
    df_null <- df_all %>%
      filter(true_es == 0.2) %>%   # keep your original null convention
      select(method, ADNE, TN, FP) %>%
      pivot_longer(
        cols = c(ADNE, TN, FP),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        theta = factor(
          "theta == 0",
          levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
        ),
        metric_label = recode(
          metric,
          ADNE = "AD",
          TN   = "TN",
          FP   = "FP"
        )
      )
    
    # theta = 0.2
    df_02 <- df_all %>%
      filter(true_es == 0.2) %>%
      select(method, ADTE, TP, FN) %>%
      pivot_longer(
        cols = c(ADTE, TP, FN),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        theta = factor(
          "theta == 0.2",
          levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
        ),
        metric_label = recode(
          metric,
          ADTE = "AD",
          TP   = "TP",
          FN   = "FN"
        )
      )
    
    # theta = 0.5
    df_05 <- df_all %>%
      filter(true_es == 0.5) %>%
      select(method, ADTE, TP, FN) %>%
      pivot_longer(
        cols = c(ADTE, TP, FN),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        theta = factor(
          "theta == 0.5",
          levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
        ),
        metric_label = recode(
          metric,
          ADTE = "AD",
          TP   = "TP",
          FN   = "FN"
        )
      )
    
    bind_rows(df_null, df_02, df_05) %>%
      mutate(
        metric_label = factor(metric_label, levels = c("AD", "TN", "FP", "TP", "FN"))
      ) %>%
      group_by(theta, method, metric_label) %>%
      summarise(value = sum(value), .groups = "drop") %>%
      group_by(theta, method) %>%
      mutate(proportion = value / sum(value)) %>%
      ungroup()
  })
  
  plot_obj <- reactive({
    plot_data <- filtered_data()
    
    # ggplot stack order is reverse of factor levels
    df_labels <- plot_data %>%
      group_by(theta, method) %>%
      arrange(desc(as.numeric(metric_label)), .by_group = TRUE) %>%
      mutate(
        ymax = cumsum(proportion),
        ymin = ymax - proportion,
        ymid = (ymin + ymax) / 2,
        label = percent(proportion, accuracy = 0.01)
      ) %>%
      ungroup() %>%
      filter(proportion >= 0.01)
    
    ggplot(plot_data, aes(x = method, y = proportion, fill = metric_label)) +
      geom_bar(
        stat = "identity",
        position = "stack",
        width = 0.7,
        color = "black",
        linewidth = 0.3
      ) +
      geom_text(
        data = df_labels,
        aes(x = method, y = ymid, label = label),
        inherit.aes = FALSE,
        color = "black",
        size = 3
      ) +
      facet_wrap(
        ~ theta,
        nrow = 1,
        labeller = labeller(theta = label_parsed)
      ) +
      scale_fill_manual(
        values = c(
          "AD" = "grey",
          "FP" = "#8ec1da",
          "TN" = "#cde1ec",
          "TP" = "#f6d6c2",
          "FN" = "salmon2"
        ),
        name = "Classification",
        labels = c(
          "AD" = "Anecdotal evidence",
          "FP" = "False positive",
          "TN" = "True negative",
          "TP" = "True positive",
          "FN" = "False negative"
        )
      ) +
      scale_y_continuous(
        limits = c(0, 1.02),
        expand = c(0, 0)
      ) +
      labs(
        x = NULL,
        y = "Proportion"
      ) +
      theme_gray(base_size = 14) +
      theme(
        panel.background = element_rect(fill = "grey92", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        
        panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        
        strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.8),
        strip.text = element_text(face = "bold", size = 12),
        
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 15),
        
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        
        plot.margin = margin(10, 10, 10, 10)
      )
  })
  
  output$piePlot <- renderPlot({
    plot_obj()
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      cutoff_lbl <- ifelse(input$bf_cutoff == "c1", "3", "10")
      paste0("collapsed_stacked_bar_BFcutoff_", cutoff_lbl, ".png")
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = plot_obj(),
        width = 12,
        height = 10,
        dpi = 600
      )
    }
  )
}

shinyApp(ui = ui, server = server)

######################
#This version include random effects meta-analysis
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
library(shiny)

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
    folder <- paste0(
      "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=",
      pcutoff, ", EScutoff_o=0/", cutoff_folder, "/"
    )
    
    load_group(folder)
    
    method_list <- list(
      BFbMA  = get(paste0("rates_BFbMA_", true_es, "null_", cutoff_folder)),
      EUBF   = get(paste0("rates_EUBF_", true_es, "null_", cutoff_folder)),
      FEMABF = get(paste0("rates_FEMABF_", true_es, "null_", cutoff_folder)),
      iBF    = get(paste0("rates_iBF_", true_es, "null_", cutoff_folder)),
      FEMA   = get(paste0("rates_FEMA_", true_es, "null_", cutoff_folder))
    )
    
    bind_rows(lapply(names(method_list), function(method) {
      method_list[[method]] %>%
        mutate(method = method) %>%
        relocate(method, .before = everything())
    })) %>%
      mutate(
        p_cutoff = as.numeric(pcutoff),
        true_es  = as.numeric(true_es)
      ) %>%
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
      selectInput(
        "bf_cutoff", "Bayes Factor Cutoff",
        choices = c("3" = "c1", "10" = "c2"),
        selected = "c1"
      ),
      br(),
      downloadButton("download_plot", "Download Plot (PNG)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Stacked Plot", plotOutput("piePlot", height = "650px"))
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
    df_all <- rates_all() %>%
      mutate(
        method = factor(
          method,
          levels = c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA"),
          labels = c("BFbMA", "EUBF", "FEMABF", "iBF", "REMA")
        )
      )
    
    # theta = 0
    # FEMA has no ADNE, so assign AD = 0 for FEMA
    df_null <- df_all %>%
      filter(true_es == 0.2) %>%   # keep your original null convention
      mutate(
        ADNE_plot = ifelse(method == "FEMA", 0, ADNE)
      ) %>%
      select(method, ADNE_plot, TN, FP) %>%
      pivot_longer(
        cols = c(ADNE_plot, TN, FP),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        theta = factor(
          "theta == 0",
          levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
        ),
        metric_label = recode(
          metric,
          ADNE_plot = "AD",
          TN        = "TN",
          FP        = "FP"
        )
      )
    
    # theta = 0.2
    # FEMA has no ADTE, so assign AD = 0 for FEMA
    df_02 <- df_all %>%
      filter(true_es == 0.2) %>%
      mutate(
        ADTE_plot = ifelse(method == "FEMA", 0, ADTE)
      ) %>%
      select(method, ADTE_plot, TP, FN) %>%
      pivot_longer(
        cols = c(ADTE_plot, TP, FN),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        theta = factor(
          "theta == 0.2",
          levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
        ),
        metric_label = recode(
          metric,
          ADTE_plot = "AD",
          TP        = "TP",
          FN        = "FN"
        )
      )
    
    # theta = 0.5
    df_05 <- df_all %>%
      filter(true_es == 0.5) %>%
      mutate(
        ADTE_plot = ifelse(method == "FEMA", 0, ADTE)
      ) %>%
      select(method, ADTE_plot, TP, FN) %>%
      pivot_longer(
        cols = c(ADTE_plot, TP, FN),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        theta = factor(
          "theta == 0.5",
          levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
        ),
        metric_label = recode(
          metric,
          ADTE_plot = "AD",
          TP        = "TP",
          FN        = "FN"
        )
      )
    
    bind_rows(df_null, df_02, df_05) %>%
      mutate(
        metric_label = factor(metric_label, levels = c("AD", "TN", "FP", "TP", "FN"))
      ) %>%
      group_by(theta, method, metric_label) %>%
      summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
      group_by(theta, method) %>%
      mutate(proportion = value / sum(value, na.rm = TRUE)) %>%
      ungroup()
  })
  
  plot_obj <- reactive({
    plot_data <- filtered_data()
    
    # ggplot stack order is reverse of factor levels
    df_labels <- plot_data %>%
      group_by(theta, method) %>%
      arrange(desc(as.numeric(metric_label)), .by_group = TRUE) %>%
      mutate(
        ymax = cumsum(proportion),
        ymin = ymax - proportion,
        ymid = (ymin + ymax) / 2,
        label = percent(proportion, accuracy = 0.01)
      ) %>%
      ungroup() %>%
      filter(proportion >= 0.01)
    
    ggplot(plot_data, aes(x = method, y = proportion, fill = metric_label)) +
      geom_bar(
        stat = "identity",
        position = "stack",
        width = 0.7,
        color = "black",
        linewidth = 0.3
      ) +
      geom_text(
        data = df_labels,
        aes(x = method, y = ymid, label = label),
        inherit.aes = FALSE,
        color = "black",
        size = 3
      ) +
      facet_wrap(
        ~ theta,
        nrow = 1,
        labeller = labeller(theta = label_parsed)
      ) +
      scale_fill_manual(
        values = c(
          "AD" = "grey",
          "FP" = "#8ec1da",
          "TN" = "#cde1ec",
          "TP" = "#f6d6c2",
          "FN" = "salmon2"
        ),
        name = "Classification",
        labels = c(
          "AD" = "Anecdotal evidence",
          "FP" = "False positive",
          "TN" = "True negative",
          "TP" = "True positive",
          "FN" = "False negative"
        )
      ) +
      scale_y_continuous(
        limits = c(0, 1.02),
        expand = c(0, 0)
      ) +
      labs(
        x = NULL,
        y = "Proportion"
      ) +
      theme_gray(base_size = 14) +
      theme(
        panel.background = element_rect(fill = "grey92", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.8),
        strip.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 15),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        plot.margin = margin(10, 10, 10, 10)
      )
  })
  
  output$piePlot <- renderPlot({
    plot_obj()
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      cutoff_lbl <- ifelse(input$bf_cutoff == "c1", "3", "10")
      paste0("classic_stackedbar_BFcutoff", cutoff_lbl, "_2sided", ".png")
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = plot_obj(),
        width = 12,
        height = 10,
        dpi = 600
      )
    }
  )
}

shinyApp(ui = ui, server = server)



######################
# This version can choose which method(s) to show
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
library(shiny)

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
    folder <- paste0(
      "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=",
      pcutoff, ", EScutoff_o=0/", cutoff_folder, "/"
    )
    
    load_group(folder)
    
    method_list <- list(
      BFbMA  = get(paste0("rates_BFbMA_", true_es, "null_", cutoff_folder)),
      EUBF   = get(paste0("rates_EUBF_", true_es, "null_", cutoff_folder)),
      FEMABF = get(paste0("rates_FEMABF_", true_es, "null_", cutoff_folder)),
      iBF    = get(paste0("rates_iBF_", true_es, "null_", cutoff_folder)),
      FEMA   = get(paste0("rates_FEMA_", true_es, "null_", cutoff_folder))
    )
    
    bind_rows(lapply(names(method_list), function(method) {
      method_list[[method]] %>%
        mutate(method = method) %>%
        relocate(method, .before = everything())
    })) %>%
      mutate(
        p_cutoff = as.numeric(pcutoff),
        true_es  = as.numeric(true_es)
      ) %>%
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
      selectInput(
        "bf_cutoff", "Bayes Factor Cutoff",
        choices = c("3" = "c1", "10" = "c2"),
        selected = "c1"
      ),
      
      checkboxGroupInput(
        "methods_to_show", "Methods to show",
        choices = c(
          "BFbMA" = "BFbMA",
          "EUBF" = "EUBF",
          "FEMABF" = "FEMABF",
          "iBF" = "iBF",
          "REMA" = "REMA"
        ),
        selected = c("BFbMA", "EUBF", "FEMABF", "iBF", "REMA")
      ),
      
      br(),
      downloadButton("download_plot", "Download Plot (PNG)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Stacked Plot", plotOutput("piePlot", height = "650px"))
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
    req(input$methods_to_show)
    
    df_all <- rates_all() %>%
      mutate(
        method = factor(
          method,
          levels = c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA"),
          labels = c("BFbMA", "EUBF", "FEMABF", "iBF", "REMA")
        )
      ) %>%
      filter(as.character(method) %in% input$methods_to_show)
    
    # theta = 0
    # FEMA has no ADNE, so assign AD = 0 for FEMA
    df_null <- df_all %>%
      filter(true_es == 0.2) %>%   # keep your original null convention
      mutate(
        ADNE_plot = ifelse(as.character(method) == "REMA", 0, ADNE)
      ) %>%
      select(method, ADNE_plot, TN, FP) %>%
      pivot_longer(
        cols = c(ADNE_plot, TN, FP),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        theta = factor(
          "theta == 0",
          levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
        ),
        metric_label = recode(
          metric,
          ADNE_plot = "AD",
          TN        = "TN",
          FP        = "FP"
        )
      )
    
    # theta = 0.2
    # FEMA has no ADTE, so assign AD = 0 for FEMA
    df_02 <- df_all %>%
      filter(true_es == 0.2) %>%
      mutate(
        ADTE_plot = ifelse(as.character(method) == "REMA", 0, ADTE)
      ) %>%
      select(method, ADTE_plot, TP, FN) %>%
      pivot_longer(
        cols = c(ADTE_plot, TP, FN),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        theta = factor(
          "theta == 0.2",
          levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
        ),
        metric_label = recode(
          metric,
          ADTE_plot = "AD",
          TP        = "TP",
          FN        = "FN"
        )
      )
    
    # theta = 0.5
    df_05 <- df_all %>%
      filter(true_es == 0.5) %>%
      mutate(
        ADTE_plot = ifelse(as.character(method) == "REMA", 0, ADTE)
      ) %>%
      select(method, ADTE_plot, TP, FN) %>%
      pivot_longer(
        cols = c(ADTE_plot, TP, FN),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        theta = factor(
          "theta == 0.5",
          levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
        ),
        metric_label = recode(
          metric,
          ADTE_plot = "AD",
          TP        = "TP",
          FN        = "FN"
        )
      )
    
    bind_rows(df_null, df_02, df_05) %>%
      mutate(
        method = factor(as.character(method), levels = input$methods_to_show),
        metric_label = factor(metric_label, levels = c("AD", "TN", "FP", "TP", "FN"))
      ) %>%
      group_by(theta, method, metric_label) %>%
      summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
      group_by(theta, method) %>%
      mutate(proportion = value / sum(value, na.rm = TRUE)) %>%
      ungroup()
  })
  
  plot_obj <- reactive({
    plot_data <- filtered_data()
    
    req(nrow(plot_data) > 0)
    
    # ggplot stack order is reverse of factor levels
    df_labels <- plot_data %>%
      group_by(theta, method) %>%
      arrange(desc(as.numeric(metric_label)), .by_group = TRUE) %>%
      mutate(
        ymax = cumsum(proportion),
        ymin = ymax - proportion,
        ymid = (ymin + ymax) / 2,
        label = percent(proportion, accuracy = 0.01)
      ) %>%
      ungroup() %>%
      filter(proportion >= 0.01)
    
    ggplot(plot_data, aes(x = method, y = proportion, fill = metric_label)) +
      geom_bar(
        stat = "identity",
        position = "stack",
        width = 0.7,
        color = "black",
        linewidth = 0.3
      ) +
      geom_text(
        data = df_labels,
        aes(x = method, y = ymid, label = label),
        inherit.aes = FALSE,
        color = "black",
        size = 3
      ) +
      facet_wrap(
        ~ theta,
        nrow = 1,
        labeller = labeller(theta = label_parsed)
      ) +
      scale_fill_manual(
        values = c(
          "AD" = "grey",
          "FP" = "#8ec1da",
          "TN" = "#cde1ec",
          "TP" = "#f6d6c2",
          "FN" = "salmon2"
        ),
        name = "Classification",
        labels = c(
          "AD" = "Anecdotal evidence",
          "FP" = "False positive",
          "TN" = "True negative",
          "TP" = "True positive",
          "FN" = "False negative"
        )
      ) +
      scale_y_continuous(
        limits = c(0, 1.02),
        expand = c(0, 0)
      ) +
      labs(
        x = NULL,
        y = "Proportion"
      ) +
      theme_gray(base_size = 14) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.8),
        strip.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 15),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        plot.margin = margin(10, 10, 10, 10)
      )
  })
  
  output$piePlot <- renderPlot({
    plot_obj()
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      cutoff_lbl <- ifelse(input$bf_cutoff == "c1", "3", "10")
      methods_lbl <- paste(input$methods_to_show, collapse = "_")
      paste0("classic_stackedbar_BFcutoff", cutoff_lbl, "_", methods_lbl, "_2sided.png")
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = plot_obj(),
        width = 12,
        height = 10,
        dpi = 900
      )
    }
  )
}

shinyApp(ui = ui, server = server)


#################
#In this version, I added FEMA2 (fixed-effects model MA)
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
library(shiny)

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
    folder <- paste0(
      "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=",
      pcutoff, ", EScutoff_o=0/", cutoff_folder, "/"
    )
    
    load_group(folder)
    
    method_list <- list(
      BFbMA  = get(paste0("rates_BFbMA_", true_es, "null_", cutoff_folder)),
      EUBF   = get(paste0("rates_EUBF_", true_es, "null_", cutoff_folder)),
      FEMABF = get(paste0("rates_FEMABF_", true_es, "null_", cutoff_folder)),
      iBF    = get(paste0("rates_iBF_", true_es, "null_", cutoff_folder)),
      FEMA   = get(paste0("rates_FEMA_", true_es, "null_", cutoff_folder)),
      FEMA2  = get(paste0("rates_FEMA2_", true_es, "null_", cutoff_folder))
    )
    
    bind_rows(lapply(names(method_list), function(method) {
      method_list[[method]] %>%
        mutate(method = method) %>%
        relocate(method, .before = everything())
    })) %>%
      mutate(
        p_cutoff = as.numeric(pcutoff),
        true_es  = as.numeric(true_es)
      ) %>%
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
      selectInput(
        "bf_cutoff", "Bayes Factor Cutoff",
        choices = c("3" = "c1", "10" = "c2"),
        selected = "c1"
      ),
      
      checkboxGroupInput(
        "methods_to_show", "Methods to show",
        choices = c(
          "BFbMA" = "BFbMA",
          "EUBF" = "EUBF",
          "FEMABF" = "FEMABF",
          "iBF" = "iBF",
          "REMA" = "REMA",
          "FEMA" = "FEMA"
        ),
        selected = c("BFbMA", "EUBF", "FEMABF", "iBF", "REMA", "FEMA")
      ),
      
      br(),
      downloadButton("download_plot", "Download Plot (PNG)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Stacked Plot", plotOutput("piePlot", height = "650px"))
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
    req(input$methods_to_show)
    
    df_all <- rates_all() %>%
      mutate(
        method = factor(
          method,
          levels = c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA", "FEMA2"),
          labels = c("BFbMA", "EUBF", "FEMABF", "iBF", "REMA", "FEMA")
        )
      ) %>%
      filter(as.character(method) %in% input$methods_to_show)
    
    # theta = 0
    # REMA and FEMA2 have no ADNE, so assign AD = 0
    df_null <- df_all %>%
      filter(true_es == 0.2) %>%   # keep your original null convention
      mutate(
        ADNE_plot = ifelse(as.character(method) %in% c("REMA", "FEMA2"), 0, ADNE)
      ) %>%
      select(method, ADNE_plot, TN, FP) %>%
      pivot_longer(
        cols = c(ADNE_plot, TN, FP),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        theta = factor(
          "theta == 0",
          levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
        ),
        metric_label = recode(
          metric,
          ADNE_plot = "AD",
          TN        = "TN",
          FP        = "FP"
        )
      )
    
    # theta = 0.2
    # REMA and FEMA2 have no ADTE, so assign AD = 0
    df_02 <- df_all %>%
      filter(true_es == 0.2) %>%
      mutate(
        ADTE_plot = ifelse(as.character(method) %in% c("REMA", "FEMA2"), 0, ADTE)
      ) %>%
      select(method, ADTE_plot, TP, FN) %>%
      pivot_longer(
        cols = c(ADTE_plot, TP, FN),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        theta = factor(
          "theta == 0.2",
          levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
        ),
        metric_label = recode(
          metric,
          ADTE_plot = "AD",
          TP        = "TP",
          FN        = "FN"
        )
      )
    
    # theta = 0.5
    # REMA and FEMA2 have no ADTE, so assign AD = 0
    df_05 <- df_all %>%
      filter(true_es == 0.5) %>%
      mutate(
        ADTE_plot = ifelse(as.character(method) %in% c("REMA", "FEMA2"), 0, ADTE)
      ) %>%
      select(method, ADTE_plot, TP, FN) %>%
      pivot_longer(
        cols = c(ADTE_plot, TP, FN),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        theta = factor(
          "theta == 0.5",
          levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
        ),
        metric_label = recode(
          metric,
          ADTE_plot = "AD",
          TP        = "TP",
          FN        = "FN"
        )
      )
    
    bind_rows(df_null, df_02, df_05) %>%
      mutate(
        method = factor(as.character(method), levels = input$methods_to_show),
        metric_label = factor(metric_label, levels = c("AD", "TN", "FP", "TP", "FN"))
      ) %>%
      group_by(theta, method, metric_label) %>%
      summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
      group_by(theta, method) %>%
      mutate(proportion = value / sum(value, na.rm = TRUE)) %>%
      ungroup()
  })
  
  plot_obj <- reactive({
    plot_data <- filtered_data()
    
    req(nrow(plot_data) > 0)
    
    # ggplot stack order is reverse of factor levels
    df_labels <- plot_data %>%
      group_by(theta, method) %>%
      arrange(desc(as.numeric(metric_label)), .by_group = TRUE) %>%
      mutate(
        ymax = cumsum(proportion),
        ymin = ymax - proportion,
        ymid = (ymin + ymax) / 2,
        label = percent(proportion, accuracy = 0.01)
      ) %>%
      ungroup() %>%
      filter(proportion >= 0.01)
    
    ggplot(plot_data, aes(x = method, y = proportion, fill = metric_label)) +
      geom_bar(
        stat = "identity",
        position = "stack",
        width = 0.7,
        color = "black",
        linewidth = 0.3
      ) +
      geom_text(
        data = df_labels,
        aes(x = method, y = ymid, label = label),
        inherit.aes = FALSE,
        color = "black",
        size = 3
      ) +
      facet_wrap(
        ~ theta,
        nrow = 1,
        labeller = labeller(theta = label_parsed)
      ) +
      scale_fill_manual(
        values = c(
          "AD" = "grey",
          "FP" = "#8ec1da",
          "TN" = "#cde1ec",
          "TP" = "#f6d6c2",
          "FN" = "salmon2"
        ),
        name = "Classification",
        labels = c(
          "AD" = "Anecdotal evidence",
          "FP" = "False positive",
          "TN" = "True negative",
          "TP" = "True positive",
          "FN" = "False negative"
        )
      ) +
      scale_y_continuous(
        limits = c(0, 1.02),
        expand = c(0, 0)
      ) +
      labs(
        x = NULL,
        y = "Proportion"
      ) +
      theme_gray(base_size = 14) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.8),
        strip.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 15),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        plot.margin = margin(10, 10, 10, 10)
      )
  })
  
  output$piePlot <- renderPlot({
    plot_obj()
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      cutoff_lbl <- ifelse(input$bf_cutoff == "c1", "3", "10")
      methods_lbl <- paste(input$methods_to_show, collapse = "_")
      paste0("classic_stackedbar_BFcutoff", cutoff_lbl, "_", methods_lbl, "_2sided.png")
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = plot_obj(),
        width = 12,
        height = 10,
        dpi = 900
      )
    }
  )
}

shinyApp(ui = ui, server = server)



##################################
#In this version, I used TS / FF / TF / FS to calculate the TN / FP / TP / FN so Fig 6 and Fig 7 in the paper match
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
library(shiny)

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
    folder <- paste0(
      "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=",
      pcutoff, ", EScutoff_o=0/", cutoff_folder, "/"
    )
    
    load_group(folder)
    
    method_list <- list(
      BFbMA  = get(paste0("rates_BFbMA_", true_es, "null_", cutoff_folder)),
      EUBF   = get(paste0("rates_EUBF_", true_es, "null_", cutoff_folder)),
      FEMABF = get(paste0("rates_FEMABF_", true_es, "null_", cutoff_folder)),
      iBF    = get(paste0("rates_iBF_", true_es, "null_", cutoff_folder)),
      FEMA   = get(paste0("rates_FEMA_", true_es, "null_", cutoff_folder)),
      FEMA2  = get(paste0("rates_FEMA2_", true_es, "null_", cutoff_folder))
    )
    
    bind_rows(lapply(names(method_list), function(method) {
      method_list[[method]] %>%
        mutate(method = method) %>%
        relocate(method, .before = everything())
    })) %>%
      mutate(
        p_cutoff = as.numeric(pcutoff),
        true_es  = as.numeric(true_es)
      ) %>%
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
      selectInput(
        "bf_cutoff", "Bayes Factor Cutoff",
        choices = c("3" = "c1", "10" = "c2"),
        selected = "c1"
      ),
      
      checkboxGroupInput(
        "methods_to_show", "Methods to show",
        choices = c(
          "BFbMA" = "BFbMA",
          "EUBF" = "EUBF",
          "FEMABF" = "FEMABF",
          "iBF" = "iBF",
          "REMA" = "REMA",
          "FEMA" = "FEMA"
        ),
        selected = c("BFbMA", "EUBF", "FEMABF", "iBF", "REMA", "FEMA")
      ),
      
      br(),
      downloadButton("download_plot", "Download Plot (PNG)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Stacked Plot", plotOutput("piePlot", height = "650px"))
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
    req(input$methods_to_show)
    
    df_all <- rates_all() %>%
      mutate(
        method = factor(
          method,
          levels = c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA", "FEMA2"),
          labels = c("BFbMA", "EUBF", "FEMABF", "iBF", "REMA", "FEMA")
        )
      ) %>%
      filter(as.character(method) %in% input$methods_to_show)
    
    # theta = 0
    # REMA and FEMA2 have no ADNE, so assign AD = 0
    # theta = 0
    # derive TN and FP from replication-success categories so the denominator matches
    df_null <- df_all %>%
      filter(true_es == 0.2) %>%   # keep your original null convention
      mutate(
        AD_plot = ifelse(as.character(method) %in% c("REMA", "FEMA"), 0, ADNE),
        TN_plot = ifelse(as.character(method) %in% c("REMA", "FEMA"), TN, TS2 + FF2),
        FP_plot = ifelse(as.character(method) %in% c("REMA", "FEMA"), FP, TF2 + FS2)
      ) %>%
      select(method, AD_plot, TN_plot, FP_plot) %>%
      pivot_longer(
        cols = c(AD_plot, TN_plot, FP_plot),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        theta = factor(
          "theta == 0",
          levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
        ),
        metric_label = recode(
          metric,
          AD_plot = "AD",
          TN_plot = "TN",
          FP_plot = "FP"
        )
      )
    
    # theta = 0.2
    # REMA and FEMA2 have no ADTE, so assign AD = 0
    # theta = 0.2
    # derive TP and FN from replication-success categories so the denominator matches
    df_02 <- df_all %>%
      filter(true_es == 0.2) %>%
      mutate(
        AD_plot = ifelse(as.character(method) %in% c("REMA", "FEMA"), 0, ADTE),
        TP_plot = ifelse(as.character(method) %in% c("REMA", "FEMA"), TP, TS1 + FF1),
        FN_plot = ifelse(as.character(method) %in% c("REMA", "FEMA"), FN, TF1 + FS1)
      ) %>%
      select(method, AD_plot, TP_plot, FN_plot) %>%
      pivot_longer(
        cols = c(AD_plot, TP_plot, FN_plot),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        theta = factor(
          "theta == 0.2",
          levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
        ),
        metric_label = recode(
          metric,
          AD_plot = "AD",
          TP_plot = "TP",
          FN_plot = "FN"
        )
      )
    
    # theta = 0.5
    # REMA and FEMA2 have no ADTE, so assign AD = 0
    # theta = 0.5
    # derive TP and FN from replication-success categories so the denominator matches
    df_05 <- df_all %>%
      filter(true_es == 0.5) %>%
      mutate(
        AD_plot = ifelse(as.character(method) %in% c("REMA", "FEMA"), 0, ADTE),
        TP_plot = ifelse(as.character(method) %in% c("REMA", "FEMA"), TP, TS1 + FF1),
        FN_plot = ifelse(as.character(method) %in% c("REMA", "FEMA"), FN, TF1 + FS1)
      ) %>%
      select(method, AD_plot, TP_plot, FN_plot) %>%
      pivot_longer(
        cols = c(AD_plot, TP_plot, FN_plot),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        theta = factor(
          "theta == 0.5",
          levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
        ),
        metric_label = recode(
          metric,
          AD_plot = "AD",
          TP_plot = "TP",
          FN_plot = "FN"
        )
      )
    
    bind_rows(df_null, df_02, df_05) %>%
      mutate(
        method = factor(as.character(method), levels = input$methods_to_show),
        metric_label = factor(metric_label, levels = c("AD", "TN", "FP", "TP", "FN"))
      ) %>%
      group_by(theta, method, metric_label) %>%
      summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
      group_by(theta, method) %>%
      mutate(proportion = value / sum(value, na.rm = TRUE)) %>%
      ungroup()
  })
  
  plot_obj <- reactive({
    plot_data <- filtered_data()
    
    req(nrow(plot_data) > 0)
    
    # ggplot stack order is reverse of factor levels
    df_labels <- plot_data %>%
      group_by(theta, method) %>%
      arrange(desc(as.numeric(metric_label)), .by_group = TRUE) %>%
      mutate(
        ymax = cumsum(proportion),
        ymin = ymax - proportion,
        ymid = (ymin + ymax) / 2,
        label = percent(proportion, accuracy = 0.01)
      ) %>%
      ungroup() %>%
      filter(proportion >= 0.01)
    
    ggplot(plot_data, aes(x = method, y = proportion, fill = metric_label)) +
      geom_bar(
        stat = "identity",
        position = "stack",
        width = 0.7,
        color = "black",
        linewidth = 0.3
      ) +
      geom_text(
        data = df_labels,
        aes(x = method, y = ymid, label = label),
        inherit.aes = FALSE,
        color = "black",
        size = 3
      ) +
      facet_wrap(
        ~ theta,
        nrow = 1,
        labeller = labeller(theta = label_parsed)
      ) +
      scale_fill_manual(
        values = c(
          "AD" = "grey",
          "FP" = "#8ec1da",
          "TN" = "#cde1ec",
          "TP" = "#f6d6c2",
          "FN" = "salmon2"
        ),
        name = "Classification",
        labels = c(
          "AD" = "Anecdotal evidence",
          "FP" = "False positive",
          "TN" = "True negative",
          "TP" = "True positive",
          "FN" = "False negative"
        )
      ) +
      scale_y_continuous(
        limits = c(0, 1.02),
        expand = c(0, 0)
      ) +
      labs(
        x = NULL,
        y = "Proportion"
      ) +
      theme_gray(base_size = 14) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.8),
        strip.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 15),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        plot.margin = margin(10, 10, 10, 10)
      )
  })
  
  output$piePlot <- renderPlot({
    plot_obj()
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      cutoff_lbl <- ifelse(input$bf_cutoff == "c1", "3", "10")
      methods_lbl <- paste(input$methods_to_show, collapse = "_")
      paste0("classic_stackedbar_BFcutoff", cutoff_lbl, "_", methods_lbl, "_1sided.png")
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = plot_obj(),
        width = 12,
        height = 10,
        dpi = 300
      )
    }
  )
}

shinyApp(ui = ui, server = server)