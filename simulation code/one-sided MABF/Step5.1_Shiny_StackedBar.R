# This script contains four Shiny apps using stacked bar plot
# Stacked bar plot showing relationship between AD, TS, FS, FF, TF
# =======================
# Packages
# =======================
library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

# =======================
# Data prep
# =======================
load_and_prepare_data <- function(cutoff_folder = "c1") {
  load_group <- function(folder_path) {
    rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
    lapply(rds_files, function(file_path) {
      dat <- readRDS(file_path)
      nm  <- tools::file_path_sans_ext(basename(file_path))
      assign(nm, dat, envir = .GlobalEnv)
    })
  }
  
  get_data_combined <- function(pcutoff, true_es) {
    folder <- paste0("./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=",
                     pcutoff, ", EScutoff_o=0/", cutoff_folder, "/")
    load_group(folder)
    
    method_list <- list(
      BFbMA  = get(paste0("rates_BFbMA_",  true_es, "null_", cutoff_folder)),
      EUBF   = get(paste0("rates_EUBF_",   true_es, "null_", cutoff_folder)),
      FEMABF = get(paste0("rates_FEMABF_", true_es, "null_", cutoff_folder)),
      iBF    = get(paste0("rates_iBF_",    true_es, "null_", cutoff_folder))
    )
    
    bind_rows(lapply(names(method_list), function(m) {
      method_list[[m]] %>% mutate(method = m) %>% relocate(method, .before = 1)
    })) %>%
      mutate(p_cutoff = as.numeric(pcutoff),
             true_es  = as.numeric(true_es)) %>%
      relocate(true_es, p_cutoff)
  }
  
  bind_rows(
    get_data_combined("0.01", "0.2"),
    get_data_combined("0.01", "0.5"),
    get_data_combined("0.05", "0.2"),
    get_data_combined("0.05", "0.5")
  )
}

# =======================
# Build analysis-ready data
# =======================
build_long <- function(df_raw, dataset_choice, p_cutoff_sel, orig_n_sel, bias_sel) {
  use_TSn     <- if (dataset_choice == "null") "TS2" else "TS1"
  true_es_num <- if (dataset_choice == "null") 0.2 else as.numeric(dataset_choice)
  
  metric_cols <- if (use_TSn == "TS1") {
    c("ADTE", "TS1", "FF1", "TF1", "FS1")
  } else {
    c("ADNE", "TS2", "FF2", "TF2", "FS2")
  }
  
  df <- df_raw %>%
    filter(true_es == true_es_num) %>%
    { if (p_cutoff_sel != "All") filter(., p_cutoff == as.numeric(p_cutoff_sel)) else . } %>%
    { if (orig_n_sel   != "All") filter(., orig.n   == as.numeric(orig_n_sel))   else . } %>%
    { if (bias_sel     != "All") filter(., censorFunc == bias_sel)               else . }
  
  df <- df %>%
    mutate(
      repN_num  = rep.number,
      nrep_num  = rep.n,
      rep.number = factor(rep.number, levels = c(2, 5, 10),
                          labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10")),
      rep.n      = factor(rep.n, levels = c(40, 100, 400),
                          labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400")),
      method     = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF"))
    ) %>%
    select(censorFunc, p_cutoff, orig.n, method, rep.n, rep.number, 
           repN_num, nrep_num, all_of(metric_cols)) %>%
    pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "value") %>%
    mutate(
      metric_label = recode(metric,
                            ADTE = "AD",  TS1 = "TS", FF1 = "FF", TF1 = "TF", FS1 = "FS",
                            ADNE = "AD",  TS2 = "TS", FF2 = "FF", TF2 = "TF", FS2 = "FS"),
      metric_label = factor(metric_label, levels = c("AD","TS","FF","TF","FS"))
    ) %>%
    group_by(censorFunc, p_cutoff, orig.n, method, rep.n, rep.number, 
             repN_num, nrep_num, metric_label) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    group_by(censorFunc, p_cutoff, orig.n, method, rep.n, rep.number, 
             repN_num, nrep_num) %>%
    mutate(proportion = value / sum(value)) %>%
    ungroup()
  
  df
}

# =======================
# UI
# =======================
ui <- fluidPage(
  titlePanel("MABF Stacked Bar Plot Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("bf_cutoff", "Bayes Factor Cutoff",
                  choices = c("3" = "c1", "10" = "c2"), selected = "c1"),
      selectInput("dataset_choice", "Data Scenario",
                  choices = c("Null" = "null", "0.2", "0.5"), selected = "null"),
      selectInput("p_cutoff", "Original p-value Cutoff",
                  choices = c("All", "0.01", "0.05"), selected = "All"),
      selectInput("orig_n", "Original Sample Size",
                  choices = c("All", "20", "50", "200"), selected = "All"),
      selectInput("bias_level", "Bias Level",
                  choices = c("All", "low", "medium", "high"), selected = "All"),
      
      checkboxGroupInput("facet_vars", "Facet By (choose one or more)",
                         choices = c("Bias Level" = "censorFunc",
                                     "Original p-value Cutoff" = "p_cutoff",
                                     "Original Sample Size" = "orig.n"),
                         selected = c("censorFunc")),
      
      checkboxGroupInput("filter_vars", "Filter by:",
                         choices = c("Replication Sample Size" = "rep.n",
                                     "Number of Replications" = "rep.number"),
                         selected = c("rep.n", "rep.number"))
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Stacked Plot", plotOutput("piePlot", height = "900px")),
        tabPanel("Summary Table", DTOutput("summaryTable"))
      )
    )
  )
)

# =======================
# Server
# =======================
server <- function(input, output, session) {
  
  rates_all <- reactive({
    load_and_prepare_data(input$bf_cutoff)
  })
  
  long_data <- reactive({
    build_long(rates_all(), input$dataset_choice,
               input$p_cutoff, input$orig_n, input$bias_level) %>%
      mutate(
        p_cutoff = factor(p_cutoff,
                          levels = c(0.01, 0.05),
                          labels = c("alpha == 0.01", "alpha == 0.05")),
        orig.n   = factor(orig.n,
                          levels = c(20, 50, 200),
                          labels = c("n[orig] == 20", "n[orig] == 50", "n[orig] == 200")),
        censorFunc = factor(censorFunc,
                            levels = c("low", "medium", "high"),
                            labels = c("Bias~level == 'low'",
                                       "Bias~level == 'medium'",
                                       "Bias~level == 'high'"))
      )
  })
  
  output$piePlot <- renderPlot({
    df <- long_data()
    group_vars <- c("method", "metric_label", input$filter_vars, input$facet_vars)
    df_panel <- df %>%
      group_by(across(any_of(group_vars))) %>%
      summarise(proportion = mean(proportion), .groups = "drop")
    
    ggplot(df_panel, aes(x = method, y = proportion, fill = metric_label)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = c(
        "AD" = "#E69F00", "TS" = "#56B4E9", "FF" = "#F0E442",
        "TF" = "#0072B2", "FS" = "#009E73"
      )) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      labs(x = "MABF Method", y = "Proportion", fill = "Evidence Category") +
      theme_gray(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text  = element_text(size = 12),
            legend.position = "bottom") +
      { if (length(input$facet_vars) == 0) 
        facet_grid(. ~ rep.n + rep.number, labeller = labeller(
          rep.number = label_parsed, rep.n = label_parsed))
        else {
          lhs <- paste(input$facet_vars, collapse = " + ")
          facet_formula <- as.formula(paste(lhs, "~ rep.n + rep.number"))
          facet_grid(facet_formula, labeller = labeller(
            rep.number = label_parsed, rep.n = label_parsed,
            p_cutoff = label_parsed, orig.n = label_parsed, censorFunc = label_parsed))
        }
      }
  })
  
  
  output$summaryTable <- DT::renderDT({
    df <- long_data()
    
    # Only facet by what the user actually checked
    facet_vars <- intersect(input$facet_vars, c("censorFunc", "p_cutoff", "orig.n"))
    
    # Clean labels -> plain values, and set factor order for sorting
    df_clean <- df %>%
      mutate(
        rep.n      = as.integer(gsub("n\\[rep\\] == ", "", as.character(rep.n))),
        rep.number = as.integer(gsub("N\\[rep\\] == ", "", as.character(rep.number))),
        p_cutoff   = gsub("alpha == ", "", as.character(p_cutoff)),
        orig.n     = as.integer(gsub("n\\[orig\\] == ", "", as.character(orig.n))),
        censorFunc = dplyr::recode(as.character(censorFunc),
                                   "Bias~level == 'low'"    = "low",
                                   "Bias~level == 'medium'" = "medium",
                                   "Bias~level == 'high'"   = "high"),
        
        rep.n      = factor(rep.n,      levels = c(40, 100, 400), ordered = TRUE),
        rep.number = factor(rep.number, levels = c(2, 5, 10),     ordered = TRUE),
        p_cutoff   = factor(p_cutoff,   levels = c("0.01","0.05"), ordered = TRUE),
        orig.n     = factor(orig.n,     levels = c(20, 50, 200),   ordered = TRUE),
        censorFunc = factor(censorFunc, levels = c("low","medium","high"), ordered = TRUE)
      )
    
    # Build grouping dynamically
    group_vars <- c("rep.n", "rep.number", facet_vars, "method", "metric_label")
    
    tab <- df_clean %>%
      group_by(across(all_of(group_vars))) %>%
      summarise(Proportion = round(mean(proportion) * 100, 2), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = metric_label, values_from = Proportion)
    
    # Sort rows
    order_vars <- c("rep.n", "rep.number", facet_vars, "method")
    tab <- tab %>% arrange(across(all_of(order_vars)))
    
    # Put Method first, then the other grouping variables, then metrics
    front_cols <- c("method", setdiff(order_vars, "method"))
    front_cols <- intersect(front_cols, names(tab))
    tab <- tab %>% dplyr::relocate(all_of(front_cols))
    
    # Rename columns
    rename_map <- c(
      rep.n       = "Sample size",
      rep.number  = "Number of replications",
      censorFunc  = "Bias mechanism",
      p_cutoff    = "Original alpha",
      orig.n      = "Original sample size",
      method      = "Method"
    )
    names(tab) <- ifelse(names(tab) %in% names(rename_map),
                         rename_map[names(tab)], names(tab))
    
    DT::datatable(
      tab,
      options = list(pageLength = 20, autoWidth = TRUE)
    )
  })
  
  
  
}


# =======================
# Run
# =======================
shinyApp(ui, server)


##########################################################
# Stacked bar plot showing relationship between AD, TP, FP
library(ggplot2)
library(dplyr)
library(tidyr)
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
                  choices = c("Null" = "null", "0.2", "0.5"), selected = "null")
    ),
    mainPanel(
      plotOutput("piePlot", height = "800px")
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
                                     ADTE = "anecdotal evidence", TP = "true positive", FN = "false negative",
                                     ADNE = "anecdotal evidence", TN = "true negative", FP = "false positive"),
        metric_label = factor(metric_label, levels = c("anecdotal evidence", "true positive", "false negative", "true negative", "false positive"))
      ) %>%
      group_by(method, rep.n, rep.number, metric_label) %>%
      summarise(value = sum(value), .groups = "drop") %>%
      group_by(method, rep.n, rep.number) %>%
      mutate(proportion = value / sum(value))
  })
  
  output$piePlot <- renderPlot({
    plot_data <- filtered_data()
    
    ggplot(plot_data, aes(x = method, y = proportion, fill = metric_label)) +
      geom_bar(stat = "identity", position = "stack") +
      facet_grid(rep.n ~ rep.number, 
                 labeller = labeller(
                   rep.number = label_parsed,
                   rep.n = label_parsed
                 )) +
      scale_fill_manual(
        values = c(
          "anecdotal evidence" = "red",
          "false positive" = "green",
          "true negative" = "blue",
          "true positive" = "blue",
          "false negative" = "green"
        ),
        name = "Metric",
        labels = c(
          "anecdotal evidence" = "Anecdotal Evidence",
          "false positive" = "False Positive",
          "true negative" = "True Negative",
          "true positive" = "True Positive",
          "false negative" = "False Negative"
        )
      ) +
      labs(
        x = "MABF Method",
        y = "Proportion",
        title = paste0("Proportions of anecdotal evidence cases, ",
                       if (input$dataset_choice == "null") "false positives, and true negatives"
                       else "true positives and false negatives",
                       "\nwhen the cutoff for anecdotal evidence cases is ",
                       ifelse(input$bf_cutoff == "c1", "3", "10"))
      ) +
      theme_gray(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 12),
        legend.position = "bottom"
      )
  })
  
}

# === Run the App ===
shinyApp(ui = ui, server = server)


############################################
#Stack bar plot comparing FEMA and MABF methods (anecdotal evidence removed)
# =======================
# Packages
# =======================
library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

# =======================
# Data prep helper
# =======================
load_and_prepare_data <- function(cutoff_folder = "c1") {
  load_group <- function(folder_path) {
    rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
    lapply(rds_files, function(file_path) {
      dat <- readRDS(file_path)
      nm  <- tools::file_path_sans_ext(basename(file_path))
      assign(nm, dat, envir = .GlobalEnv)
    })
  }
  
  get_data_combined <- function(pcutoff, true_es) {
    folder <- paste0("./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=",
                     pcutoff, ", EScutoff_o=0/", cutoff_folder, "/")
    load_group(folder)
    
    method_list <- list(
      FEMA   = get(paste0("rates_FEMA_",   true_es, "null_", cutoff_folder)),
      BFbMA  = get(paste0("rates_BFbMA_",  true_es, "null_", cutoff_folder)),
      EUBF   = get(paste0("rates_EUBF_",   true_es, "null_", cutoff_folder)),
      FEMABF = get(paste0("rates_FEMABF_", true_es, "null_", cutoff_folder)),
      iBF    = get(paste0("rates_iBF_",    true_es, "null_", cutoff_folder))
    )
    
    bind_rows(lapply(names(method_list), function(m) {
      method_list[[m]] %>% mutate(method = m) %>% relocate(method, .before = 1)
    })) %>%
      mutate(p_cutoff = as.numeric(pcutoff),
             true_es  = as.numeric(true_es)) %>%
      relocate(true_es, p_cutoff)
  }
  
  bind_rows(
    get_data_combined("0.01", "0.2"),
    get_data_combined("0.01", "0.5"),
    get_data_combined("0.05", "0.2"),
    get_data_combined("0.05", "0.5")
  )
}

# =======================
# Build analysis-ready data (no AD, FEMA included)
# =======================
build_long <- function(df_raw, dataset_choice, p_cutoff_sel, orig_n_sel, bias_sel) {
  use_TSn     <- if (dataset_choice == "null") "TS2" else "TS1"
  true_es_num <- if (dataset_choice == "null") 0.2 else as.numeric(dataset_choice)
  
  # No AD included
  metric_cols <- if (use_TSn == "TS1") {
    c("TS1", "FF1", "TF1", "FS1")
  } else {
    c("TS2", "FF2", "TF2", "FS2")
  }
  
  df <- df_raw %>%
    filter(true_es == true_es_num) %>%
    { if (p_cutoff_sel != "All") filter(., p_cutoff == as.numeric(p_cutoff_sel)) else . } %>%
    { if (orig_n_sel   != "All") filter(., orig.n   == as.numeric(orig_n_sel))   else . } %>%
    { if (bias_sel     != "All") filter(., censorFunc == bias_sel)               else . }
  
  df <- df %>%
    mutate(
      repN_num  = rep.number,
      nrep_num  = rep.n,
      rep.number = factor(rep.number, levels = c(2, 5, 10),
                          labels = c("Nrep = 2", "Nrep = 5", "Nrep = 10")),
      rep.n      = factor(rep.n, levels = c(40, 100, 400),
                          labels = c("nrep = 40", "nrep = 100", "nrep = 400")),
      method     = factor(method, levels = c("FEMA", "BFbMA", "EUBF", "FEMABF", "iBF"))
    ) %>%
    select(censorFunc, p_cutoff, orig.n, method, rep.n, rep.number, 
           repN_num, nrep_num, all_of(metric_cols)) %>%
    pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "value") %>%
    mutate(
      metric_label = recode(metric,
                            TS1 = "TS", FF1 = "FF", TF1 = "TF", FS1 = "FS",
                            TS2 = "TS", FF2 = "FF", TF2 = "TF", FS2 = "FS"),
      metric_label = factor(metric_label, levels = c("TS","FF","TF","FS"))
    ) %>%
    group_by(censorFunc, p_cutoff, orig.n, method, rep.n, rep.number, 
             repN_num, nrep_num, metric_label) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    group_by(censorFunc, p_cutoff, orig.n, method, rep.n, rep.number, 
             repN_num, nrep_num) %>%
    mutate(proportion = value / sum(value)) %>%
    ungroup()
  
  df
}

# =======================
# UI
# =======================
ui <- fluidPage(
  titlePanel("MABF + FEMA Stacked Bar Plot Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("bf_cutoff", "Bayes Factor Cutoff",
                  choices = c("3" = "c1", "10" = "c2"), selected = "c1"),
      selectInput("dataset_choice", "Data Scenario",
                  choices = c("Null" = "null", "0.2", "0.5"), selected = "null"),
      selectInput("p_cutoff", "Original p-value Cutoff",
                  choices = c("All", "0.01", "0.05"), selected = "All"),
      selectInput("orig_n", "Original Sample Size",
                  choices = c("All", "20", "50", "200"), selected = "All"),
      selectInput("bias_level", "Bias Level",
                  choices = c("All", "low", "medium", "high"), selected = "All"),
      
      checkboxGroupInput("facet_vars", "Facet By (choose one or more)",
                         choices = c("Bias Level" = "censorFunc",
                                     "Original p-value Cutoff" = "p_cutoff",
                                     "Original Sample Size" = "orig.n"),
                         selected = c("censorFunc")),
      
      checkboxGroupInput("filter_vars", "Filter by:",
                         choices = c("Replication Sample Size" = "rep.n",
                                     "Number of Replications" = "rep.number"),
                         selected = c("rep.n", "rep.number"))
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Stacked Plot", plotOutput("barPlot", height = "900px")),
        tabPanel("Summary Table", DTOutput("summaryTable"))
      )
    )
  )
)

# =======================
# Server
# =======================
server <- function(input, output, session) {
  
  rates_all <- reactive({
    load_and_prepare_data(input$bf_cutoff)
  })
  
  long_data <- reactive({
    build_long(rates_all(), input$dataset_choice,
               input$p_cutoff, input$orig_n, input$bias_level) %>%
      mutate(
        p_cutoff = factor(p_cutoff,
                          levels = c(0.01, 0.05),
                          labels = c("alpha = 0.01", "alpha = 0.05")),
        orig.n   = factor(orig.n,
                          levels = c(20, 50, 200),
                          labels = c("n_orig = 20", "n_orig = 50", "n_orig = 200")),
        censorFunc = factor(censorFunc,
                            levels = c("low", "medium", "high"),
                            labels = c("bias = low",
                                       "bias = medium",
                                       "bias = high"))
      )
  })
  
  output$barPlot <- renderPlot({
    df <- long_data()
    group_vars <- c("method", "metric_label", input$filter_vars, input$facet_vars)
    df_panel <- df %>%
      group_by(across(any_of(group_vars))) %>%
      summarise(proportion = mean(proportion), .groups = "drop")
    
    ggplot(df_panel, aes(x = method, y = proportion, fill = metric_label)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = c(
        "TS" = "#56B4E9", "FF" = "#F0E442",
        "TF" = "#0072B2", "FS" = "#009E73"
      )) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      labs(x = "Method", y = "Proportion", fill = "Evidence Category") +
      theme_gray(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text  = element_text(size = 12),
            legend.position = "bottom") +
      { if (length(input$facet_vars) == 0) 
        facet_grid(. ~ rep.number + rep.n)
        else {
          lhs <- paste(input$facet_vars, collapse = " + ")
          facet_formula <- as.formula(paste(lhs, "~ rep.number + rep.n"))
          facet_grid(facet_formula)
        }
      }
  })
  
  output$summaryTable <- DT::renderDT({
    df <- long_data()
    
    # Facet vars in the exact order the user selected
    facet_vars <- input$facet_vars
    facet_vars <- facet_vars[facet_vars %in% c("censorFunc", "p_cutoff", "orig.n")]
    
    # Clean labels to plain values; handle both math-style and plain-style labels
    df_clean <- df %>%
      mutate(
        # rep labels can be "N[rep] == 2"/"n[rep] == 40" or "Nrep = 2"/"nrep = 40"
        rep.n      = as.integer(gsub("^(n\\[rep\\] == |nrep = )", "", as.character(rep.n))),
        rep.number = as.integer(gsub("^(N\\[rep\\] == |Nrep = )", "", as.character(rep.number))),
        # alpha labels can be "alpha == 0.01" or "alpha = 0.01"
        p_cutoff   = gsub("^alpha( ==| =) ", "", as.character(p_cutoff)),
        # orig labels can be "n[orig] == 20" or "n_orig = 20"
        orig.n     = as.integer(gsub("^(n\\[orig\\] == |n_orig = )", "", as.character(orig.n))),
        # bias labels can be "Bias~level == 'low'" or "bias = low"
        censorFunc = dplyr::recode(as.character(censorFunc),
                                   "Bias~level == 'low'" = "low",
                                   "Bias~level == 'medium'" = "medium",
                                   "Bias~level == 'high'" = "high",
                                   "bias = low" = "low",
                                   "bias = medium" = "medium",
                                   "bias = high" = "high",
                                   .default = as.character(censorFunc))
      ) %>%
      mutate(
        rep.n      = factor(rep.n,      levels = c(40, 100, 400), ordered = TRUE),
        rep.number = factor(rep.number, levels = c(2, 5, 10),     ordered = TRUE),
        p_cutoff   = factor(p_cutoff,   levels = c("0.01","0.05"), ordered = TRUE),
        orig.n     = factor(orig.n,     levels = c(20, 50, 200),   ordered = TRUE),
        censorFunc = factor(censorFunc, levels = c("low","medium","high"), ordered = TRUE)
      )
    
    # Group by rep sizes + selected facets + method + metric
    group_vars <- c("rep.n", "rep.number", facet_vars, "method", "metric_label")
    
    tab <- df_clean %>%
      group_by(across(all_of(group_vars))) %>%
      summarise(
        Count = sum(value),
        Proportion = round(mean(proportion) * 100, 2),
        .groups = "drop"
      ) %>%
      tidyr::pivot_wider(
        names_from = metric_label,
        values_from = c(Proportion, Count),
        names_glue = "{metric_label}_{.value}"
      ) %>%
      mutate(
        Total = coalesce(TS_Count, 0) + coalesce(FF_Count, 0) +
          coalesce(TF_Count, 0) + coalesce(FS_Count, 0)
      )
    
    # Sort rows by rep sizes, selected facets, then method
    order_vars <- c("rep.n", "rep.number", facet_vars, "method")
    tab <- tab %>% arrange(across(all_of(order_vars)))
    
    # Column order: Method, rep sizes, selected facets (in selected order), metric pairs, Total
    front_cols <- c(
      "method", "rep.n", "rep.number", facet_vars,
      "TS_Proportion","TS_Count",
      "FF_Proportion","FF_Count",
      "TF_Proportion","TF_Count",
      "FS_Proportion","FS_Count",
      "Total"
    )
    front_cols <- intersect(front_cols, names(tab))
    tab <- tab %>% dplyr::select(all_of(front_cols))
    
    # Friendly headers
    hdr <- names(tab)
    hdr[hdr == "method"]      <- "Method"
    hdr[hdr == "rep.n"]       <- "Sample size"
    hdr[hdr == "rep.number"]  <- "Number of replications"
    hdr[hdr == "censorFunc"]  <- "Bias mechanism"
    hdr[hdr == "p_cutoff"]    <- "Original alpha"
    hdr[hdr == "orig.n"]      <- "Original sample size"
    hdr <- gsub("TS_Proportion", "TS (%)", hdr)
    hdr <- gsub("FF_Proportion", "FF (%)", hdr)
    hdr <- gsub("TF_Proportion", "TF (%)", hdr)
    hdr <- gsub("FS_Proportion", "FS (%)", hdr)
    hdr <- gsub("TS_Count", "TS", hdr)
    hdr <- gsub("FF_Count", "FF", hdr)
    hdr <- gsub("TF_Count", "TF", hdr)
    hdr <- gsub("FS_Count", "FS", hdr)
    names(tab) <- hdr
    
    DT::datatable(
      tab,
      options = list(pageLength = 20, autoWidth = TRUE)
    )
  })
  
  
  
}

# =======================
# Run
# =======================
shinyApp(ui, server)


##########################################################
# Stacked bar plot showing relationship between AD, medium evidence, strong evidence
######
# plan(multisession, workers = parallel::detectCores() - 1)
# 
# ## Import all datasets
# # Set the path to the directory containing RDS files
# folder_path <- "./MABFanalyses/row-wise/"
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
# ## Combine all datasets into a list
# all_names <- ls(pattern = "^MABF_rates_.*_rowise_BFcutoff[0-9]+$")
# all_dfs <- mget(all_names)
# 
# ## Define the summarizing function
# summarize_sim_chunk <- function(df) {
#   df %>%
#     filter(method != "iBF2") %>% 
#     mutate(group_id = (row_number() - 1) %/% 500 + 1) %>%
#     group_by(group_id) %>%
#     summarise(
#       # Keep simulation conditions (e.g., columns 1 to 12)
#       across(1:10, first),
#       # Use column positions 13 to 17 for metrics
#       across(13:17, ~mean(.x, na.rm = TRUE), .names = "{.col}"),
#       across(23:24, ~mean(.x, na.rm = TRUE), .names = "{.col}"),
#       across(29:32, ~mean(.x, na.rm = TRUE), .names = "{.col}"),
#       .groups = "drop"
#     )
# }
# # Process the data in parallel
# summary_list <- future_map(all_dfs, summarize_sim_chunk)
# 
# 
# # Assign the datasets in summary_list to the global environment:
# list2env(summary_list, envir = .GlobalEnv)
# 
# # Get the names of the summary datasets
# summary_names <- names(summary_list)
# 
# # Remove all other objects
# rm(list = setdiff(ls(), summary_names))
# 
# # Stop parallel mode
# plan(sequential)
# 
# # Save the workspace
# save.image(file = "./MABFanalyses/row-wise/row2matrix.RData")
######
# Load the workspace (data converted from row-wise datasets into matrix-wise datasets)
# The imported datasets include 12 datasets, with 3 categories (null, 0.2, and 0.5), within each category, 4 BF cut offs (1,3,10,30).
# But actually, the datasets within the same category are identical -_- |||, so we can pick anyone of the datasets within a category.
#load(file = "./MABFanalyses/row-wise/row2matrix.RData")
######
# ##Underlying effect is null
# # Prepare the data
# plot_data_null <- MABF_rates_null_rowise_BFcutoff3 %>%
#   mutate(orig.alpha = factor (orig.alpha, levels = c(0.01, 0.05))) %>% 
#   mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
#   mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
#   mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
#   mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>%
#   select(method, rep.n, rep.number, ADR_null_anecdotal, ADR_null_moderate, ADR_null_strong, ADR_null_vstrong, ADR_null_estrong) %>%
#   gather(key = "metric", value = "value", ADR_null_anecdotal, ADR_null_moderate, ADR_null_strong, ADR_null_vstrong, ADR_null_estrong) %>%
#   mutate(metric = factor(metric, levels = c("ADR_null_anecdotal", "ADR_null_moderate", "ADR_null_strong", "ADR_null_vstrong", "ADR_null_estrong"))) %>%
#   mutate(metric_label = factor(metric, 
#                                levels = c("ADR_null_anecdotal", "ADR_null_moderate", "ADR_null_strong", "ADR_null_vstrong", "ADR_null_estrong"),
#                                labels = c("anecdotal", "moderate", "strong", "very strong", "extremely strong"))) %>%
#   group_by(method, rep.n, rep.number) %>%
#   mutate(proportion = value / sum(value))
# 
# # Create custom gradient colors for the AD_null metrics
# ADR_null_colors <- c("#B0E0E6","#87CEEB", "#1E90FF",  "#0000CD","#00008B" )
# 
# # Create a named vector for scale_fill_manual
# color_mapping <- c(anecdotal = ADR_null_colors[1],
#                    moderate = ADR_null_colors[2],
#                    strong = ADR_null_colors[3],
#                    "very strong" = ADR_null_colors[4],
#                    "extremely strong" = ADR_null_colors[5])
# 
# # Custom labeller for metrics
# metric_labeller <- c(
#   ADR_null_anecdotal = "anecdotal",
#   ADR_null_moderate = "moderate",
#   ADR_null_strong = "strong",
#   ADR_null_vstrong = "very strong",
#   ADR_null_estrong = "extremely strong"
# )
# 
# # Create the plot
# plot <- ggplot(plot_data_null, aes(x = factor(method), y = proportion, fill = metric_label)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = color_mapping) +  # Custom colors
#   labs(x = "MABF Method", y = "Proportion", fill = "Evidence Strength", title = "Proportion of Bayes Factors Categoried into Varying Levels of Evidence When the Underlying Effect is Null") +
#   facet_grid(rep.n ~ rep.number, scales = "free_x", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")
# print(plot)
# 
# ##Underlying effect is 0.2
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# ##Underlying effect is true
# # Prepare the data
# plot_data_true <- MABF_rates_0.2_rowise_BFcutoff3 %>%
#   mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
#   mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
#   mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
#   mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>%
#   select(method, rep.n, rep.number, ADR_true_anecdotal, ADR_true_moderate, ADR_true_strong, ADR_true_vstrong, ADR_true_estrong) %>%
#   gather(key = "metric", value = "value", ADR_true_anecdotal, ADR_true_moderate, ADR_true_strong, ADR_true_vstrong, ADR_true_estrong) %>%
#   mutate(metric = factor(metric, levels = c("ADR_true_anecdotal", "ADR_true_moderate", "ADR_true_strong", "ADR_true_vstrong", "ADR_true_estrong"))) %>%
#   mutate(metric_label = factor(metric, 
#                                levels = c("ADR_true_anecdotal", "ADR_true_moderate", "ADR_true_strong", "ADR_true_vstrong", "ADR_true_estrong"),
#                                labels = c("anecdotal", "moderate", "strong", "very strong", "extremely strong"))) %>%
#   group_by(method, rep.n, rep.number) %>%
#   mutate(proportion = value / sum(value))
# 
# # Create custom gradient colors for the AD_true metrics
# ADR_true_colors <- c("#B0E0E6","#87CEEB", "#1E90FF",  "#0000CD","#00008B" )
# 
# # Create a named vector for scale_fill_manual
# color_mapping <- c(anecdotal = ADR_true_colors[1],
#                    moderate = ADR_true_colors[2],
#                    strong = ADR_true_colors[3],
#                    "very strong" = ADR_true_colors[4],
#                    "extremely strong" = ADR_true_colors[5])
# 
# # Custom labeller for metrics
# metric_labeller <- c(
#   ADR_true_anecdotal = "anecdotal",
#   ADR_true_moderate = "moderate",
#   ADR_true_strong = "strong",
#   ADR_true_vstrong = "very strong",
#   ADR_true_estrong = "extremely strong"
# )
# 
# # Create the plot
# plot <- ggplot(plot_data_true, aes(x = factor(method), y = proportion, fill = metric_label)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = color_mapping) +  # Custom colors
#   labs(x = "MABF Method", y = "Proportion", fill = "Evidence Strength", title = "Proportion of Bayes Factors Categoried into Varying Levels of Evidence When the Underlying Effect is true") +
#   facet_grid(rep.n ~ rep.number, scales = "free_x", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")
# 
# print(plot)
# 
# ##Underlying effect is 0.5
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# ##Underlying effect is true
# # Prepare the data
# plot_data_true <- MABF_rates_0.5_rowise_BFcutoff3 %>%
#   mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
#   mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
#   mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
#   mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>%
#   select(method, rep.n, rep.number, ADR_true_anecdotal, ADR_true_moderate, ADR_true_strong, ADR_true_vstrong, ADR_true_estrong) %>%
#   gather(key = "metric", value = "value", ADR_true_anecdotal, ADR_true_moderate, ADR_true_strong, ADR_true_vstrong, ADR_true_estrong) %>%
#   mutate(metric = factor(metric, levels = c("ADR_true_anecdotal", "ADR_true_moderate", "ADR_true_strong", "ADR_true_vstrong", "ADR_true_estrong"))) %>%
#   mutate(metric_label = factor(metric, 
#                                levels = c("ADR_true_anecdotal", "ADR_true_moderate", "ADR_true_strong", "ADR_true_vstrong", "ADR_true_estrong"),
#                                labels = c("anecdotal", "moderate", "strong", "very strong", "extremely strong"))) %>%
#   group_by(method, rep.n, rep.number) %>%
#   mutate(proportion = value / sum(value))
# 
# # Create custom gradient colors for the AD_true metrics
# ADR_true_colors <- c("#B0E0E6","#87CEEB", "#1E90FF",  "#0000CD","#00008B" )
# 
# # Create a named vector for scale_fill_manual
# color_mapping <- c(anecdotal = ADR_true_colors[1],
#                    moderate = ADR_true_colors[2],
#                    strong = ADR_true_colors[3],
#                    "very strong" = ADR_true_colors[4],
#                    "extremely strong" = ADR_true_colors[5])
# 
# # Custom labeller for metrics
# metric_labeller <- c(
#   ADR_true_anecdotal = "anecdotal",
#   ADR_true_moderate = "moderate",
#   ADR_true_strong = "strong",
#   ADR_true_vstrong = "very strong",
#   ADR_true_estrong = "extremely strong"
# )
# 
# # Create the plot
# plot <- ggplot(plot_data_true, aes(x = factor(method), y = proportion, fill = metric_label)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = color_mapping) +  # Custom colors
#   labs(x = "MABF Method", y = "Proportion", fill = "Evidence Strength", title = "Proportion of Bayes Factors Categoried into Varying Levels of Evidence When the Underlying Effect is true") +
#   facet_grid(rep.n ~ rep.number, scales = "free_x", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")
# 
# print(plot)
######
# Stacked bar plot showing relationship between levels of evidence strength
library(ggplot2)
library(dplyr)
library(tidyr)

ui <- fluidPage(
  titlePanel("MABF Evidence Strength Proportions"),
  sidebarLayout(
    sidebarPanel(
      selectInput("effect_size", "Underlying Effect Size:",
                  choices = c("null", "0.2", "0.5"), selected = "null")
    ),
    mainPanel(
      plotOutput("bar_plot", width = "1200px", height = "800px")
    )
  )
)

server <- function(input, output, session) {
  
  # Always load the RData once at startup
  env <- new.env()
  try(load("./MABFanalyses/row-wise/row2matrix.RData", envir = env), silent = TRUE)
  loaded_env <- reactiveVal(env)
  
  output$bar_plot <- renderPlot({
    req(loaded_env())
    
    # Build dataset name from selected effect size
    data_name <- paste0(
      "MABF_rates_",
      ifelse(input$effect_size == "null", "null", input$effect_size),
      "_rowise_BFcutoff3"  # fixed cutoff of 3
    )
    
    # Retrieve dataset from the loaded environment
    dataset <- tryCatch({
      get(data_name, envir = loaded_env())
    }, error = function(e) NULL)
    
    if (is.null(dataset)) return(NULL)
    
    # Choose columns based on effect size
    if (input$effect_size == "null") {
      metric_cols <- c("ADR_null_anecdotal", "ADR_null_moderate", "ADR_null_strong", 
                       "ADR_null_vstrong", "ADR_null_estrong")
      labels <- c("anecdotal", "moderate", "strong", "very strong", "extremely strong")
    } else {
      metric_cols <- c("ADR_true_anecdotal", "ADR_true_moderate", "ADR_true_strong", 
                       "ADR_true_vstrong", "ADR_true_estrong")
      labels <- c("anecdotal", "moderate", "strong", "very strong", "extremely strong")
    }
    
    # Prepare the data
    plot_data <- dataset %>%
      mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
      mutate(rep.number = factor(rep.number, levels = c(2, 5, 10),
                                 labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
      mutate(rep.n = factor(rep.n, levels = c(40, 100, 400),
                            labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
      mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>%
      select(method, rep.n, rep.number, all_of(metric_cols)) %>%
      pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "value") %>%
      mutate(metric_label = factor(metric, levels = metric_cols, labels = labels)) %>%
      group_by(method, rep.n, rep.number) %>%
      mutate(proportion = value / sum(value)) %>%
      ungroup()
    
    # Color palette
    color_mapping <- setNames(
      c("#B0E0E6", "#87CEEB", "#1E90FF", "#0000CD", "#00008B"),
      labels
    )
    
    # Plot
    ggplot(plot_data, aes(x = factor(method), y = proportion, fill = metric_label)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = color_mapping) +
      labs(
        x = "MABF Method", y = "Proportion", fill = "Evidence Strength",
        title = paste("Bayes Factor Evidence Strength by Method\nEffect Size:", input$effect_size)
      ) +
      facet_grid(rep.n ~ rep.number, scales = "free_x",
                 labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")
  })
}

shinyApp(ui = ui, server = server)
