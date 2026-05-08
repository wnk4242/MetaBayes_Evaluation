# This Shiny generate tables for replication success measures (AD, TS, FF, TF, FS) for the MABF methods

# =======================
# Packages
# =======================
library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
library(kableExtra)
library(DT)
library(stringr)
library(purrr)

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
# Build long “analysis-ready” data
# =======================
build_long <- function(df_raw, dataset_choice, p_cutoff, orig_n, bias_level) {
  use_TSn     <- if (dataset_choice == "null") "TS2" else "TS1"
  true_es_num <- if (dataset_choice == "null") 0.2 else as.numeric(dataset_choice)
  
  metric_cols <- if (use_TSn == "TS1") {
    c("ADTE", "TS1", "FF1", "TF1", "FS1")
  } else {
    c("ADNE", "TS2", "FF2", "TF2", "FS2")
  }
  

  cutoff_val <- if (p_cutoff == "All") NULL else as.numeric(p_cutoff)
  orig_val   <- if (orig_n   == "All") NULL else as.numeric(orig_n)
  bias_val   <- if (bias_level == "All") NULL else bias_level
  
  df <- df_raw %>%
    filter(true_es == true_es_num) %>%
    { if (!is.null(cutoff_val)) filter(., p_cutoff == cutoff_val) else . } %>%
    { if (!is.null(orig_val))   filter(., orig.n   == orig_val)   else . } %>%
    { if (!is.null(bias_val))   filter(., censorFunc == bias_val) else . }
  
  df <- df %>%
    mutate(
      repN_num   = rep.number,
      nrep_num   = rep.n,
      censorFunc = factor(censorFunc, levels = c("low", "medium", "high"))
    ) %>%
    mutate(
      rep.number = factor(rep.number, levels = c(2, 5, 10),
                          labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10")),
      rep.n      = factor(rep.n, levels = c(40, 100, 400),
                          labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400")),
      method     = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF"))
    ) %>%
    select(censorFunc, p_cutoff, orig.n, method, rep.n, rep.number, repN_num, nrep_num, all_of(metric_cols)) %>%
    pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "value") %>%
    mutate(
      metric_label = recode(metric,
                            ADTE = "AD", TS1 = "TS", FF1 = "FF", TF1 = "TF", FS1 = "FS",
                            ADNE = "AD", TS2 = "TS", FF2 = "FF", TF2 = "TF", FS2 = "FS"),
      metric_label = factor(metric_label, levels = c("AD","TS","FF","TF","FS"))
    ) %>%
    group_by(censorFunc, p_cutoff, orig.n, method, rep.n, rep.number, repN_num, nrep_num, metric_label) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    group_by(censorFunc, p_cutoff, orig.n, method, rep.n, rep.number, repN_num, nrep_num) %>%
    mutate(proportion = value / sum(value)) %>%
    ungroup()
  
  df
}

# =======================
# Build APA-style table
# =======================
make_apa_table <- function(df_long, facet_var) {
  method_order  <- c("BFbMA","EUBF","FEMABF","iBF")
  metric_order  <- c("AD","TS","FF","TF","FS")
  
  facet_label_prefix <- switch(facet_var,
                               censorFunc = "Bias level",
                               p_cutoff   = "alpha",
                               orig.n     = "n_orig")
  
  base <- df_long %>%
    mutate(percent = round(proportion * 100, 2),
           method = factor(method, levels = method_order),
           metric_label = factor(metric_label, levels = metric_order)) %>%
    select(all_of(facet_var), repN_num, nrep_num, method, metric_label, percent) %>%
    unite(".col", method, metric_label, sep = "_") %>%
    pivot_wider(names_from = ".col", values_from = percent,
                values_fn = mean, values_fill = 0) %>%
    arrange(.data[[facet_var]], repN_num, nrep_num)
  

  method_major_cols <- unlist(lapply(method_order, function(m)
    paste(m, metric_order, sep = "_")
  ))
  
  wanted_cols <- c(facet_var, "repN_num", "nrep_num", method_major_cols)
  for (cn in setdiff(wanted_cols, names(base))) base[[cn]] <- NA_real_
  base <- base[, wanted_cols]
  
  metric_cols <- setdiff(names(base), c(facet_var, "repN_num", "nrep_num"))
  base[metric_cols] <- lapply(base[metric_cols], function(x) {
    ifelse(is.na(x), "", format(round(x, 2), nsmall = 2))
  })
  

  clean_labels <- function(var, val) {
    if (var == "censorFunc") {
      recode(val,
             "Bias~level == 'low'"    = "Bias level: low",
             "Bias~level == 'medium'" = "Bias level: medium",
             "Bias~level == 'high'"   = "Bias level: high")
    } else if (var == "p_cutoff") {
      recode(val,
             "alpha == 0.01" = "$\\alpha = 0.01$",
             "alpha == 0.05" = "$\\alpha = 0.05$")
    } else if (var == "orig.n") {
      recode(val,
             "n[orig] == 20"  = "$n_{orig} = 20$",
             "n[orig] == 50"  = "$n_{orig} = 50$",
             "n[orig] == 200" = "$n_{orig} = 200$")
    } else {
      val
    }
  }
  
  # Split by facet level and add section headers
  pieces <- split(base, base[[facet_var]], drop = TRUE)
  
  out <- purrr::imap_dfr(pieces, function(df_piece, level_val) {
    df_piece <- df_piece %>%
      group_by(repN_num) %>%
      mutate(repN_show = if_else(row_number() == 1, as.character(repN_num), "")) %>%
      ungroup() %>%
      transmute(
        Nrep = repN_show,
        nrep = as.character(nrep_num),
        !!!rlang::syms(metric_cols)
      )
    
    header <- tibble(
      Nrep = paste0("\\textbf{", clean_labels(facet_var, level_val), "}"),
      nrep = "",
      !!!setNames(rep("", length(metric_cols)), metric_cols)
    )
    
    bind_rows(header, df_piece)
  })
  
  out
}

# =======================
# UI
# =======================
ui <- fluidPage(
  titlePanel("MABF Stacked Bar Plot + APA Tables"),
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
      
      # NEW: optional third filter (can leave empty)
      checkboxGroupInput("facet_row", "Optional Third Filter",
                         choices = c("Bias Level" = "censorFunc",
                                     "Original p-value Cutoff" = "p_cutoff",
                                     "Original Sample Size" = "orig.n"),
                         selected = NULL),
      
      checkboxGroupInput("filter_vars", "Filter by:",
                         choices = c("Replication Sample Size" = "rep.n",
                                     "Number of Replications" = "rep.number"),
                         selected = c("rep.n", "rep.number")),
      br(),
      downloadButton("downloadPlot", "Download Plot (PNG)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Stacked Plot",
                 plotOutput("piePlot", height = "900px")
        ),
        tabPanel("APA Table",
                 htmlOutput("apa_html"),
                 br(),
                 downloadButton("download_apa_csv", "Download APA Table (CSV)"),
                 br(),
                 downloadButton("download_apa_tex", "Download APA Table (LaTeX)"),
                 br(),
                 verbatimTextOutput("apa_latex"))
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
    df <- build_long(rates_all(), input$dataset_choice, input$p_cutoff, input$orig_n, input$bias_level)
    
    # Recode facet variables into math expressions for label_parsed
    df <- df %>%
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
    df
  })
  
  # -------- Plot --------
  output$piePlot <- renderPlot({
    df <- long_data()
    
    # Group only by selected filters + method + metric_label + optional facet_row
    group_vars <- c("method", "metric_label", input$filter_vars, input$facet_row)
    df_panel <- df %>%
      group_by(across(any_of(group_vars))) %>%
      summarise(proportion = mean(proportion), .groups = "drop")
    
    p <- ggplot(df_panel, aes(x = method, y = proportion, fill = metric_label)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = c(
        "AD" = "#E69F00", "TS" = "#56B4E9", "FF" = "#F0E442",
        "TF" = "#0072B2", "FS" = "#009E73"
      )) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      labs(x = "MABF Method", y = "Proportion", fill = "Evidence Category") +
      theme_gray(base_size = 14) +
      theme(
        axis.text.x     = element_text(angle = 45, hjust = 1),
        strip.text      = element_text(size = 12),
        legend.position = "bottom"
      )
    
    # Add facetting
    if (length(input$facet_row) == 1) {
      p <- p + facet_grid(
        rows = vars(!! rlang::sym(input$facet_row)),
        cols = vars(!!!rlang::syms(input$filter_vars)),
        labeller = labeller(
          rep.number = label_parsed,
          rep.n      = label_parsed,
          p_cutoff   = label_parsed,
          orig.n     = label_parsed,
          censorFunc = label_parsed
        )
      )
    } else {
      p <- p + facet_grid(
        cols = vars(!!!rlang::syms(input$filter_vars)),
        labeller = labeller(
          rep.number = label_parsed,
          rep.n      = label_parsed,
          p_cutoff   = label_parsed,
          orig.n     = label_parsed,
          censorFunc = label_parsed
        )
      )
    }
    
    p
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      # Build filename parts
      bf_cutoff   <- ifelse(input$bf_cutoff == "c1", "BFcutoff3", "BFcutoff10")
      dataset     <- paste0("theta", input$dataset_choice)
      filters     <- paste(c(input$filter_vars, input$facet_row), collapse = "_")
      if (filters == "") filters <- "nofilters"
      
      paste0("StackedPlot_", bf_cutoff, "_", dataset, "_", filters, ".png")
    },
    content = function(file) {
      # Use same plot code as output$piePlot
      df <- long_data()
      group_vars <- c("method", "metric_label", input$filter_vars, input$facet_row)
      df_panel <- df %>%
        group_by(across(any_of(group_vars))) %>%
        summarise(proportion = mean(proportion), .groups = "drop")
      
      p <- ggplot(df_panel, aes(x = method, y = proportion, fill = metric_label)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_manual(values = c(
          "AD" = "#E69F00", "TS" = "#56B4E9", "FF" = "#F0E442",
          "TF" = "#0072B2", "FS" = "#009E73"
        )) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        labs(x = "MABF Method", y = "Proportion", fill = "Evidence Category") +
        theme_gray(base_size = 14) +
        theme(
          axis.text.x     = element_text(angle = 45, hjust = 1),
          strip.text      = element_text(size = 12),
          legend.position = "bottom"
        )
      
      if (length(input$facet_row) == 1) {
        p <- p + facet_grid(
          rows = vars(!! rlang::sym(input$facet_row)),
          cols = vars(!!!rlang::syms(input$filter_vars)),
          labeller = labeller(
            rep.number = label_parsed,
            rep.n      = label_parsed,
            p_cutoff   = label_parsed,
            orig.n     = label_parsed,
            censorFunc = label_parsed
          )
        )
      } else {
        p <- p + facet_grid(
          cols = vars(!!!rlang::syms(input$filter_vars)),
          labeller = labeller(
            rep.number = label_parsed,
            rep.n      = label_parsed,
            p_cutoff   = label_parsed,
            orig.n     = label_parsed,
            censorFunc = label_parsed
          )
        )
      }
      
      ggsave(file, plot = p, width = 12, height = 8, dpi = 300)
    }
  )
  
  # -------- APA Table --------
  apa_tbl <- reactive({
    if (length(input$facet_row) == 1) {
      make_apa_table(long_data(), input$facet_row)
    } else {
      long_data() %>%
        group_by(method, repN_num, nrep_num, metric_label) %>%
        summarise(percent = round(mean(proportion)*100, 2), .groups="drop") %>%
        pivot_wider(names_from = c(method, metric_label), values_from = percent)
    }
  })
  
  output$apa_html <- renderUI({
    tab <- apa_tbl()
    if (length(input$facet_row) == 1) {
      header_methods <- c(" " = 2, "BFbMA" = 5, "EUBF" = 5, "FEMABF" = 5, "iBF" = 5)
      header_rates <- c(" " = 2,
                        rep(c("AD (%)","TS (%)","FF (%)","TF (%)","FS (%)"), times = 4))
      
      kable(tab, format = "html", escape = FALSE,
            align = "lrrrrrrrrrrrrrrrrrrrr",
            col.names = c("Nrep", "nrep", rep("", 20))) %>%
        kable_styling(full_width = FALSE, bootstrap_options = c("condensed", "hover")) %>%
        add_header_above(header_rates, bold = FALSE) %>%
        add_header_above(header_methods, bold = FALSE, line = FALSE) %>%
        HTML()
    } else {
      DT::datatable(tab)
    }
  })
  
  apa_tbl_latex <- reactive({
    if (length(input$facet_row) == 1) {
      tab <- apa_tbl()
      header_methods <- c(" " = 2, "BFbMA" = 5, "EUBF" = 5, "FEMABF" = 5, "iBF" = 5)
      header_rates   <- c("$N_{\\\\text{rep}}$", "$n_{\\\\text{rep}}$",
                          rep(c("AD (\\\\%)","TS (\\\\%)","FF (\\\\%)",
                                "TF (\\\\%)","FS (\\\\%)"), times = 4))
      
      theta_text <- switch(input$dataset_choice,
                           "null" = " for true effect size $\\theta = 0$",
                           "0.2"  = " for true effect size $\\theta = 0.2$",
                           "0.5"  = " for true effect size $\\theta = 0.5$")
      
      facet_nice <- c(
        censorFunc = "bias level",
        p_cutoff   = "original p-value cutoff",
        orig.n     = "original sample size"
      )[input$facet_row]
      
      caption_text <- paste0(
        "Replication success measures (AD, TS, FF, TF, FS) across MABF methods by ",
        "replication sample size ($n_{\\text{rep}}$) and number of replications ($N_{\\text{rep}}$)",
        theta_text, ", stratified by ", facet_nice, "."
      )
      
      kable(tab, format = "latex", booktabs = TRUE, escape = FALSE,
            align = c('c','c', rep('r', 20)),
            col.names = rep("", 22),
            caption = caption_text) %>%
        kable_styling(latex_options = c("hold_position")) %>%
        add_header_above(header_rates,   bold = FALSE, line = FALSE, escape = FALSE) %>%
        add_header_above(header_methods, bold = FALSE,  line = TRUE,  escape = FALSE)
    } else {
      NULL
    }
  })
  
  output$apa_latex <- renderText({
    if (!is.null(apa_tbl_latex())) {
      paste(capture.output(print(apa_tbl_latex())), collapse = "\n")
    } else {
      "APA LaTeX table not available (no third filter selected)."
    }
  })
  
  output$download_apa_tex <- downloadHandler(
    filename = function() paste0("APA_table_", ifelse(length(input$facet_row)==1, input$facet_row, "collapsed"), ".tex"),
    content = function(file) {
      if (!is.null(apa_tbl_latex())) {
        tex_code <- paste(capture.output(print(apa_tbl_latex())), collapse = "\n")
        writeLines(tex_code, file)
      } else {
        writeLines("No APA table available without third filter.", file)
      }
    }
  )
  
  output$download_apa_csv <- downloadHandler(
    filename = function() paste0("APA_table_", ifelse(length(input$facet_row)==1, input$facet_row, "collapsed"), ".csv"),
    content  = function(file) write.csv(apa_tbl(), file, row.names = FALSE)
  )
}

# =======================
# Run
# =======================
shinyApp(ui, server)

###########################################################
##### This version collapses over replication designs and focuses on methodological factor influence
library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# =========================================================
# 1. Load and combine data
# =========================================================
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
    folder <- paste0(
      "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=",
      pcutoff, ", EScutoff_o=0/", cutoff_folder, "/"
    )
    
    load_group(folder)
    
    method_list <- list(
      BFbMA  = get(paste0("rates_BFbMA_",  true_es, "null_", cutoff_folder)),
      EUBF   = get(paste0("rates_EUBF_",   true_es, "null_", cutoff_folder)),
      FEMABF = get(paste0("rates_FEMABF_", true_es, "null_", cutoff_folder)),
      iBF    = get(paste0("rates_iBF_",    true_es, "null_", cutoff_folder)),
      FEMA   = get(paste0("rates_FEMA_",   true_es, "null_", cutoff_folder))
    )
    
    bind_rows(lapply(names(method_list), function(m) {
      method_list[[m]] %>%
        mutate(method = m) %>%
        relocate(method, .before = 1)
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

# =========================================================
# 2. Build long-format proportions
# =========================================================
build_long <- function(df_raw,
                       dataset_choice = "null",
                       p_cutoff = "All",
                       orig_n = "All",
                       bias_level = "All") {
  
  use_TSn     <- if (dataset_choice == "null") "TS2" else "TS1"
  true_es_num <- if (dataset_choice == "null") 0.2 else as.numeric(dataset_choice)
  
  cutoff_val <- if (p_cutoff == "All") NULL else as.numeric(p_cutoff)
  orig_val   <- if (orig_n == "All") NULL else as.numeric(orig_n)
  bias_val   <- if (bias_level == "All") NULL else bias_level
  
  df <- df_raw %>%
    filter(true_es == true_es_num) %>%
    { if (!is.null(cutoff_val)) filter(., p_cutoff == cutoff_val) else . } %>%
    { if (!is.null(orig_val))   filter(., orig.n == orig_val) else . } %>%
    { if (!is.null(bias_val))   filter(., censorFunc == bias_val) else . }
  
  # FEMA/REMA has no anecdotal category, so set AD = 0 for plotting
  if (use_TSn == "TS1") {
    if (!"ADTE" %in% names(df)) df$ADTE <- NA_real_
    df <- df %>%
      mutate(AD_plot = ifelse(method == "FEMA", 0, ADTE))
    metric_cols_plot <- c("AD_plot", "TS1", "FF1", "TF1", "FS1")
  } else {
    if (!"ADNE" %in% names(df)) df$ADNE <- NA_real_
    df <- df %>%
      mutate(AD_plot = ifelse(method == "FEMA", 0, ADNE))
    metric_cols_plot <- c("AD_plot", "TS2", "FF2", "TF2", "FS2")
  }
  
  df_long <- df %>%
    mutate(
      repN_num   = rep.number,
      nrep_num   = rep.n,
      censorFunc = factor(censorFunc, levels = c("low", "medium", "high"))
    ) %>%
    mutate(
      rep.number = factor(
        rep.number,
        levels = c(2, 5, 10),
        labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10")
      ),
      rep.n = factor(
        rep.n,
        levels = c(40, 100, 400),
        labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400")
      ),
      method = factor(
        method,
        levels = c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA"),
        labels = c("BFbMA", "EUBF", "FEMABF", "iBF", "REMA")
      )
    ) %>%
    select(
      true_es, censorFunc, p_cutoff, orig.n, method, rep.n, rep.number,
      repN_num, nrep_num, all_of(metric_cols_plot)
    ) %>%
    pivot_longer(
      cols = all_of(metric_cols_plot),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric_label = recode(
        metric,
        AD_plot = "AD",
        TS1 = "TS", FF1 = "FF", TF1 = "TF", FS1 = "FS",
        TS2 = "TS", FF2 = "FF", TF2 = "TF", FS2 = "FS"
      ),
      metric_label = factor(metric_label, levels = c("AD", "TS", "FF", "TF", "FS"))
    ) %>%
    group_by(
      true_es, censorFunc, p_cutoff, orig.n, method, rep.n, rep.number,
      repN_num, nrep_num, metric_label
    ) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    group_by(
      true_es, censorFunc, p_cutoff, orig.n, method, rep.n, rep.number,
      repN_num, nrep_num
    ) %>%
    mutate(proportion = value / sum(value, na.rm = TRUE)) %>%
    ungroup()
  
  df_long
}

# =========================================================
# 3. Build all three datasets together
# =========================================================
build_long_all_datasets <- function(df_raw,
                                    p_cutoff = "All",
                                    orig_n = "All",
                                    bias_level = "All") {
  
  bind_rows(
    build_long(df_raw, dataset_choice = "null", p_cutoff = p_cutoff, orig_n = orig_n, bias_level = bias_level) %>%
      mutate(dataset_label = "theta == 0"),
    
    build_long(df_raw, dataset_choice = "0.2", p_cutoff = p_cutoff, orig_n = orig_n, bias_level = bias_level) %>%
      mutate(dataset_label = "theta == 0.2"),
    
    build_long(df_raw, dataset_choice = "0.5", p_cutoff = p_cutoff, orig_n = orig_n, bias_level = bias_level) %>%
      mutate(dataset_label = "theta == 0.5")
  ) %>%
    mutate(
      dataset_label = factor(
        dataset_label,
        levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
      )
    )
}

# =========================================================
# 4. Plotting function
# =========================================================
plot_combined_factors_all_datasets <- function(df_long_all,
                                               method_keep = NULL,
                                               rep_n_keep = c(40, 100, 400),
                                               rep_number_keep = c(2, 5, 10)) {
  
  df_plot <- df_long_all %>%
    { if (!is.null(method_keep))      filter(., method %in% method_keep) else . } %>%
    { if (!is.null(rep_n_keep))       filter(., nrep_num %in% rep_n_keep) else . } %>%
    { if (!is.null(rep_number_keep))  filter(., repN_num %in% rep_number_keep) else . } %>%
    mutate(
      orig.n     = factor(orig.n, levels = c(20, 50, 200),
                          labels = c("n == 20", "n == 50", "n == 200")),
      p_cutoff   = factor(p_cutoff, levels = c(0.01, 0.05),
                          labels = c("alpha == .01", "alpha == .05")),
      censorFunc = factor(censorFunc, levels = c("low", "medium", "high"),
                          labels = c("'Low'", "'Medium'", "'High'"))
    ) %>%
    group_by(dataset_label, orig.n, p_cutoff, censorFunc, metric_label) %>%
    summarise(proportion = mean(proportion), .groups = "drop")
  
  df_orig <- df_plot %>%
    transmute(
      dataset_label,
      factor_type = "Original sample size",
      factor_level = as.character(orig.n),
      metric_label,
      proportion
    ) %>%
    group_by(dataset_label, factor_type, factor_level, metric_label) %>%
    summarise(proportion = mean(proportion), .groups = "drop")
  
  df_alpha <- df_plot %>%
    transmute(
      dataset_label,
      factor_type = "Original α level",
      factor_level = as.character(p_cutoff),
      metric_label,
      proportion
    ) %>%
    group_by(dataset_label, factor_type, factor_level, metric_label) %>%
    summarise(proportion = mean(proportion), .groups = "drop")
  
  df_bias <- df_plot %>%
    transmute(
      dataset_label,
      factor_type = "Bias level",
      factor_level = as.character(censorFunc),
      metric_label,
      proportion
    ) %>%
    group_by(dataset_label, factor_type, factor_level, metric_label) %>%
    summarise(proportion = mean(proportion), .groups = "drop")
  
  df_all <- bind_rows(df_orig, df_alpha, df_bias) %>%
    mutate(
      factor_type = factor(
        factor_type,
        levels = c("Original sample size", "Original α level", "Bias level")
      ),
      factor_level = factor(
        factor_level,
        levels = c("n == 20", "n == 50", "n == 200",
                   "alpha == .01", "alpha == .05",
                   "'Low'", "'Medium'", "'High'")
      ),
      metric_label = factor(metric_label, levels = c("AD", "TS", "FF", "TF", "FS"))
    )
  
  df_labels <- df_all %>%
    group_by(dataset_label, factor_type, factor_level) %>%
    arrange(desc(as.numeric(metric_label)), .by_group = TRUE) %>%
    mutate(
      ymax = cumsum(proportion),
      ymin = ymax - proportion,
      ymid = (ymin + ymax) / 2
    ) %>%
    ungroup() %>%
    mutate(
      show_label = proportion >= 0.01
    ) %>%
    filter(show_label) %>%
    mutate(label = percent(proportion, accuracy = 0.01))
  
  ggplot(df_all, aes(x = factor_level, y = proportion, fill = metric_label)) +
    geom_bar(
      stat = "identity",
      position = "stack",
      width = 0.7,
      color = "black",
      linewidth = 0.3
    ) +
    geom_text(
      data = df_labels,
      aes(x = factor_level, y = ymid, label = label),
      inherit.aes = FALSE,
      size = 3,
      color = "black"
    ) +
    facet_grid(
      dataset_label ~ factor_type,
      scales = "free_x",
      space = "free_x",
      labeller = labeller(dataset_label = label_parsed)
    ) +
    scale_x_discrete(
      labels = function(x) parse(text = x)
    ) +
    scale_y_continuous(
      limits = c(0, 1.02),
      breaks = seq(0, 1, 0.2),
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      values = c(
        "AD" = "grey",
        "TS" = "lightskyblue",
        "FF" = "khaki1",
        "TF" = "bisque1",
        "FS" = "cornsilk2"
      ),
      labels = c(
        "AD" = "Anecdotal",
        "TS" = "True success",
        "FF" = "False failure",
        "TF" = "True failure",
        "FS" = "False success"
      )
    ) +
    labs(
      x = NULL,
      y = "Proportion",
      fill = "Classification"
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
}

# =========================================================
# 5. UI
# =========================================================
ui <- fluidPage(
  titlePanel("MABF Stacked Bar Plot Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "bf_cutoff",
        "Bayes factor cutoff",
        choices = c("3" = "c1", "10" = "c2"),
        selected = "c1"
      ),
      
      selectInput(
        "method_keep",
        "Method",
        choices = c("FEMABF", "iBF", "EUBF", "BFbMA", "REMA"),
        selected = "iBF"
      ),
      
      downloadButton("download_plot", "Download Plot (PNG)")
    ),
    
    mainPanel(
      plotOutput("combined_plot", height = "1200px", width = "800px")
    )
  )
)

# =========================================================
# 6. Server
# =========================================================
server <- function(input, output, session) {
  
  df_raw_reactive <- reactive({
    load_and_prepare_data(cutoff_folder = input$bf_cutoff)
  })
  
  df_long_all_reactive <- reactive({
    build_long_all_datasets(
      df_raw = df_raw_reactive(),
      p_cutoff = "All",
      orig_n = "All",
      bias_level = "All"
    )
  })
  
  plot_reactive <- reactive({
    plot_combined_factors_all_datasets(
      df_long_all = df_long_all_reactive(),
      method_keep = input$method_keep,
      rep_n_keep = c(40, 100, 400),
      rep_number_keep = c(2, 5, 10)
    )
  })
  
  output$combined_plot <- renderPlot({
    plot_reactive()
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0(
        "repsuccess_stackedbar_",
        "BFcutoff", ifelse(input$bf_cutoff == "c1", "3", "10"),
        "_method_", input$method_keep,
        "_1sided",
        ".png"
      )
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = plot_reactive(),
        width = 10,
        height = 14,
        dpi = 900
      )
    }
  )
}

# =========================================================
# 7. Run app
# =========================================================
shinyApp(ui = ui, server = server)


##########################################################
##### This version collapses over replication designs and focuses on methodological factor influence
#This version added FEMA2 (fixed-effects MA)
library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# =========================================================
# 1. Load and combine data
# =========================================================
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
    folder <- paste0(
      "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=",
      pcutoff, ", EScutoff_o=0/", cutoff_folder, "/"
    )
    
    load_group(folder)
    
    method_list <- list(
      BFbMA  = get(paste0("rates_BFbMA_",  true_es, "null_", cutoff_folder)),
      EUBF   = get(paste0("rates_EUBF_",   true_es, "null_", cutoff_folder)),
      FEMABF = get(paste0("rates_FEMABF_", true_es, "null_", cutoff_folder)),
      iBF    = get(paste0("rates_iBF_",    true_es, "null_", cutoff_folder)),
      FEMA   = get(paste0("rates_FEMA_",   true_es, "null_", cutoff_folder)),
      FEMA2  = get(paste0("rates_FEMA2_",  true_es, "null_", cutoff_folder))
    )
    
    bind_rows(lapply(names(method_list), function(m) {
      method_list[[m]] %>%
        mutate(method = m) %>%
        relocate(method, .before = 1)
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

# =========================================================
# 2. Build long-format proportions
# =========================================================
build_long <- function(df_raw,
                       dataset_choice = "null",
                       p_cutoff = "All",
                       orig_n = "All",
                       bias_level = "All") {
  
  use_TSn     <- if (dataset_choice == "null") "TS2" else "TS1"
  true_es_num <- if (dataset_choice == "null") 0.2 else as.numeric(dataset_choice)
  
  cutoff_val <- if (p_cutoff == "All") NULL else as.numeric(p_cutoff)
  orig_val   <- if (orig_n == "All") NULL else as.numeric(orig_n)
  bias_val   <- if (bias_level == "All") NULL else bias_level
  
  df <- df_raw %>%
    filter(true_es == true_es_num) %>%
    { if (!is.null(cutoff_val)) filter(., p_cutoff == cutoff_val) else . } %>%
    { if (!is.null(orig_val))   filter(., orig.n == orig_val) else . } %>%
    { if (!is.null(bias_val))   filter(., censorFunc == bias_val) else . }
  
  # FEMA/REMA and FEMA2/FEMA have no anecdotal category, so set AD = 0 for plotting
  if (use_TSn == "TS1") {
    if (!"ADTE" %in% names(df)) df$ADTE <- NA_real_
    df <- df %>%
      mutate(AD_plot = ifelse(method %in% c("FEMA", "FEMA2"), 0, ADTE))
    metric_cols_plot <- c("AD_plot", "TS1", "FF1", "TF1", "FS1")
  } else {
    if (!"ADNE" %in% names(df)) df$ADNE <- NA_real_
    df <- df %>%
      mutate(AD_plot = ifelse(method %in% c("FEMA", "FEMA2"), 0, ADNE))
    metric_cols_plot <- c("AD_plot", "TS2", "FF2", "TF2", "FS2")
  }
  
  df_long <- df %>%
    mutate(
      repN_num   = rep.number,
      nrep_num   = rep.n,
      censorFunc = factor(censorFunc, levels = c("low", "medium", "high"))
    ) %>%
    mutate(
      rep.number = factor(
        rep.number,
        levels = c(2, 5, 10),
        labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10")
      ),
      rep.n = factor(
        rep.n,
        levels = c(40, 100, 400),
        labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400")
      ),
      method = factor(
        method,
        levels = c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA", "FEMA2"),
        labels = c("BFbMA", "EUBF", "FEMABF", "iBF", "REMA", "FEMA")
      )
    ) %>%
    select(
      true_es, censorFunc, p_cutoff, orig.n, method, rep.n, rep.number,
      repN_num, nrep_num, all_of(metric_cols_plot)
    ) %>%
    pivot_longer(
      cols = all_of(metric_cols_plot),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric_label = recode(
        metric,
        AD_plot = "AD",
        TS1 = "TS", FF1 = "FF", TF1 = "TF", FS1 = "FS",
        TS2 = "TS", FF2 = "FF", TF2 = "TF", FS2 = "FS"
      ),
      metric_label = factor(metric_label, levels = c("AD", "TS", "FF", "TF", "FS"))
    ) %>%
    group_by(
      true_es, censorFunc, p_cutoff, orig.n, method, rep.n, rep.number,
      repN_num, nrep_num, metric_label
    ) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    group_by(
      true_es, censorFunc, p_cutoff, orig.n, method, rep.n, rep.number,
      repN_num, nrep_num
    ) %>%
    mutate(proportion = value / sum(value, na.rm = TRUE)) %>%
    ungroup()
  
  df_long
}

# =========================================================
# 3. Build all three datasets together
# =========================================================
build_long_all_datasets <- function(df_raw,
                                    p_cutoff = "All",
                                    orig_n = "All",
                                    bias_level = "All") {
  
  bind_rows(
    build_long(df_raw, dataset_choice = "null", p_cutoff = p_cutoff, orig_n = orig_n, bias_level = bias_level) %>%
      mutate(dataset_label = "theta == 0"),
    
    build_long(df_raw, dataset_choice = "0.2", p_cutoff = p_cutoff, orig_n = orig_n, bias_level = bias_level) %>%
      mutate(dataset_label = "theta == 0.2"),
    
    build_long(df_raw, dataset_choice = "0.5", p_cutoff = p_cutoff, orig_n = orig_n, bias_level = bias_level) %>%
      mutate(dataset_label = "theta == 0.5")
  ) %>%
    mutate(
      dataset_label = factor(
        dataset_label,
        levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
      )
    )
}

# =========================================================
# 4. Plotting function
# =========================================================
plot_combined_factors_all_datasets <- function(df_long_all,
                                               method_keep = NULL,
                                               rep_n_keep = c(40, 100, 400),
                                               rep_number_keep = c(2, 5, 10)) {
  
  df_plot <- df_long_all %>%
    { if (!is.null(method_keep))      filter(., method %in% method_keep) else . } %>%
    { if (!is.null(rep_n_keep))       filter(., nrep_num %in% rep_n_keep) else . } %>%
    { if (!is.null(rep_number_keep))  filter(., repN_num %in% rep_number_keep) else . } %>%
    mutate(
      orig.n     = factor(orig.n, levels = c(20, 50, 200),
                          labels = c("n == 20", "n == 50", "n == 200")),
      p_cutoff   = factor(p_cutoff, levels = c(0.01, 0.05),
                          labels = c("alpha == .01", "alpha == .05")),
      censorFunc = factor(censorFunc, levels = c("low", "medium", "high"),
                          labels = c("'Low'", "'Medium'", "'High'"))
    ) %>%
    group_by(dataset_label, orig.n, p_cutoff, censorFunc, metric_label) %>%
    summarise(proportion = mean(proportion), .groups = "drop")
  
  df_orig <- df_plot %>%
    transmute(
      dataset_label,
      factor_type = "Original sample size",
      factor_level = as.character(orig.n),
      metric_label,
      proportion
    ) %>%
    group_by(dataset_label, factor_type, factor_level, metric_label) %>%
    summarise(proportion = mean(proportion), .groups = "drop")
  
  df_alpha <- df_plot %>%
    transmute(
      dataset_label,
      factor_type = "Original α level",
      factor_level = as.character(p_cutoff),
      metric_label,
      proportion
    ) %>%
    group_by(dataset_label, factor_type, factor_level, metric_label) %>%
    summarise(proportion = mean(proportion), .groups = "drop")
  
  df_bias <- df_plot %>%
    transmute(
      dataset_label,
      factor_type = "Bias level",
      factor_level = as.character(censorFunc),
      metric_label,
      proportion
    ) %>%
    group_by(dataset_label, factor_type, factor_level, metric_label) %>%
    summarise(proportion = mean(proportion), .groups = "drop")
  
  df_all <- bind_rows(df_orig, df_alpha, df_bias) %>%
    mutate(
      factor_type = factor(
        factor_type,
        levels = c("Original sample size", "Original α level", "Bias level")
      ),
      factor_level = factor(
        factor_level,
        levels = c("n == 20", "n == 50", "n == 200",
                   "alpha == .01", "alpha == .05",
                   "'Low'", "'Medium'", "'High'")
      ),
      metric_label = factor(metric_label, levels = c("AD", "TS", "FF", "TF", "FS"))
    )
  
  df_labels <- df_all %>%
    group_by(dataset_label, factor_type, factor_level) %>%
    arrange(desc(as.numeric(metric_label)), .by_group = TRUE) %>%
    mutate(
      ymax = cumsum(proportion),
      ymin = ymax - proportion,
      ymid = (ymin + ymax) / 2
    ) %>%
    ungroup() %>%
    mutate(
      show_label = proportion >= 0.01
    ) %>%
    filter(show_label) %>%
    mutate(label = percent(proportion, accuracy = 0.01))
  
  ggplot(df_all, aes(x = factor_level, y = proportion, fill = metric_label)) +
    geom_bar(
      stat = "identity",
      position = "stack",
      width = 0.7,
      color = "black",
      linewidth = 0.3
    ) +
    geom_text(
      data = df_labels,
      aes(x = factor_level, y = ymid, label = label),
      inherit.aes = FALSE,
      size = 3,
      color = "black"
    ) +
    facet_grid(
      dataset_label ~ factor_type,
      scales = "free_x",
      space = "free_x",
      labeller = labeller(dataset_label = label_parsed)
    ) +
    scale_x_discrete(
      labels = function(x) parse(text = x)
    ) +
    scale_y_continuous(
      limits = c(0, 1.02),
      breaks = seq(0, 1, 0.2),
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      values = c(
        "AD" = "grey",
        "TS" = "lightskyblue",
        "FF" = "khaki1",
        "TF" = "bisque1",
        "FS" = "cornsilk2"
      ),
      labels = c(
        "AD" = "Anecdotal",
        "TS" = "True success",
        "FF" = "False failure",
        "TF" = "True failure",
        "FS" = "False success"
      )
    ) +
    labs(
      x = NULL,
      y = "Proportion",
      fill = "Classification"
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
}

# =========================================================
# 5. UI
# =========================================================
ui <- fluidPage(
  titlePanel("MABF Stacked Bar Plot Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "bf_cutoff",
        "Bayes factor cutoff",
        choices = c("3" = "c1", "10" = "c2"),
        selected = "c1"
      ),
      
      selectInput(
        "method_keep",
        "Method",
        choices = c(
          "FEMABF" = "FEMABF",
          "iBF"    = "iBF",
          "EUBF"   = "EUBF",
          "BFbMA"  = "BFbMA",
          "REMA"   = "REMA",
          "FEMA"   = "FEMA"
        ),
        selected = "iBF"
      ),
      
      downloadButton("download_plot", "Download Plot (PNG)")
    ),
    
    mainPanel(
      plotOutput("combined_plot", height = "1200px", width = "800px")
    )
  )
)

# =========================================================
# 6. Server
# =========================================================
server <- function(input, output, session) {
  
  df_raw_reactive <- reactive({
    load_and_prepare_data(cutoff_folder = input$bf_cutoff)
  })
  
  df_long_all_reactive <- reactive({
    build_long_all_datasets(
      df_raw = df_raw_reactive(),
      p_cutoff = "All",
      orig_n = "All",
      bias_level = "All"
    )
  })
  
  plot_reactive <- reactive({
    plot_combined_factors_all_datasets(
      df_long_all = df_long_all_reactive(),
      method_keep = input$method_keep,
      rep_n_keep = c(40, 100, 400),
      rep_number_keep = c(2, 5, 10)
    )
  })
  
  output$combined_plot <- renderPlot({
    plot_reactive()
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0(
        "repsuccess_stackedbar_",
        "BFcutoff", ifelse(input$bf_cutoff == "c1", "3", "10"),
        "_method_", input$method_keep,
        "_1sided",
        ".png"
      )
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = plot_reactive(),
        width = 10,
        height = 14,
        dpi = 300
      )
    }
  )
}

# =========================================================
# 7. Run app
# =========================================================
shinyApp(ui = ui, server = server)

###############################################################
# This Shiny generates tables for replication success measures (AD, TS, FF, TF, FS)
# for Bayes-factor methods plus meta-analysis methods (FEMA -> REMA, FEMA2 -> FEMA)
# This version is based on the Shiny app simulation outcome stackedbar (D:\wnk\PhD\PhD research\side projects\Shiny App 4-color replication stackedbar)
# =======================
# Packages
# =======================
library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
library(kableExtra)
library(DT)
library(stringr)
library(purrr)

# =======================
# Plot settings
# =======================
class_colors <- c(
  "AD" = "grey",
  "TS" = "lightskyblue",
  "FF" = "khaki1",
  "TF" = "bisque1",
  "FS" = "cornsilk2"
)

class_labels <- c(
  "AD" = "Anecdotal evidence",
  "TS" = "True success",
  "FF" = "False failure",
  "TF" = "True failure",
  "FS" = "False success"
)

plot_labeller <- labeller(
  rep.number = label_parsed,
  rep.n      = label_parsed,
  p_cutoff   = label_parsed,
  orig.n     = label_parsed,
  censorFunc = label_parsed,
  theta_label = label_parsed
)

# =======================
# Data prep helper
# =======================
load_and_prepare_data <- function(cutoff_folder = "c1") {
  
  get_one_method_data <- function(folder, method, true_es) {
    rds_files <- list.files(folder, pattern = "\\.RDS$", full.names = TRUE)
    file_base <- tools::file_path_sans_ext(basename(rds_files))
    
    # Avoid FEMA matching FEMA2
    if (method == "FEMA2") {
      hit_idx <- which(
        grepl("FEMA2", file_base, ignore.case = FALSE) &
          grepl(paste0(true_es, "null"), file_base, fixed = TRUE)
      )
    } else if (method == "FEMA") {
      hit_idx <- which(
        grepl("FEMA", file_base, ignore.case = FALSE) &
          !grepl("FEMA2", file_base, ignore.case = FALSE) &
          grepl(paste0(true_es, "null"), file_base, fixed = TRUE)
      )
    } else {
      hit_idx <- which(
        grepl(method, file_base, ignore.case = FALSE) &
          grepl(paste0(true_es, "null"), file_base, fixed = TRUE)
      )
    }
    
    if (length(hit_idx) == 0) {
      stop(
        paste0(
          "Could not find an RDS file for method '", method,
          "' and true_es '", true_es, "null' in folder:\n", folder,
          "\n\nAvailable files were:\n",
          paste(basename(rds_files), collapse = "\n")
        )
      )
    }
    
    if (length(hit_idx) > 1) {
      warning(
        paste0(
          "Multiple files matched method '", method,
          "' and true_es '", true_es,
          "null'. Using the first match:\n",
          basename(rds_files[hit_idx[1]])
        )
      )
    }
    
    readRDS(rds_files[hit_idx[1]])
  }
  
  get_data_combined <- function(pcutoff, true_es) {
    folder <- paste0(
      "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=",
      pcutoff, ", EScutoff_o=0/", cutoff_folder, "/"
    )
    
    method_list <- list(
      BFbMA  = get_one_method_data(folder, "BFbMA",  true_es),
      EUBF   = get_one_method_data(folder, "EUBF",   true_es),
      FEMABF = get_one_method_data(folder, "FEMABF", true_es),
      iBF    = get_one_method_data(folder, "iBF",    true_es),
      FEMA   = get_one_method_data(folder, "FEMA",   true_es),
      FEMA2  = get_one_method_data(folder, "FEMA2",  true_es)
    )
    
    bind_rows(lapply(names(method_list), function(m) {
      method_list[[m]] %>%
        mutate(method = m) %>%
        relocate(method, .before = 1)
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

# =======================
# Build long “analysis-ready” data
# =======================
build_long <- function(df_raw, dataset_choice, p_cutoff, orig_n, bias_level) {
  
  cutoff_val <- if (p_cutoff == "All") NULL else as.numeric(p_cutoff)
  orig_val   <- if (orig_n   == "All") NULL else as.numeric(orig_n)
  bias_val   <- if (bias_level == "All") NULL else bias_level
  
  mabf_methods <- c("BFbMA", "EUBF", "FEMABF", "iBF")
  meta_methods <- c("FEMA", "FEMA2")
  
  df_filtered <- df_raw %>%
    { if (!is.null(cutoff_val)) filter(., p_cutoff == cutoff_val) else . } %>%
    { if (!is.null(orig_val))   filter(., orig.n   == orig_val)   else . } %>%
    { if (!is.null(bias_val))   filter(., censorFunc == bias_val) else . }
  
  make_one_theta <- function(df_source, use_TSn, theta_value) {
    
    df_source <- df_source %>%
      mutate(
        repN_num   = rep.number,
        nrep_num   = rep.n,
        censorFunc = factor(censorFunc, levels = c("low", "medium", "high")),
        rep.number = factor(
          rep.number,
          levels = c(2, 5, 10),
          labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10")
        ),
        rep.n = factor(
          rep.n,
          levels = c(40, 100, 400),
          labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400")
        ),
        method = factor(
          method,
          levels = c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA", "FEMA2")
        )
      )
    
    if (use_TSn == "TS1") {
      df_mabf <- df_source %>%
        filter(method %in% mabf_methods) %>%
        select(
          true_es, censorFunc, p_cutoff, orig.n, method, rep.n, rep.number, repN_num, nrep_num,
          ADTE, TS1, FF1, TF1, FS1
        ) %>%
        pivot_longer(
          cols = c(ADTE, TS1, FF1, TF1, FS1),
          names_to = "metric",
          values_to = "value"
        ) %>%
        mutate(
          metric_label = recode(
            metric,
            ADTE = "AD", TS1 = "TS", FF1 = "FF", TF1 = "TF", FS1 = "FS"
          )
        )
      
      df_meta <- df_source %>%
        filter(method %in% meta_methods) %>%
        select(
          true_es, censorFunc, p_cutoff, orig.n, method, rep.n, rep.number, repN_num, nrep_num,
          TS1, FF1, TF1, FS1
        ) %>%
        pivot_longer(
          cols = c(TS1, FF1, TF1, FS1),
          names_to = "metric",
          values_to = "value"
        ) %>%
        mutate(
          metric_label = recode(
            metric,
            TS1 = "TS", FF1 = "FF", TF1 = "TF", FS1 = "FS"
          )
        )
    } else {
      df_mabf <- df_source %>%
        filter(method %in% mabf_methods) %>%
        select(
          true_es, censorFunc, p_cutoff, orig.n, method, rep.n, rep.number, repN_num, nrep_num,
          ADNE, TS2, FF2, TF2, FS2
        ) %>%
        pivot_longer(
          cols = c(ADNE, TS2, FF2, TF2, FS2),
          names_to = "metric",
          values_to = "value"
        ) %>%
        mutate(
          metric_label = recode(
            metric,
            ADNE = "AD", TS2 = "TS", FF2 = "FF", TF2 = "TF", FS2 = "FS"
          )
        )
      
      df_meta <- df_source %>%
        filter(method %in% meta_methods) %>%
        select(
          true_es, censorFunc, p_cutoff, orig.n, method, rep.n, rep.number, repN_num, nrep_num,
          TS2, FF2, TF2, FS2
        ) %>%
        pivot_longer(
          cols = c(TS2, FF2, TF2, FS2),
          names_to = "metric",
          values_to = "value"
        ) %>%
        mutate(
          metric_label = recode(
            metric,
            TS2 = "TS", FF2 = "FF", TF2 = "TF", FS2 = "FS"
          )
        )
    }
    
    bind_rows(df_mabf, df_meta) %>%
      mutate(
        metric_label = factor(metric_label, levels = c("AD", "TS", "FF", "TF", "FS")),
        theta_label = factor(
          paste0("theta == ", theta_value),
          levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
        ),
        method_display = recode(
          as.character(method),
          "FEMA" = "REMA",
          "FEMA2" = "FEMA",
          .default = as.character(method)
        ),
        method_display = factor(
          method_display,
          levels = c("BFbMA", "EUBF", "FEMABF", "iBF", "REMA", "FEMA")
        )
      ) %>%
      group_by(
        theta_label, censorFunc, p_cutoff, orig.n, method, method_display,
        rep.n, rep.number, repN_num, nrep_num, metric_label
      ) %>%
      summarise(value = sum(value), .groups = "drop") %>%
      group_by(
        theta_label, censorFunc, p_cutoff, orig.n, method, method_display,
        rep.n, rep.number, repN_num, nrep_num
      ) %>%
      mutate(proportion = value / sum(value)) %>%
      ungroup()
  }
  
  if (dataset_choice == "all") {
    bind_rows(
      make_one_theta(df_filtered %>% filter(true_es == 0.2), "TS2", 0),
      make_one_theta(df_filtered %>% filter(true_es == 0.2), "TS1", 0.2),
      make_one_theta(df_filtered %>% filter(true_es == 0.5), "TS1", 0.5)
    )
  } else if (dataset_choice == "null") {
    make_one_theta(df_filtered %>% filter(true_es == 0.2), "TS2", 0)
  } else if (dataset_choice == "0.2") {
    make_one_theta(df_filtered %>% filter(true_es == 0.2), "TS1", 0.2)
  } else {
    make_one_theta(df_filtered %>% filter(true_es == 0.5), "TS1", 0.5)
  }
}

# =======================
# Build stacked plot
# =======================
make_stacked_plot <- function(df, filter_vars, facet_row, methods_to_show) {
  
  df <- df %>%
    filter(method_display %in% methods_to_show) %>%
    mutate(method_display = factor(method_display, levels = methods_to_show))
  
  group_vars <- c("theta_label", "method_display", "metric_label", filter_vars, facet_row)
  
  df_panel <- df %>%
    group_by(across(any_of(group_vars))) %>%
    summarise(proportion = mean(proportion), .groups = "drop") %>%
    mutate(
      label = ifelse(
        proportion > 0.01,
        paste0(sprintf("%.2f", proportion * 100), "%"),
        ""
      )
    )
  
  p <- ggplot(df_panel, aes(x = method_display, y = proportion, fill = metric_label)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3, width = 0.7) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 4
    ) +
    scale_fill_manual(
      values = class_colors,
      labels = class_labels,
      breaks = c("AD", "TS", "FF", "TF", "FS")
    ) +
    scale_y_continuous(
      limits = c(0, 1.02),
      breaks = seq(0, 1, 0.25),
      expand = c(0, 0)
    ) +
    labs(
      x = NULL,
      y = "Proportion",
      fill = "Classification"
    ) +
    theme_bw(base_size = 16) +
    theme(
      panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.8),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = "grey80", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey85", colour = "black", linewidth = 0.8),
      strip.text = element_text(size = 16),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.line = element_blank(),
      legend.position = "bottom",
      legend.key = element_rect(colour = "black", linewidth = 0.2),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
  
  facet_cols <- c("theta_label", filter_vars)
  
  if (length(facet_row) == 1) {
    p <- p + facet_grid(
      rows = vars(!!rlang::sym(facet_row)),
      cols = vars(!!!rlang::syms(facet_cols)),
      labeller = plot_labeller
    )
  } else {
    p <- p + facet_grid(
      cols = vars(!!!rlang::syms(facet_cols)),
      labeller = plot_labeller
    )
  }
  
  p
}

# =======================
# Build APA-style table
# =======================
make_apa_table <- function(df_long, facet_var) {
  method_order <- c("BFbMA", "EUBF", "FEMABF", "iBF", "REMA", "FEMA")
  metric_order <- c("AD", "TS", "FF", "TF", "FS")
  
  base <- df_long %>%
    mutate(
      percent = round(proportion * 100, 2),
      method_display = factor(method_display, levels = method_order),
      metric_label = factor(metric_label, levels = metric_order)
    ) %>%
    select(all_of(facet_var), repN_num, nrep_num, method_display, metric_label, percent) %>%
    unite(".col", method_display, metric_label, sep = "_") %>%
    pivot_wider(
      names_from = ".col",
      values_from = percent,
      values_fn = mean,
      values_fill = NA_real_
    ) %>%
    arrange(.data[[facet_var]], repN_num, nrep_num)
  
  method_major_cols <- unlist(lapply(method_order, function(m) {
    paste(m, metric_order, sep = "_")
  }))
  
  wanted_cols <- c(facet_var, "repN_num", "nrep_num", method_major_cols)
  for (cn in setdiff(wanted_cols, names(base))) base[[cn]] <- NA_real_
  base <- base[, wanted_cols]
  
  metric_cols <- setdiff(names(base), c(facet_var, "repN_num", "nrep_num"))
  base[metric_cols] <- lapply(base[metric_cols], function(x) {
    ifelse(is.na(x), "", format(round(x, 2), nsmall = 2))
  })
  
  clean_labels <- function(var, val) {
    if (var == "censorFunc") {
      recode(
        val,
        "Bias~level == 'low'"    = "Bias level: low",
        "Bias~level == 'medium'" = "Bias level: medium",
        "Bias~level == 'high'"   = "Bias level: high"
      )
    } else if (var == "p_cutoff") {
      recode(
        val,
        "alpha == 0.01" = "$\\alpha = 0.01$",
        "alpha == 0.05" = "$\\alpha = 0.05$"
      )
    } else if (var == "orig.n") {
      recode(
        val,
        "n[orig] == 20"  = "$n_{orig} = 20$",
        "n[orig] == 50"  = "$n_{orig} = 50$",
        "n[orig] == 200" = "$n_{orig} = 200$"
      )
    } else if (var == "theta_label") {
      recode(
        val,
        "theta == 0"   = "$\\theta = 0$",
        "theta == 0.2" = "$\\theta = 0.2$",
        "theta == 0.5" = "$\\theta = 0.5$"
      )
    } else {
      val
    }
  }
  
  pieces <- split(base, base[[facet_var]], drop = TRUE)
  
  out <- purrr::imap_dfr(pieces, function(df_piece, level_val) {
    df_piece <- df_piece %>%
      group_by(repN_num) %>%
      mutate(repN_show = if_else(row_number() == 1, as.character(repN_num), "")) %>%
      ungroup() %>%
      transmute(
        Nrep = repN_show,
        nrep = as.character(nrep_num),
        !!!rlang::syms(metric_cols)
      )
    
    header <- tibble(
      Nrep = paste0("\\textbf{", clean_labels(facet_var, level_val), "}"),
      nrep = "",
      !!!setNames(rep("", length(metric_cols)), metric_cols)
    )
    
    bind_rows(header, df_piece)
  })
  
  out
}

# =======================
# UI
# =======================
ui <- fluidPage(
  titlePanel("Replication Success Stacked Bar Plot + APA Tables"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "bf_cutoff", "Bayes Factor Cutoff",
        choices = c("3" = "c1", "10" = "c2"), selected = "c1"
      ),
      selectInput(
        "dataset_choice", "Data Scenario",
        choices = c("All" = "all", "Null" = "null", "0.2" = "0.2", "0.5" = "0.5"),
        selected = "all"
      ),
      selectInput(
        "p_cutoff", "Original p-value Cutoff",
        choices = c("All", "0.01", "0.05"), selected = "All"
      ),
      selectInput(
        "orig_n", "Original Sample Size",
        choices = c("All", "20", "50", "200"), selected = "All"
      ),
      selectInput(
        "bias_level", "Bias Level",
        choices = c("All", "low", "medium", "high"), selected = "All"
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
      checkboxGroupInput(
        "facet_row", "Optional Third Filter",
        choices = c(
          "Bias Level" = "censorFunc",
          "Original p-value Cutoff" = "p_cutoff",
          "Original Sample Size" = "orig.n"
        ),
        selected = NULL
      ),
      checkboxGroupInput(
        "filter_vars", "Additional Facet Columns",
        choices = c(
          "Replication Sample Size" = "rep.n",
          "Number of Replications" = "rep.number"
        ),
        selected = NULL
      ),
      br(),
      downloadButton("downloadPlot", "Download Plot (PNG)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Stacked Plot",
          plotOutput("piePlot", height = "900px")
        ),
        tabPanel(
          "APA Table",
          htmlOutput("apa_html"),
          br(),
          downloadButton("download_apa_csv", "Download APA Table (CSV)"),
          br(),
          downloadButton("download_apa_tex", "Download APA Table (LaTeX)"),
          br(),
          verbatimTextOutput("apa_latex")
        )
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
    df <- build_long(
      rates_all(),
      input$dataset_choice,
      input$p_cutoff,
      input$orig_n,
      input$bias_level
    )
    
    df %>%
      mutate(
        p_cutoff = factor(
          p_cutoff,
          levels = c(0.01, 0.05),
          labels = c("alpha == 0.01", "alpha == 0.05")
        ),
        orig.n = factor(
          orig.n,
          levels = c(20, 50, 200),
          labels = c("n[orig] == 20", "n[orig] == 50", "n[orig] == 200")
        ),
        censorFunc = factor(
          censorFunc,
          levels = c("low", "medium", "high"),
          labels = c(
            "Bias~level == 'low'",
            "Bias~level == 'medium'",
            "Bias~level == 'high'"
          )
        ),
        method_display = factor(
          method_display,
          levels = c("BFbMA", "EUBF", "FEMABF", "iBF", "REMA", "FEMA")
        ),
        theta_label = factor(
          theta_label,
          levels = c("theta == 0", "theta == 0.2", "theta == 0.5")
        )
      )
  })
  
  output$piePlot <- renderPlot({
    req(length(input$methods_to_show) > 0)
    make_stacked_plot(
      long_data(),
      input$filter_vars,
      input$facet_row,
      input$methods_to_show
    )
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      bf_cutoff <- ifelse(input$bf_cutoff == "c1", "BFcutoff3", "BFcutoff10")
      dataset   <- paste0("theta_", input$dataset_choice)
      filters   <- paste(c(input$filter_vars, input$facet_row), collapse = "_")
      methods   <- paste(input$methods_to_show, collapse = "-")
      
      if (filters == "") filters <- "nofilters"
      if (methods == "") methods <- "nomethods"
      
      paste0("StackedPlot_", bf_cutoff, "_", dataset, "_", filters, "_", methods, ".png")
    },
    content = function(file) {
      req(length(input$methods_to_show) > 0)
      p <- make_stacked_plot(
        long_data(),
        input$filter_vars,
        input$facet_row,
        input$methods_to_show
      )
      ggsave(file, plot = p, width = 16, height = 12, dpi = 300)
    }
  )
  
  apa_tbl <- reactive({
    if (input$dataset_choice == "all") {
      make_apa_table(long_data(), "theta_label")
    } else if (length(input$facet_row) == 1) {
      make_apa_table(long_data(), input$facet_row)
    } else {
      long_data() %>%
        group_by(theta_label, method_display, repN_num, nrep_num, metric_label) %>%
        summarise(percent = round(mean(proportion) * 100, 2), .groups = "drop") %>%
        pivot_wider(names_from = c(method_display, metric_label), values_from = percent)
    }
  })
  
  output$apa_html <- renderUI({
    tab <- apa_tbl()
    
    if (input$dataset_choice == "all" || length(input$facet_row) == 1) {
      header_methods <- c(
        " " = 2,
        "BFbMA" = 5,
        "EUBF" = 5,
        "FEMABF" = 5,
        "iBF" = 5,
        "REMA" = 5,
        "FEMA" = 5
      )
      
      header_rates <- c(
        " " = 2,
        rep(c("AD (%)", "TS (%)", "FF (%)", "TF (%)", "FS (%)"), times = 6)
      )
      
      kable(
        tab,
        format = "html",
        escape = FALSE,
        align = c("l", "l", rep("r", 30)),
        col.names = c("Nrep", "nrep", rep("", 30))
      ) %>%
        kable_styling(
          full_width = FALSE,
          bootstrap_options = c("condensed", "hover")
        ) %>%
        add_header_above(header_rates, bold = FALSE) %>%
        add_header_above(header_methods, bold = FALSE, line = FALSE) %>%
        HTML()
    } else {
      DT::datatable(tab)
    }
  })
  
  apa_tbl_latex <- reactive({
    if (input$dataset_choice == "all" || length(input$facet_row) == 1) {
      tab <- apa_tbl()
      
      header_methods <- c(
        " " = 2,
        "BFbMA" = 5,
        "EUBF" = 5,
        "FEMABF" = 5,
        "iBF" = 5,
        "REMA" = 5,
        "FEMA" = 5
      )
      
      header_rates <- c(
        "$N_{\\\\text{rep}}$", "$n_{\\\\text{rep}}$",
        rep(c("AD (\\\\%)", "TS (\\\\%)", "FF (\\\\%)", "TF (\\\\%)", "FS (\\\\%)"), times = 6)
      )
      
      if (input$dataset_choice == "all") {
        theta_text <- " for true effect sizes $\\theta = 0$, $\\theta = 0.2$, and $\\theta = 0.5$"
        facet_nice <- "true effect size"
      } else {
        theta_text <- switch(
          input$dataset_choice,
          "null" = " for true effect size $\\theta = 0$",
          "0.2"  = " for true effect size $\\theta = 0.2$",
          "0.5"  = " for true effect size $\\theta = 0.5$"
        )
        
        facet_nice <- c(
          censorFunc = "bias level",
          p_cutoff   = "original p-value cutoff",
          orig.n     = "original sample size"
        )[input$facet_row]
      }
      
      caption_text <- paste0(
        "Replication success measures (AD, TS, FF, TF, FS) across Bayes-factor and meta-analytic methods by ",
        "replication sample size ($n_{\\text{rep}}$) and number of replications ($N_{\\text{rep}}$)",
        theta_text, ", stratified by ", facet_nice, "."
      )
      
      kable(
        tab,
        format = "latex",
        booktabs = TRUE,
        escape = FALSE,
        align = c("c", "c", rep("r", 30)),
        col.names = rep("", 32),
        caption = caption_text
      ) %>%
        kable_styling(latex_options = c("hold_position")) %>%
        add_header_above(header_rates, bold = FALSE, line = FALSE, escape = FALSE) %>%
        add_header_above(header_methods, bold = FALSE, line = TRUE, escape = FALSE)
    } else {
      NULL
    }
  })
  
  output$apa_latex <- renderText({
    if (!is.null(apa_tbl_latex())) {
      paste(capture.output(print(apa_tbl_latex())), collapse = "\n")
    } else {
      "APA LaTeX table not available."
    }
  })
  
  output$download_apa_tex <- downloadHandler(
    filename = function() {
      paste0(
        "APA_table_",
        ifelse(input$dataset_choice == "all", "theta", ifelse(length(input$facet_row) == 1, input$facet_row, "collapsed")),
        ".tex"
      )
    },
    content = function(file) {
      if (!is.null(apa_tbl_latex())) {
        tex_code <- paste(capture.output(print(apa_tbl_latex())), collapse = "\n")
        writeLines(tex_code, file)
      } else {
        writeLines("No APA table available.", file)
      }
    }
  )
  
  output$download_apa_csv <- downloadHandler(
    filename = function() {
      paste0(
        "APA_table_",
        ifelse(input$dataset_choice == "all", "theta", ifelse(length(input$facet_row) == 1, input$facet_row, "collapsed")),
        ".csv"
      )
    },
    content = function(file) {
      write.csv(apa_tbl(), file, row.names = FALSE)
    }
  )
}

# =======================
# Run
# =======================
shinyApp(ui, server)