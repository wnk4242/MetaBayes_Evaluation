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
