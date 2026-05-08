# This Shiny generates tables for TPR, FPR for the MABF methods (one-sided test)
library(shiny)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(tibble)
library(knitr)
library(kableExtra)
library(scales)
library(htmltools)

# ====== Load data =======

load('./MABF4ROC/4TSFS/MABF_lists_0.20.5null_regrouped_deltap.RData')

# ====== Table =======
make_rate_table <- function(effect_size = "0",
                            rate_type   = c("FPR", "TPR"),
                            cutoffs     = c(1, 3, 10)) {
  rate_type <- match.arg(rate_type)
  
  
  if (rate_type == "FPR") {
    summarize_fun <- function(method_name, cutoff, effect_size) {
      effect_key <- if (effect_size == "0") "0.2" else effect_size
      data_list <- get(paste0(method_name, "_lists_", effect_key, "null_regrouped_deltap"))
      purrr::map_dfr(names(data_list), function(dataset_name) {
        dataset <- data_list[[dataset_name]]
        purrr::map_dfr(names(dataset), function(setting_name) {
          mat <- dataset[[setting_name]][1:500, 4:503]   # null studies
          count <- sum(mat > cutoff)
          tibble(method = method_name, setting = setting_name,
                 Rate = count / 250000)
        })
      }) %>%
        separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
        group_by(method, rep_num, rep_n) %>%
        summarise(Rate = mean(Rate), .groups = "drop")
    }
  } else {
    summarize_fun <- function(method_name, cutoff, effect_size) {
      data_list <- get(paste0(method_name, "_lists_", effect_size, "null_regrouped_deltap"))
      purrr::map_dfr(names(data_list), function(dataset_name) {
        dataset <- data_list[[dataset_name]]
        purrr::map_dfr(names(dataset), function(setting_name) {
          mat <- dataset[[setting_name]][501:1000, 4:503]  # true-effect studies
          count <- sum(mat > cutoff)
          tibble(method = method_name, setting = setting_name,
                 Rate = count / 250000)
        })
      }) %>%
        separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
        group_by(method, rep_num, rep_n) %>%
        summarise(Rate = mean(Rate), .groups = "drop")
    }
  }
  
  
  summarize_cutoff <- function(cutoff, effect_size) {
    methods <- c("FEMABF", "BFbMA", "iBF", "EUBF")
    dat <- purrr::map_dfr(methods, ~ summarize_fun(.x, cutoff, effect_size))
    dat %>%
      pivot_wider(names_from = method, values_from = Rate) %>%
      mutate(cutoff = cutoff) %>%
      # Column order: Nrep, nrep, FEMABF, iBF, EUBF, BFbMA
      select(cutoff, rep_num, rep_n, FEMABF, iBF, EUBF, BFbMA)
  }
  
  all_dat <- purrr::map_dfr(cutoffs, ~ summarize_cutoff(.x, effect_size))
  
  
  all_fmt <- all_dat %>%
    mutate(
      cutoff  = factor(cutoff,  levels = cutoffs, ordered = TRUE),
      rep_num = factor(rep_num, levels = c(2, 5, 10), ordered = TRUE),
      rep_n   = factor(rep_n,   levels = c(40, 100, 400), ordered = TRUE)
    ) %>%
    mutate(across(FEMABF:BFbMA,
                  ~ gsub("%", "\\\\%", scales::percent(.x, accuracy = 0.01))))
  
  
  blocks <- split(all_fmt, all_fmt$cutoff)
  
  nested <- purrr::imap_dfr(blocks, function(df, key) {
    df <- arrange(df, rep_num, rep_n)
    rep_blocks <- split(df, df$rep_num)
    
    rep_nested <- purrr::imap_dfr(rep_blocks, function(df2, repkey) {
      df2 %>%
        select(-cutoff) %>%
        mutate(rep_num = as.character(rep_num),
               rep_n   = as.character(rep_n)) %>%
        # keep rep_num only on the first row; blank for remaining rows
        mutate(rep_num = ifelse(row_number() == 1, rep_num, ""))
    })
    
    hdr1 <- tibble(
      rep_num = paste0("\\textbf{BF cutoff = ", key, "}"),
      rep_n   = "",
      FEMABF = "", iBF = "", EUBF = "", BFbMA = ""
    )
    bind_rows(hdr1, rep_nested)
  })
  
  
  cap <- if (rate_type == "FPR") {
    "False Positive Rates (FPR)"
  } else {
    paste0("True Positive Rates (TPR) for true effect size $\\theta = ", effect_size, "$")
  }
  
  
  list(
    latex = kable(
      nested,
      format   = "latex",
      booktabs = TRUE,
      caption  = paste(
        cap,
        "across MABF methods (one-sided test) by replication sample size ($n_{\\text{rep}}$) and number of replications ($N_{\\text{rep}}$), for Bayes factor cutoffs",
        paste(cutoffs, collapse = ", "), "."
      ),
      col.names = c("$N_{\\text{rep}}$", "$n_{\\text{rep}}$", "FEMABF", "iBF", "EUBF", "BFbMA"),
      escape   = FALSE,
      linesep  = "",
      align    = c("l","r","r","r","r","r")
    ) %>% kable_styling(latex_options = "hold_position", full_width = FALSE),
    
    html = kable(
      nested,
      format   = "html",
      caption  = paste(
        cap,
        "across MABF methods (one-sided test) by replication sample size ($n_{\\text{rep}}$) and number of replications ($N_{\\text{rep}}$), for Bayes factor cutoffs",
        paste(cutoffs, collapse = ", "), "."
      ),
      col.names = c("Nrep", "nrep", "FEMABF", "iBF", "EUBF", "BFbMA"),
      escape   = FALSE
    ) %>%
      kable_styling(full_width = FALSE, bootstrap_options = c("striped", "condensed"))
  )
}

# ======= UI =======
ui <- fluidPage(
  titlePanel("MABF Rates Table Builder (APA LaTeX)"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons(
        "rate_type", "Rate Type",
        choices = c("FPR", "TPR"), selected = "FPR", inline = TRUE
      ),
      conditionalPanel(
        condition = "input.rate_type == 'TPR'",
        selectInput(
          "effect_size", "True Effect Size (θ)",
          choices = c("0.2", "0.5"), selected = "0.2"
        )
      ),
      conditionalPanel(
        condition = "input.rate_type == 'FPR'",
        helpText("FPR uses θ = 0 internally (null studies).")
      ),
      selectizeInput(
        "cutoffs", "Bayes Factor Cutoffs",
        choices = c(1, 3, 10), selected = c(1, 3, 10),
        multiple = TRUE, options = list(create = TRUE, placeholder = "Type to add more...")
      ),
      tags$hr(),
      downloadButton("download_tex", "Download LaTeX (.tex)"),
      width = 3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Preview (HTML)",
                 htmlOutput("table_html")
        ),
        tabPanel("LaTeX",
                 tags$small("Copy and paste this into Overleaf:"),
                 tags$br(),
                 verbatimTextOutput("table_latex", placeholder = TRUE)
        )
      ),
      width = 9
    )
  )
)

# ======= Server =======
server <- function(input, output, session) {
  
  # reactive parameters
  params <- reactive({
    list(
      rate_type   = input$rate_type,
      effect_size = if (input$rate_type == "TPR") input$effect_size else "0",
      cutoffs     = as.numeric(input$cutoffs)
    )
  })
  
  # build tables
  tables <- reactive({
    p <- params()
    validate(
      need(length(p$cutoffs) > 0, "Please select at least one BF cutoff.")
    )
    make_rate_table(
      effect_size = p$effect_size,
      rate_type   = p$rate_type,
      cutoffs     = sort(unique(p$cutoffs))
    )
  })
  
  # HTML preview
  output$table_html <- renderUI({
    as.character(tables()$html) |> HTML()
  })
  
  # LaTeX
  output$table_latex <- renderText({
    as.character(tables()$latex)
  })
  
  # Download handler for .tex
  output$download_tex <- downloadHandler(
    filename = function() {
      p <- params()
      if (p$rate_type == "FPR") {
        "FPR_table.tex"
      } else {
        paste0("TPR_theta_", gsub("\\.", "_", p$effect_size), "_table.tex")
      }
    },
    content = function(file) {
      writeLines(as.character(tables()$latex), con = file, useBytes = TRUE)
    }
  )
}

shinyApp(ui, server)
