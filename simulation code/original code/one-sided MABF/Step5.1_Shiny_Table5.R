library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(knitr)
library(kableExtra)
library(htmltools)
library(tibble)

ui <- fluidPage(
  titlePanel("MABF Evidence Strength Proportions"),
  sidebarLayout(
    sidebarPanel(
      selectInput("effect_size", "Underlying Effect Size:",
                  choices = c("null", "0.2", "0.5"), selected = "null")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Stacked Bar Plot",
                 plotOutput("bar_plot", width = "1200px", height = "800px")),
        tabPanel("Summary Table",
                 DTOutput("summary_table")),
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

server <- function(input, output, session) {
  # Load once
  env <- new.env()
  try(load("./MABFanalyses/row-wise/row2matrix.RData", envir = env), silent = TRUE)
  loaded_env <- reactiveVal(env)
  
  # metric definitions
  metric_cols <- reactive({
    if (input$effect_size == "null") {
      c("ADR_null_anecdotal","ADR_null_moderate","ADR_null_strong",
        "ADR_null_vstrong","ADR_null_estrong")
    } else {
      c("ADR_true_anecdotal","ADR_true_moderate","ADR_true_strong",
        "ADR_true_vstrong","ADR_true_estrong")
    }
  })
  # labels with APA style
  metric_labels <- reactive(c("anecdotal","moderate","strong","vx strong","ex strong"))
  
  raw_ds <- reactive({
    data_name <- paste0(
      "MABF_rates_",
      ifelse(input$effect_size == "null","null",input$effect_size),
      "_rowise_BFcutoff3"
    )
    ds <- tryCatch(get(data_name, envir = loaded_env()), error = function(e) NULL)
    validate(need(!is.null(ds), "Dataset not found."))
    ds
  })
  
  # -------- Plot data --------
  plot_data <- reactive({
    raw_ds() %>%
      mutate(
        rep.number = factor(rep.number, levels = c(2,5,10),
                            labels = c("N[rep]==2","N[rep]==5","N[rep]==10")),
        rep.n      = factor(rep.n, levels = c(40,100,400),
                            labels = c("n[rep]==40","n[rep]==100","n[rep]==400"))
      ) %>%
      select(method, rep.number, rep.n, all_of(metric_cols())) %>%
      tidyr::pivot_longer(all_of(metric_cols()), names_to = "metric", values_to = "value") %>%
      mutate(metric_label = factor(metric, levels = metric_cols(), labels = metric_labels())) %>%
      group_by(method, rep.number, rep.n, metric_label) %>%
      summarise(proportion = mean(value, na.rm = TRUE), .groups = "drop")
  })
  
  
  output$bar_plot <- renderPlot({
    pd <- plot_data()
    cols <- setNames(c("#B0E0E6","#87CEEB","#1E90FF","#0000CD","#00008B"), levels(pd$metric_label))
    ggplot(pd, aes(x = factor(method), y = proportion, fill = metric_label)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = cols, name = "Evidence Strength") +
      labs(x = "MABF Method", y = "Proportion") +
      facet_grid(rep.n ~ rep.number, scales = "free_x",
                 labeller = labeller(rep.number = label_parsed,
                                     rep.n = label_parsed)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "bottom")
  })
  
  
  # -------- Summary table --------
  output$summary_table <- DT::renderDT({
    raw_ds() %>%
      mutate(
        rep.number = factor(rep.number, levels = c(2,5,10)),
        rep.n      = factor(rep.n,      levels = c(40,100,400))
      ) %>%
      select(method, rep.number, rep.n, all_of(metric_cols())) %>%
      group_by(method, rep.number, rep.n) %>%
      summarise(across(all_of(metric_cols()), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
      arrange(method, rep.number, rep.n) %>%
      mutate(across(all_of(metric_cols()), ~ formatC(., format = "f", digits = 6))) %>%
      DT::datatable(options = list(pageLength = 20, autoWidth = TRUE))
  })
  
  # -------- APA table --------
  apa_tbl <- reactive({
    df <- raw_ds() %>%
      mutate(
        rep.number = factor(rep.number, levels = c(2,5,10)),   # Nrep
        rep.n      = factor(rep.n,      levels = c(40,100,400))# nrep
      ) %>%
      select(method, rep.number, rep.n, all_of(metric_cols())) %>%
      group_by(method, rep.number, rep.n) %>%
      summarise(across(all_of(metric_cols()), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
      rename_with(~ metric_labels(), all_of(metric_cols())) %>%
      mutate(method = factor(method, levels = c("BFbMA","EUBF","FEMABF","iBF"))) %>%
      arrange(method, rep.number, rep.n)
    
    #  convert to percentages
    df <- df %>%
      mutate(across(all_of(metric_labels()), ~ formatC(. * 100, format = "f", digits = 2)))
    
    #  blank duplicates for APA style
    df <- df %>%
      group_by(method) %>%
      mutate(Method_disp = if_else(row_number() == 1, as.character(method), "")) %>%
      ungroup() %>%
      group_by(method, rep.number) %>%
      mutate(Nrep_disp   = if_else(row_number() == 1, as.character(rep.number), "")) %>%
      ungroup() %>%
      mutate(nrep_disp   = as.character(rep.n)) %>%
      select(Method_disp, Nrep_disp, nrep_disp, all_of(metric_labels()))
    
    df
  })
  
  output$apa_html <- renderUI({
    df <- apa_tbl()
    col_labels <- c("Method","Nrep","nrep",
                    paste0(metric_labels(), " (%)"))
    knitr::kable(df, format = "html", booktabs = TRUE,
                 col.names = col_labels, escape = FALSE) %>%
      kableExtra::kable_styling(full_width = FALSE,
                                bootstrap_options = c("condensed","hover")) %>%
      as.character() %>% htmltools::HTML()
  })
  
  output$apa_latex <- renderText({
    df <- apa_tbl()
    col_labels <- c("Method","$N_{rep}$","$n_{rep}$",
                    paste0(metric_labels(), " (\\%)"))
    as.character(
      knitr::kable(df, format = "latex", booktabs = TRUE,
                   col.names = col_labels, escape = FALSE,
                   caption = paste("Proportions of evidence strength categories across MABF methods (one-sided test) when effect size =",
                                   input$effect_size, "and cutoff = 3.")) %>%
        kableExtra::kable_styling(latex_options = "hold_position", full_width = FALSE)
    )
  })
  
  output$download_apa_tex <- downloadHandler(
    filename = function() paste0("APA_table_effectsize_",input$effect_size,".tex"),
    content = function(file) {
      writeLines(output$apa_latex(), file)
    }
  )
}


shinyApp(ui = ui, server = server)
