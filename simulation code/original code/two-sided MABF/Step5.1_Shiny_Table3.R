# This Shiny app generates the AUC table for underlying effect = 0.2 and 0.5 (two-sided test)
library(shiny)
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)

ui <- fluidPage(
  titlePanel("AUC Summary Tables for MABF Methods"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("theta", "Select underlying effect (θ):",
                  choices = c("0.2", "0.5"),
                  selected = "0.2")
    ),
    
    mainPanel(
      h4("LaTeX Table Output"),
      verbatimTextOutput("latex_table")
    )
  )
)

server <- function(input, output, session) {
  
  # Load once
  load("./MABFanalyses/matrix-wise/AUC/MABF_0.2null_AUC.RData")
  load("./MABFanalyses/matrix-wise/AUC/MABF_0.5null_AUC.RData")
  
  AUC0.2 <- bind_rows(BFbMA_0.2null_AUC,
                      EUBF_0.2null_AUC,
                      FEMABF_0.2null_AUC,
                      iBF_0.2null_AUC,
                      FEMA_0.2null_AUC)   # FEMA moved to the end
  
  AUC0.5 <- bind_rows(BFbMA_0.5null_AUC,
                      EUBF_0.5null_AUC,
                      FEMABF_0.5null_AUC,
                      iBF_0.5null_AUC,
                      FEMA_0.5null_AUC)   # FEMA moved to the end
  
  # Reactive dataset
  selected_data <- reactive({
    if (input$theta == "0.2") {
      AUC0.2
    } else {
      AUC0.5
    }
  })
  
  output$latex_table <- renderText({
    dat <- selected_data() %>%
      mutate(
        ci95 = paste0(formatC(mc_ci_lo, format = "f", digits = 3),
                      " -- ",
                      formatC(mc_ci_hi, format = "f", digits = 3)),
        mean_AUC = formatC(mean_AUC, format = "f", digits = 3),
        sd_AUC   = formatC(sd_AUC, format = "f", digits = 3)
      ) %>%
      select(method, rep.number, rep.n, mean_AUC, sd_AUC, ci95)
    
    dat_wide <- dat %>%
      pivot_wider(
        id_cols = c(rep.number, rep.n),
        names_from = method,
        values_from = c(mean_AUC, sd_AUC, ci95),
        names_glue = "{method}_{.value}"
      ) %>%
      arrange(rep.number, rep.n) %>%
      select(rep.number, rep.n,
             BFbMA_mean_AUC, BFbMA_sd_AUC, BFbMA_ci95,
             EUBF_mean_AUC,  EUBF_sd_AUC,  EUBF_ci95,
             FEMABF_mean_AUC,FEMABF_sd_AUC,FEMABF_ci95,
             iBF_mean_AUC,   iBF_sd_AUC,   iBF_ci95,
             FEMA_mean_AUC,  FEMA_sd_AUC,  FEMA_ci95)   # FEMA placed last
    
    # Column labels
    col_names <- c("$N_{\\text{rep}}$", "$n_{\\text{rep}}$",
                   "AUC", "SD", "95\\% CI",
                   "AUC", "SD", "95\\% CI",
                   "AUC", "SD", "95\\% CI",
                   "AUC", "SD", "95\\% CI",
                   "AUC", "SD", "95\\% CI")
    
    kable(
      dat_wide,
      format = "latex",
      booktabs = TRUE,
      align = c("c","c", rep("c", 15)),
      col.names = col_names,
      caption = paste0("Area under the curve (AUC), standard deviation (SD), and 95\\% confidence intervals (CI) across MABF methods (two-sided test) by replication sample size ($n_{\\text{rep}}$) and number of replications ($N_{\\text{rep}}$), for true effect size $\\theta = ", input$theta, "$."),
      escape = FALSE
    ) %>%
      add_header_above(c(" " = 2, "BFbMA" = 3, "EUBF" = 3, "FEMABF" = 3, "iBF" = 3, "FEMA" = 3)) %>%
      as.character()
  })
}

shinyApp(ui, server)
