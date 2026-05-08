#####
# library(shiny)
# library(ggplot2)
# library(dplyr)
# 
# # Define scenario lookup table
# scenario_df <- data.frame(
#   scenario = 1:9,
#   true_es = 0.2,
#   rep.number = rep(c(2, 5, 10), times = 3),
#   rep.n = rep(c(40, 100, 400), each = 3)
# )
# 
# # === UI ===
# ui <- fluidPage(
#   tags$head(
#     tags$style(HTML("
#       .sidebar {
#         background-image: url('https://www.transparenttextures.com/patterns/cubes.png');
#       }
#     "))
#   ),
#   
#   tags$div(style = "max-width: 1400px; margin: auto;",
#            titlePanel("Evaluating Meta-Analytic Bayes Factors with ROC Curves"),
#            fluidRow(
#              column(width = 4,
#                     tags$div(class = "sidebar",
#                              selectInput("es_choice", "True Effect Size:", choices = c("0.2", "0.5")),
#                              selectInput("rep_number", "Number of Replications:", choices = c(2, 5, 10)),
#                              selectInput("rep_n", "Replication Sample Size:", choices = c(40, 100, 400)),
#                              actionButton("plot_btn", "Generate ROC Plot"),
#                              downloadButton("downloadPlot", "Download ROC Plot (PNG)"),
#                              
#                              tags$div(
#                                style = "margin-top: 30px; padding: 15px; background-color: #f7f7f7; border: 1px solid #ddd; border-radius: 5px;",
#                                HTML("<strong>Meta-Analytic Bayes Factors (MABFs):</strong><br/><br/>
#             <ul style='margin-bottom: 0; padding-left: 20px;'>
#               <li><strong>FEMABF</strong>: Fixed-Effect Meta-Analytic Bayes Factor</li>
#               <li><strong>BFbMA</strong>: Bayes Factor Based on Meta-Analysis</li>
#               <li><strong>EUBF</strong>: Evidence-Updating Bayes Factor</li>
#               <li><strong>iBF</strong>: Inclusion Bayes Factor</li>
#               <li><strong>REMA</strong>: Random-Effects Meta-Analysis</li>
#             </ul>")
#                              )
#                     )
#              ),
#              
#              column(width = 8,
#                     tags$div(class = "main",
#                              tabsetPanel(
#                                tabPanel("Study Overview",
#                                         tags$div(style = "margin-top: 20px;",
#                                                  HTML("<p>Recent concerns about the replication crisis in psychology have underscored the need for better frameworks to assess whether original findings can be reliably replicated. This study aims to compare the performance of several meta-analytic Bayes factors (FEMABF, BFbMA, EUBF, and iBF) both relative to one another and against the traditional fixed-effect meta-analysis in evaluating replication success.</p>
#                 <p>To do so, we use a large-scale simulation to generate original studies across varied scenarios by manipulating true effect size, research environment (levels of p-hacking and publication bias), and original sample size. Replication studies are then simulated by varying replication sample size and number of replications, yielding 243 distinct scenarios and over 60 million data points. This Shiny app visualizes ROC curves to compare how well each method distinguishes between true and null effects across these conditions.</p>")
#                                         )
#                                ),
#                                tabPanel("ROC Curves",
#                                         uiOutput("plotTitle"),   # Title above plot
#                                         plotOutput("rocPlot", height = "600px")
#                                )
#                              )
#                     )
#              )
#            )
#   )
# )
# 
# # === Server ===
# server <- function(input, output) {
#   scenario_selected <- eventReactive(input$plot_btn, {
#     scenario_df %>%
#       filter(rep.number == input$rep_number,
#              rep.n == input$rep_n) %>%
#       pull(scenario)
#   })
#   
#   plot_event <- eventReactive(input$plot_btn, {
#     # Use isolate to avoid reactivity to es_choice dropdown
#     selected_es <- isolate(input$es_choice)
#     selected_scenario <- isolate(scenario_selected())
#     
#     # Load the corresponding RDS file based on true effect size (precomputed data)
#     file_path <- switch(selected_es,
#                         "0.2" = "./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.2null.RDS",
#                         "0.5" = "./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.5null.RDS")
#     df <- readRDS(file_path)
#     
#     condition <- scenario_df %>% filter(scenario == selected_scenario)
#     
#     rep_number <- condition$rep.number[1]
#     rep_n <- condition$rep.n[1]
#     
#     title_text <- paste0(
#       "ROC curves of the MABF methods (two-sided test) and random-effects meta-analysis when the number of replications is ",
#       rep_number,
#       ", replication sample size is ",
#       rep_n,
#       ", true effect is ",
#       selected_es
#     )
#     
#     list(
#       plot = df %>%
#         mutate(method = recode(method,
#                                "FEMA"   = "REMA")) %>% 
#         mutate(method = factor(method,
#                                levels = c("BFbMA", "EUBF", "REMA", "FEMABF","iBF"))) %>% 
#         filter(scenario == as.character(selected_scenario)) %>%
#         group_by(method) %>%
#         ggplot(aes(x = FPR, y = TPR, color = method)) +
#         geom_step() +
#         labs(x = "False Positive Rate", 
#              y = "True Positive Rate",
#              color = "Method") +  # removed title
#         scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
#         scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
#         geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
#         theme_minimal(),
#       title = title_text
#     )
#   })
#   
#   output$rocPlot <- renderPlot({
#     plot_event()$plot
#   })
#   
#   output$plotTitle <- renderUI({
#     req(plot_event())
#     h4(plot_event()$title)   # Title displayed above plot
#   })
#   
#   output$downloadPlot <- downloadHandler(
#     filename = function() {
#       paste0("ROC_plot_", input$rep_number, "_repN", input$rep_n, "_ES", input$es_choice, "_onesided", ".png")
#     },
#     content = function(file) {
#       ggsave(file, plot = plot_event()$plot, width = 8, height = 6, dpi = 300, bg = "white")
#     }
#   )
# }
# 
# # Run the app
# shinyApp(ui = ui, server = server)


# #######All scenario version###############
# # This version displays all scenarios together in a 3x3 grid
# # === Packages ===
# library(shiny)
# library(ggplot2)
# library(dplyr)
# library(patchwork)
# library(cowplot)   # for get_legend()
# 
# # === Scenario Lookup Table ===
# scenario_df <- data.frame(
#   scenario = 1:9,
#   true_es = 0.2,
#   rep.number = rep(c(2, 5, 10), times = 3),
#   rep.n = rep(c(40, 100, 400), each = 3)
# )
# 
# # === UI ===
# ui <- fluidPage(
#   tags$head(
#     tags$style(HTML("
#       .sidebar {
#         background-image: url('https://www.transparenttextures.com/patterns/cubes.png');
#       }
#     "))
#   ),
#   
#   tags$div(style = "max-width: 1400px; margin: auto;",
#            titlePanel("Evaluating Meta-Analytic Bayes Factors with ROC Curves"),
#            
#            fluidRow(
#              # Sidebar
#              column(width = 4,
#                     tags$div(class = "sidebar",
#                              selectInput("es_choice", "True Effect Size:", choices = c("0.2", "0.5")),
#                              actionButton("plot_btn", "Generate ROC Plot"),
#                              downloadButton("downloadPlot", "Download ROC Plot (PNG)"),
#                              
#                              tags$div(
#                                style = "margin-top: 30px; padding: 15px; background-color: #f7f7f7; border: 1px solid #ddd; border-radius: 5px;",
#                                HTML("<strong>Meta-Analytic Bayes Factors (MABFs):</strong><br/><br/>
#             <ul style='margin-bottom: 0; padding-left: 20px;'>
#               <li><strong>FEMABF</strong>: Fixed-Effect Meta-Analytic Bayes Factor</li>
#               <li><strong>BFbMA</strong>: Bayes Factor Based on Meta-Analysis</li>
#               <li><strong>EUBF</strong>: Evidence-Updating Bayes Factor</li>
#               <li><strong>iBF</strong>: Inclusion Bayes Factor</li>
#               <li><strong>REMA</strong>: Random-Effects Meta-Analysis</li>
#             </ul>")
#                              )
#                     )
#              ),
#              
#              # Main panel
#              column(width = 8,
#                     tags$div(class = "main",
#                              tabsetPanel(
#                                tabPanel("Study Overview",
#                                         tags$div(style = "margin-top: 20px;",
#                                                  HTML("<p>Recent concerns about the replication crisis in psychology have underscored the need for better frameworks to assess whether original findings can be reliably replicated. This study aims to compare the performance of several meta-analytic Bayes factors (FEMABF, BFbMA, EUBF, and iBF) both relative to one another and against the traditional fixed-effect meta-analysis in evaluating replication success.</p>
#                 <p>To do so, we use a large-scale simulation to generate original studies across varied scenarios by manipulating true effect size, research environment (levels of p-hacking and publication bias), and original sample size. Replication studies are then simulated by varying replication sample size and number of replications, yielding 243 distinct scenarios and over 60 million data points. This Shiny app visualizes ROC curves to compare how well each method distinguishes between true and null effects across these conditions.</p>")
#                                         )
#                                ),
#                                tabPanel("ROC Curves",
#                                         uiOutput("plotTitle"),
#                                         plotOutput("rocPlot", height = "900px")
#                                )
#                              )
#                     )
#              )
#            )
#   )
# )
# 
# # === Server ===
# server <- function(input, output) {
#   
#   plot_event <- eventReactive(input$plot_btn, {
#     selected_es <- isolate(input$es_choice)
#     
#     # Load corresponding RDS file
#     file_path <- switch(selected_es,
#                         "0.2" = "./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.2null.RDS",
#                         "0.5" = "./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.5null.RDS")
#     df <- readRDS(file_path)
#     
#     # Prepare data
#     df <- df %>%
#       mutate(method = ifelse(method == "FEMA", "REMA", method)) %>%
#       mutate(method = factor(method, levels = c("BFbMA", "EUBF", "REMA", "FEMABF", "iBF")),
#              scenario = as.numeric(scenario))
#     
#     # Generate base ROC plots without legend
#     base_plots <- lapply(1:9, function(scn) {
#       cond <- scenario_df %>% filter(scenario == scn)
#       ggplot(df %>% filter(scenario == scn),
#              aes(x = FPR, y = TPR, color = method)) +
#         geom_step(size = 0.7) +
#         geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
#         labs(
#           title = bquote(N[rep] == .(cond$rep.number) * "," ~ n[rep] == .(cond$rep.n)),
#           x = "False Positive Rate",
#           y = "True Positive Rate"
#         ) +
#         scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
#         scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
#         theme_minimal(base_size = 11) +
#         theme(
#           plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
#           legend.position = "none"
#         )
#     })
#     
#     # Extract shared legend (cowplot)
#     legend_plot <- ggplot(df, aes(x = FPR, y = TPR, color = method)) +
#       geom_step() +
#       theme_minimal() +
#       theme(
#         legend.position = "right",
#         legend.title = element_text(size = 10),
#         legend.text = element_text(size = 9)
#       ) +
#       guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
#       labs(color = "Method")
#     
#     legend <- cowplot::get_legend(legend_plot)
#     
#     # Combine 3x3 grid and legend
#     combined_plot <- wrap_plots(base_plots, ncol = 3) +
#       plot_annotation(
#         title = paste0("ROC Curves of MABF Methods (True Effect = ", selected_es, ")"),
#         theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
#       )
#     
#     final_plot <- combined_plot | patchwork::wrap_elements(legend)
#     final_plot <- final_plot + plot_layout(widths = c(3, 0.5))
#     
#     list(plot = final_plot,
#          title = paste0("ROC Curves for All Nine Scenarios (True Effect = ", selected_es, ")"))
#   })
#   
#   # Render title
#   output$plotTitle <- renderUI({
#     req(plot_event())
#     h4(plot_event()$title)
#   })
#   
#   # Render plot
#   output$rocPlot <- renderPlot({
#     plot_event()$plot
#   })
#   
#   # Download handler
#   output$downloadPlot <- downloadHandler(
#     filename = function() {
#       paste0("ROC_GridPlot_ES", input$es_choice, "_AllScenarios.png")
#     },
#     content = function(file) {
#       ggsave(file, plot = plot_event()$plot, width = 14, height = 10, dpi = 300, bg = "white")
#     }
#   )
# }
# 
# # === Run ===
# shinyApp(ui = ui, server = server)
# 


######
#Old version, no option to select ShinyROC_0.2null or ShinyROC_0.5null datasets
# library(shiny)
# library(ggplot2)
# library(dplyr)
# 
# # Load your precomputed simulation data
# rates_0.2null <- readRDS(file = "./MABFanalyses/matrix-wise/rates4ROC/shinyROC.RDS")
# 
# # Define scenario lookup table
# scenario_df <- data.frame(
#   scenario = 1:9,
#   true_es = 0.2,
#   rep.number = rep(c(2, 5, 10), times = 3),
#   rep.n = rep(c(40, 100, 400), each = 3)
# )
# 
# # UI
# ui <- fluidPage(
#   tags$head(
#     tags$style(HTML("
#       .sidebar {
#         background-image: url('https://www.transparenttextures.com/patterns/cubes.png');
#       }
#     "))
#   ),
#   
#   tags$div(style = "max-width: 1400px; margin: auto;",
#            titlePanel("Evaluating Meta-Analytic Bayes Factors with ROC Curves"),
#            fluidRow(
#              column(width = 4,
#                     tags$div(class = "sidebar",
#                              selectInput("rep_number", "Number of Replications:", choices = c(2, 5, 10)),
#                              selectInput("rep_n", "Replication Sample Size:", choices = c(40, 100, 400)),
#                              actionButton("plot_btn", "Generate ROC Plot"),
#                              
#                              tags$div(
#                                style = "margin-top: 30px; padding: 15px; background-color: #f7f7f7; border: 1px solid #ddd; border-radius: 5px;",
#                                HTML("<strong>Meta-Analytic Bayes Factors (MABFs):</strong><br/><br/>
#             <ul style='margin-bottom: 0; padding-left: 20px;'>
#               <li><strong>FEMABF</strong>: Fixed-Effect Meta-Analytic Bayes Factor</li>
#               <li><strong>BFbMA</strong>: Bayes Factor Based on Meta-Analysis</li>
#               <li><strong>EUBF</strong>: Evidence-Updating Bayes Factor</li>
#               <li><strong>iBF</strong>: Inclusion Bayes Factor</li>
#               <li><strong>FEMA</strong>: Fixed-Effect Meta-Analysis</li>
#             </ul>")
#                              )
#                     )
#              ),
#              
#              column(width = 8,
#                     tags$div(class = "main",
#                              tabsetPanel(
#                                tabPanel("Study Overview",
#                                         tags$div(style = "margin-top: 20px;",
#                                                  HTML("<p>Recent concerns about the replication crisis in psychology have underscored the need for better frameworks to assess whether original findings can be reliably replicated. This study aims to compare the performance of several meta-analytic Bayes factors (FEMABF, BFbMA, EUBF, and iBF) both relative to one another and against the traditional fixed-effect meta-analysis in evaluating replication success.</p>
#                 <p>To do so, we use a large-scale simulation to generate original studies across varied scenarios by manipulating true effect size, research environment (levels of p-hacking and publication bias), and original sample size. Replication studies are then simulated by varying replication sample size and number of replications, yielding 243 distinct scenarios and over 60 million data points. This Shiny app visualizes ROC curves to compare how well each method distinguishes between true and null effects across these conditions.</p>")
#                                         )
#                                ),
#                                tabPanel("ROC Curves",
#                                         plotOutput("rocPlot", height = "600px")
#                                )
#                              )
#                     )
#              )
#            )
#   )
# )
# 
# # Server
# server <- function(input, output) {
#   scenario_selected <- eventReactive(input$plot_btn, {
#     scenario_df %>%
#       filter(rep.number == input$rep_number,
#              rep.n == input$rep_n) %>%
#       pull(scenario)
#   })
#   
#   output$rocPlot <- renderPlot({
#     req(scenario_selected())
#     
#     scenario_number <- scenario_selected()
#     condition <- scenario_df %>% filter(scenario == scenario_number)
#     
#     true_es <- condition$true_es[1]
#     rep_number <- condition$rep.number[1]
#     rep_n <- condition$rep.n[1]
#     
#     title_text <- paste0(
#       "ROC Curves when the number of replications is ",
#       rep_number,
#       ", replication sample size is ",
#       rep_n,
#       ", true effect is ",
#       true_es
#     )
#     
#     rates_0.2null %>%
#       filter(scenario == as.character(scenario_number)) %>%
#       group_by(method) %>%
#       ggplot(aes(x = FPR, y = TPR, color = method)) +
#       geom_step() +
#       labs(title = title_text, x = "False Positive Rate", y = "True Positive Rate") +
#       scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
#       scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
#       geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
#       theme_minimal()
#   })
# }
# 
# # Run the app
# shinyApp(ui = ui, server = server)

#####
# Update 3/30/2026: TPR and FPR only appear on bottom and very left, each graph has x and y axes
# This version shows truncated ROC curves in all scenarios together in a 3 x 3 grid
# === Packages ===
# library(shiny)
# library(ggplot2)
# library(dplyr)
# library(patchwork)
# library(cowplot)   # for get_legend()
# 
# # === Scenario Lookup Table ===
# scenario_df <- data.frame(
#   scenario = 1:9,
#   true_es = 0.2,
#   rep.number = rep(c(2, 5, 10), times = 3),
#   rep.n = rep(c(40, 100, 400), each = 3)
# )
# 
# # === UI ===
# ui <- fluidPage(
#   tags$head(
#     tags$style(HTML("
#       .sidebar {
#         background-image: url('https://www.transparenttextures.com/patterns/cubes.png');
#       }
#     "))
#   ),
#   
#   tags$div(
#     style = "max-width: 1400px; margin: auto;",
#     titlePanel("Evaluating Meta-Analytic Bayes Factors with ROC Curves"),
#     
#     fluidRow(
#       # Sidebar
#       column(
#         width = 4,
#         tags$div(
#           class = "sidebar",
#           
#           selectInput("es_choice", "True Effect Size:", choices = c("0.2", "0.5")),
#           
#           checkboxGroupInput(
#             "methods_choice",
#             "Methods to show:",
#             choices = c(
#               "BFbMA"  = "BFbMA",
#               "EUBF"   = "EUBF",
#               "FEMABF" = "FEMABF",
#               "iBF"    = "iBF",
#               "REMA"   = "FEMA"   # display REMA, keep FEMA internally
#             ),
#             selected = c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA")
#           ),
#           
#           actionButton("plot_btn", "Generate ROC Plot"),
#           downloadButton("downloadPlot", "Download ROC Plot (PNG)"),
#           
#           tags$div(
#             style = "margin-top: 30px; padding: 15px; background-color: #f7f7f7; border: 1px solid #ddd; border-radius: 5px;",
#             HTML("<strong>Meta-Analytic Bayes Factors (MABFs):</strong><br/><br/>
#               <ul style='margin-bottom: 0; padding-left: 20px;'>
#                 <li><strong>FEMABF</strong>: Fixed-Effect Meta-Analytic Bayes Factor</li>
#                 <li><strong>BFbMA</strong>: Bayes Factor Based on Meta-Analysis</li>
#                 <li><strong>EUBF</strong>: Evidence-Updating Bayes Factor</li>
#                 <li><strong>iBF</strong>: Inclusion Bayes Factor</li>
#                 <li><strong>REMA</strong>: Random-Effects Meta-Analysis</li>
#               </ul>")
#           )
#         )
#       ),
#       
#       # Main panel
#       column(
#         width = 8,
#         tags$div(
#           class = "main",
#           tabsetPanel(
#             tabPanel(
#               "Study Overview",
#               tags$div(
#                 style = "margin-top: 20px;",
#                 HTML("<p>Recent concerns about the replication crisis in psychology have underscored the need for better frameworks to assess whether original findings can be reliably replicated. This study aims to compare the performance of several meta-analytic Bayes factors (FEMABF, BFbMA, EUBF, and iBF) both relative to one another and against the traditional random-effects meta-analysis in evaluating replication success.</p>
#                 <p>To do so, we use a large-scale simulation to generate original studies across varied scenarios by manipulating true effect size, research environment (levels of p-hacking and publication bias), and original sample size. Replication studies are then simulated by varying replication sample size and number of replications, yielding 243 distinct scenarios and over 60 million data points. This Shiny app visualizes ROC curves to compare how well each method distinguishes between true and null effects across these conditions.</p>")
#               )
#             ),
#             tabPanel(
#               "ROC Curves",
#               uiOutput("plotTitle"),
#               plotOutput("rocPlot", height = "900px")
#             )
#           )
#         )
#       )
#     )
#   )
# )
# 
# # === Server ===
# server <- function(input, output) {
#   
#   plot_event <- eventReactive(input$plot_btn, {
#     selected_es <- isolate(input$es_choice)
#     selected_methods <- isolate(input$methods_choice)
#     req(length(selected_methods) > 0)
#     
#     # fixed method colors no matter which methods are selected
#     method_colors <- c(
#       "BFbMA"  = "#CDBA96",
#       "EUBF"   = "#1E90FF",
#       "FEMABF" = "black",
#       "iBF"    = "#FF4040",
#       "FEMA"   = "#00BA38"
#     )
#     
#     # legend labels: change FEMA to REMA only at the end
#     method_labels <- c(
#       "BFbMA"  = "BFbMA",
#       "EUBF"   = "EUBF",
#       "FEMABF" = "FEMABF",
#       "iBF"    = "iBF",
#       "FEMA"   = "REMA"
#     )
#     
#     # Load corresponding RDS file
#     file_path <- switch(
#       selected_es,
#       "0.2" = "./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.2null.RDS",
#       "0.5" = "./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.5null.RDS"
#     )
#     df <- readRDS(file_path)
#     
#     # Prepare data
#     df <- df %>%
#       filter(method %in% selected_methods) %>%
#       mutate(
#         method = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA")),
#         scenario = as.numeric(scenario)
#       )
#     
#     # === Base 3x3 plots ===
#     base_plots <- lapply(1:9, function(scn) {
#       cond <- scenario_df %>% filter(scenario == scn)
#       
#       
#       mabf_methods <- intersect(c("EUBF", "FEMABF", "iBF", "BFbMA"), selected_methods)
#       all_methods  <- intersect(c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA"), selected_methods)
#       row_id <- ceiling(scn / 3)
#       col_id <- ((scn - 1) %% 3) + 1
#       x_lab <- if (row_id == 3) "False Positive Rate" else NULL
#       y_lab <- if (col_id == 1) "True Positive Rate" else NULL
#       
#       roc_mabf <- df %>%
#         filter(
#           scenario == scn,
#           method %in% mabf_methods,
#           threshold >= 3,
#           threshold <= 30
#         )
#       
#       roc_fema <- df %>%
#         filter(
#           scenario == scn,
#           method == "FEMA",
#           threshold >= 0.001,
#           threshold <= 0.05
#         )
#       
#       roc_data <- bind_rows(roc_mabf, roc_fema) %>%
#         mutate(
#           method = factor(method, levels = all_methods),
#           line_type = ifelse(method == "BFbMA", "solid", "dotted")
#         ) %>%
#         arrange(method, threshold)
#       
#       # closest available BF cutoffs
#       mabf_targets <- c(3, 10, 30)
#       cutoff_points_mabf <- roc_mabf %>%
#         group_by(method) %>%
#         group_modify(~{
#           purrr::map_dfr(mabf_targets, function(tar) {
#             .x %>%
#               slice_min(order_by = abs(threshold - tar), n = 1, with_ties = FALSE) %>%
#               mutate(cutoff_label = paste0("BF > ", tar))
#           })
#         }) %>%
#         ungroup()
#       
#       # closest available p cutoffs for FEMA/REMA
#       fema_targets <- c(0.001, 0.01, 0.05)
#       cutoff_points_fema <- roc_fema %>%
#         group_by(method) %>%
#         group_modify(~{
#           purrr::map_dfr(fema_targets, function(tar) {
#             .x %>%
#               slice_min(order_by = abs(threshold - tar), n = 1, with_ties = FALSE) %>%
#               mutate(cutoff_label = paste0("p < ", format(tar, nsmall = 2)))
#           })
#         }) %>%
#         ungroup()
#       
#       cutoff_points <- bind_rows(cutoff_points_mabf, cutoff_points_fema) %>%
#         mutate(
#           method = factor(method, levels = all_methods),
#           cutoff_label = factor(
#             cutoff_label,
#             levels = c("BF > 3", "BF > 10", "BF > 30", "p < 0.001", "p < 0.01", "p < 0.05")
#           )
#         )
#       
#       x_upper <- max(roc_data$FPR, na.rm = TRUE) * 1.05
#       x_breaks <- pretty(c(0, x_upper), n = 4)
#       x_breaks <- x_breaks[x_breaks >= 0 & x_breaks <= x_upper]
#       
#       y_lower <- max(0, min(roc_data$TPR, na.rm = TRUE) - 0.02)
#       y_upper <- min(1, max(roc_data$TPR, na.rm = TRUE) + 0.02)
#       y_breaks <- pretty(c(y_lower, y_upper), n = 4)
#       y_breaks <- y_breaks[y_breaks >= y_lower & y_breaks <= y_upper]
#       
#       ggplot(
#         roc_data,
#         aes(
#           x = FPR,
#           y = TPR,
#           color = method,
#           linetype = line_type,
#           group = method
#         )
#       ) +
#         geom_step(linewidth = 0.7) +
#         geom_point(
#           data = cutoff_points,
#           mapping = aes(x = FPR, y = TPR, color = method, shape = cutoff_label),
#           inherit.aes = FALSE,
#           size = 2.5,
#           stroke = 1
#         ) +
#         geom_abline(
#           slope = 1,
#           intercept = 0,
#           linetype = "dotted",
#           color = "black"
#         ) +
#         labs(
#           title = bquote(N[rep] == .(cond$rep.number) * "," ~ n[rep] == .(cond$rep.n)),
#           x = x_lab,
#           y = y_lab
#         ) +
#         scale_linetype_identity() +
#         scale_shape_manual(
#           values = c(
#             "BF > 3"   = 15,
#             "BF > 10"  = 17,
#             "BF > 30"  = 18,
#             "p < 0.001" = 5,
#             "p < 0.01" = 2,
#             "p < 0.05" = 0
#           )
#         ) +
#         scale_color_manual(
#           values = method_colors,
#           limits = names(method_colors),
#           breaks = intersect(names(method_colors), selected_methods),
#           labels = method_labels[intersect(names(method_colors), selected_methods)]
#         ) +
#         scale_x_continuous(
#           limits = c(0, x_upper),
#           breaks = x_breaks,
#           labels = function(x) formatC(x, format = "f", digits = 3)
#         ) +
#         scale_y_continuous(
#           limits = c(y_lower, y_upper),
#           breaks = y_breaks,
#           labels = function(y) formatC(y, format = "f", digits = 2)
#         ) +
#         theme_minimal(base_size = 11) +
#         theme(
#           plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
#           legend.position = "none",
#           axis.text.x = element_text(),
#           axis.ticks.x = element_line(),
#           axis.text.y = element_text(),
#           axis.ticks.y = element_line(),
#           axis.line = element_line(color = "black")
#         )
#     })
#     
#     # === Shared legend ===
#     legend_df_mabf <- df %>%
#       filter(method %in% c("BFbMA", "EUBF", "FEMABF", "iBF"), threshold >= 3, threshold <= 30)
#     
#     legend_df_fema <- df %>%
#       filter(method == "FEMA", threshold >= 0.001, threshold <= 0.05)
#     
#     legend_cutoff_points_mabf <- legend_df_mabf %>%
#       group_by(method) %>%
#       group_modify(~{
#         purrr::map_dfr(c(3, 10, 30), function(tar) {
#           .x %>%
#             slice_min(order_by = abs(threshold - tar), n = 1, with_ties = FALSE) %>%
#             mutate(cutoff_label = paste0("BF > ", tar))
#         })
#       }) %>%
#       ungroup()
#     
#     legend_cutoff_points_fema <- legend_df_fema %>%
#       group_by(method) %>%
#       group_modify(~{
#         purrr::map_dfr(c(0.001, 0.01, 0.05), function(tar) {
#           .x %>%
#             slice_min(order_by = abs(threshold - tar), n = 1, with_ties = FALSE) %>%
#             mutate(cutoff_label = paste0("p < ", format(tar, nsmall = 2)))
#         })
#       }) %>%
#       ungroup()
#     
#     legend_cutoff_points <- bind_rows(legend_cutoff_points_mabf, legend_cutoff_points_fema) %>%
#       mutate(
#         method = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA")),
#         cutoff_label = factor(
#           cutoff_label,
#           levels = c("BF > 3", "BF > 10", "BF > 30", "p < 0.001", "p < 0.01", "p < 0.05")
#         )
#       )
#     
#     legend_plot <- ggplot(
#       bind_rows(legend_df_mabf, legend_df_fema) %>%
#         mutate(
#           method = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA")),
#           line_type = ifelse(method == "FEMA", "dotted", "solid")
#         ),
#       aes(x = FPR, y = TPR, color = method, linetype = line_type, group = method)
#     ) +
#       geom_step(linewidth = 0.7) +
#       geom_point(
#         data = legend_cutoff_points,
#         mapping = aes(x = FPR, y = TPR, color = method, shape = cutoff_label),
#         inherit.aes = FALSE,
#         size = 2.5,
#         stroke = 1
#       ) +
#       scale_linetype_identity() +
#       scale_shape_manual(
#         values = c(
#           "BF > 3"   = 15,
#           "BF > 10"  = 17,
#           "BF > 30"  = 18,
#           "p < 0.001" = 5,
#           "p < 0.01" = 2,
#           "p < 0.05" = 0
#         )
#       ) +
#       scale_color_manual(
#         values = method_colors,
#         limits = names(method_colors),
#         breaks = intersect(names(method_colors), selected_methods),
#         labels = method_labels[intersect(names(method_colors), selected_methods)]
#       ) +
#       guides(
#         color = guide_legend(
#           order = 1,
#           override.aes = list(
#             shape = 15,
#             size = 5,
#             linetype = 0,
#             linewidth = 0
#           )
#         ),
#         shape = guide_legend(order = 2)
#       ) +
#       labs(color = "Method", shape = "Cutoff") +
#       theme_minimal() +
#       theme(
#         legend.position = "right",
#         legend.title = element_text(size = 10),
#         legend.text = element_text(size = 9)
#       )
#     
#     legend <- cowplot::get_legend(legend_plot)
#     
#     # Combine 3x3 grid and legend
#     combined_plot <- wrap_plots(base_plots, ncol = 3) +
#       plot_annotation(
#         title = paste0("ROC Curves of MABF Methods (True Effect = ", selected_es, ")"),
#         theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
#       )
#     
#     final_plot <- combined_plot | patchwork::wrap_elements(legend)
#     final_plot <- final_plot + plot_layout(widths = c(3, 0.5))
#     
#     list(
#       plot = final_plot,
#       title = paste0("ROC Curves for All Nine Scenarios (True Effect = ", selected_es, ")")
#     )
#   })
#   
#   # Render title
#   output$plotTitle <- renderUI({
#     req(plot_event())
#     h4(plot_event()$title)
#   })
#   
#   # Render plot
#   output$rocPlot <- renderPlot({
#     plot_event()$plot
#   })
#   
#   # Download handler
#   output$downloadPlot <- downloadHandler(
#     filename = function() {
#       paste0("ROC_ES", input$es_choice, "_3x3_2sided.png")
#     },
#     content = function(file) {
#       ggsave(file, plot = plot_event()$plot, width = 14, height = 10, dpi = 900, bg = "white")
#     }
#   )
# }
# 
# # === Run App ===
# shinyApp(ui = ui, server = server)


########################
# 4/20 Update: Added fixed-effect MA (FEMA2)
# === Packages ===
library(shiny)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)   # for get_legend()

# === Scenario Lookup Table ===
scenario_df <- data.frame(
  scenario = 1:9,
  true_es = 0.2,
  rep.number = rep(c(2, 5, 10), times = 3),
  rep.n = rep(c(40, 100, 400), each = 3)
)

# === UI ===
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .sidebar {
        background-image: url('https://www.transparenttextures.com/patterns/cubes.png');
      }
    "))
  ),
  
  tags$div(
    style = "max-width: 1400px; margin: auto;",
    titlePanel("Evaluating Meta-Analytic Bayes Factors with ROC Curves"),
    
    fluidRow(
      # Sidebar
      column(
        width = 4,
        tags$div(
          class = "sidebar",
          
          selectInput("es_choice", "True Effect Size:", choices = c("0.2", "0.5")),
          
          checkboxGroupInput(
            "methods_choice",
            "Methods to show:",
            choices = c(
              "BFbMA"  = "BFbMA",
              "EUBF"   = "EUBF",
              "FEMABF" = "FEMABF",
              "iBF"    = "iBF",
              "REMA"   = "FEMA",   # display REMA, keep FEMA internally
              "FEMA"  = "FEMA2"
            ),
            selected = c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA", "FEMA2")
          ),
          
          actionButton("plot_btn", "Generate ROC Plot"),
          downloadButton("downloadPlot", "Download ROC Plot (PNG)"),
          
          tags$div(
            style = "margin-top: 30px; padding: 15px; background-color: #f7f7f7; border: 1px solid #ddd; border-radius: 5px;",
            HTML("<strong>Meta-Analytic Bayes Factors (MABFs):</strong><br/><br/>
              <ul style='margin-bottom: 0; padding-left: 20px;'>
                <li><strong>FEMABF</strong>: Fixed-Effect Meta-Analytic Bayes Factor</li>
                <li><strong>BFbMA</strong>: Bayes Factor Based on Meta-Analysis</li>
                <li><strong>EUBF</strong>: Evidence-Updating Bayes Factor</li>
                <li><strong>iBF</strong>: Inclusion Bayes Factor</li>
                <li><strong>REMA</strong>: Random-Effects Meta-Analysis</li>
                <li><strong>FEMA</strong>: Fixed-Effects Meta-Analysis</li>
              </ul>")
          )
        )
      ),
      
      # Main panel
      column(
        width = 8,
        tags$div(
          class = "main",
          tabsetPanel(
            tabPanel(
              "Study Overview",
              tags$div(
                style = "margin-top: 20px;",
                HTML("<p>Recent concerns about the replication crisis in psychology have underscored the need for better frameworks to assess whether original findings can be reliably replicated. This study aims to compare the performance of several meta-analytic Bayes factors (FEMABF, BFbMA, EUBF, and iBF) both relative to one another and against the traditional random-effects meta-analysis in evaluating replication success.</p>
                <p>To do so, we use a large-scale simulation to generate original studies across varied scenarios by manipulating true effect size, research environment (levels of p-hacking and publication bias), and original sample size. Replication studies are then simulated by varying replication sample size and number of replications, yielding 243 distinct scenarios and over 60 million data points. This Shiny app visualizes ROC curves to compare how well each method distinguishes between true and null effects across these conditions.</p>")
              )
            ),
            tabPanel(
              "ROC Curves",
              uiOutput("plotTitle"),
              plotOutput("rocPlot", height = "900px")
            )
          )
        )
      )
    )
  )
)

# === Server ===
server <- function(input, output) {
  
  plot_event <- eventReactive(input$plot_btn, {
    selected_es <- isolate(input$es_choice)
    selected_methods <- isolate(input$methods_choice)
    req(length(selected_methods) > 0)
    
    # fixed method colors no matter which methods are selected
    method_colors <- c(
      "BFbMA"  = "khaki3",
      "EUBF"   = "#1E90FF",
      "FEMABF" = "black",
      "iBF"    = "#FF4040",
      "FEMA"   = "#00BA38",
      "FEMA2"  = "burlywood3"
    )
    
    # legend labels: change FEMA to REMA only at the end
    method_labels <- c(
      "BFbMA"  = "BFbMA",
      "EUBF"   = "EUBF",
      "FEMABF" = "FEMABF",
      "iBF"    = "iBF",
      "FEMA"   = "REMA",
      "FEMA2"  = "FEMA"
    )
    
    # Load corresponding RDS file
    file_path <- switch(
      selected_es,
      "0.2" = "./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.2null.RDS",
      "0.5" = "./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.5null.RDS"
    )
    df <- readRDS(file_path)
    
    # Prepare data
    df <- df %>%
      filter(method %in% selected_methods) %>%
      mutate(
        method = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA", "FEMA2")),
        scenario = as.numeric(scenario)
      )
    
    # === Base 3x3 plots ===
    base_plots <- lapply(1:9, function(scn) {
      cond <- scenario_df %>% filter(scenario == scn)
      
      
      mabf_methods <- intersect(c("EUBF", "FEMABF", "iBF", "BFbMA"), selected_methods)
      all_methods  <- intersect(c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA", "FEMA2"), selected_methods)
      row_id <- ceiling(scn / 3)
      col_id <- ((scn - 1) %% 3) + 1
      x_lab <- if (row_id == 3) "False Positive Rate" else NULL
      y_lab <- if (col_id == 1) "True Positive Rate" else NULL
      
      roc_mabf <- df %>%
        filter(
          scenario == scn,
          method %in% mabf_methods,
          threshold >= 3,
          threshold <= 30
        )
      
      roc_fema <- df %>%
        filter(
          scenario == scn,
          method %in% c("FEMA", "FEMA2"),
          threshold >= 0.001,
          threshold <= 0.05
        )
      
      roc_data <- bind_rows(roc_mabf, roc_fema) %>%
        mutate(
          method = factor(method, levels = all_methods),
          line_type = ifelse(method == "BFbMA", "solid", "dotted")
        ) %>%
        arrange(method, threshold)
      
      # closest available BF cutoffs
      mabf_targets <- c(3, 10, 30)
      cutoff_points_mabf <- roc_mabf %>%
        group_by(method) %>%
        group_modify(~{
          purrr::map_dfr(mabf_targets, function(tar) {
            .x %>%
              slice_min(order_by = abs(threshold - tar), n = 1, with_ties = FALSE) %>%
              mutate(cutoff_label = paste0("BF > ", tar))
          })
        }) %>%
        ungroup()
      
      # closest available p cutoffs for FEMA/REMA
      fema_targets <- c(0.001, 0.01, 0.05)
      cutoff_points_fema <- roc_fema %>%
        group_by(method) %>%
        group_modify(~{
          purrr::map_dfr(fema_targets, function(tar) {
            .x %>%
              slice_min(order_by = abs(threshold - tar), n = 1, with_ties = FALSE) %>%
              mutate(cutoff_label = paste0("p < ", format(tar, nsmall = 2)))
          })
        }) %>%
        ungroup()
      
      cutoff_points <- bind_rows(cutoff_points_mabf, cutoff_points_fema) %>%
        mutate(
          method = factor(method, levels = all_methods),
          cutoff_label = factor(
            cutoff_label,
            levels = c("BF > 3", "BF > 10", "BF > 30", "p < 0.001", "p < 0.01", "p < 0.05")
          )
        )
      
      x_upper <- max(roc_data$FPR, na.rm = TRUE) * 1.05
      x_breaks <- pretty(c(0, x_upper), n = 4)
      x_breaks <- x_breaks[x_breaks >= 0 & x_breaks <= x_upper]
      
      y_lower <- max(0, min(roc_data$TPR, na.rm = TRUE) - 0.02)
      y_upper <- min(1, max(roc_data$TPR, na.rm = TRUE) + 0.02)
      y_breaks <- pretty(c(y_lower, y_upper), n = 4)
      y_breaks <- y_breaks[y_breaks >= y_lower & y_breaks <= y_upper]
      
      ggplot(
        roc_data,
        aes(
          x = FPR,
          y = TPR,
          color = method,
          linetype = line_type,
          group = method
        )
      ) +
        geom_step(linewidth = 0.7) +
        geom_point(
          data = cutoff_points,
          mapping = aes(x = FPR, y = TPR, color = method, shape = cutoff_label),
          inherit.aes = FALSE,
          size = 2.5,
          stroke = 1
        ) +
        geom_abline(
          slope = 1,
          intercept = 0,
          linetype = "dotted",
          color = "black"
        ) +
        labs(
          title = bquote(N[rep] == .(cond$rep.number) * "," ~ n[rep] == .(cond$rep.n)),
          x = x_lab,
          y = y_lab
        ) +
        scale_linetype_identity() +
        scale_shape_manual(
          values = c(
            "BF > 3"   = 15,
            "BF > 10"  = 17,
            "BF > 30"  = 18,
            "p < 0.001" = 5,
            "p < 0.01" = 2,
            "p < 0.05" = 0
          )
        ) +
        scale_color_manual(
          values = method_colors,
          limits = names(method_colors),
          breaks = intersect(names(method_colors), selected_methods),
          labels = method_labels[intersect(names(method_colors), selected_methods)]
        ) +
        scale_x_continuous(
          limits = c(0, x_upper),
          breaks = x_breaks,
          labels = function(x) formatC(x, format = "f", digits = 3)
        ) +
        scale_y_continuous(
          limits = c(y_lower, y_upper),
          breaks = y_breaks,
          labels = function(y) formatC(y, format = "f", digits = 2)
        ) +
        theme_minimal(base_size = 11) +
        theme(
          plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
          legend.position = "none",
          axis.text.x = element_text(),
          axis.ticks.x = element_line(),
          axis.text.y = element_text(),
          axis.ticks.y = element_line(),
          axis.line = element_line(color = "black")
        )
    })
    
    # === Shared legend ===
    legend_df_mabf <- df %>%
      filter(method %in% c("BFbMA", "EUBF", "FEMABF", "iBF"), threshold >= 3, threshold <= 30)
    
    legend_df_fema <- df %>%
      filter(method %in% c("FEMA", "FEMA2"), threshold >= 0.001, threshold <= 0.05)
    
    legend_cutoff_points_mabf <- legend_df_mabf %>%
      group_by(method) %>%
      group_modify(~{
        purrr::map_dfr(c(3, 10, 30), function(tar) {
          .x %>%
            slice_min(order_by = abs(threshold - tar), n = 1, with_ties = FALSE) %>%
            mutate(cutoff_label = paste0("BF > ", tar))
        })
      }) %>%
      ungroup()
    
    legend_cutoff_points_fema <- legend_df_fema %>%
      group_by(method) %>%
      group_modify(~{
        purrr::map_dfr(c(0.001, 0.01, 0.05), function(tar) {
          .x %>%
            slice_min(order_by = abs(threshold - tar), n = 1, with_ties = FALSE) %>%
            mutate(cutoff_label = paste0("p < ", format(tar, nsmall = 2)))
        })
      }) %>%
      ungroup()
    
    legend_cutoff_points <- bind_rows(legend_cutoff_points_mabf, legend_cutoff_points_fema) %>%
      mutate(
        method = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA", "FEMA2")),
        cutoff_label = factor(
          cutoff_label,
          levels = c("BF > 3", "BF > 10", "BF > 30", "p < 0.001", "p < 0.01", "p < 0.05")
        )
      )
    
    legend_plot <- ggplot(
      bind_rows(legend_df_mabf, legend_df_fema) %>%
        mutate(
          method = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA", "FEMA2")),
          line_type = ifelse(method %in% c("FEMA", "FEMA2"), "dotted", "solid")
        ),
      aes(x = FPR, y = TPR, color = method, linetype = line_type, group = method)
    ) +
      geom_step(linewidth = 0.7) +
      geom_point(
        data = legend_cutoff_points,
        mapping = aes(x = FPR, y = TPR, color = method, shape = cutoff_label),
        inherit.aes = FALSE,
        size = 2.5,
        stroke = 1
      ) +
      scale_linetype_identity() +
      scale_shape_manual(
        values = c(
          "BF > 3"   = 15,
          "BF > 10"  = 17,
          "BF > 30"  = 18,
          "p < 0.001" = 5,
          "p < 0.01" = 2,
          "p < 0.05" = 0
        )
      ) +
      scale_color_manual(
        values = method_colors,
        limits = names(method_colors),
        breaks = intersect(names(method_colors), selected_methods),
        labels = method_labels[intersect(names(method_colors), selected_methods)]
      ) +
      guides(
        color = guide_legend(
          order = 1,
          override.aes = list(
            shape = 15,
            size = 5,
            linetype = 0,
            linewidth = 0
          )
        ),
        shape = guide_legend(order = 2)
      ) +
      labs(color = "Method", shape = "Cutoff") +
      theme_minimal() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )
    
    legend <- cowplot::get_legend(legend_plot)
    
    # Combine 3x3 grid and legend
    combined_plot <- wrap_plots(base_plots, ncol = 3) +
      plot_annotation(
        title = paste0("ROC Curves of MABF Methods (True Effect = ", selected_es, ")"),
        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
      )
    
    final_plot <- combined_plot | patchwork::wrap_elements(legend)
    final_plot <- final_plot + plot_layout(widths = c(3, 0.5))
    
    list(
      plot = final_plot,
      title = paste0("ROC Curves for All Nine Scenarios (True Effect = ", selected_es, ")")
    )
  })
  
  # Render title
  output$plotTitle <- renderUI({
    req(plot_event())
    h4(plot_event()$title)
  })
  
  # Render plot
  output$rocPlot <- renderPlot({
    plot_event()$plot
  })
  
  # Download handler
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("ROC_ES", input$es_choice, "_3x3_2sided.png")
    },
    content = function(file) {
      ggsave(file, plot = plot_event()$plot, width = 14, height = 10, dpi = 300, bg = "white")
    }
  )
}

# === Run App ===
shinyApp(ui = ui, server = server)
