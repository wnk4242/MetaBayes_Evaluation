library(shiny)
library(shinyjs)
library(tidyverse)
library(glmmTMB)
library(lmtest)
library(car)
library(emmeans)
library(httr)
library(jsonlite)
library(ggplot2)
library(ggeffects)

# Load RDS data once
folder_path <- "./MABFanalyses/row-wise"
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
for (file_path in rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

# Modeling function
analyze_MABF_effects <- function(.method, BFcutoff, true.effect, dv, ivs, interactions = NULL) {
  tpr_rds <- sprintf("MABF_rates_%.1f_rowise_BFcutoff%d", true.effect, BFcutoff)
  fpr_rds <- sprintf("MABF_rates_null_rowise_BFcutoff%d", BFcutoff)
  data_tpr <- get(tpr_rds, envir = .GlobalEnv)
  data_fpr <- get(fpr_rds, envir = .GlobalEnv)
  
  data_tpr <- data_tpr %>% select(method = 1, true.effect = 3, orig.n = 4, orig.alpha = 5,
                                  bias.level = 6, rep.number = 9, rep.n = 10,
                                  ADRt = 11, TPR = 23, TSR1 = 29)
  data_fpr <- data_fpr %>% select(method = 1, true.effect = 3, orig.n = 4, orig.alpha = 5,
                                  bias.level = 6, rep.number = 9, rep.n = 10,
                                  ADRn = 11, FPR = 24, TSR2 = 29, FSR2 = 30)
  
  to_factor <- c("method", "true.effect", "orig.n", "orig.alpha", "bias.level", "rep.number", "rep.n")
  data_tpr <- data_tpr %>% mutate(across(all_of(to_factor), as.factor))
  data_fpr <- data_fpr %>% mutate(across(all_of(to_factor), as.factor))
  
  data_tpr <- data_tpr %>% filter(method == .method)
  data_fpr <- data_fpr %>% filter(method == .method)
  
  data_tpr <- data_tpr %>% mutate(scenario_id = rep(1:162, each = 500))
  data_fpr <- data_fpr %>% mutate(scenario_id = rep(1:162, each = 500))
  
  factor_levels <- list(rep.n = c("40", "100", "400"),
                        rep.number = c("2", "5", "10"),
                        orig.n = c("20", "50", "200"),
                        orig.alpha = c("0.01", "0.05"),
                        bias.level = c("low", "medium", "high"))
  for (var in names(factor_levels)) {
    data_tpr[[var]] <- factor(data_tpr[[var]], levels = factor_levels[[var]])
    data_fpr[[var]] <- factor(data_fpr[[var]], levels = factor_levels[[var]])
  }
  
  data_combined <- data_tpr %>%
    mutate(FPR = data_fpr$FPR,
           ADRn = data_fpr$ADRn,
           TSR2 = data_fpr$TSR2,
           FSR2 = data_fpr$FSR2)
  
  epsilon <- 1e-6
  data_combined <- data_combined %>%
    mutate(across(c(TPR, FPR, ADRt, ADRn, TSR1, TSR2, FSR2),
                  ~ ifelse(. <= 0, epsilon, ifelse(. >= 1, 1 - epsilon, .))))
  
  main_formula <- paste(ivs, collapse = " + ")
  if (!is.null(interactions) && length(interactions) > 0) {
    inter_formula <- paste(interactions, collapse = " + ")
    full_rhs <- paste(main_formula, inter_formula, sep = " + ")
  } else {
    full_rhs <- main_formula
  }
  full_formula <- reformulate(full_rhs, response = dv)
  
  full_model <- glmmTMB(full_formula, data = data_combined, family = beta_family(link = "logit"))
  wald_test <- Anova(full_model, type = 3)
  drop1_result <- drop1(full_model, test = "Chisq")
  
  terms_to_summarize <- unique(c(ivs, interactions))
  emmeans_results <- list()
  for (term in terms_to_summarize) {
    tryCatch({
      emmeans_results[[term]] <- emmeans(full_model, specs = as.formula(paste("~", term)))
    }, error = function(e) {
      emmeans_results[[term]] <- paste("Error in emmeans for:", term)
    })
  }
  
  list(
    Wald_Type3 = wald_test,
    Drop1_Table = drop1_result,
    EMMeans = emmeans_results,
    Full_Model = summary(full_model),
    Model_Object = full_model
  )
}

# --- UI ---
ui <- fluidPage(
  useShinyjs(),
  tags$head(
    tags$style(HTML("

/* Sidebar styling */
.sidebar {
  background-image: url('https://www.transparenttextures.com/patterns/cubes.png');
  background-color: #ffffff;
  padding: 15px;
  border-radius: 5px;
  position: fixed;
  top: 100px;
  left: 0;
  width: 400px;
  max-height: 85vh;
  overflow-y: auto;
  z-index: 999;
}

/* Invisible scrollbar */
.sidebar::-webkit-scrollbar {
  width: 10px;
}
.sidebar::-webkit-scrollbar-track {
  background: transparent;
}
.sidebar::-webkit-scrollbar-thumb {
  background-color: transparent;
  border: none;
}

/* Main panel layout */
.main-panel {
  margin-left: 420px;
  padding: 20px;
}
.tab-content {
  background-color: white;
  padding: 20px;
  border-radius: 5px;
}

/* Output containers */
#gpt_output,
#beta_output,
#lr_output,
#full_model_output,
#wald_output,
#drop1_output,
#emmeans_output {
  white-space: pre-wrap;
  overflow-y: auto;         /* enables vertical scroll */
  max-height: 800px;        /* height limit for scroll to appear */
  font-family: monospace;
  padding-right: 10px;      /* leaves room for scrollbar */
  border: 1px solid #ddd;   /* make it visible even when empty */
  min-height: 150px;        /* always show a box with height */
}

/* Transparent scrollbars for output areas */
    #gpt_output::-webkit-scrollbar,
    #beta_output::-webkit-scrollbar,
    #lr_output::-webkit-scrollbar,
    #full_model_output::-webkit-scrollbar,
    #wald_output::-webkit-scrollbar,
    #drop1_output::-webkit-scrollbar,
    #emmeans_output::-webkit-scrollbar {
      width: 10px;
    }

    #gpt_output::-webkit-scrollbar-track,
    #beta_output::-webkit-scrollbar-track,
    #lr_output::-webkit-scrollbar-track,
    #full_model_output::-webkit-scrollbar-track,
    #wald_output::-webkit-scrollbar-track,
    #drop1_output::-webkit-scrollbar-track,
    #emmeans_output::-webkit-scrollbar-track {
      background: transparent;
    }

    #gpt_output::-webkit-scrollbar-thumb,
    #beta_output::-webkit-scrollbar-thumb,
    #lr_output::-webkit-scrollbar-thumb,
    #full_model_output::-webkit-scrollbar-thumb,
    #wald_output::-webkit-scrollbar-thumb,
    #drop1_output::-webkit-scrollbar-thumb,
    #emmeans_output::-webkit-scrollbar-thumb {
      background-color: transparent;
      border: none;
    }

html, body {
  overflow: hidden;
  height: 100%;
}

/* --- Tab styling with texture and raised active tab --- */
.nav-tabs {
border-bottom: 1px solid transparent;    /* Light gray line (now is transparent) */
  margin-bottom: 0px;
  padding-bottom: 0px;              /* Add space between tabs and line */
}


.nav-tabs > li > a {
  background-image: url('https://www.transparenttextures.com/patterns/cubes.png');
  border: 1px solid #dcd7d7;          /* consistent border width */
  background-color: rgba(255, 255, 255, 0.6);
  background-blend-mode: multiply;
  
  color: black;
  font-weight: bold;               /* consistent font weight */
  padding: 10px 20px;
  margin-right: 0px;               /* tighter spacing */
  border-radius: 5px;
  box-shadow: none;
  position: relative;
  top: 0;                          /* fix base position */
  transition: all 0.15s ease-in-out;
}

#plot_box {
  background-color: white;
  border: 1px solid #ccc;
  padding: 10px;
  border-radius: 4px;
  min-height: 400px;
  box-sizing: border-box;
  overflow: hidden;
  display: flex;
  align-items: center;
  justify-content: center;
}


/* Active tab - visually raised upward */
.nav-tabs > li.active > a,
.nav-tabs > li.active > a:focus,
.nav-tabs > li.active > a:hover {
  background-image: url('https://www.transparenttextures.com/patterns/cubes.png');
  background-color: rgba(255, 255, 255, 1);
  background-blend-mode: multiply;
  border: 1px solid #999;
  color: black;
  box-shadow: none;
  
  z-index: 10;
}
"))
  ),
  
  titlePanel("Beta-Regression Model and visualization"),
  
  sidebarLayout(
    sidebarPanel(
      class = "sidebar",
      selectInput("method", "MABF Method:", choices = c("FEMABF", "EUBF", "BFbMA", "iBF")),
      selectInput("bfcutoff", "Bayes Factor Cutoff:", choices = c(1, 3, 10, 30)),
      selectInput("effectsize", "True Effect Size:", choices = c(0.2, 0.5)),
      selectInput("dv", "Outcome Measure:", choices = c("TPR", "FPR", "ADRt", "ADRn", "TSR1", "TSR2", "FSR2"), selected = "TPR"),
      checkboxGroupInput("main_effects", "Explanatory Variables:",
                         choices = c(
                           "Replication Sample Size" = "rep.n",
                           "Number of Replications" = "rep.number"
                         ),
                         selected = c("rep.n", "rep.number")
      ),
      
      checkboxGroupInput("covariates", "Covariates:",
                         choices = c(
                           "Level of Bias" = "bias.level",
                           "Original Study Sample Size" = "orig.n",
                           "Original Study α Level" = "orig.alpha"
                         )
      ),
      uiOutput("interaction_ui"),
      actionButton("run", "Run Analysis"),
      hr(),
      fileInput("apikey", "Upload OpenAI API Key (txt)"),
      textAreaInput("prompt", "ChatGPT Prompt", 
                    value = "TPR is true positive rate. FPR is false positive rate. rep.n is the sample size of replication studies. rep.number is the number of replication studies. 
                             Please interpret these statistical results (including main and interaction effects) in APA style."),
      selectInput("gpt_model", "ChatGPT Model:", choices = c("gpt-3.5-turbo", "gpt-4")),
      actionButton("askgpt", "Ask ChatGPT")
    ),
    
    mainPanel(
      class = "main-panel",
      tabsetPanel(
        tabPanel("Study Overview",
                 tags$div(style = "margin-top: 20px;",
                          HTML("
      <p>Recent concerns about the replication crisis in psychology have underscored the need for better frameworks to assess whether original findings can be reliably replicated. This study aims to compare the performance of several meta-analytic Bayes factors (FEMABF, BFbMA, EUBF, and iBF) both relative to one another and against the traditional fixed-effect meta-analysis in evaluating replication success.</p>
      <p>To do so, we use a large-scale simulation to generate original studies across varied scenarios by manipulating true effect size, research environment (levels of p-hacking and publication bias), and original sample size. Replication studies are then simulated by varying replication sample size and number of replications, yielding 243 distinct scenarios and over 60 million data points.</p>
      <p>This Shiny app explores how design features of replication studies and original studies impact key outcome measures of meta-analytic Bayes factor methods (e.g., True Positive Rate and False Positive Rate). </p>
    ")
                 )
        ),
        
        tabPanel("Full Model Summary",
                 fluidRow(
                   column(9,
                          div(style = "padding: 10px;",
                              p("Displays the complete output of the beta regression model, including estimated coefficients, standard errors, z-values, and significance levels for both the mean and precision (phi) submodels."),
                              verbatimTextOutput("full_model_output")
                          )
                   ),
                   column(3,
                          div(style = "padding: 10px; border-left: 1px solid transparent;",  #transparent vertical line
                              textAreaInput("fullmodel_notes", "Notes:", value = "", width = "100%", height = "250px")
                          )
                   )
                 )
        ),
        
        tabPanel("Wald Type III Tests",
                 fluidRow(
                   column(9,
                          div(style = "padding: 10px;",
                              p("Presents Type III Wald chi-square tests for each main effect and interaction term, allowing you to assess the significance of each predictor while accounting for all other terms in the model."),
                              verbatimTextOutput("wald_output")
                          )
                   ),
                   column(3,
                          div(style = "padding: 10px; border-left: 1px solid transparent;",
                              textAreaInput("wald_notes", "Notes:", value = "", width = "100%", height = "250px")
                          )
                   )
                 )
        ),
        
        tabPanel("Drop1 LRT",
                 fluidRow(
                   column(9,
                          div(style = "padding: 10px;",
                              p("Provides likelihood ratio tests by sequentially dropping one term at a time from the full model, useful for evaluating how much each term contributes to model fit."),
                              verbatimTextOutput("drop1_output")
                          )
                   ),
                   column(3,
                          div(style = "padding: 10px; border-left: 1px solid transparent;",
                              textAreaInput("drop1_notes", "Notes:", value = "", width = "100%", height = "250px")
                          )
                   )
                 )
        ),
        
        tabPanel("EMMeans",
                 fluidRow(
                   column(9,
                          div(style = "padding: 10px;",
                              p("Shows estimated marginal means (also known as least-squares means), which summarize the predicted values of the outcome variable for each level of the predictors, averaging over other factors."),
                              verbatimTextOutput("emmeans_output")
                          )
                   ),
                   column(3,
                          div(style = "padding: 10px; border-left: 1px solid transparent;",
                              textAreaInput("emmeans_notes", "Notes:", value = "", width = "100%", height = "250px")
                          )
                   )
                 )
        ),
        
        tabPanel("Visualization",
                 fluidRow(
                   column(9,
                          div(style = "padding: 10px;",
                              p("Displays a visual summary of the predicted values from the beta regression model. The plot shows how the outcome variable varies across levels of the main predictors, with optional faceting by covariates or interactions to help identify patterns and effect differences."),
                              div(id = "plot_box",
                                  plotOutput("plot_output", height = "600px")
                              ),
                              div(style = "text-align: right;",
                                  downloadButton("download_plot", "Download Plot")
                              )
                          )
                   ),
                   column(3,
                          div(style = "padding: 10px; border-left: 1px solid transparent;",
                              textAreaInput("visual_notes", "Notes:", value = "", width = "100%", height = "250px")
                          )
                   )
                 )
        ),
        
        tabPanel("ChatGPT Interpretation",
                 fluidRow(
                   column(9,
                          div(style = "padding: 10px;",
                              p("Generates a plain-language interpretation of the statistical results using ChatGPT. It summarizes the significance of effects, interactions, and estimated marginal means in APA style, based on the full model output and your custom prompt."),
                              verbatimTextOutput("gpt_output")
                          )
                   ),
                   column(3,
                          div(style = "padding: 10px; border-left: 1px solid transparent;",
                              textAreaInput("gpt_notes", "Notes:", value = "", width = "100%", height = "250px")
                          )
                   )
                 )
        )
        
      )
    )
    
  )
)


# --- Server ---
server <- function(input, output, session) {
  # Track the most recent plot
  last_plot <- reactiveVal(NULL)
  
  # Run model only when "Run Analysis" is clicked
  results <- eventReactive(input$run, {
    analyze_MABF_effects(
      .method = input$method,
      BFcutoff = as.numeric(input$bfcutoff),
      true.effect = as.numeric(input$effectsize),
      dv = input$dv,
      ivs = c(input$main_effects, input$covariates),
      interactions = input$interactions
    )
  })
  
  # Full model summary
  output$full_model_output <- renderPrint({
    tryCatch({
      req(results())
      print(results()$Full_Model)
    }, error = function(e) {
      cat("Click Run Analysis to show results.")
    })
  })
  
  # Wald Type III
  output$wald_output <- renderPrint({
    tryCatch({
      req(results())
      print(results()$Wald_Type3)
    }, error = function(e) {
      cat("Click Run Analysis to show results.")
    })
  })
  
  # Drop1 LRT
  output$drop1_output <- renderPrint({
    tryCatch({
      req(results())
      print(results()$Drop1_Table)
    }, error = function(e) {
      cat("Click Run Analysis to show results.")
    })
  })
  
  # EMMeans
  output$emmeans_output <- renderPrint({
    tryCatch({
      req(results())
      all_emm <- lapply(names(results()$EMMeans), function(term) {
        paste0("Estimated Marginal Means for ", term, ":\n",
               paste(capture.output(print(results()$EMMeans[[term]])), collapse = "\n"),
               "\n==============================\n")
      })
      cat(paste(all_emm, collapse = "\n\n"))
    }, error = function(e) {
      cat("Click Run Analysis to show results.")
    })
  })
  
  
  # Dynamically generate all 2-way interaction terms
  output$interaction_ui <- renderUI({
    ivs_selected <- c(input$main_effects, input$covariates)
    
    if (length(ivs_selected) >= 2) {
      label_lookup <- c(
        "rep.n" = "Replication Sample Size",
        "rep.number" = "Number of Replications",
        "bias.level" = "Level of Bias",
        "orig.n" = "Original Study Sample Size",
        "orig.alpha" = "Original Study α Level"
      )
      
      interaction_choices <- combn(ivs_selected, 2, FUN = function(x) {
        value <- paste(x, collapse = ":")
        label <- paste(label_lookup[x], collapse = " × ")
        setNames(value, label)
      }, simplify = FALSE)
      
      choices_named <- setNames(
        sapply(interaction_choices, identity),
        sapply(interaction_choices, names)
      )
    } else {
      choices_named <- character(0)
    }
    
    checkboxGroupInput(
      inputId = "interactions",
      label = "Interaction Terms (optional):",
      choices = choices_named
    )
  })
  
  # Prepare plot input data only after Run button
  plot_data <- eventReactive(input$run, {
    req(results())
    list(
      model = results()$Model_Object,
      data = results()$data_used,
      dv = isolate(input$dv),
      covars = isolate(input$covariates)
    )
  })
  
  output$plot_output <- renderPlot({
    tryCatch({
      pd <- plot_data()
      req(pd)
      
      model <- pd$model
      dv <- pd$dv
      covars <- pd$covars
      
      selected_interaction <- isolate(input$interactions)
      covar_interactions <- selected_interaction[grepl(":", selected_interaction)]
      valid_covars <- c("bias.level", "orig.n", "orig.alpha")
      
      var_labels <- c(
        "rep.n" = "Replication Sample Size",
        "rep.number" = "Number of Replications",
        "bias.level" = "Level of Bias",
        "orig.n" = "Original Sample Size",
        "orig.alpha" = "Original α Level"
      )
      
      if (length(covar_interactions) == 1) {
        interaction_vars <- strsplit(covar_interactions, ":")[[1]]
      } else {
        interaction_vars <- character(0)
      }
      use_grid <- length(interaction_vars) == 2 && all(interaction_vars %in% valid_covars)
      
      label_map <- list(
        bias.level = c("low" = "Bias = Low", "medium" = "Bias = Medium", "high" = "Bias = High"),
        orig.n = c("20" = "Original Sample Size = 20", "50" = "Original Sample Size = 50", "200" = "Original Sample Size = 200"),
        orig.alpha = c("0.01" = "α = .01", "0.05" = "α = .05")
      )
      
      plot_obj <- if (use_grid) {
        row_var <- interaction_vars[1]
        col_var <- interaction_vars[2]
        
        preds <- ggpredict(model, terms = c("rep.n", "rep.number", col_var, row_var))
        
        ggplot(preds, aes(x = x, y = predicted, color = group)) +
          geom_point(size = 2) +
          geom_line(aes(group = group), linewidth = 1) +
          facet_grid(facet ~ panel,
                     labeller = labeller(
                       facet = function(x) paste(var_labels[[col_var]], "=", x),
                       panel = function(x) paste(var_labels[[row_var]], "=", x)
                     )) +
          labs(
            x = "Replication Sample Size",
            y = paste("Predicted", dv),
            color = "Number of Replications",
            title = paste("Predicted", dv, "by",
                          paste0(var_labels["rep.n"], ","),
                          paste0(var_labels["rep.number"], ", and"),
                          var_labels[[row_var]], "×", var_labels[[col_var]])
          ) +
          theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            strip.background = element_rect(fill = "white", color = "black"),
            plot.margin = margin(10, 10, 10, 10),
            text = element_text(size = 10),
            legend.key.size = unit(0.7, "lines")
          )
        
      } else if (length(covars) == 0) {
        preds <- ggpredict(model, terms = c("rep.n", "rep.number"))
        
        ggplot(preds, aes(x = x, y = predicted, color = group)) +
          geom_point(size = 2) +
          geom_line(aes(group = group), linewidth = 1) +
          labs(
            x = "Replication Sample Size",
            y = paste("Predicted", dv),
            color = "Number of Replications",
            title = paste("Predicted", dv, "by",
                          paste0(var_labels["rep.n"], " and"),
                          var_labels["rep.number"])
          ) +
          theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            strip.background = element_rect(fill = "white", color = "black"),
            plot.margin = margin(10, 10, 10, 10),
            text = element_text(size = 10),
            legend.key.size = unit(0.7, "lines")
          )
        
      } else {
        plot_list <- lapply(covars, function(facet_var) {
          preds <- ggpredict(model, terms = c("rep.n", "rep.number", facet_var))
          
          ggplot(preds, aes(x = x, y = predicted, color = group)) +
            geom_point(size = 2) +
            geom_line(aes(group = group), linewidth = 1) +
            facet_wrap(~ facet,
                       labeller = labeller(facet = function(x) paste(var_labels[[facet_var]], "=", x))) +
            labs(
              x = "Replication Sample Size",
              y = paste("Predicted", dv),
              color = "Number of Replications",
              title = paste("Predicted", dv, "by",
                            paste0(var_labels["rep.n"], ","),
                            paste0(var_labels["rep.number"], ", and"),
                            var_labels[[facet_var]])
            ) +
            theme(
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              strip.background = element_rect(fill = "white", color = "black"),
              plot.margin = margin(10, 10, 10, 10),
              text = element_text(size = 10),
              legend.key.size = unit(0.7, "lines")
            )
        })
        
        patchwork::wrap_plots(plot_list, ncol = 1)
      }
      
      last_plot(plot_obj)
      plot_obj
    }, error = function(e) {
      # Placeholder before "Run Analysis"
      ggplot() +
        theme_void() +
        labs(title = "Click Run Analysis to show plot.") +
        theme(
          plot.title = element_text(hjust = 0, size = 14, color = "gray50"),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
        )
    })
  })
  
  
  # Download plot
  output$download_plot <- downloadHandler(
    filename = function() {
      method <- input$method
      bf <- input$bfcutoff
      effect <- input$effectsize
      dv <- input$dv
      
      covars <- if (length(input$covariates) > 0) {
        paste(input$covariates, collapse = "&")
      } else {
        "NoCovars"
      }
      
      inters <- if (length(input$interactions) > 0) {
        paste(gsub(":", "x", input$interactions), collapse = "&")
      } else {
        "NoInter"
      }
      
      sprintf("Predicted%s_%s_BF%s_Eff%s_Covars-%s_Inter-%s.png",
              dv, method, bf, effect, covars, inters)
    },
    content = function(file) {
      ggsave(file, plot = last_plot(), device = "png", width = 10, height = 7, dpi = 300)
    }
  )
  
  # ChatGPT Response
  # Show placeholder before any interaction
  output$gpt_output <- renderPrint({
    cat("First upload API then click Ask ChatGPT.")
  })
  observeEvent(input$askgpt, {
    req(input$apikey)
    api_key <- tryCatch({
      readLines(input$apikey$datapath)[1]
    }, error = function(e) return(NULL))
    
    if (is.null(api_key)) {
      output$gpt_output <- renderPrint({ cat("Failed to read API key.") })
      return()
    }
    
    gpt_text <- paste(
      "Wald Type III Test:\n", paste(capture.output(print(results()$Wald_Type3)), collapse = "\n"), "\n\n",
      "Drop1 Likelihood Ratio Test:\n", paste(capture.output(print(results()$Drop1_Table)), collapse = "\n"), "\n\n",
      "Estimated Marginal Means:\n",
      paste(
        lapply(names(results()$EMMeans), function(term) {
          paste0("EMMeans for ", term, ":\n",
                 paste(capture.output(print(results()$EMMeans[[term]])), collapse = "\n"))
        }), collapse = "\n\n"
      )
    )
    
    response <- POST(
      url = "https://api.openai.com/v1/chat/completions",
      add_headers(
        Authorization = paste("Bearer", api_key),
        `Content-Type` = "application/json"
      ),
      body = toJSON(list(
        model = input$gpt_model,
        messages = list(
          list(role = "user", content = gpt_text),
          list(role = "user", content = input$prompt)
        )
      ), auto_unbox = TRUE)
    )
    
    raw_response <- content(response, as = "text", encoding = "UTF-8")
    parsed <- fromJSON(raw_response, simplifyVector = FALSE)
    
    output$gpt_output <- renderPrint({
      if (!is.null(parsed$choices)) {
        cat(parsed$choices[[1]]$message$content)
      } else {
        cat("No interpretation returned. Check API key or prompt.")
      }
    })
  })
}

# Run app
shinyApp(ui, server)
