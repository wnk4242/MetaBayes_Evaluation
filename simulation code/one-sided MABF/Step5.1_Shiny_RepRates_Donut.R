# This is the Shiny App for the Step4.0_Data_Visualization_RepRates_Donut
# We use the donut chart to showcase the relationship between FPR/TPR and FSR/TSR
library(shiny)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(scales)

# Load RDS files once on startup
folder_path <- "./MABF4ROC/4TSFS"
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
for (file_path in rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

# ---- UI ----
ui <- fluidPage(
  titlePanel("Bayesian Replication Classification Donut"),
  sidebarLayout(
    sidebarPanel(
      selectInput("method_name", "MABF Method",
                  choices = c("FEMABF", "BFbMA", "EUBF", "iBF")),
      selectInput("effect_size", "Effect Size",
                  choices = c("0", "0.2", "0.5")),
      selectInput("orig_n", "Original Sample Size",
                  choices = c(20, 50, 200), selected = 20),
      selectInput("bias_level", "Bias Level",
                  choices = c("low", "medium", "high")),
      selectInput("rep_number", "Number of Replications",
                  choices = c(2, 5, 10), selected = 2),
      selectInput("rep_n", "Replication Sample Size",
                  choices = c(40, 100, 400), selected = 40),
      numericInput("cutoff", "BF Cutoff", value = 1, min = 1),
      selectInput("pcutoff", "Original Alpha Level", choices = c(0.01, 0.05, 0.1), selected = 0.05),
      checkboxInput("show_adr", "Display Anecdotal Cases", value = FALSE),
      radioButtons("label_type", "Label Type", choices = c("count", "percent"), inline = TRUE)
    ),
    mainPanel(
      plotOutput("donutPlot", height = "700px")
    )
  )
)


# ---- Server ----
server <- function(input, output) {
  output$donutPlot <- renderPlot({
    bias_map <- c("low" = "none_low", "medium" = "medium_medium", "high" = "high_high")
    bias_str <- bias_map[[input$bias_level]]
    
    effect_key <- if (input$effect_size == "0") "0.2" else input$effect_size
    data_list <- get(paste0(input$method_name, "_lists_", effect_key, "null_regrouped_deltap"))
    dataset_name <- paste0("0,", effect_key, "_", input$orig_n, "_", bias_str)
    setting_name <- paste0(input$rep_number, "_", input$rep_n)
    
    # Use first or second half of data depending on effect size == 0
    df <- if (input$effect_size == "0") {
      data_list[[dataset_name]][[setting_name]][1:500, ] %>% as.data.frame()
    } else {
      data_list[[dataset_name]][[setting_name]][501:1000, ] %>% as.data.frame()
    }
    
    cutoff <- input$cutoff
    pcutoff <- input$pcutoff
    
    if (input$effect_size == "0") {
      # Null true effect
      FS2 <- df %>% filter(original_p < pcutoff) %>% select(4:ncol(.)) %>% as.matrix() %>% { sum(. > cutoff) }
      TF2 <- df %>% filter(original_p > pcutoff) %>% select(4:ncol(.)) %>% as.matrix() %>% { sum(. > cutoff) }
      FF2 <- df %>% filter(original_p < pcutoff) %>% select(4:ncol(.)) %>% as.matrix() %>% { sum(. < 1 / cutoff) }
      TS2 <- df %>% filter(original_p > pcutoff) %>% select(4:ncol(.)) %>% as.matrix() %>% { sum(. < 1 / cutoff) }
      AD <- df %>% select(4:ncol(.)) %>% as.matrix() %>% { sum((. > 1 / cutoff) & (. < cutoff)) }
      
      if (input$show_adr) {
        group1_total <- FS2 + TS2
        group2_total <- TF2 + FF2
        total <- group1_total + group2_total + AD
      } else {
        group1_total <- FS2 + TS2
        group2_total <- TF2 + FF2
        total <- FS2 + TF2 + FF2 + TS2
      }
      
      FP <- FS2 + TF2
      TN <- FF2 + TS2
      FPR <- FP / total
      TNR <- TN / total
      FSR2 <- FS2 / group1_total
      TSR2 <- TS2 / group1_total
      TFR2 <- TF2 / group2_total
      FFR2 <- FF2 / group2_total
      ADR <- AD / total
      
      inner_df <- if (input$show_adr) {
        tibble(
          Category = c("FPR", "TNR", "ADR"),
          Count = c(FP, TN, AD),
          x = 2,
          Ring = "Inner",
          Label = if (input$label_type == "count") {
            c(paste0("FP: ", comma(FP)), paste0("TN: ", comma(TN)), paste0("AD: ", comma(AD)))
          } else {
            c(paste0("FPR (", percent(FPR, accuracy = 0.1), ")"), paste0("TNR (", percent(TNR, accuracy = 0.1), ")"), paste0("ADR (", percent(ADR, accuracy = 0.1), ")"))
          }
        )
      } else {
        tibble(
          Category = c("FPR", "TNR"),
          Count = c(FP, TN),
          x = 2,
          Ring = "Inner",
          Label = if (input$label_type == "count") {
            c(paste0("FP: ", comma(FP)), paste0("TN: ", comma(TN)))
          } else {
            c(paste0("FPR (", percent(FPR, accuracy = 0.1), ")"), paste0("TNR (", percent(TNR, accuracy = 0.1), ")"))
          }
        )
      }
      
      outer_df <- if (input$show_adr) {
        tibble(
          Category = c("FSR2", "TFR2", "FFR2", "TSR2", "ADR"),
          Count = c(FS2, TF2, FF2, TS2, AD),
          x = 3,
          Ring = "Outer",
          Label = if (input$label_type == "count") {
            c(paste0("FS2: ", comma(FS2)), paste0("TF2: ", comma(TF2)), paste0("FF2: ", comma(FF2)), paste0("TS2: ", comma(TS2)), "")
          } else {
            c(paste0("FSR2: ", percent(FSR2, accuracy = 0.1)), paste0("TFR2: ", percent(TFR2, accuracy = 0.1)),
              paste0("FFR2: ", percent(FFR2, accuracy = 0.1)), paste0("TSR2: ", percent(TSR2, accuracy = 0.1)), "")
          }
        )
      } else {
        tibble(
          Category = c("FSR2", "TFR2", "FFR2", "TSR2"),
          Count = c(FS2, TF2, FF2, TS2),
          x = 3,
          Ring = "Outer",
          Label = if (input$label_type == "count") {
            c(paste0("FS2: ", comma(FS2)), paste0("TF2: ", comma(TF2)), paste0("FF2: ", comma(FF2)), paste0("TS2: ", comma(TS2)))
          } else {
            c(paste0("FSR2: ", percent(FSR2, accuracy = 0.1)), paste0("TFR2: ", percent(TFR2, accuracy = 0.1)),
              paste0("FFR2: ", percent(FFR2, accuracy = 0.1)), paste0("TSR2: ", percent(TSR2, accuracy = 0.1)))
          }
        )
      }
      
    } else {
      # True effect present
      FS1 <- df %>% filter(original_p > pcutoff) %>% select(4:ncol(.)) %>% as.matrix() %>% { sum(. < 1 / cutoff) }
      TF1 <- df %>% filter(original_p < pcutoff ) %>% select(4:ncol(.)) %>% as.matrix() %>% { sum(. < 1 / cutoff) }
      FF1 <- df %>% filter(original_p > pcutoff) %>% select(4:ncol(.)) %>% as.matrix() %>% { sum(. > cutoff) }
      TS1 <- df %>% filter(original_p < pcutoff ) %>% select(4:ncol(.)) %>% as.matrix() %>% { sum(. > cutoff) }
      AD <- df %>% select(4:ncol(.)) %>% as.matrix() %>% { sum((. > 1 / cutoff) & (. < cutoff)) }
      
      if (input$show_adr) {
        group1_total <- FS1 + TS1
        group2_total <- TF1 + FF1
        total <- group1_total + group2_total + AD
      } else {
        group1_total <- FS1 + TS1
        group2_total <- TF1 + FF1
        total <- FS1 + TF1 + FF1 + TS1
      }
      
      FN <- FS1 + TF1
      TP <- FF1 + TS1
      FNR <- FN / total
      TPR <- TP / total
      FSR1 <- FS1 / group1_total
      TSR1 <- TS1 / group1_total
      TFR1 <- TF1 / group2_total
      FFR1 <- FF1 / group2_total
      ADR <- AD / total
      
      inner_df <- if (input$show_adr) {
        tibble(
          Category = c("FNR", "TPR", "ADR"),
          Count = c(FN, TP, AD),
          x = 2,
          Ring = "Inner",
          Label = if (input$label_type == "count") {
            c(paste0("FN: ", comma(FN)), paste0("TP: ", comma(TP)), paste0("AD: ", comma(AD)))
          } else {
            c(paste0("FNR (", percent(FNR, accuracy = 0.1), ")"), paste0("TPR (", percent(TPR, accuracy = 0.1), ")"), paste0("ADR (", percent(ADR, accuracy = 0.1), ")"))
          }
        )
      } else {
        tibble(
          Category = c("FNR", "TPR"),
          Count = c(FN, TP),
          x = 2,
          Ring = "Inner",
          Label = if (input$label_type == "count") {
            c(paste0("FN: ", comma(FN)), paste0("TP: ", comma(TP)))
          } else {
            c(paste0("FNR (", percent(FNR, accuracy = 0.1), ")"), paste0("TPR (", percent(TPR, accuracy = 0.1), ")"))
          }
        )
      }
      
      outer_df <- if (input$show_adr) {
        tibble(
          Category = c("FSR1", "TFR1", "FFR1", "TSR1", "ADR"),
          Count = c(FS1, TF1, FF1, TS1, AD),
          x = 3,
          Ring = "Outer",
          Label = if (input$label_type == "count") {
            c(paste0("FS1: ", comma(FS1)), paste0("TF1: ", comma(TF1)), paste0("FF1: ", comma(FF1)), paste0("TS1: ", comma(TS1)), "")
          } else {
            c(paste0("FSR1: ", percent(FSR1, accuracy = 0.1)), paste0("TFR1: ", percent(TFR1, accuracy = 0.1)),
              paste0("FFR1: ", percent(FFR1, accuracy = 0.1)), paste0("TSR1: ", percent(TSR1, accuracy = 0.1)), "")
          }
        )
      } else {
        tibble(
          Category = c("FSR1", "TFR1", "FFR1", "TSR1"),
          Count = c(FS1, TF1, FF1, TS1),
          x = 3,
          Ring = "Outer",
          Label = if (input$label_type == "count") {
            c(paste0("FS1: ", comma(FS1)), paste0("TF1: ", comma(TF1)), paste0("FF1: ", comma(FF1)), paste0("TS1: ", comma(TS1)))
          } else {
            c(paste0("FSR1: ", percent(FSR1, accuracy = 0.1)), paste0("TFR1: ", percent(TFR1, accuracy = 0.1)),
              paste0("FFR1: ", percent(FFR1, accuracy = 0.1)), paste0("TSR1: ", percent(TSR1, accuracy = 0.1)))
          }
        )
      }
    }
    
    all_df <- bind_rows(inner_df, outer_df) %>%
      mutate(fill_color = case_when(
        Category %in% c("FPR", "FNR")  ~ "#f9938b",
        Category %in% c("TNR", "TPR")  ~ "#5ec9cc",
        grepl("FSR|FS1", Category) ~ "#fbd4d4",
        grepl("TFR|TF1", Category) ~ "#f5b7b1",
        grepl("FFR", Category)     ~ "#b6e3e3",
        grepl("TSR", Category)     ~ "#d4f3f5",
        Category == "ADR"          ~ "lightyellow"
      ))
    
    # Set consistent stacking order to align inner and outer rings
    all_df$Category <- factor(
      all_df$Category,
      levels = if (input$effect_size == "0") {
        if (input$show_adr) {
          c("FSR2", "TFR2", "FPR", "FFR2", "TSR2", "TNR", "ADR")
        } else {
          c("FSR2", "TFR2", "FPR", "FFR2", "TSR2", "TNR")
        }
      } else {
        if (input$show_adr) {
          c("FSR1", "TFR1", "FNR", "FFR1", "TSR1", "TPR", "ADR")
        } else {
          c("FSR1", "TFR1", "FNR", "FFR1", "TSR1", "TPR")
        }
      }
    )
    
    # Arrange the rows in matching stacking order
    all_df <- all_df %>% arrange(match(Category, levels(Category)))
    
    
    
    
    
    ggplot(all_df, aes(x = x, y = Count, fill = Category)) +
      geom_col(width = 1, color = "gray95", show.legend = FALSE) +
      coord_polar(theta = "y") +
      geom_text(
        data = all_df %>% filter(Ring == "Inner"),
        aes(label = Label),
        position = position_stack(vjust = 0.5),
        size = 4,
        fontface = "bold"
      ) +
      geom_label_repel(
        data = all_df %>% filter(Ring == "Outer") %>% mutate(x = 3.5),
        aes(x = x, label = Label),
        position = position_stack(vjust = 0.1),
        size = 3,
        box.padding = 0.5,
        segment.size = 0.7,
        direction = "y",
        force = 9,
        show.legend = FALSE
      ) +
      scale_fill_manual(values = setNames(all_df$fill_color, all_df$Category)) +
      xlim(0.5, 3.5) +
      theme_void(base_size = 20) +
      annotate("text", x = 0.5, y = 15, label = "Classification", size = 4, fontface = "bold")
  })
}


# ---- Run ----
shinyApp(ui = ui, server = server)
