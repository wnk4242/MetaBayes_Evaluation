#This script uses donut chart to showcase the relationship between FPR/TPR and FSR/TSR

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(scales)
library(gridExtra)

# Set the path to the directory containing RDS files
folder_path <- "./MABF4ROC/4TSFS" 
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}




rates_donut <- function(method_name = "FEMABF",
                        effect_size = "0",
                        orig_n = 20,
                        bias_level = "low",
                        rep_number = 2,
                        rep_n = 40,
                        cutoff = 1,
                        pcutoff = 0.05,
                        label_type = c("count", "percent")) {
  
  label_type <- match.arg(label_type)
  
  # Bias level mapping
  bias_map <- c("low" = "none_low", "medium" = "medium_medium", "high" = "high_high")
  bias_str <- bias_map[[bias_level]]
  
  # Load data
  dataset_name <- paste0("0,0.2_", orig_n, "_", bias_str)
  setting_name <- paste0(rep_number, "_", rep_n)
  effect_key <- if (effect_size == "0") "0.2" else effect_size
  data_list <- get(paste0(method_name, "_lists_", effect_key, "null_regrouped_deltap"))
  df <- data_list[[dataset_name]][[setting_name]][1:500, ] %>% as.data.frame()
  
  # Classify
  FS2 <- df %>% filter(original_p < pcutoff, observed_es > 0) %>% select(4:ncol(.)) %>% as.matrix() %>% { sum(. > cutoff) }
  TF2 <- df %>% filter(original_p > pcutoff) %>% select(4:ncol(.)) %>% as.matrix() %>% { sum(. > cutoff) }
  FF2 <- df %>% filter(original_p < pcutoff, observed_es > 0) %>% select(4:ncol(.)) %>% as.matrix() %>% { sum(. < 1 / cutoff) }
  TS2 <- df %>% filter(original_p > pcutoff) %>% select(4:ncol(.)) %>% as.matrix() %>% { sum(. < 1 / cutoff) }
  
  group1_total <- FS2 + TS2
  group2_total <- TF2 + FF2
  total <- FS2 + TF2 + FF2 + TS2
  FP <- FS2 + TF2
  TN <- FF2 + TS2
  FPR <- FP / total
  TNR <- TN / total
  FSR2 <- FS2 / group1_total
  TSR2 <- TS2 / group1_total
  TFR2 <- TF2 / group2_total
  FFR2 <- FF2 / group2_total
  
  # Inner ring (FPR, TNR)
  # Inner ring (FPR, TNR)
  inner_df <- tibble(
    Category = c("FPR", "TNR"),
    Count = c(FP, TN),
    x = 2,
    Ring = "Inner",
    Label = case_when(
      label_type == "count" ~ c(
        paste0("FP: ", comma(FP)),
        paste0("TN: ", comma(TN))
      ),
      label_type == "percent" ~ c(
        paste0("FPR (", percent(FPR, accuracy = 0.1), ")"),
        paste0("TNR (", percent(TNR, accuracy = 0.1), ")")
      )
    )
  )
  
  
  # Outer ring (FSR2 etc.)
  outer_df <- tibble(
    Category = c("FSR2", "TFR2", "FFR2", "TSR2"),
    Count = c(FS2, TF2, FF2, TS2),
    x = 3,
    Ring = "Outer",
    Label = case_when(
      label_type == "count" ~ c(
        paste0("FS2: ", comma(FS2)),
        paste0("TF2: ", comma(TF2)),
        paste0("FF2: ", comma(FF2)),
        paste0("TS2: ", comma(TS2))
      ),
      label_type == "percent" ~ c(
        paste0("FSR2: ", percent(FSR2, accuracy = 0.1)),
        paste0("TFR2: ", percent(TFR2, accuracy = 0.1)),
        paste0("FFR2: ", percent(FFR2, accuracy = 0.1)),
        paste0("TSR2: ", percent(TSR2, accuracy = 0.1))
      )
    )
  )
  
  # Combine
  all_df <- bind_rows(inner_df, outer_df) %>%
    mutate(
      fill_color = case_when(
        Category == "FPR" ~ "#da7166",
        Category == "TNR" ~ "#5ec9cc",
        Category == "FSR2" ~ "#f1948a",
        Category == "TSR2" ~ "#d4f3f5",
        Category == "TFR2" ~ "#f5b7b1",
        Category == "FFR2" ~ "#b6e3e3"
      )
    )
  
  all_df$Category <- factor(all_df$Category,
                            levels = c("FPR", "TNR", "FSR2", "TFR2", "FFR2", "TSR2"))
  
  # Plot
  p <- ggplot(all_df, aes(x = x, y = Count, fill = Category)) +
    geom_col(width = 1, color = "white", show.legend = FALSE) +
    coord_polar(theta = "y") +
    
    # FPR & TNR inner ring centered with position_stack
    geom_text(
      data = all_df %>% filter(Ring == "Inner"),
      aes(label = Label),
      position = position_stack(vjust = 0.5),
      size = 3,
      fontface = "bold"
    ) +
    
    # Outer ring labels repelled
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
    annotate("text", x = 0.5, y = 15, label = "Classification", size = 3.5, fontface = "bold")
  
  print(p)
}


rates_donut(method_name = "FEMABF",
            effect_size = "0",
            orig_n = 20,
            bias_level = "low",
            rep_number = 2,
            rep_n = 40,
            cutoff = 1,
            pcutoff = 0.05,
            label_type = "percent")

