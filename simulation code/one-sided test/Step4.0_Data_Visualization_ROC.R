# NOTE: You need to take a look at the ROC Shiny App you built. You cut data from groups 10 to 81.
# Load necessary packages
library(ggplot2)
library(tidyverse)
library(purrr)
rm(list = ls()) # housekeeping

# Set the path to the directory containing RDS files
folder_path <- "./MABFanalyses/matrix-wise/rates4ROC"
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

# Combine rate outcomes of different methods in order to compare them in ROC plot
rates_0.2null <- rbind(EUBF_TPFP4ROC_0.2null, FEMABF_TPFP4ROC_0.2null, iBF_TPFP4ROC_0.2null, BFbMA_TPFP4ROC_0.2null, FEMA_TPFP4ROC_0.2null, FEMA2_TPFP4ROC_0.2null) 
# Convert character to numeric
rates_0.2null <- rates_0.2null %>%
  mutate(across(c(`rep number`,rep.n), as.numeric)) %>% 
  mutate(true_es = 0.2) %>% 
  relocate(true_es, .after = `true effect`) %>% 
  relocate(group, .before = `method`) %>% 
  rename(scenario = group)

saveRDS(rates_0.2null, file = "./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.2null.RDS")

# Combine rate outcomes of different methods in order to compare them in ROC plot
rates_0.5null <- rbind(EUBF_TPFP4ROC_0.5null, FEMABF_TPFP4ROC_0.5null, iBF_TPFP4ROC_0.5null, BFbMA_TPFP4ROC_0.5null, FEMA_TPFP4ROC_0.5null, FEMA2_TPFP4ROC_0.5null) 
# Convert character to numeric
rates_0.5null <- rates_0.5null %>%
  mutate(across(c(`rep number`,rep.n), as.numeric)) %>% 
  mutate(true_es = 0.5) %>% 
  relocate(true_es, .after = `true effect`) %>% 
  relocate(group, .before = `method`) %>% 
  rename(scenario = group)

saveRDS(rates_0.5null, file = "./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.5null.RDS")



# Summarize simulation scenarios showing group number and associated variables
scenario_df <- data.frame(
  scenario = c(1:81),
  true_es = 0.2,
  orig.n = rep(c(20, 50, 200), each=9, times=3),
  rep.number = rep(c(2, 5, 10), times = 27),
  rep.n = rep(c(40,100,400), each = 3, times = 9),
  bias = rep(c('none','medium','high'), each=27)
)
scenario_df

# Function used to create ROC curves
plot_ROC <- function(scenario_number) {
  condition <- scenario_df %>% filter(scenario == scenario_number)
  
  scenario_number <- condition$scenario[1]
  true_es <- condition$true_es[1]
  orig_n <- condition$orig.n[1]
  rep_number <- condition$rep.number[1]
  rep_n <- condition$rep.n[1]
  bias <- condition$bias[1]
  
  print(paste("Scenario:", scenario_number,
              "Replications:", rep_number,
              "Sample Size:", rep_n))
  
  mabf_methods <- c("EUBF", "FEMABF", "iBF", "BFbMA")
  all_methods  <- c("BFbMA", "EUBF", "FEMABF", "iBF", "FEMA")
  
  # ----------------------------
  # Method-specific threshold ranges
  # ----------------------------
  roc_mabf <- rates_0.2null %>%
    filter(
      scenario == as.character(scenario_number),
      method %in% mabf_methods,
      threshold >= 3,
      threshold <= 30
    )
  
  roc_fema <- rates_0.2null %>%
    filter(
      scenario == as.character(scenario_number),
      method == "FEMA",
      threshold >= 0.001,
      threshold <= 0.05
    )
  
  roc_data <- bind_rows(roc_mabf, roc_fema) %>%
    mutate(
      method = factor(method, levels = all_methods),
      line_type = ifelse(method == "FEMA", "dotted", "solid")
    ) %>%
    arrange(method, threshold)
  
  # ----------------------------
  # Closest available cutoffs for MABF
  # ----------------------------
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
  
  # ----------------------------
  # Closest available cutoffs for FEMA
  # ----------------------------
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
  
  # ----------------------------
  # Dynamic axis ranges
  # ----------------------------
  x_upper <- max(roc_data$FPR, na.rm = TRUE) * 1.05
  x_breaks <- pretty(c(0, x_upper), n = 5)
  x_breaks <- x_breaks[x_breaks >= 0 & x_breaks <= x_upper]
  
  y_lower <- max(0, min(roc_data$TPR, na.rm = TRUE) - 0.02)
  y_upper <- min(1, max(roc_data$TPR, na.rm = TRUE) + 0.02)
  y_breaks <- pretty(c(y_lower, y_upper), n = 5)
  y_breaks <- y_breaks[y_breaks >= y_lower & y_breaks <= y_upper]
  
  title_text <- paste0(
    "ROC Curves with method-specific cutoffs when the number of replications is ",
    rep_number,
    ", replication sample size is ",
    rep_n,
    ", true effect is ",
    true_es,
    ", and bias level is ",
    bias
  )
  
  roc_plot <- ggplot(
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
      size = 3,
      stroke = 1
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dotted",
      color = "black"
    ) +
    labs(
      title = title_text,
      x = "False Positive Rate",
      y = "True Positive Rate",
      shape = "Cutoff"
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
      values = c(
        "BFbMA"  = "#EEA9B8",
        "EUBF"   = "#1E90FF",
        "FEMABF" = "black",
        "iBF"    = "#F8766D",
        "FEMA"   = "#00BA38"
      ),
      labels = c(
        "BFbMA"  = "BFbMA",
        "EUBF"   = "EUBF",
        "FEMABF" = "FEMABF",
        "iBF"    = "iBF",
        "FEMA"   = "REMA"
      )
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          shape = 15,
          size = 5,
          linetype = 0,
          linewidth = 0
        )
      )
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
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  file_name <- paste0(
    "./plots/ROC/0.2_ROC_group/",
    scenario_number, "_", orig_n, "_", rep_number, "_", rep_n, "_", bias,
    "_method_specific_cutoffs.jpg"
  )
  
  ggsave(file_name, plot = roc_plot, width = 12, height = 8, dpi = 300)
  print(roc_plot)
}

# Automatically generate ROC plots and save to designated folder
# No need to run all 81 conditions because QRP and PB do not affect FPR and TPR
# Just run the conditions under no QPR or PB scenarios
walk(1, plot_ROC)
# You can also use a for loop to do the same thing of course
# for (group_number in 1:2) {
#   plot_ROC(group_number)
# }

