#In this script, we will visualize replication TPR and FPR of MABF methods using heatmaps
#Each heat map visualize all 81 scenarios for every MABF method

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(scales)


# Set the path to the directory containing RDS files
folder_path <- "./MABF4ROC/4TPFP" 
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

# False Positive Rate
# Manually create heatmaps
summary_counts <- map_dfr(names(FEMABF_lists_0.2null_regrouped), function(dataset_name) {
  dataset <- FEMABF_lists_0.2null_regrouped[[dataset_name]]
  
  map_dfr(names(dataset), function(setting_name) {
    mat <- dataset[[setting_name]][1:500, ]
    count_above_1 <- sum(mat > 1)
    
    tibble(
      dataset = dataset_name,
      setting = setting_name,
      count = count_above_1
    )
  })
}) %>%
  separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
  mutate(
    rep_num = factor(rep_num, levels = c(10, 5, 2)),
    rep_n = factor(rep_n, levels = c(40, 100, 400)),
    
    # Extract bias and original sample size from dataset name
    orig_n = case_when(
      str_detect(dataset, "_20_") ~ "20",
      str_detect(dataset, "_50_") ~ "50",
      str_detect(dataset, "_200_") ~ "200"
    ),
    bias_level = case_when(
      str_detect(dataset, "none_low") ~ "Low",
      str_detect(dataset, "medium_medium") ~ "Medium",
      str_detect(dataset, "high_high") ~ "High"
    ),
    
    orig_n = factor(orig_n, levels = c("20", "50", "200")),
    bias_level = factor(bias_level, levels = c("Low", "Medium", "High"))
  )

summary_counts <- summary_counts %>%
  mutate(
    false_positive_rate = count / 250000
  )

ggplot(summary_counts, aes(x = rep_n, y = rep_num, fill = false_positive_rate)) +
  geom_tile(color = "white") +
  facet_grid(bias_level ~ orig_n, labeller = label_both) +
  scale_fill_gradient(
    low = "#fee5d9", high = "#a50f15",
    name = "% FEMABF > 1",
    labels = scales::percent_format(accuracy = 0.1),
    breaks = pretty(summary_counts$false_positive_rate, n = 5)
  ) + geom_text(aes(label = count), color = "black", size = 3) +
  labs(
    title = "False Positive Rate (FEMABF > 1) in Null Rows",
    subtitle = "Each tile = a scenario; 250,000 values per scenario",
    x = "Replication Sample Size (rep_n)",
    y = "# of Replications (rep_num)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 11, face = "bold")
  )



# Helper function creates heatmaps
FPR_heatmap <- function(method_name = "FEMABF", cutoff = 1) {
  
  # Dynamically get the object by name
  data_list <- get(paste0(method_name, "_lists_0.2null_regrouped"))
  
  summary_counts <- map_dfr(names(data_list), function(dataset_name) {
    dataset <- data_list[[dataset_name]]
    
    map_dfr(names(dataset), function(setting_name) {
      mat <- dataset[[setting_name]][1:500, ]
      count_above_cutoff <- sum(mat > cutoff)
      
      tibble(
        method = method_name,
        dataset = dataset_name,
        setting = setting_name,
        count = count_above_cutoff
      )
    })
  }) %>%
    separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
    mutate(
      rep_num = factor(rep_num, levels = c(10, 5, 2)),
      rep_n = factor(rep_n, levels = c(40, 100, 400)),
      orig_n = case_when(
        str_detect(dataset, "_20_") ~ "20",
        str_detect(dataset, "_50_") ~ "50",
        str_detect(dataset, "_200_") ~ "200"
      ),
      bias_level = case_when(
        str_detect(dataset, "none_low") ~ "Low",
        str_detect(dataset, "medium_medium") ~ "Medium",
        str_detect(dataset, "high_high") ~ "High"
      ),
      orig_n = factor(orig_n, levels = c("20", "50", "200")),
      bias_level = factor(bias_level, levels = c("Low", "Medium", "High")),
      false_positive_rate = count / 250000
    )
  
  ggplot(summary_counts, aes(x = rep_n, y = rep_num, fill = false_positive_rate)) +
    geom_tile(color = "white") +
    #geom_text(aes(label = count), color = "black", size = 3) + #absolute count
    geom_text(aes(label = scales::percent(false_positive_rate, accuracy = 0.1)), color = "black", size = 3) +
    facet_grid(bias_level ~ orig_n, labeller = label_both) +
    scale_fill_gradient(
      low = "#fee5d9", high = "#a50f15",
      name = paste0("% ", method_name, " > ", cutoff),
      labels = percent_format(accuracy = 0.1),
      breaks = pretty(summary_counts$false_positive_rate, n = 5)
    ) +
    labs(
      title = paste0("False Positive Rate (", method_name, " > ", cutoff, ") in Null Rows"),
      subtitle = "Each tile = 250,000 FEMABF values per scenario (500 rows × 500 columns)",
      x = "Replication Sample Size (rep_n)",
      y = "# of Replications (rep_num)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(size = 11, face = "bold")
    )
}

FPR_heatmap("FEMABF", cutoff = 3)


# Collapsed FPR heatmap
FPR_heatmap_collapsed <- function(method_name = "FEMABF", cutoff = 1) {
  
  data_list <- get(paste0(method_name, "_lists_0.2null_regrouped"))
  
  summary_counts <- map_dfr(names(data_list), function(dataset_name) {
    dataset <- data_list[[dataset_name]]
    
    map_dfr(names(dataset), function(setting_name) {
      mat <- dataset[[setting_name]][1:500, ]
      count_above_cutoff <- sum(mat > cutoff)
      
      tibble(
        method = method_name,
        dataset = dataset_name,
        setting = setting_name,
        count = count_above_cutoff
      )
    })
  }) %>%
    separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
    mutate(
      rep_num = factor(rep_num, levels = c(10, 5, 2)),
      rep_n = factor(rep_n, levels = c(40, 100, 400))
    ) %>%
    group_by(rep_num, rep_n) %>%
    summarize(
      total_count = sum(count),
      n_scenarios = n(),  # number of datasets contributing
      total_possible = n_scenarios * 250000,
      proportion = total_count / total_possible,
      .groups = "drop"
    )
  
  ggplot(summary_counts, aes(x = rep_n, y = rep_num, fill = proportion)) +
    geom_tile(color = "white") +
    geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), color = "black", size = 4) +
    scale_fill_gradient(
      low = "#fee5d9", high = "#a50f09",
      name = paste0("% ", method_name, " > ", cutoff),
      labels = percent_format(accuracy = 0.1),
      breaks = pretty(summary_counts$proportion, n = 5)
    ) +
    labs(
      title = paste0("Collapsed False Positive Rate (", method_name, " > ", cutoff, ")"),
      subtitle = "Proportion across all datasets with same rep_num and rep_n",
      x = "Replication Sample Size (rep_n)",
      y = "# of Replications (rep_num)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(size = 11, face = "bold")
    )
}
FPR_heatmap_collapsed("FEMABF", cutoff = 3)




# True Positive Rate
# Manually create heatmaps
summary_counts <- map_dfr(names(FEMABF_lists_0.2null_regrouped), function(dataset_name) {
  dataset <- FEMABF_lists_0.2null_regrouped[[dataset_name]]
  
  map_dfr(names(dataset), function(setting_name) {
    mat <- dataset[[setting_name]][501:1000, ]  # TRUE EFFECT rows
    count_above_1 <- sum(mat > 1)
    
    tibble(
      dataset = dataset_name,
      setting = setting_name,
      count = count_above_1
    )
  })
}) %>%
  separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
  mutate(
    rep_num = factor(rep_num, levels = c(10, 5, 2)),
    rep_n = factor(rep_n, levels = c(40, 100, 400)),
    orig_n = case_when(
      str_detect(dataset, "_20_") ~ "20",
      str_detect(dataset, "_50_") ~ "50",
      str_detect(dataset, "_200_") ~ "200"
    ),
    bias_level = case_when(
      str_detect(dataset, "none_low") ~ "Low",
      str_detect(dataset, "medium_medium") ~ "Medium",
      str_detect(dataset, "high_high") ~ "High"
    ),
    orig_n = factor(orig_n, levels = c("20", "50", "200")),
    bias_level = factor(bias_level, levels = c("Low", "Medium", "High")),
    true_positive_rate = count / 250000
  )

ggplot(summary_counts, aes(x = rep_n, y = rep_num, fill = true_positive_rate)) +
  geom_tile(color = "white") +
  #geom_text(aes(label = count), color = "black", size = 3) + #absolute count
  geom_text(aes(label = scales::percent(true_positive_rate, accuracy = 0.1)), color = "black", size = 3)+
  facet_grid(bias_level ~ orig_n, labeller = label_both) +
  scale_fill_gradient(
    low = "#edf8e9", high = "#006d2c",
    name = "% FEMABF > 1",
    labels = scales::percent_format(accuracy = 0.1),
    breaks = pretty(summary_counts$true_positive_rate, n = 5)
  ) +
  labs(
    title = "True Positive Rate (FEMABF > 1) in Rows with True Effect",
    subtitle = "Each tile = 250,000 values (500 true-effect rows × 500 columns)",
    x = "Replication Sample Size (rep_n)",
    y = "# of Replications (rep_num)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 11, face = "bold")
  )
# Helper function to visualize TPR
  TPR_heatmap <- function(method_name = "FEMABF", cutoff = 1) {
    
    # Dynamically get the object by name
    data_list <- get(paste0(method_name, "_lists_0.2null_regrouped"))
    
    summary_counts <- purrr::map_dfr(names(data_list), function(dataset_name) {
      dataset <- data_list[[dataset_name]]
      
      purrr::map_dfr(names(dataset), function(setting_name) {
        mat <- dataset[[setting_name]][501:1000, ]  # rows with true effects
        count_above_cutoff <- sum(mat > cutoff)
        
        tibble::tibble(
          method = method_name,
          dataset = dataset_name,
          setting = setting_name,
          count = count_above_cutoff
        )
      })
    }) %>%
      tidyr::separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
      dplyr::mutate(
        rep_num = factor(rep_num, levels = c(10, 5, 2)),
        rep_n = factor(rep_n, levels = c(40, 100, 400)),
        orig_n = dplyr::case_when(
          stringr::str_detect(dataset, "_20_") ~ "20",
          stringr::str_detect(dataset, "_50_") ~ "50",
          stringr::str_detect(dataset, "_200_") ~ "200"
        ),
        bias_level = dplyr::case_when(
          stringr::str_detect(dataset, "none_low") ~ "Low",
          stringr::str_detect(dataset, "medium_medium") ~ "Medium",
          stringr::str_detect(dataset, "high_high") ~ "High"
        ),
        orig_n = factor(orig_n, levels = c("20", "50", "200")),
        bias_level = factor(bias_level, levels = c("Low", "Medium", "High")),
        true_positive_rate = count / 250000
      )
    
    ggplot2::ggplot(summary_counts, ggplot2::aes(x = rep_n, y = rep_num, fill = true_positive_rate)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(
        aes(label = scales::percent(true_positive_rate, accuracy = 0.1)),
        color = "black",
        size = 3
      ) +
      ggplot2::facet_grid(bias_level ~ orig_n, labeller = ggplot2::label_both) +
      ggplot2::scale_fill_gradient(
        low = "#edf8e9", high = "#006d2c",
        name = paste0("% ", method_name, " > ", cutoff),
        labels = scales::percent_format(accuracy = 0.1),
        breaks = pretty(summary_counts$true_positive_rate, n = 5)
      ) +
      ggplot2::labs(
        title = paste0("True Positive Rate (", method_name, " > ", cutoff, ") in Rows with True Effect"),
        subtitle = "Each tile = 250,000 values (500 true-effect rows × 500 columns)",
        x = "Replication Sample Size (rep_n)",
        y = "# of Replications (rep_num)"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(size = 11, face = "bold")
      )
  }
TPR_heatmap("FEMABF", 1)

# Collapsed TPR heatmap
TPR_heatmap_collapsed <- function(method_name = "FEMABF", cutoff = 1) {
  data_list <- get(paste0(method_name, "_lists_0.2null_regrouped"))
  
  summary_counts <- purrr::map_dfr(names(data_list), function(dataset_name) {
    dataset <- data_list[[dataset_name]]
    
    purrr::map_dfr(names(dataset), function(setting_name) {
      mat <- dataset[[setting_name]][501:1000, ]  # rows with true effects
      count_above_cutoff <- sum(mat > cutoff)
      
      tibble::tibble(
        method = method_name,
        dataset = dataset_name,
        setting = setting_name,
        count = count_above_cutoff
      )
    })
  }) %>%
    tidyr::separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
    dplyr::mutate(
      rep_num = factor(rep_num, levels = c(10, 5, 2)),
      rep_n = factor(rep_n, levels = c(40, 100, 400))
    ) %>%
    group_by(rep_num, rep_n) %>%
    summarise(
      total_count = sum(count),
      n_scenarios = n(),
      total_possible = n_scenarios * 250000,
      proportion = total_count / total_possible,
      .groups = "drop"
    )
  
  ggplot(summary_counts, aes(x = rep_n, y = rep_num, fill = proportion)) +
    geom_tile(color = "white") +
    geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), color = "black", size = 4) +
    scale_fill_gradient(
      low = "#edf8e9", high = "#006d2c",
      name = paste0("% ", method_name, " > ", cutoff),
      labels = percent_format(accuracy = 0.1),
      breaks = pretty(summary_counts$proportion, n = 5)
    ) +
    labs(
      title = paste0("Collapsed True Positive Rate (", method_name, " > ", cutoff, ")"),
      subtitle = "Proportion across all datasets with same rep_num and rep_n",
      x = "Replication Sample Size (rep_n)",
      y = "# of Replications (rep_num)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(size = 11, face = "bold")
    )
}

TPR_heatmap_collapsed("FEMABF", cutoff = 1)
