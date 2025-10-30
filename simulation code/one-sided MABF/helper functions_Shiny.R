# These helper functions are used in Step5.1_Shiny_HeatMap
####### False positive rate heatmap ######
#cutoff is BFcutoff, which is Bayes factor cutoff
# 9 x 9 grid
FPR_heatmap <- function(method_name = "FEMABF", cutoff = 1, label_type = "percent", effect_size= "0") {
  # Use 0.2 data for FPR by default, even if effect_size == "0"
  effect_key <- if (effect_size == "0") "0.2" else effect_size
  data_list <- get(paste0(method_name, "_lists_", effect_key, "null_regrouped_deltap"))
  
  summary_counts <- purrr::map_dfr(names(data_list), function(dataset_name) {
    dataset <- data_list[[dataset_name]]
    purrr::map_dfr(names(dataset), function(setting_name) {
      mat <- dataset[[setting_name]][1:500, 4:503]
      count <- sum(mat > cutoff)
      tibble(method = method_name, dataset = dataset_name, setting = setting_name, count = count)
    })
  }) %>%
    tidyr::separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
    mutate(
      rep_num = factor(rep_num, levels = c(10, 5, 2)),
      rep_n   = factor(rep_n,   levels = c(40, 100, 400)),
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
      rate = count / 250000
    )
  
  p <- ggplot(summary_counts, aes(x = rep_n, y = rep_num, fill = rate)) +
    geom_tile(color = "white") +
    facet_grid(bias_level ~ orig_n, labeller = label_both) +
    scale_fill_gradient(
      low = "#fee5d9", high = "#a50f15",
      name = paste0("% ", method_name, " > ", cutoff),
      labels = percent_format(accuracy = 0.1),
      breaks = pretty(summary_counts$rate, n = 5)
    ) +
    labs(
      title = paste0("False Positive Rate (", method_name, " > ", cutoff, ")"),
      subtitle = "Each tile = 250,000 null values",
      x = "Replication Sample Size (rep_n)",
      y = "# of Replications (rep_num)"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  
  if (label_type == "percent") {
    p <- p + geom_text(aes(label = scales::percent(rate, accuracy = 0.1)), size = 3)
  } else if (label_type == "count") {
    p <- p + geom_text(aes(label = count), size = 3)
  }
  
  return(p)
}
#FPR_heatmap(method_name = "FEMABF", cutoff = 1, label_type = "percent", effect_size= "0")
# 3 x 3 grid
FPR_heatmap_collapsed <- function(method_name = "FEMABF", cutoff = 1, label_type = "percent", effect_size= "0") {
  # Use 0.2 data for FPR by default, even if effect_size == "0"
  effect_key <- if (effect_size == "0") "0.2" else effect_size
  data_list <- get(paste0(method_name, "_lists_", effect_key, "null_regrouped_deltap"))
  
  summary <- purrr::map_dfr(names(data_list), function(dataset_name) {
    dataset <- data_list[[dataset_name]]
    purrr::map_dfr(names(dataset), function(setting_name) {
      mat <- dataset[[setting_name]][1:500, 4:503]
      count <- sum(mat > cutoff)
      tibble(dataset = dataset_name, setting = setting_name, count = count)
    })
  }) %>%
    tidyr::separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
    group_by(rep_num, rep_n) %>%
    summarise(
      total_count = sum(count),
      n = n(),
      total_possible = n * 250000,
      rate = total_count / total_possible,
      .groups = "drop"
    ) %>%
    mutate(
      rep_num = factor(rep_num, levels = c(10, 5, 2)),
      rep_n = factor(rep_n, levels = c(40, 100, 400))
    )
  
  p <- ggplot(summary, aes(x = rep_n, y = rep_num, fill = rate)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      low = "#fee5d9", high = "#a50f15",
      name = paste0("% ", method_name, " > ", cutoff),
      labels = percent_format(accuracy = 0.1),
      breaks = pretty(summary$rate, n = 5)
    ) +
    labs(
      title = paste0("Collapsed FPR (", method_name, " > ", cutoff, ")"),
      subtitle = "Collapsed across all null-effect datasets",
      x = "Replication Sample Size",
      y = "# of Replications"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  
  if (label_type == "percent") {
    p <- p + geom_text(aes(label = scales::percent(rate, accuracy = 0.1)), size = 4)
  } else if (label_type == "count") {
    p <- p + geom_text(aes(label = total_count), size = 4)
  }
  
  return(p)
}
#FPR_heatmap_collapsed(method_name = "FEMABF", cutoff = 1, label_type = "percent", effect_size= "0")

####### True positive rate heatmap ######
# 9 x 9 grid
TPR_heatmap <- function(method_name = "FEMABF", cutoff = 1, label_type = "percent", effect_size= "0") {
  # Use 0.2 data for FPR by default, even if effect_size == "0"
  effect_key <- if (effect_size == "0") "0.2" else effect_size
  data_list <- get(paste0(method_name, "_lists_", effect_key, "null_regrouped_deltap"))
  
  summary_counts <- purrr::map_dfr(names(data_list), function(dataset_name) {
    dataset <- data_list[[dataset_name]]
    purrr::map_dfr(names(dataset), function(setting_name) {
      mat <- dataset[[setting_name]][501:1000, 4:503]
      count <- sum(mat > cutoff)
      tibble(method = method_name, dataset = dataset_name, setting = setting_name, count = count)
    })
  }) %>%
    tidyr::separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
    mutate(
      rep_num = factor(rep_num, levels = c(10, 5, 2)),
      rep_n   = factor(rep_n,   levels = c(40, 100, 400)),
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
      rate = count / 250000
    )
  
  p <- ggplot(summary_counts, aes(x = rep_n, y = rep_num, fill = rate)) +
    geom_tile(color = "white") +
    facet_grid(bias_level ~ orig_n, labeller = label_both) +
    scale_fill_gradient(
      low = "#edf8e9", high = "#006d2c",
      name = paste0("% ", method_name, " > ", cutoff),
      labels = percent_format(accuracy = 0.1),
      breaks = pretty(summary_counts$rate, n = 5)
    ) +
    labs(
      title = paste0("True Positive Rate (", method_name, " > ", cutoff, ")"),
      subtitle = "Each tile = 250,000 true-effect values",
      x = "Replication Sample Size (rep_n)",
      y = "# of Replications (rep_num)"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  
  if (label_type == "percent") {
    p <- p + geom_text(aes(label = scales::percent(rate, accuracy = 0.1)), size = 3)
  } else if (label_type == "count") {
    p <- p + geom_text(aes(label = count), size = 3)
  }
  
  return(p)
}
#TPR_heatmap(method_name = "FEMABF", cutoff = 1, label_type = "percent", effect_size= "0")

# 3 x 3 grid
TPR_heatmap_collapsed <- function(method_name = "FEMABF", cutoff = 1, label_type = "percent", effect_size= "0") {
  # Use 0.2 data for FPR by default, even if effect_size == "0"
  effect_key <- if (effect_size == "0") "0.2" else effect_size
  data_list <- get(paste0(method_name, "_lists_", effect_key, "null_regrouped_deltap"))
  
  summary <- purrr::map_dfr(names(data_list), function(dataset_name) {
    dataset <- data_list[[dataset_name]]
    purrr::map_dfr(names(dataset), function(setting_name) {
      mat <- dataset[[setting_name]][501:1000, 4:503]
      count <- sum(mat > cutoff)
      tibble(dataset = dataset_name, setting = setting_name, count = count)
    })
  }) %>%
    tidyr::separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
    group_by(rep_num, rep_n) %>%
    summarise(
      total_count = sum(count),
      n = n(),
      total_possible = n * 250000,
      rate = total_count / total_possible,
      .groups = "drop"
    ) %>%
    mutate(
      rep_num = factor(rep_num, levels = c(10, 5, 2)),
      rep_n = factor(rep_n, levels = c(40, 100, 400))
    )
  
  p <- ggplot(summary, aes(x = rep_n, y = rep_num, fill = rate)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      low = "#edf8e9", high = "#006d2c",
      name = paste0("% ", method_name, " > ", cutoff),
      labels = percent_format(accuracy = 0.1),
      breaks = pretty(summary$rate, n = 5)
    ) +
    labs(
      title = paste0("Collapsed TPR (", method_name, " > ", cutoff, ")"),
      subtitle = "Collapsed across all true-effect datasets",
      x = "Replication Sample Size",
      y = "# of Replications"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  
  if (label_type == "percent") {
    p <- p + geom_text(aes(label = scales::percent(rate, accuracy = 0.1)), size = 4)
  } else if (label_type == "count") {
    p <- p + geom_text(aes(label = total_count), size = 4)
  }
  
  return(p)
}
#TPR_heatmap_collapsed(method_name = "FEMABF", cutoff = 1, label_type = "percent", effect_size= "0")
####### False Success Rate (type 2) heatmap######
# 9 x 9 grid
FSR2_heatmap <- function(method_name = "FEMABF", cutoff = 1, pcutoff = 0.05, label_type = "percent", effect_size= "0") {
  effect_key <- if (effect_size == "0") "0.2" else effect_size
  data_list <- get(paste0(method_name, "_lists_", effect_key, "null_regrouped_deltap"))
  
  summary_counts <- purrr::map_dfr(names(data_list), function(dataset_name) {
    dataset <- data_list[[dataset_name]]
    purrr::map_dfr(names(dataset), function(setting_name) {
      mat <- dataset[[setting_name]][1:500, 4:503]
      orig.p <- dataset[[setting_name]][1:500, 2]
      orig.es <- dataset[[setting_name]][1:500, 3]
      count_FS2 <- sum(mat > cutoff & orig.p < pcutoff & orig.es > 0)
      count_TS2 <- sum(mat < round(1/cutoff, 2) & orig.p > pcutoff)
      tibble(method = method_name, dataset = dataset_name, setting = setting_name, count_FS2 = count_FS2, count_TS2 = count_TS2)
    })
  }) %>%
    tidyr::separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
    mutate(
      rep_num = factor(rep_num, levels = c(10, 5, 2)),
      rep_n   = factor(rep_n,   levels = c(40, 100, 400)),
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
      rate = count_FS2 / (count_FS2 + count_TS2)
    )
  
  p <- ggplot(summary_counts, aes(x = rep_n, y = rep_num, fill = rate)) +
    geom_tile(color = "white") +
    facet_grid(bias_level ~ orig_n, labeller = label_both) +
    scale_fill_gradient(
      low = "#fee5d9", high = "#a50f15",
      name = paste0("% ", method_name, " > ", cutoff),
      labels = percent_format(accuracy = 0.1),
      breaks = pretty(summary_counts$rate, n = 5)
    ) +
    labs(
      title = paste0("False Success Rate II (", method_name, " > ", cutoff, ")"),
      subtitle = "Each tile = 250,000 null values",
      x = "Replication Sample Size (rep_n)",
      y = "# of Replications (rep_num)"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  
  if (label_type == "percent") {
    p <- p + geom_text(aes(label = scales::percent(rate, accuracy = 0.1)), size = 3)
  } else if (label_type == "count") {
    p <- p + geom_text(aes(label = count_FS2), size = 3)
  }
  
  return(p)
}

#FSR2_heatmap(method_name = "FEMABF", cutoff = 1, pcutoff = 0.05, label_type = "percent", effect_size= "0")


# 3 x 3 grid
FSR2_heatmap_collapsed <- function(method_name = "FEMABF", cutoff = 1, pcutoff = 0.05, label_type = "percent", effect_size= "0") {
  effect_key <- if (effect_size == "0") "0.2" else effect_size
  data_list <- get(paste0(method_name, "_lists_", effect_key, "null_regrouped_deltap"))
  
  summary <- purrr::map_dfr(names(data_list), function(dataset_name) {
    dataset <- data_list[[dataset_name]]
    purrr::map_dfr(names(dataset), function(setting_name) {
      mat <- dataset[[setting_name]][1:500, 4:503]
      orig.p <- dataset[[setting_name]][1:500, 2]
      orig.es <- dataset[[setting_name]][1:500, 3]
      count_FS2 <- sum(mat > cutoff & orig.p < pcutoff & orig.es > 0)
      count_TS2 <- sum(mat < round(1/cutoff, 2) & orig.p > pcutoff)
      tibble(dataset = dataset_name, setting = setting_name, count_FS2 = count_FS2, count_TS2 = count_TS2)
    })
  }) %>%
    tidyr::separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
    group_by(rep_num, rep_n) %>%
    summarise(
      total_FS2 = sum(count_FS2),
      total_TS2 = sum(count_TS2),
      rate = total_FS2 / (total_FS2 + total_TS2),
      .groups = "drop"
    ) %>%
    mutate(
      rep_num = factor(rep_num, levels = c(10, 5, 2)),
      rep_n = factor(rep_n, levels = c(40, 100, 400))
    )
  
  p <- ggplot(summary, aes(x = rep_n, y = rep_num, fill = rate)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      low = "#fee5d9", high = "#a50f15",
      name = paste0("% ", method_name, " > ", cutoff),
      labels = percent_format(accuracy = 0.1),
      breaks = pretty(summary$rate, n = 5)
    ) +
    labs(
      title = paste0("Collapsed FSR II (", method_name, " > ", cutoff, ")"),
      subtitle = "Collapsed across all null-effect datasets",
      x = "Replication Sample Size",
      y = "# of Replications"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  
  if (label_type == "percent") {
    p <- p + geom_text(aes(label = scales::percent(rate, accuracy = 0.1)), size = 4)
  } else if (label_type == "count") {
    p <- p + geom_text(aes(label = total_count), size = 4)
  }
  
  return(p)
}

#FSR2_heatmap_collapsed(method_name = "FEMABF", cutoff = 1, pcutoff = 0.05, label_type = "percent", effect_size= "0")


###### True Success Rate (type 1) ######
# 9 x 9 grid
TSR1_heatmap <- function(method_name = "FEMABF", cutoff = 1, pcutoff = 0.05, label_type = "percent", effect_size= "0") {
  effect_key <- if (effect_size == "0") "0.2" else effect_size
  data_list <- get(paste0(method_name, "_lists_", effect_key, "null_regrouped_deltap"))
  
  summary_counts <- purrr::map_dfr(names(data_list), function(dataset_name) {
    dataset <- data_list[[dataset_name]]
    purrr::map_dfr(names(dataset), function(setting_name) {
      mat <- dataset[[setting_name]][501:1000, 4:503]
      orig.p <- dataset[[setting_name]][501:1000, 2]
      orig.es <- dataset[[setting_name]][501:1000, 3]
      count_TS1 <- sum(mat > cutoff & orig.p < pcutoff & orig.es > 0)
      count_FS1 <- sum(mat < round(1/cutoff, 2) & orig.p > pcutoff)
      tibble(method = method_name, dataset = dataset_name, setting = setting_name, count_TS1 = count_TS1, count_FS1 = count_FS1)
    })
  }) %>%
    tidyr::separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
    mutate(
      rep_num = factor(rep_num, levels = c(10, 5, 2)),
      rep_n   = factor(rep_n,   levels = c(40, 100, 400)),
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
      rate = count_TS1 / (count_TS1 + count_FS1)
    )
  
  p <- ggplot(summary_counts, aes(x = rep_n, y = rep_num, fill = rate)) +
    geom_tile(color = "white") +
    facet_grid(bias_level ~ orig_n, labeller = label_both) +
    scale_fill_gradient(
      low = "#edf8e9", high = "#006d2c",
      name = paste0("% ", method_name, " > ", cutoff),
      labels = percent_format(accuracy = 0.1),
      breaks = pretty(summary_counts$rate, n = 5)
    ) +
    labs(
      title = paste0("True Success Rate I (", method_name, " > ", cutoff, ")"),
      subtitle = "Each tile = 250,000 null values",
      x = "Replication Sample Size (rep_n)",
      y = "# of Replications (rep_num)"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  
  if (label_type == "percent") {
    p <- p + geom_text(aes(label = scales::percent(rate, accuracy = 0.1)), size = 3)
  } else if (label_type == "count") {
    p <- p + geom_text(aes(label = count_TS1), size = 3)
  }
  
  return(p)
}

#TSR1_heatmap(method_name = "FEMABF", cutoff = 1, pcutoff = 0.05, label_type = "percent", effect_size= "0")


# 3 x 3 grid
TSR1_heatmap_collapsed <- function(method_name = "FEMABF", cutoff = 1, pcutoff = 0.05, label_type = "percent", effect_size= "0") {
  effect_key <- if (effect_size == "0") "0.2" else effect_size
  data_list <- get(paste0(method_name, "_lists_", effect_key, "null_regrouped_deltap"))
  
  summary <- purrr::map_dfr(names(data_list), function(dataset_name) {
    dataset <- data_list[[dataset_name]]
    purrr::map_dfr(names(dataset), function(setting_name) {
      mat <- dataset[[setting_name]][501:1000, 4:503]
      orig.p <- dataset[[setting_name]][501:1000, 2]
      orig.es <- dataset[[setting_name]][501:1000, 3]
      count_TS1 <- sum(mat > cutoff & orig.p < pcutoff & orig.es > 0)
      count_FS1 <- sum(mat < round(1/cutoff, 2) & orig.p > pcutoff)
      tibble(dataset = dataset_name, setting = setting_name, count_TS1 = count_TS1, count_FS1 = count_FS1)
    })
  }) %>%
    tidyr::separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
    group_by(rep_num, rep_n) %>%
    summarise(
      total_TS1 = sum(count_TS1),
      total_FS1 = sum(count_FS1),
      rate = total_TS1 / (total_TS1 + total_FS1),
      .groups = "drop"
    ) %>%
    mutate(
      rep_num = factor(rep_num, levels = c(10, 5, 2)),
      rep_n = factor(rep_n, levels = c(40, 100, 400))
    )
  
  p <- ggplot(summary, aes(x = rep_n, y = rep_num, fill = rate)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      low = "#edf8e9", high = "#006d2c",
      name = paste0("% ", method_name, " > ", cutoff),
      labels = percent_format(accuracy = 0.1),
      breaks = pretty(summary$rate, n = 5)
    ) +
    labs(
      title = paste0("Collapsed TSR I (", method_name, " > ", cutoff, ")"),
      subtitle = "Collapsed across all true-effect datasets",
      x = "Replication Sample Size",
      y = "# of Replications"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  
  if (label_type == "percent") {
    p <- p + geom_text(aes(label = scales::percent(rate, accuracy = 0.1)), size = 4)
  } else if (label_type == "count") {
    p <- p + geom_text(aes(label = total_count), size = 4)
  }
  
  return(p)
}

#TSR1_heatmap_collapsed(method_name = "FEMABF", cutoff = 3, pcutoff = 0.05, label_type = "percent", effect_size= "0")



######Successful Correction Rates 1######
# 9 x 9 grid
SCR1_heatmap <- function(method_name = "FEMABF", cutoff = 1, pcutoff = 0.05, label_type = "percent", effect_size= "0") {
  effect_key <- if (effect_size == "0") "0.2" else effect_size
  data_list <- get(paste0(method_name, "_lists_", effect_key, "null_regrouped_deltap"))
  
  summary_counts <- purrr::map_dfr(names(data_list), function(dataset_name) {
    dataset <- data_list[[dataset_name]]
    purrr::map_dfr(names(dataset), function(setting_name) {
      mat <- dataset[[setting_name]][501:1000, 4:503]
      orig.p <- dataset[[setting_name]][501:1000, 2]
      orig.es <- dataset[[setting_name]][501:1000, 3]
      count_FF1 <- sum(mat > cutoff & orig.p > pcutoff)
      count_FS1 <- sum(mat < round(1/cutoff, 2) & orig.p > pcutoff)
      tibble(method = method_name, dataset = dataset_name, setting = setting_name, count_FF1 = count_FF1, count_FS1 = count_FS1)
    })
  }) %>%
    tidyr::separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
    mutate(
      rep_num = factor(rep_num, levels = c(10, 5, 2)),
      rep_n   = factor(rep_n,   levels = c(40, 100, 400)),
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
      rate = count_FF1 / (count_FF1 + count_FS1)
    )
  
  p <- ggplot(summary_counts, aes(x = rep_n, y = rep_num, fill = rate)) +
    geom_tile(color = "white") +
    facet_grid(bias_level ~ orig_n, labeller = label_both) +
    scale_fill_gradient(
      low = "#edf8e9", high = "#006d2c",
      name = paste0("% ", method_name, " > ", cutoff),
      labels = percent_format(accuracy = 0.1),
      breaks = pretty(summary_counts$rate, n = 5)
    ) +
    labs(
      title = paste0("Successful Correction Rate I (", method_name, " > ", cutoff, ")"),
      subtitle = "Each tile = 250,000 null values",
      x = "Replication Sample Size (rep_n)",
      y = "# of Replications (rep_num)"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  
  if (label_type == "percent") {
    p <- p + geom_text(aes(label = scales::percent(rate, accuracy = 0.1)), size = 3)
  } else if (label_type == "count") {
    p <- p + geom_text(aes(label = count_FF1), size = 3)
  }
  
  return(p)
}

#SCR1_heatmap(method_name = "FEMABF", cutoff = 1, pcutoff = 0.05, label_type = "percent", effect_size= "0")


# 3 x 3 grid
SCR1_heatmap_collapsed <- function(method_name = "FEMABF", cutoff = 1, pcutoff = 0.05, label_type = "percent", effect_size= "0") {
  effect_key <- if (effect_size == "0") "0.2" else effect_size
  data_list <- get(paste0(method_name, "_lists_", effect_key, "null_regrouped_deltap"))
  
  summary <- purrr::map_dfr(names(data_list), function(dataset_name) {
    dataset <- data_list[[dataset_name]]
    purrr::map_dfr(names(dataset), function(setting_name) {
      mat <- dataset[[setting_name]][501:1000, 4:503]
      orig.p <- dataset[[setting_name]][501:1000, 2]
      orig.es <- dataset[[setting_name]][501:1000, 3]
      count_FF1 <- sum(mat > cutoff & orig.p > pcutoff)
      count_FS1 <- sum(mat < round(1/cutoff, 2) & orig.p > pcutoff)
      tibble(dataset = dataset_name, setting = setting_name, count_FF1 = count_FF1, count_FS1 = count_FS1)
    })
  }) %>%
    tidyr::separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
    group_by(rep_num, rep_n) %>%
    summarise(
      total_FF1 = sum(count_FF1),
      total_FS1 = sum(count_FS1),
      rate = total_FF1 / (total_FF1 + total_FS1),
      .groups = "drop"
    ) %>%
    mutate(
      rep_num = factor(rep_num, levels = c(10, 5, 2)),
      rep_n = factor(rep_n, levels = c(40, 100, 400))
    )
  
  p <- ggplot(summary, aes(x = rep_n, y = rep_num, fill = rate)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      low = "#edf8e9", high = "#006d2c",
      name = paste0("% ", method_name, " > ", cutoff),
      labels = percent_format(accuracy = 0.1),
      breaks = pretty(summary$rate, n = 5)
    ) +
    labs(
      title = paste0("Collapsed SCR I (", method_name, " > ", cutoff, ")"),
      subtitle = "Collapsed across all true-effect datasets",
      x = "Replication Sample Size",
      y = "# of Replications"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  
  if (label_type == "percent") {
    p <- p + geom_text(aes(label = scales::percent(rate, accuracy = 0.1)), size = 4)
  } else if (label_type == "count") {
    p <- p + geom_text(aes(label = total_count), size = 4)
  }
  
  return(p)
}


######Successful Correction Rates 2######
# 9 x 9 grid
SCR2_heatmap <- function(method_name = "FEMABF", cutoff = 1, pcutoff = 0.05, label_type = "percent", effect_size= "0") {
  effect_key <- if (effect_size == "0") "0.2" else effect_size
  data_list <- get(paste0(method_name, "_lists_", effect_key, "null_regrouped_deltap"))
  
  summary_counts <- purrr::map_dfr(names(data_list), function(dataset_name) {
    dataset <- data_list[[dataset_name]]
    purrr::map_dfr(names(dataset), function(setting_name) {
      mat <- dataset[[setting_name]][1:500, 4:503]
      orig.p <- dataset[[setting_name]][1:500, 2]
      orig.es <- dataset[[setting_name]][1:500, 3]
      count_FF2 <- sum(mat < round(1/cutoff, 2) & orig.p < pcutoff & orig.es > 0)
      count_FS2 <- sum(mat > cutoff & orig.p < pcutoff & orig.es > 0)
      tibble(method = method_name, dataset = dataset_name, setting = setting_name, count_FF2 = count_FF2, count_FS2 = count_FS2)
    })
  }) %>%
    tidyr::separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
    mutate(
      rep_num = factor(rep_num, levels = c(10, 5, 2)),
      rep_n   = factor(rep_n,   levels = c(40, 100, 400)),
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
      rate = count_FF2 / (count_FF2 + count_FS2)
    )
  
  p <- ggplot(summary_counts, aes(x = rep_n, y = rep_num, fill = rate)) +
    geom_tile(color = "white") +
    facet_grid(bias_level ~ orig_n, labeller = label_both) +
    scale_fill_gradient(
      low = "#edf8e9", high = "#006d2c",
      name = paste0("% ", method_name, " > ", cutoff),
      labels = percent_format(accuracy = 0.1),
      breaks = pretty(summary_counts$rate, n = 5)
    ) +
    labs(
      title = paste0("Successful Correction Rate II (", method_name, " > ", cutoff, ")"),
      subtitle = "Each tile = 250,000 null values",
      x = "Replication Sample Size (rep_n)",
      y = "# of Replications (rep_num)"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  
  if (label_type == "percent") {
    p <- p + geom_text(aes(label = scales::percent(rate, accuracy = 0.1)), size = 3)
  } else if (label_type == "count") {
    p <- p + geom_text(aes(label = count_FF2), size = 3)
  }
  
  return(p)
}

#SCR2_heatmap(method_name = "FEMABF", cutoff = 1, pcutoff = 0.05, label_type = "percent", effect_size= "0")


# 3 x 3 grid
SCR2_heatmap_collapsed <- function(method_name = "FEMABF", cutoff = 1, pcutoff = 0.05, label_type = "percent", effect_size= "0") {
  effect_key <- if (effect_size == "0") "0.2" else effect_size
  data_list <- get(paste0(method_name, "_lists_", effect_key, "null_regrouped_deltap"))
  
  summary <- purrr::map_dfr(names(data_list), function(dataset_name) {
    dataset <- data_list[[dataset_name]]
    purrr::map_dfr(names(dataset), function(setting_name) {
      mat <- dataset[[setting_name]][1:500, 4:503]
      orig.p <- dataset[[setting_name]][1:500, 2]
      orig.es <- dataset[[setting_name]][1:500, 3]
      count_FF2 <- sum(mat < round(1/cutoff, 2) & orig.p < pcutoff & orig.es > 0)
      count_FS2 <- sum(mat > cutoff & orig.p < pcutoff & orig.es > 0)
      tibble(dataset = dataset_name, setting = setting_name, count_FF2 = count_FF2, count_FS2 = count_FS2)
    })
  }) %>%
    tidyr::separate(setting, into = c("rep_num", "rep_n"), sep = "_", convert = TRUE) %>%
    group_by(rep_num, rep_n) %>%
    summarise(
      total_FF2 = sum(count_FF2),
      total_FS2 = sum(count_FS2),
      rate = total_FF2 / (total_FF2 + total_FS2),
      .groups = "drop"
    ) %>%
    mutate(
      rep_num = factor(rep_num, levels = c(10, 5, 2)),
      rep_n = factor(rep_n, levels = c(40, 100, 400))
    )
  
  p <- ggplot(summary, aes(x = rep_n, y = rep_num, fill = rate)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      low = "#edf8e9", high = "#006d2c",
      name = paste0("% ", method_name, " > ", cutoff),
      labels = percent_format(accuracy = 0.1),
      breaks = pretty(summary$rate, n = 5)
    ) +
    labs(
      title = paste0("Collapsed SCR II (", method_name, " > ", cutoff, ")"),
      subtitle = "Collapsed across all true-effect datasets",
      x = "Replication Sample Size",
      y = "# of Replications"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  
  if (label_type == "percent") {
    p <- p + geom_text(aes(label = scales::percent(rate, accuracy = 0.1)), size = 4)
  } else if (label_type == "count") {
    p <- p + geom_text(aes(label = total_count), size = 4)
  }
  
  return(p)
}