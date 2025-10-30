#This script will visualize original studies' false positive and false negative rates under different conditions 
#in individual scattered dot plot and facet grid dot histogram, and heatmap
#A Shiny app is built based on this script in Step 5.1_Shiny_origp_GridHistogram&HeatMap&Bar
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)
# Load necessary helper functions: ratesCalc_TNFP and ratesCalc_TPFN 
source("./helper functions.R")
source("./helper functions_XXR.R")
# Set the path to the directory containing RDS files
folder_path <- "./OGDG/10000" 
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}
df_lists <- df_lists_10000 #for convenience purpose
remove(df_lists_10000)



#####
#Create individual plot showing false positives from null effect datasets using a helper function
# Scattered dot (no pattern)
# plot_false_positives_jitter <- function(condition_name, df_list) {
#   # Check condition exists
#   if (!condition_name %in% names(df_list)) {
#     stop("Condition not found in df_list.")
#   }
#   
#   # Extract and label data
#   df <- df_list[[condition_name]] %>%
#     mutate(is_false_positive = p < 0.05)
#   
#   fpr <- mean(df$is_false_positive)
#   
#   # Parse condition name
#   parts <- strsplit(condition_name, "_")[[1]]
#   effect_size <- as.numeric(parts[1])
#   orig_n <- as.integer(parts[3])
#   p_hacking <- parts[4]
#   pub_bias <- parts[5]
#   
#   # Bias description logic
#   bias_desc <- case_when(
#     p_hacking == "none" & pub_bias == "low" ~ "No p-hacking or publication bias",
#     p_hacking == "medium" & pub_bias == "medium" ~ "Medium p-hacking and publication bias",
#     p_hacking == "high" & pub_bias == "high" ~ "High p-hacking and publication bias",
#     TRUE ~ paste0("P-hacking: ", p_hacking, ", Publication bias: ", pub_bias)
#   )
#   
#   # Subtitle text
#   subtitle_text <- paste0(
#     "True Effect Size = ", effect_size,
#     ", Original Sample Size = ", orig_n,
#     ", ", bias_desc, "\n",
#     "False Positive Rate = ", round(fpr, 3)
#   )
#   
#   # Plot
#   ggplot(df, aes(x = d, color = is_false_positive)) +
#     geom_jitter(aes(y = 0), height = 0.1, alpha = 0.5) +
#     scale_color_manual(
#       values = c("FALSE" = "gray60", "TRUE" = "firebrick"),
#       labels = c("TRUE" = "False Positive", "FALSE" = "True Negative")
#     ) +
#     labs(
#       title = "False Positives When True Effect Size = 0",
#       subtitle = subtitle_text,
#       x = "Observed Effect Size (d)",
#       y = "",
#       color = "Result Type"
#     ) +
#     theme_minimal() +
#     theme(
#       axis.text.y = element_blank(),
#       axis.ticks.y = element_blank(),
#       legend.position = "bottom"
#     )
# }
# 
# 
# plot_false_positives_jitter("0_0.02_20_none_low", df_lists)
# plot_false_positives_jitter("0_0.02_20_medium_medium", df_lists)
# plot_false_positives_jitter("0_0.02_20_high_high", df_lists)

#####
# Create a grid plot showing false positives from null effect datasets under different simulation conditions
#####
# Scattered dot (no pattern)
# Step 1: Define the 9 conditions to include
# conditions_to_plot <- c(
#   "0_0.02_20_none_low",
#   "0_0.02_50_none_low",
#   "0_0.02_200_none_low",
#   "0_0.02_20_medium_medium",
#   "0_0.02_50_medium_medium",
#   "0_0.02_200_medium_medium",
#   "0_0.02_20_high_high",
#   "0_0.02_50_high_high",
#   "0_0.02_200_high_high"
# )
# 
# # Step 2: Extract and combine into one data frame
# plot_df <- map_dfr(conditions_to_plot, function(cond_name) {
#   parts <- str_split(cond_name, "_")[[1]]
#   effect_size <- as.numeric(parts[1])
#   orig_n <- as.integer(parts[3])
#   p_hacking <- parts[4]
#   pub_bias <- parts[5]
#   
#   bias_level <- case_when(
#     p_hacking == "none" & pub_bias == "low" ~ "none",
#     p_hacking == "medium" & pub_bias == "medium" ~ "medium",
#     p_hacking == "high" & pub_bias == "high" ~ "high",
#     TRUE ~ "other"
#   )
#   
#   df_lists[[cond_name]] %>%
#     mutate(
#       is_false_positive = p < 0.05,
#       orig_n = orig_n,
#       bias_level = bias_level
#     )
# })
# 
# # Make bias_level and orig_n into factors for facet ordering
# plot_df <- plot_df %>%
#   mutate(
#     orig_n = factor(orig_n, levels = c(20, 50, 200)),
#     bias_level = factor(bias_level, levels = c("none", "medium", "high"))
#   )
# 
# # Step 3: Plot
# ggplot(plot_df, aes(x = d, color = is_false_positive)) +
#   geom_jitter(aes(y = 0), height = 0.1, alpha = 0.4) +
#   scale_color_manual(
#     values = c("FALSE" = "gray60", "TRUE" = "firebrick"),
#     labels = c("TRUE" = "False Positive", "FALSE" = "True Negative")
#   ) +
#   labs(
#     title = "False Positives When True Effect Size = 0",
#     subtitle = "Faceted by Bias Level (rows) and Original Sample Size (columns)",
#     x = "Observed Effect Size (d)",
#     y = "",
#     color = "Result Type"
#   ) +
#   facet_grid(bias_level ~ orig_n) +
#   theme_minimal() +
#   theme(
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     legend.position = "bottom"
#   )
#####
# Histogram
# Step 1: list of 9 null-effect simulation conditions
conditions_to_plot <- c(
  "0_0.02_20_none_low",
  "0_0.02_50_none_low",
  "0_0.02_200_none_low",
  "0_0.02_20_medium_medium",
  "0_0.02_50_medium_medium",
  "0_0.02_200_medium_medium",
  "0_0.02_20_high_high",
  "0_0.02_50_high_high",
  "0_0.02_200_high_high"
)

# Step 2: extract and label
plot_df <- map_dfr(conditions_to_plot, function(cond_name) {
  parts <- str_split(cond_name, "_")[[1]]
  effect_size <- as.numeric(parts[1])
  orig_n <- as.integer(parts[3])
  p_hacking <- parts[4]
  pub_bias <- parts[5]
  
  bias_level <- case_when(
    p_hacking == "none" & pub_bias == "low" ~ "none",
    p_hacking == "medium" & pub_bias == "medium" ~ "medium",
    p_hacking == "high" & pub_bias == "high" ~ "high",
    TRUE ~ "other"
  )
  
  df_lists[[cond_name]] %>%
    mutate(
      is_false_positive = p < 0.05,
      orig_n = orig_n,
      bias_level = bias_level
    )
})

# Step 3: format variables
plot_df <- plot_df %>%
  mutate(
    orig_n = factor(orig_n, levels = c(20, 50, 200)),
    bias_level = factor(bias_level, levels = c("none", "medium", "high")),
    dot_color = ifelse(is_false_positive, "False Positive", "True Negative")
  )

# Step 4: calculate FPR per condition
fpr_df <- plot_df %>%
  group_by(bias_level, orig_n) %>%
  summarise(FPR = mean(is_false_positive), .groups = "drop") %>%
  mutate(
    label = paste0("FPR = ", round(FPR, 3)),
    x = -0.5,
    y = 3000  # Move text into the upper white space (adjust as needed)
  )

# Step 5: plot with density, dots, and annotation
ggplot(plot_df, aes(x = d)) +
  # Histogram with fill by false positive status
  geom_histogram(
    aes(fill = is_false_positive),
    bins = 50,
    position = "stack",  # stack false positive and true negative cases
    alpha = 0.4,
    color = "gray30"
  ) +
  # FPR annotation
  geom_text(
    data = fpr_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3,
    fontface = "italic"
  ) +
  # Fill colors
  scale_fill_manual(
    values = c(`TRUE` = "firebrick", `FALSE` = "lightblue"),
    labels = c("True Negative","False Positive"),
    name = "Result Type"
  ) +
  facet_grid(bias_level ~ orig_n, labeller = label_both) +
  labs(
    title = "Observed Effect Size (d) under the Null",
    subtitle = "Red = False Positives (p < .05); Blue = True Negatives",
    x = "Observed Effect Size (d)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

# Heatmap of FPR
fpr_df <- fpr_df %>%
  mutate(bias_level = factor(bias_level, levels = c("high", "medium", "none")))

ggplot(fpr_df, aes(x = orig_n, y = bias_level, fill = FPR)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(FPR, 3)), color = "black", size = 4) +
  scale_fill_gradient(
    low = "#fef6d6", high = "#fb972c",
    name = "False Positive Rate"
  ) +
  labs(
    title = "False Positive Rate (p < .05) by Bias Level and Sample Size",
    x = "Original Sample Size",
    y = "Bias Level"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold")
  )

#####
# Create a grid plot showing true positives from 0.2 effect datasets under different simulation conditions
#####
# Step 1: list of 9 effect_size = 0.2 simulation conditions
conditions_to_plot <- c(
  "0.2_0.05_20_none_low",
  "0.2_0.05_50_none_low",
  "0.2_0.05_200_none_low",
  "0.2_0.05_20_medium_medium",
  "0.2_0.05_50_medium_medium",
  "0.2_0.05_200_medium_medium",
  "0.2_0.05_20_high_high",
  "0.2_0.05_50_high_high",
  "0.2_0.05_200_high_high"
)

# Step 2: extract and label
plot_df <- map_dfr(conditions_to_plot, function(cond_name) {
  parts <- str_split(cond_name, "_")[[1]]
  effect_size <- as.numeric(parts[1])
  orig_n <- as.integer(parts[3])
  p_hacking <- parts[4]
  pub_bias <- parts[5]
  
  bias_level <- case_when(
    p_hacking == "none" & pub_bias == "low" ~ "none",
    p_hacking == "medium" & pub_bias == "medium" ~ "medium",
    p_hacking == "high" & pub_bias == "high" ~ "high",
    TRUE ~ "other"
  )
  
  df_lists[[cond_name]] %>%
    mutate(
      is_true_positive = p < 0.05,
      is_false_negative = p >=0.05,
      orig_n = orig_n,
      bias_level = bias_level
    )
})

# Step 3: format for plotting
plot_df <- plot_df %>%
  mutate(
    orig_n = factor(orig_n, levels = c(20, 50, 200)),
    bias_level = factor(bias_level, levels = c("none", "medium", "high")),
    dot_color = ifelse(is_true_positive, "True Positive", "False Negative")
  )

# Step 4: calculate TPR (true positive rate)
tpr_df <- plot_df %>%
  group_by(bias_level, orig_n) %>%
  summarise(TPR = mean(is_true_positive), .groups = "drop") %>%
  mutate(
    label = paste0("TPR = ", round(TPR, 3)),
    x = -0.5,
    y = 2000
  )
# Step 4.1: calculate FNR
fnr_df <- plot_df %>%
  group_by(bias_level, orig_n) %>%
  summarise(FNR = mean(is_false_negative), .groups = "drop") %>%
  mutate(
    label = paste0("FNR = ", round(FNR, 3)),
    x = -0.5,
    y = 4
  )
# Step 5: plot
ggplot(plot_df, aes(x = d)) +
  # Histogram with fill by false positive status
  geom_histogram(
    aes(fill = is_true_positive),
    bins = 50,
    position = "stack", 
    alpha = 0.4,
    color = "gray30"
  ) +
  # TPR annotation
  geom_text(
    data = tpr_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3,
    fontface = "italic"
  ) +
  # Fill colors
  scale_fill_manual(
    values = c(`TRUE` = "green", `FALSE` = "lightblue"),
    labels = c("False Negative","True Positive"),
    name = "Result Type"
  ) +
  facet_grid(bias_level ~ orig_n, labeller = label_both) +
  labs(
    title = "Observed Effect Size (d) under the Null",
    subtitle = "Green = True Positives (p < .05); Blue = False Negatives",
    x = "Observed Effect Size (d)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

# Heatmap of FNR
fnr_df <- fnr_df %>%
  mutate(bias_level = factor(bias_level, levels = c("high", "medium", "none")))

ggplot(fnr_df, aes(x = orig_n, y = bias_level, fill = FNR)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(FNR, 3)), color = "black", size = 4) +
  scale_fill_gradient(
    low = "#c9dffc", high = "#2c87fb",
    name = "False Negative Rate"
  ) +
  labs(
    title = "False Negative Rate (p > .05) by Bias Level and Sample Size",
    x = "Original Sample Size",
    y = "Bias Level"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold")
  )


# Create a grid plot showing true positives from 0.5 effect datasets under different simulation conditions
# Step 1: list of 9 effect_size = 0.2 simulation conditions
conditions_to_plot <- c(
  "0.5_0.125_20_none_low",
  "0.5_0.125_50_none_low",
  "0.5_0.125_200_none_low",
  "0.5_0.125_20_medium_medium",
  "0.5_0.125_50_medium_medium",
  "0.5_0.125_200_medium_medium",
  "0.5_0.125_20_high_high",
  "0.5_0.125_50_high_high",
  "0.5_0.125_200_high_high"
)

# Step 2: extract and label
plot_df <- map_dfr(conditions_to_plot, function(cond_name) {
  parts <- str_split(cond_name, "_")[[1]]
  effect_size <- as.numeric(parts[1])
  orig_n <- as.integer(parts[3])
  p_hacking <- parts[4]
  pub_bias <- parts[5]
  
  bias_level <- case_when(
    p_hacking == "none" & pub_bias == "low" ~ "none",
    p_hacking == "medium" & pub_bias == "medium" ~ "medium",
    p_hacking == "high" & pub_bias == "high" ~ "high",
    TRUE ~ "other"
  )
  
  df_lists[[cond_name]] %>%
    mutate(
      is_true_positive = p < 0.05,
      is_false_negative = p >= 0.05,
      orig_n = orig_n,
      bias_level = bias_level
    )
})

# Step 3: format for plotting
plot_df <- plot_df %>%
  mutate(
    orig_n = factor(orig_n, levels = c(20, 50, 200)),
    bias_level = factor(bias_level, levels = c("none", "medium", "high")),
    dot_color = ifelse(is_true_positive, "True Positive", "False Negative")
  )

# Step 4: calculate TPR (true positive rate)
tpr_df <- plot_df %>%
  group_by(bias_level, orig_n) %>%
  summarise(TPR = mean(is_true_positive), .groups = "drop") %>%
  mutate(
    label = paste0("TPR = ", round(TPR, 3)),
    x = 0.5,
    y = 2000
  )
# Step 4.1: calculate FNR
fnr_df <- plot_df %>%
  group_by(bias_level, orig_n) %>%
  summarise(FNR = mean(is_false_negative), .groups = "drop") %>%
  mutate(
    label = paste0("FNR = ", round(FNR, 3)),
    x = 0.5,
    y = 4
  )
# Step 5: plot
ggplot(plot_df, aes(x = d)) +
  # Histogram with fill by false positive status
  geom_histogram(
    aes(fill = is_true_positive),
    bins = 50,
    position = "stack", 
    alpha = 0.4,
    color = "gray30"
  ) +
  # TPR annotation
  geom_text(
    data = tpr_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3,
    fontface = "italic"
  ) +
  # Fill colors
  scale_fill_manual(
    values = c(`TRUE` = "green", `FALSE` = "lightblue"),
    labels = c("False Negative","True Positive"),
    name = "Result Type"
  ) +
  facet_grid(bias_level ~ orig_n, labeller = label_both) +
  labs(
    title = "Observed Effect Size (d) under the Null",
    subtitle = "Green = True Positives (p < .05); Blue = False Negatives",
    x = "Observed Effect Size (d)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )


# Heatmap of FNR
fnr_df <- fnr_df %>%
  mutate(bias_level = factor(bias_level, levels = c("high", "medium", "none")))

ggplot(fnr_df, aes(x = orig_n, y = bias_level, fill = FNR)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(FNR, 3)), color = "black", size = 4) +
  scale_fill_gradient(
    low = "#c9dffc", high = "#2c87fb",
    name = "False Negative Rate"
  ) +
  labs(
    title = "False Negative Rate (p > .05) by Bias Level and Sample Size",
    x = "Original Sample Size",
    y = "Bias Level"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold")
  )
