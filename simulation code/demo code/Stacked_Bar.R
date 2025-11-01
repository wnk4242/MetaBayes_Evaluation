library(ggplot2)
library(tidyverse)
library(scales)
# Use row-wise rates to create stacked bar plots
# Set the path to the directory containing RDS files
folder_path <- "./data files/data for stacked bar" 
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

#############################################################################
#Compare AD, FP, TN rates across MABF methods when underlying effect is null#
#############################################################################
# rates when underlying ES is null
# Save the original working directory
original_directory <- getwd()
# List of dataset names
dataset_names <- c("MABF_rates_null_BFcutoff3")
# Custom labels for metrics
metric_labels <- c("AD_null" = "Anecdotal Evidence", "FP" = "False Positive", "TN" = "True Negative")
# Save the plot as JPG
plot_dir <- "./plots/stacked bar"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}
setwd(plot_dir)
# Loop through the datasets
for (dataset_name in dataset_names) {
  # Extract the BF cutoff value from the dataset name
  BFcutoff_value <- sub(".*_BFcutoff", "", dataset_name)
  # Modify the dataset
  plot_data_null <- get(dataset_name) %>%
    mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
    mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N = 2", "N = 5", "N = 10"))) %>%
    mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n = 40", "n = 100", "n = 400"))) %>%
    mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>%
    select(method, rep.n, rep.number, AD_null, FP, TN) %>%
    gather(key = "metric", value = "value", AD_null, FP, TN) %>%
    group_by(method, rep.n, rep.number) %>%
    mutate(proportion = value / sum(value))
  # Create the plot
  plot <- ggplot(plot_data_null, aes(x = factor(method), y = proportion, fill = metric)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("AD_null" = "red", "FP" = "green", "TN" = "blue"),
                      labels = as_labeller(metric_labels)) +  # Custom labels using as_labeller()
    labs(x = "Replication Sample Size", y = "Proportion", fill = "Metric", title = paste("Proportions of anecdotal evidence cases, false positives, and true negatives when the cutoff for anecdotal evidence cases is", BFcutoff_value, "when θ = 0")) +
    facet_grid(rep.n ~ rep.number, scales = "free_x", labeller = label_value) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")
  ggsave(filename = paste0("null_stackedFPTNAD_BFcutoff", BFcutoff_value, ".jpg"), plot = plot, device = "jpg", width = 14, height = 10, bg = "white")
}
# Reset working directory to the original directory
setwd(original_directory)

############################################################################
#Compare AD, TP, FN rates across MABF methods when underlying effect is 0.2#
############################################################################
# rates when underlying ES is true
# Save the original working directory
original_directory <- getwd()
# List of dataset names
dataset_names <- c("MABF_rates_0.2_BFcutoff3")
# Custom labels for metrics
metric_labels <- c("AD_true" = "Anecdotal Evidence", "TP" = "True Positive", "FN" = "False Negative")
# Save the plot as JPG
plot_dir <- "./plots/stacked bar"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}
setwd(plot_dir)

# Loop through the datasets
for (dataset_name in dataset_names) {
  # Extract the number from the dataset name
  BFcutoff_value <- sub(".*_BFcutoff", "", dataset_name)
  # Modify the dataset
  plot_data_0.2 <- get(dataset_name)  %>%
    mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>% 
    mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N = 2", "N = 5", "N = 10"))) %>%
    mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n = 40", "n = 100", "n = 400"))) %>%
    mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>% 
    select(method, rep.n, rep.number, AD_true, TP, FN) %>%
    gather(key = "metric", value = "value", AD_true, TP, FN) %>%
    group_by(method, rep.n, rep.number) %>%
    mutate(proportion = value / sum(value))
  # Create stacked bar chart showing proportions
  plot <- ggplot(plot_data_0.2, aes(x = factor(method), y = proportion, fill = metric)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("AD_true" = "red", "TP" = "blue", "FN" = "green"),
                      labels = as_labeller(metric_labels)) +  # Custom labels using as_labeller()
    labs(x = "Replication Sample Size", y = "Proportion", fill = "Metric", title = paste("Proportions of anecdotal evidence cases, true positives, and false negatives when the cutoff for anecdotal evidence cases is", BFcutoff_value, "when θ = 0.2")) +
    facet_grid(rep.n ~ rep.number, scales = "free_x", labeller = label_value) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")
  ggsave(filename = paste0("0.2_stackedTPFNAD_BFcutoff", BFcutoff_value, ".jpg"), plot = plot, device = "jpg", width = 14, height = 10, bg = "white")
}
# Reset working directory to the original directory
setwd(original_directory)


############################################################################
#Compare AD, TP, FN rates across MABF methods when underlying effect is 0.5#
############################################################################
# rates when underlying ES is true
# Save the original working directory
original_directory <- getwd()
# List of dataset names
dataset_names <- c("MABF_rates_0.5_BFcutoff3")
# Custom labels for metrics
metric_labels <- c("AD_true" = "Anecdotal Evidence", "TP" = "True Positive", "FN" = "False Negative")
# Save the plot as JPG
plot_dir <- "./plots/stacked bar"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}
setwd(plot_dir)

# Loop through the datasets
for (dataset_name in dataset_names) {
  # Extract the number from the dataset name
  BFcutoff_value <- sub(".*_BFcutoff", "", dataset_name)
  # Modify the dataset
  plot_data_0.5 <- get(dataset_name)  %>%
    mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>% 
    mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N = 2", "N = 5", "N = 10"))) %>%
    mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n = 40", "n = 100", "n = 400"))) %>%
    mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>% 
    select(method, rep.n, rep.number, AD_true, TP, FN) %>%
    gather(key = "metric", value = "value", AD_true, TP, FN) %>%
    group_by(method, rep.n, rep.number) %>%
    mutate(proportion = value / sum(value))
  # Create stacked bar chart showing proportions
  plot <- ggplot(plot_data_0.5, aes(x = factor(method), y = proportion, fill = metric)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("AD_true" = "red", "TP" = "blue", "FN" = "green"),
                      labels = as_labeller(metric_labels)) +  # Custom labels using as_labeller()
    labs(x = "Replication Sample Size", y = "Proportion", fill = "Metric", title = paste("Proportions of anecdotal evidence cases, true positives, and false negatives when the cutoff for anecdotal evidence cases is", BFcutoff_value, "when θ = 0.5")) +
    facet_grid(rep.n ~ rep.number, scales = "free_x", labeller = label_value) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")
  ggsave(filename = paste0("0.5_stackedTPFNAD_BFcutoff", BFcutoff_value, ".jpg"), plot = plot, device = "jpg", width = 14, height = 10, bg = "white")
}
# Reset working directory to the original directory
setwd(original_directory)


############################################################################################################
#Compare AD, TS, FS, TF, FF rates filtered by bias level across MABF methods when underlying effect is null#
############################################################################################################
# Load the data
BiasFilter_null <- MABF_rates_null_BFcutoff3 %>%
  mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
  mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
  mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
  mutate(bias.level = factor(QRP.level, levels = c("none", "medium", "high"), labels = c("bias = low", "bias = medium", "bias = high"))) %>%
  select(method, rep.n, rep.number, AD_null, TS2, FS2, FF2, TF2, bias.level) %>%
  gather(key = "metric", value = "value", AD_null, TS2, FS2, FF2, TF2) %>%
  mutate(metric = factor(metric, levels = c("AD_null", "TS2", "FS2", "FF2", "TF2"))) %>%
  mutate(metric_label = factor(metric, 
                               levels = c("AD_null", "TS2", "FS2", "FF2", "TF2"),
                               labels = c("anecdotal", "true success", "false success", "false failure", "true failure"))) %>%
  group_by(method, rep.n, rep.number, bias.level) %>%
  mutate(proportion = value / sum(value))

# Create modern high contrast colors for the AD_null metrics
modern_high_contrast_colors <- c("#E69F00", "#56B4E9", "coral2", "#F0E442", "darkorchid")

# Create a named vector for scale_fill_manual
color_mapping <- c(anecdotal = modern_high_contrast_colors[1],
                   "true success" = modern_high_contrast_colors[2],
                   "false success" = modern_high_contrast_colors[3],
                   "false failure" = modern_high_contrast_colors[4],
                   "true failure" = modern_high_contrast_colors[5])

# Custom labeller for metrics
metric_labeller <- c(
  AD_null = "anecdotal",
  TS2 = "true success",
  FS2 = "false success",
  FF2 = "false failure",
  TF2 = "true failure"
)

# Create the plot
stackedplot_BiasFilter_null <- ggplot(BiasFilter_null, aes(x = factor(method), y = proportion, fill = metric_label)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_mapping) +  # Custom colors
  labs(x = "MABF Method", y = "Proportion", fill = "Evidence Category", title = "Proportion of Bayes Factors Categorized into Success and Failure Influenced by Bias Levels When the Underlying Effect is Null") +
  facet_grid(bias.level ~ rep.number + rep.n, scales = "free_x", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")
print(stackedplot_BiasFilter_null)
# Save the plot
ggsave("./plots/stacked bar/null_stackedADTSFS_BiasFilter_pcutoff=0.05_BFcutoff=3.jpg", plot = stackedplot_BiasFilter_null, width = 15, height = 10)


###########################################################################################################
#Compare AD, TS, FS, TF, FF rates filtered by bias level across MABF methods when underlying effect is 0.2#
###########################################################################################################
# Load the data
BiasFilter_0.2 <- MABF_rates_0.2_BFcutoff3 %>%
  mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
  mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
  mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
  mutate(bias.level = factor(QRP.level, levels = c("none", "medium", "high"), labels = c("bias = low", "bias = medium", "bias = high"))) %>%
  select(method, rep.n, rep.number, AD_true, TS1, FS1, FF1, TF1, bias.level) %>%
  gather(key = "metric", value = "value", AD_true, TS1, FS1, FF1, TF1) %>%
  mutate(metric = factor(metric, levels = c("AD_true", "TS1", "FS1", "FF1", "TF1"))) %>%
  mutate(metric_label = factor(metric, 
                               levels = c("AD_true", "TS1", "FS1", "FF1", "TF1"),
                               labels = c("anecdotal", "true success", "false success", "false failure", "true failure"))) %>%
  group_by(method, rep.n, rep.number, bias.level) %>%
  mutate(proportion = value / sum(value))

# Create modern high contrast colors for the AD_true metrics
modern_high_contrast_colors <- c("#E69F00", "#56B4E9", "coral2", "#F0E442", "darkorchid")

# Create a named vector for scale_fill_manual
color_mapping <- c(anecdotal = modern_high_contrast_colors[1],
                   "true success" = modern_high_contrast_colors[2],
                   "false success" = modern_high_contrast_colors[3],
                   "false failure" = modern_high_contrast_colors[4],
                   "true failure" = modern_high_contrast_colors[5])

# Custom labeller for metrics
metric_labeller <- c(
  AD_true = "anecdotal",
  TS1 = "true success",
  FS1 = "false success",
  FF1 = "false failure",
  TF1 = "true failure"
)

# Create the plot
stackedplot_BiasFilter_0.2 <- ggplot(BiasFilter_0.2, aes(x = factor(method), y = proportion, fill = metric_label)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_mapping) +  # Custom colors
  labs(x = "MABF Method", y = "Proportion", fill = "Evidence Category", title = "Proportion of Bayes Factors Categorized into Success and Failure Influenced by Bias Levels When the Underlying Effect is 0.2") +
  facet_grid(bias.level ~ rep.number + rep.n, scales = "free_x", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")
print(stackedplot_BiasFilter_0.2)
# Save the plot
ggsave("./plots/stacked bar/stackedADTSFS_0.2_BiasFilter_pcutoff=0.05_BFcutoff=3.jpg", plot = stackedplot_BiasFilter_0.2, width = 15, height = 10)


###########################################################################################################
#Compare AD, TS, FS, TF, FF rates filtered by bias level across MABF methods when underlying effect is 0.5#
###########################################################################################################
# Load the data
BiasFilter_0.5 <- MABF_rates_0.5_BFcutoff3 %>%
  mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
  mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
  mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
  mutate(bias.level = factor(QRP.level, levels = c("none", "medium", "high"), labels = c("bias = low", "bias = medium", "bias = high"))) %>%
  select(method, rep.n, rep.number, AD_true, TS1, FS1, FF1, TF1, bias.level) %>%
  gather(key = "metric", value = "value", AD_true, TS1, FS1, FF1, TF1) %>%
  mutate(metric = factor(metric, levels = c("AD_true", "TS1", "FS1", "FF1", "TF1"))) %>%
  mutate(metric_label = factor(metric, 
                               levels = c("AD_true", "TS1", "FS1", "FF1", "TF1"),
                               labels = c("anecdotal", "true success", "false success", "false failure", "true failure"))) %>%
  group_by(method, rep.n, rep.number, bias.level) %>%
  mutate(proportion = value / sum(value))

# Create modern high contrast colors for the AD_true metrics
modern_high_contrast_colors <- c("#E69F00", "#56B4E9", "coral2", "#F0E442", "darkorchid")

# Create a named vector for scale_fill_manual
color_mapping <- c(anecdotal = modern_high_contrast_colors[1],
                   "true success" = modern_high_contrast_colors[2],
                   "false success" = modern_high_contrast_colors[3],
                   "false failure" = modern_high_contrast_colors[4],
                   "true failure" = modern_high_contrast_colors[5])

# Custom labeller for metrics
metric_labeller <- c(
  AD_true = "anecdotal",
  TS1 = "true success",
  FS1 = "false success",
  FF1 = "false failure",
  TF1 = "true failure"
)

# Create the plot
stackedplot_BiasFilter_0.5 <- ggplot(BiasFilter_0.5, aes(x = factor(method), y = proportion, fill = metric_label)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_mapping) +  # Custom colors
  labs(x = "MABF Method", y = "Proportion", fill = "Evidence Category", title = "Proportion of Bayes Factors Categorized into Success and Failure Influenced by Bias Levels When the Underlying Effect is 0.5") +
  facet_grid(bias.level ~ rep.number + rep.n, scales = "free_x", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")
print(stackedplot_BiasFilter_0.5)
# Save the plot
ggsave("./plots/stacked bar/stackedADTSFS_0.5_BiasFilter_pcutoff=0.05_BFcutoff=3.jpg", plot = stackedplot_BiasFilter_0.5, width = 15, height = 10)
rm(list = ls())
