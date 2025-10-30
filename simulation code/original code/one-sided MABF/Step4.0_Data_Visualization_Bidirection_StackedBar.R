# The bidirection stacked bar plot demonstrate the number of (strong) BF for null hypothesis compared to number of (strong) BF for alternative hypothesis 

library(ggplot2)
library(dplyr)
library(tidyverse)
# Use row-wise rates to create stacked bar plots
# Set the path to the directory containing RDS files
folder_path <- "./MABFanalyses/row-wise/new" 
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

#####null
# Sample data transformation
plot_data_null <- get("MABF_rates_null_rowise_pcutoff0.05_BFcutoff3") %>%
  mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
  mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N = 2", "N = 5", "N = 10"))) %>%
  mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n = 40", "n = 100", "n = 400"))) %>%
  mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>%
  mutate(AD_null_allstrong_N = AD_null_strong_N + AD_null_vstrong_N + AD_null_estrong_N) %>%
  mutate(AD_null_allstrong_A = AD_null_strong_A + AD_null_vstrong_A + AD_null_estrong_A) %>%
  select(method, rep.n, rep.number, AD_null_anecdotal_A, AD_null_moderate_A, AD_null_allstrong_A,
         AD_null_anecdotal_N, AD_null_moderate_N, AD_null_allstrong_N)

# Reshape the data and calculate median
data_median <- plot_data_null %>%
  pivot_longer(cols = starts_with("AD_"),
               names_to = c("type", "null", "strength", "direction"),
               names_sep = "_",
               values_to = "count") %>%
  mutate(strength = recode_factor(paste(null, strength, sep = "_"),
                                  null_allstrong = "strong",
                                  null_moderate = "moderate",
                                  null_anecdotal = "anecdotal"),
         direction = ifelse(direction == "A", "Correct", "Wrong")) %>%
  group_by(method, rep.n, rep.number, strength, direction) %>%
  summarise(count = median(count), .groups = 'drop')

# Define colors for each method
method_colors <- c("FEMABF" = "red", "BFbMA" = "green", "EUBF" = "blue", "iBF" = "purple")

# Create the plot with reduced bar width and dashed lines
plot <- ggplot(data_median, aes(x = strength, fill = method)) +
  geom_bar(
    aes(y = ifelse(direction == "Correct", count, 0)),
    stat = "identity",
    position = position_dodge(width = 1),
    width = 0.6,  # Reduce bar width
    alpha = 1
  ) +
  geom_bar(
    aes(y = ifelse(direction == "Wrong", -count, 0)),
    stat = "identity",
    position = position_dodge(width = 1),
    width = 0.6,  # Reduce bar width
    alpha = 0.5
  ) +
  scale_fill_manual(values = method_colors) +
  coord_flip() +
  facet_grid(rep.n ~ rep.number) +
  labs(x = "Evidence Strength", y = "Count", fill = "Method", alpha = "Direction",
       title = "Count Numbers of Meta-Analytic Bayes Factors Representing Varying Degrees of Evidence Strength When the Underlying Effect is Null") +
  scale_y_continuous(labels = abs) +
  geom_vline(xintercept = c(1.5, 2.5), color = "grey", size = 0.1, linetype = "dashed") +  # Add dashed vertical lines to separate levels
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 0, hjust = 1)
  )

# Create the directory if it does not exist
# dir.create("./plots/stacked bar plot/bidirection stacked bar plot", recursive = TRUE, showWarnings = FALSE)

# Save the plot with increased width and height
ggsave(filename = "./plots/stacked bar plot/bidirection stacked bar plot/null_bidirectional_stacked.png",
       plot = plot, width = 15, height = 12)


#####0.2
# Sample data transformation
plot_data_true <- get("MABF_rates_0.2_rowise_pcutoff0.05_BFcutoff3") %>%
  mutate(rep.number = as.numeric(rep.number)) %>%
  mutate(rep.n = as.numeric(rep.n)) %>%
  mutate(AD_true_allstrong_N = AD_true_strong_N + AD_true_vstrong_N + AD_true_estrong_N) %>%
  mutate(AD_true_allstrong_A = AD_true_strong_A + AD_true_vstrong_A + AD_true_estrong_A) %>%
  select(method, rep.n, rep.number, AD_true_anecdotal_A, AD_true_moderate_A, AD_true_allstrong_A,
         AD_true_anecdotal_N, AD_true_moderate_N, AD_true_allstrong_N)

# Reshape the data and calculate median
data_median <- plot_data_true %>%
  pivot_longer(cols = starts_with("AD_"),
               names_to = c("type", "true", "strength", "direction"),
               names_sep = "_",
               values_to = "count") %>%
  mutate(strength = recode_factor(paste(true, strength, sep = "_"),
                                  true_allstrong = "strong",
                                  true_moderate = "moderate",
                                  true_anecdotal = "anecdotal"),
         direction = ifelse(direction == "A", "Correct", "Wrong")) %>%
  group_by(method, rep.n, rep.number, strength, direction) %>%
  summarise(count = median(count), .groups = 'drop')

# Define colors for each method
method_colors <- c("FEMABF" = "red", "BFbMA" = "green", "EUBF" = "blue", "iBF" = "purple")

# Create the plot with reduced bar width and dashed lines
plot <- ggplot(data_median, aes(x = strength, fill = method)) +
  geom_bar(
    aes(y = ifelse(direction == "Correct", count, 0)),
    stat = "identity",
    position = position_dodge(width = 1),
    width = 0.6,  # Reduce bar width
    alpha = 1
  ) +
  geom_bar(
    aes(y = ifelse(direction == "Wrong", -count, 0)),
    stat = "identity",
    position = position_dodge(width = 1),
    width = 0.6,  # Reduce bar width
    alpha = 0.5
  ) +
  scale_fill_manual(values = method_colors) +
  coord_flip() +
  facet_grid(rep.n ~ rep.number) +
  labs(x = "Evidence Strength", y = "Count", fill = "Method", alpha = "Direction",
       title = "Count Numbers of Meta-Analytic Bayes Factors Representing Varying Degrees of Evidence Strength When the Underlying Effect is 0.2") +
  scale_y_continuous(labels = abs) +
  geom_vline(xintercept = c(1.5, 2.5), color = "black", size = 0.1, linetype = "dashed") +  # Add dashed vertical lines to separate levels
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save the plot with increased width and height
ggsave(filename = "./plots/stacked bar plot/bidirection stacked bar plot/0.2_bidirectional_stacked.png",
       plot = plot, width = 15, height = 12)

#####0.5
# Sample data transformation
plot_data_true <- get("MABF_rates_0.5_rowise_pcutoff0.05_BFcutoff3") %>%
  mutate(rep.number = as.numeric(rep.number)) %>%
  mutate(rep.n = as.numeric(rep.n)) %>%
  mutate(AD_true_allstrong_N = AD_true_strong_N + AD_true_vstrong_N + AD_true_estrong_N) %>%
  mutate(AD_true_allstrong_A = AD_true_strong_A + AD_true_vstrong_A + AD_true_estrong_A) %>%
  select(method, rep.n, rep.number, AD_true_anecdotal_A, AD_true_moderate_A, AD_true_allstrong_A,
         AD_true_anecdotal_N, AD_true_moderate_N, AD_true_allstrong_N)

# Reshape the data and calculate median
data_median <- plot_data_true %>%
  pivot_longer(cols = starts_with("AD_"),
               names_to = c("type", "true", "strength", "direction"),
               names_sep = "_",
               values_to = "count") %>%
  mutate(strength = recode_factor(paste(true, strength, sep = "_"),
                                  true_allstrong = "strong",
                                  true_moderate = "moderate",
                                  true_anecdotal = "anecdotal"),
         direction = ifelse(direction == "A", "Correct", "Wrong")) %>%
  group_by(method, rep.n, rep.number, strength, direction) %>%
  summarise(count = median(count), .groups = 'drop')

# Define colors for each method
method_colors <- c("FEMABF" = "red", "BFbMA" = "green", "EUBF" = "blue", "iBF" = "purple")

# Create the plot with reduced bar width and dashed lines
plot <- ggplot(data_median, aes(x = strength, fill = method)) +
  geom_bar(
    aes(y = ifelse(direction == "Correct", count, 0)),
    stat = "identity",
    position = position_dodge(width = 1),
    width = 0.6,  # Reduce bar width
    alpha = 1
  ) +
  geom_bar(
    aes(y = ifelse(direction == "Wrong", -count, 0)),
    stat = "identity",
    position = position_dodge(width = 1),
    width = 0.6,  # Reduce bar width
    alpha = 0.5
  ) +
  scale_fill_manual(values = method_colors) +
  coord_flip() +
  facet_grid(rep.n ~ rep.number) +
  labs(x = "Evidence Strength", y = "Count", fill = "Method", alpha = "Direction",
       title = "Count Numbers of Meta-Analytic Bayes Factors Representing Varying Degrees of Evidence Strength When the Underlying Effect is 0.5") +
  scale_y_continuous(labels = abs) +
  geom_vline(xintercept = c(1.5, 2.5), color = "black", size = 0.1, linetype = "dashed") +  # Add dashed vertical lines to separate levels
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save the plot with increased width and height
ggsave(filename = "./plots/stacked bar plot/bidirection stacked bar plot/0.5_bidirectional_stacked.png",
       plot = plot, width = 15, height = 12)


#####Another possibility
# Sample data transformation
plot_data_null <- get("MABF_rates_null_rowise_pcutoff0.05_BFcutoff3") %>%
  mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
  mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N = 2", "N = 5", "N = 10"))) %>%
  mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n = 40", "n = 100", "n = 400"))) %>%
  mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>%
  mutate(AD_null_allstrong_N = AD_null_strong_N + AD_null_vstrong_N + AD_null_estrong_N) %>%
  mutate(AD_null_allstrong_A = AD_null_strong_A + AD_null_vstrong_A + AD_null_estrong_A) %>%
  select(method, rep.n, rep.number, AD_null_anecdotal_A, AD_null_moderate_A, AD_null_allstrong_A,
         AD_null_anecdotal_N, AD_null_moderate_N, AD_null_allstrong_N)

# Reshape the data and calculate median
data_median <- plot_data_null %>%
  pivot_longer(cols = starts_with("AD_"),
               names_to = c("type", "null", "strength", "direction"),
               names_sep = "_",
               values_to = "count") %>%
  mutate(strength = recode_factor(paste(null, strength, sep = "_"),
                                  null_allstrong = "strong",
                                  null_moderate = "moderate",
                                  null_anecdotal = "anecdotal"),
         direction = ifelse(direction == "A", "Correct", "Wrong")) %>%
  group_by(method, rep.n, rep.number, strength, direction) %>%
  summarise(count = median(count), .groups = 'drop')

# Define colors for each method
method_colors <- c("FEMABF" = "red", "BFbMA" = "green", "EUBF" = "blue", "iBF" = "purple")

# Create the plot with reduced bar width and dashed lines
ggplot(data_median, aes(x = strength, fill = method, alpha = direction)) +
  geom_bar(
    aes(y = ifelse(direction == "Correct", count, 0)),
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7  # Reduce bar width
  ) +
  geom_bar(
    aes(y = ifelse(direction == "Wrong", -count, 0)),
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7  # Reduce bar width
  ) +
  scale_fill_manual(values = method_colors) +
  scale_alpha_manual(values = c("Correct" = 1, "Wrong" = 0.5),
                     labels = c("Alternative Hypothesis", "Null Hypothesis")) +
  coord_flip() +
  facet_grid(rep.n ~ rep.number) +
  labs(
    title = "Count Numbers of Meta-Analytic Bayes Factors Representing Varying Degrees of Evidence Strength When the Underlying Effect is Null",
    x = "Evidence Strength",
    y = "Count",
    fill = "Method",
    alpha = "Direction"
  ) +
  scale_y_continuous(labels = abs) +
  geom_vline(xintercept = c(1.5, 2.5), color = "grey", size = 0.1, linetype = "dashed") +  # Add dashed vertical lines to separate levels
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 0, hjust = 1)
  )
