# This R scripts generate ROC curves for all MABF methods and the FEMA method.
# This R scripts calculate AUC for for all MABF methods and the FEMA method.
# Load necessary packages
library(ggplot2)
library(tidyverse)
library(purrr)
library(caTools)
rm(list = ls())

# Set the path to the directory containing RDS files
folder_path <- "./data files/data for ROC/"
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
# If you want to create ROC plots for underlying effect sizes that are equal to null and 0.2
rates_0.2null <- rbind(EUBF_0.2null, FEMABF_0.2null, iBF_0.2null, iBF2_0.2null, BFbMA_0.2null, FEMA_0.2null)
# Convert character to numeric
# Add a new column named true_es
rates_0.2null <- rates_0.2null %>%
  mutate(across(c(`rep number`,rep.n), as.numeric)) %>% 
  mutate(true_es = 0.2) %>% 
  relocate(true_es, .after = `true effect`) %>% 
  relocate(group, .before = `method`) %>% 
  rename(scenario = group)

# Summarize simulation scenarios showing scenario numbers and associated variables
# The first number in the title of an ROC plot corresponds to the scenario number
# For instance: 1_20_2_40_none.jpg means that scenario number is 1, original sample size is 20, 
# number of replications is 2, replication sample size is 40, and there's no bias. 
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
  # Ensure we extract single values for rep_number and rep_n
  scenario_number <- condition$scenario[1]
  true_es <- condition$true_es[1]
  orig_n <- condition$orig.n[1]
  rep_number <- condition$rep.number[1]
  rep_n <- condition$rep.n[1]
  bias <- condition$bias[1]
  # Print extracted values for debugging
  print(paste("Scenario:", scenario_number, "Replications:", rep_number, "Sample Size:", rep_n))
  title_text <- paste0(
    'ROC Curves when the number of replications is ', 
    rep_number, 
    ', replication sample size is ', 
    rep_n, 
    ', true effect is ',
    true_es,
    ', and bias level is ',
    bias)
  roc_plot <- rates_0.2null %>%
    filter(scenario == as.character(scenario_number)) %>%
    group_by(method) %>%
    ggplot(aes(x = FPR, y = TPR, color = method)) +
    geom_step() + 
    labs(title = title_text) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +  
    theme_minimal()
  # Save the plot
  file_name <- paste0("./plots/ROC/", scenario_number, "_", orig_n, "_", rep_number, "_", rep_n, "_", bias, ".jpg")
  # Save the plot
  ggsave(file_name, plot = roc_plot, width = 12, height = 8)
  print(roc_plot)
}

# Automatically generate ROC plots and save to a designated folder: ./plots/ROC/
# Number 1 to 27 are conditions under no QPR or PB scenarios
# Number 28 to 54 are conditions under medium level of QPR or PB scenarios
# Number 55 to 81 are conditions under high level of QPR or PB scenarios
walk(1:27, plot_ROC)
#walk(28:54, plot_ROC)
#walk(55:81, plot_ROC)

########
# AUC calculation
auc_results <- rates_0.2null %>%
  group_by(method) %>%
  arrange(FPR) %>%  # ensure FPR is in ascending order for correct integration
  summarise(auc = trapz(FPR, TPR))

print(auc_results)
