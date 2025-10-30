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
rates_0.2null <- rbind(EUBF_TPFP4ROC_0.2null, FEMABF_TPFP4ROC_0.2null, iBF_TPFP4ROC_0.2null, BFbMA_TPFP4ROC_0.2null, FEMA_TPFP4ROC_0.2null) 
# Convert character to numeric
rates_0.2null <- rates_0.2null %>%
  mutate(across(c(`rep number`,rep.n), as.numeric)) %>% 
  mutate(true_es = 0.2) %>% 
  relocate(true_es, .after = `true effect`) %>% 
  relocate(group, .before = `method`) %>% 
  rename(scenario = group)

saveRDS(rates_0.2null, file = "./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.2null.RDS")

# Combine rate outcomes of different methods in order to compare them in ROC plot
rates_0.5null <- rbind(EUBF_TPFP4ROC_0.5null, FEMABF_TPFP4ROC_0.5null, iBF_TPFP4ROC_0.5null, BFbMA_TPFP4ROC_0.5null, FEMA_TPFP4ROC_0.5null) 
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
  file_name <- paste0("./plots/ROC/0.2_ROC_group/", scenario_number, "_", orig_n, "_", rep_number, "_", rep_n, "_", bias, ".jpg")
  # Save the plot
  ggsave(file_name, plot = roc_plot, width = 12, height = 8)
  print(roc_plot)
}

# Automatically generate ROC plots and save to designated folder
# No need to run all 81 conditions because QRP and PB do not affect FPR and TPR
# Just run the conditions under no QPR or PB scenarios
walk(1:3, plot_ROC)
# You can also use a for loop to do the same thing of course
# for (group_number in 1:2) {
#   plot_ROC(group_number)
# }

