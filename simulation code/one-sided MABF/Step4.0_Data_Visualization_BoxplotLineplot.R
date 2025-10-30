# This is matrix-wise rates, not row-wise rates
################################################
#Step 0: Load necessary libraries and functions#
################################################
rm(list = ls()) #housekeeping
library(tidyverse)
library(scales) #for creating plots
library(ggplot2)
setwd("D:/wnk/PhD/PhD research/dissertation/code/hprc")

# Define the base folder paths
# The following for loop just iterate over the dataset in the following folders
base_folder_paths <- c(
  "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.1, EScutoff_o=0",
  "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.05, EScutoff_o=0",
  "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.01, EScutoff_o=0",
  "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.001, EScutoff_o=0"
)

# Loop through each base folder path
for (base_folder_path in base_folder_paths) {
  #####################################################
  # Step 1: Load R data files containing classification accuracy rates of MABFs
  #####################################################
  
  # List all criteria folders
  criteria_folders <- list.dirs(base_folder_path, full.names = TRUE, recursive = FALSE)
  
  # Loop through each criteria folder
  for (criteria_folder in criteria_folders) {
    # List all RDS files in the current criteria folder
    rds_files <- list.files(criteria_folder, pattern = "\\.RDS$", full.names = TRUE)
    
    # Read in all RDS files into the workspace
    for (file_path in rds_files) {
      # Extract the base name without the extension
      file_name <- tools::file_path_sans_ext(basename(file_path))
      
      # Create a variable with the name of the file and assign the data from readRDS
      assign(file_name, readRDS(file_path), envir = .GlobalEnv)
    }
  }
  
  #######################################
  # Step 2: Combine datasets for analyses
  #######################################
  
  # Define a function to combine datasets for different criteria and underlying effects
  combine_datasets <- function(effect_size, criteria_number) {
    rbind(
      get(paste0("rates_FEMABF_", effect_size, "null_c", criteria_number))[, -c(11, 12, 13, 18, 19, 20)], #exclude ADR related columns
      get(paste0("rates_BFbMA_", effect_size, "null_c", criteria_number))[, -c(11, 12, 13, 18, 19, 20)],
      get(paste0("rates_EUBF_", effect_size, "null_c", criteria_number))[, -c(11, 12, 13, 18, 19, 20)],
      get(paste0("rates_iBF_", effect_size, "null_c", criteria_number))[, -c(11, 12, 13, 18, 19, 20)],
      get(paste0("rates_FEMA_", effect_size, "null_c", criteria_number))
    ) %>%
      mutate(Method = rep(c("FEMABF", "BFbMA", "EUBF", "iBF", "FEMA"), each = 81),
             criteria = rep(paste0("Criteria ", criteria_number), 81 * 5)) %>%
      select(Method, criteria, everything()) %>%
      rownames_to_column()
  }
  
  rates_0.2null <- do.call(rbind, lapply(1:3, function(c) combine_datasets("0.2", c)))
  rates_0.5null <- do.call(rbind, lapply(1:3, function(c) combine_datasets("0.5", c)))
  
  # Specify the names of the objects you want to keep
  keep_objects <- c("rates_0.2null", "rates_0.2null_c1", "rates_0.2null_c2", "rates_0.2null_c3",
                    "rates_0.5null", "rates_0.5null_c1", "rates_0.5null_c2", "rates_0.5null_c3", "base_folder_path")
  
  # List all objects in the global environment
  all_objects <- ls()
  
  # Determine which objects to remove (all objects except the ones you want to keep)
  remove_objects <- setdiff(all_objects, keep_objects)
  
  # Remove the specified objects from the global environment
  rm(list = remove_objects)
  rm(all_objects, remove_objects)
  
  ###########################
  # Step 3: Data manipulation
  ###########################
  
  # Define a function to manipulate data
  manipulate_data <- function(data) {
    data %>%
      mutate(rep.n.level = factor(rep.n, levels = c(40, 100, 400), labels = c("40", "100", "400"))) %>%
      mutate(orig.n.level = factor(orig.n, levels = c(20, 50, 200), labels = c("20", "50", "200"))) %>%
      mutate(rep.number.level = factor(rep.number, levels = c(2, 5, 10), labels = c("2", "5", "10"))) %>%
      mutate(bias.level = factor(qrpEnv, levels = c("none", "medium", "high"), labels = c(1, 2, 3))) %>%
      mutate(criteria.level = factor(criteria, levels = c("Criteria 1", "Criteria 2", "Criteria 3"), labels = c(1, 2, 3))) %>%
      select(-rep.n.level, everything()) %>%
      relocate(rep.n.level, .after = rep.n) %>%
      select(-orig.n.level, everything()) %>%
      relocate(orig.n.level, .after = orig.n) %>%
      select(-rep.number.level, everything()) %>%
      relocate(rep.number.level, .after = rep.number) %>%
      select(-bias.level, everything()) %>%
      relocate(bias.level, .after = censorFunc) %>%
      select(-criteria.level, everything()) %>%
      relocate(criteria.level, .after = criteria) %>%
      mutate(across(c(Method, criteria), as.factor))
  }
  
  rates_0.2null <- manipulate_data(rates_0.2null) #These two datasets are what we use to create plots
  rates_0.5null <- manipulate_data(rates_0.5null)
  
  #######################
  # Step 4: Create plots
  #######################
  
  # Load plot functions
  source("./helper functions_XXR.R")
  
  # Labels to show in the plot
  rep.number_labeller <- as_labeller(c(`2` = "N[rep] == 2", `5` = "N[rep] == 5", `10` = "N[rep] == 10"), label_parsed)
  bias.level_labeller <- as_labeller(c(`1` = "Bias == low", `2` = "Bias == medium", `3` = "Bias == high"), label_parsed)
  orig.n_labeller <- as_labeller(c(`20` = "n[orig] == 20", `50` = "n[orig] == 50", `200` = "n[orig] == 200"), label_parsed)
  criteria_labeller_fxorig <- as_labeller(c(`Criteria 1` = "BF[10] > 3 *',' ~ italic(p[rep]) < .05", `Criteria 2` = "BF[10] > 10 *',' ~ italic(p[rep]) < .01", `Criteria 3` = "BF[10] > 30 *',' ~ italic(p[rep]) < .001"), label_parsed)
  #criteria_labeller_fxrep <- as_labeller(c(`Criteria 1` = "p[orig] == .05", `Criteria 2` = "p[orig] == .01", `Criteria 3` = "p[orig] == .001"), label_parsed) #fixed replication criteria
  
  # Set up parameters to create and save the plot
  # Note: The third variable (first two are rep.n and rep.number) can be criteria.level, bias.level, and orig.n.level
  # Note: TPR plots are the same across all 4 folders (base_folder_path) because TPR doesn't depend on original study p value (same goes to FPR)
  data <- rates_0.2null #can be changed to rates_0.5null
  trueES <- 0.2 #can be changed to 0.5
  ratetype <- "TPR"#can be changed to TPR, FPR, TNR, FNR, TSR, FSR, TFR, FFR, SAR, SCR
  ratetype_char = "True Positive Rate"
  use_third_var = TRUE
  third_var <- "bias.level"
  thirdvar_labeller = bias.level_labeller
  apply_filter = TRUE
  filter_var = "criteria.level" #When filter is turned off, this argument won't activate
  filter_level = 1 #When filter is turned off, this argument won't activate
  
  plot <- plot_line(data = data, trueES = trueES, ratetype = ratetype, ratetype_char = ratetype_char, use_third_var = use_third_var, third_var = third_var, thirdvar_labeller = thirdvar_labeller, apply_filter = apply_filter, filter_var = filter_var, filter_level = filter_level)
  
  # Define plot directory
  # Don't forget to change the directory to line plot if you create line plots
  plot_dir <- "./plots/line plot/fixed original cutoff"
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  # Save the plot with the new naming convention
  # Don't forget to change the name to line plot if you create line plots
  # plot_filename <- sprintf("%s/boxplot_trueES=%.1f, ratetype=%s, thirdvar=%s, fixed%s.jpg", plot_dir, trueES, ratetype, third_var, basename(base_folder_path))
  # If you use a filter, use the following line of code instead of the above line:
  plot_filename <- sprintf("%s/trendplot_trueES=%.1f, ratetype=%s, thirdvar=%s, filter_var=%s, filter_level=%s, fixed%s.jpg", plot_dir, trueES, ratetype, third_var, filter_var, filter_level, basename(base_folder_path))
  ggsave(filename = plot_filename, plot = plot, width = 10, height = 6)
}

rm(list = ls())



# ###################################The following script investigates ADR##################################
# ###################################It uses the same script as above, excluding FEMA#######################
# ################################################
# #Step 0: Load necessary libraries and functions#
# ################################################
# rm(list = ls()) #housekeeping
# library(tidyverse)
# library(scales) #for creating plots
# library(ggplot2)
# setwd("D:/wnk/PhD/PhD research/dissertation/code/hprc")
# 
# # Define the base folder paths
# # The following for loop just iterate over the dataset in the following folders
# base_folder_paths <- c(
#   "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.1, EScutoff_o=0",
#   "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.05, EScutoff_o=0",
#   "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.01, EScutoff_o=0",
#   "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.001, EScutoff_o=0"
# )
# 
# # Loop through each base folder path
# for (base_folder_path in base_folder_paths) {
#   #####################################################
#   # Step 1: Load R data files containing classification accuracy rates of MABFs
#   #####################################################
#   
#   # List all criteria folders
#   criteria_folders <- list.dirs(base_folder_path, full.names = TRUE, recursive = FALSE)
#   
#   # Loop through each criteria folder
#   for (criteria_folder in criteria_folders) {
#     # List all RDS files in the current criteria folder
#     rds_files <- list.files(criteria_folder, pattern = "\\.RDS$", full.names = TRUE)
#     
#     # Read in all RDS files into the workspace
#     for (file_path in rds_files) {
#       # Extract the base name without the extension
#       file_name <- tools::file_path_sans_ext(basename(file_path))
#       
#       # Create a variable with the name of the file and assign the data from readRDS
#       assign(file_name, readRDS(file_path), envir = .GlobalEnv)
#     }
#   }
#   
#   #######################################
#   # Step 2: Combine datasets for analyses
#   #######################################
#   
#   # Define a function to combine datasets for different criteria and underlying effects
#   combine_datasets <- function(effect_size, criteria_number) {
#     rbind(
#       get(paste0("rates_FEMABF_", effect_size, "null_c", criteria_number)), #Including ADR related columns
#       get(paste0("rates_BFbMA_", effect_size, "null_c", criteria_number)),
#       get(paste0("rates_EUBF_", effect_size, "null_c", criteria_number)),
#       get(paste0("rates_iBF_", effect_size, "null_c", criteria_number))
#     ) %>%
#       mutate(Method = rep(c("FEMABF", "BFbMA", "EUBF", "iBF"), each = 81),
#              criteria = rep(paste0("Criteria ", criteria_number), 81 * 4)) %>%
#       select(Method, criteria, everything()) %>%
#       rownames_to_column()
#   }
#   
#   rates_0.2null <- do.call(rbind, lapply(1:3, function(c) combine_datasets("0.2", c)))
#   rates_0.5null <- do.call(rbind, lapply(1:3, function(c) combine_datasets("0.5", c)))
#   
#   # Specify the names of the objects you want to keep
#   keep_objects <- c("rates_0.2null", "rates_0.2null_c1", "rates_0.2null_c2", "rates_0.2null_c3",
#                     "rates_0.5null", "rates_0.5null_c1", "rates_0.5null_c2", "rates_0.5null_c3", "base_folder_path")
#   
#   # List all objects in the global environment
#   all_objects <- ls()
#   
#   # Determine which objects to remove (all objects except the ones you want to keep)
#   remove_objects <- setdiff(all_objects, keep_objects)
#   
#   # Remove the specified objects from the global environment
#   rm(list = remove_objects)
#   rm(all_objects, remove_objects)
#   
#   ###########################
#   # Step 3: Data manipulation
#   ###########################
#   
#   # Define a function to manipulate data
#   manipulate_data <- function(data) {
#     data %>%
#       mutate(rep.n.level = factor(rep.n, levels = c(40, 100, 400), labels = c("40", "100", "400"))) %>%
#       mutate(orig.n.level = factor(orig.n, levels = c(20, 50, 200), labels = c("20", "50", "200"))) %>%
#       mutate(rep.number.level = factor(rep.number, levels = c(2, 5, 10), labels = c("2", "5", "10"))) %>%
#       mutate(bias.level = factor(qrpEnv, levels = c("none", "medium", "high"), labels = c(1, 2, 3))) %>%
#       mutate(criteria.level = factor(criteria, levels = c("Criteria 1", "Criteria 2", "Criteria 3"), labels = c(1, 2, 3))) %>%
#       select(-rep.n.level, everything()) %>%
#       relocate(rep.n.level, .after = rep.n) %>%
#       select(-orig.n.level, everything()) %>%
#       relocate(orig.n.level, .after = orig.n) %>%
#       select(-rep.number.level, everything()) %>%
#       relocate(rep.number.level, .after = rep.number) %>%
#       select(-bias.level, everything()) %>%
#       relocate(bias.level, .after = censorFunc) %>%
#       select(-criteria.level, everything()) %>%
#       relocate(criteria.level, .after = criteria) %>%
#       mutate(across(c(Method, criteria), as.factor))
#   }
#   
#   rates_0.2null <- manipulate_data(rates_0.2null)
#   rates_0.5null <- manipulate_data(rates_0.5null)
#   
#   #######################
#   # Step 4: Create plots
#   #######################
#   
#   # Load plot functions
#   source("./helper functions_XXR.R")
#   
#   # Labels to show in the plot
#   rep.number_labeller <- as_labeller(c(`2` = "N[rep] == 2", `5` = "N[rep] == 5", `10` = "N[rep] == 10"), label_parsed)
#   bias.level_labeller <- as_labeller(c(`1` = "Bias == low", `2` = "Bias == medium", `3` = "Bias == high"), label_parsed)
#   orig.n_labeller <- as_labeller(c(`20` = "n[orig] == 20", `50` = "n[orig] == 50", `200` = "n[orig] == 200"), label_parsed)
#   criteria_labeller_fxorig <- as_labeller(c(`Criteria 1` = "BF[10] > 3 *',' ~ italic(p[rep]) < .05", `Criteria 2` = "BF[10] > 10 *',' ~ italic(p[rep]) < .01", `Criteria 3` = "BF[10] > 30 *',' ~ italic(p[rep]) < .001"), label_parsed)
#   #criteria_labeller_fxrep <- as_labeller(c(`Criteria 1` = "p[orig] == .05", `Criteria 2` = "p[orig] == .01", `Criteria 3` = "p[orig] == .001"), label_parsed) #fixed replication criteria
#   
#   # Set up parameters to create and save the plot
#   # Note: The third variable (first two are rep.n and rep.number) can be criteria.level, bias.level, and orig.n.level
#   # Note: TPR plots are the same across all 4 folders (base_folder_path) because TPR doesn't depend on original study p value (same goes to FPR)
#   data <- rates_0.2null #can be changed to rates_0.5null
#   trueES <- 0.2 #can be changed to 0.5
#   ratetype <- "ADR"#can be changed to TPR, FPR, SAR, SCR,...
#   ratetype_char = "Anecdotal Evidence Rate"
#   use_third_var = TRUE
#   third_var <- "criteria"
#   thirdvar_labeller = criteria_labeller_fxorig
#   apply_filter = FALSE
#   filter_var = "bias.level" #filter is turned off, so this argument won't activate
#   filter_level = 3 #filter is turned off, so this argument won't activate
#   
#   plot <- plot_line(data = data, trueES = trueES, ratetype = ratetype, ratetype_char = ratetype_char, use_third_var = use_third_var, third_var = third_var, thirdvar_labeller = thirdvar_labeller, apply_filter = apply_filter, filter_var = filter_var, filter_level = filter_level)
#   
#   # Define plot directory
#   plot_dir <- "./plots/box plot/fixed original cutoff"
#   if (!dir.exists(plot_dir)) {
#     dir.create(plot_dir, recursive = TRUE)
#   }
#   
#   # Save the plot with the new naming convention
#   plot_filename <- sprintf("%s/boxplot_trueES=%.1f, ratetype=%s, thirdvar=%s, fixed%s.jpg", plot_dir, trueES, ratetype, third_var, basename(base_folder_path))
#   # If you use a filter, use the following line of code:
#   # plot_filename <- sprintf("%s/boxplot_trueES=%.1f, ratetype=%s, thirdvar=%s, filter_var=%s, filter_level=%s, fixed%s.jpg", plot_dir, trueES, ratetype, third_var, filter_var, filter_level, basename(base_folder_path))
#   ggsave(filename = plot_filename, plot = plot, width = 10, height = 6)
# }
# 
# rm(list = ls())
