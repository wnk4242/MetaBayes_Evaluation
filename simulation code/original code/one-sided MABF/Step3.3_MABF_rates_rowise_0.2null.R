# Load necessary packages
library(ggplot2)
library(dplyr)
library(tidyr)
source("./Rfiles/helper functions_XXR.R")

# Set the path to the directory containing RDS files
folder_path <- "./RDS/MABF4ROC/4TSFS" #These dataset contain original study p and ES for calculating TSR, FSR, etc.
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

# Function to loop through MABF datasets using row-wise ratesCalc function to deliver all types of rates
process_datasets <- function(BFcutoff, pcutoff) {
  ########
  #FEMABF#
  ########
  # Specify which RDS file is being used
  dataset <- FEMABF_lists_0.2null_regrouped_deltap
  
  # Apply ratesCalc_rowise function to obtain row-wise rates
  FEMABF_rateslists_0.2null <- ratesCalc_rowise(dataset, pcutoff = pcutoff, EScutoff = 0, BFcutoff = BFcutoff)
  
  # Save the results in their respective objects for analysis
  FEMABF_rates_null <- FEMABF_rateslists_0.2null$results_null_df %>%
    mutate(method = rep("FEMABF", n())) %>%
    select(method, everything())
  
  FEMABF_rates_0.2 <- FEMABF_rateslists_0.2null$results_true_df %>%
    mutate(method = rep("FEMABF", n())) %>%
    select(method, everything())
  
  #######
  #BFbMA#
  #######
  # Specify which RDS file is being used
  dataset <- BFbMA_lists_0.2null_regrouped_deltap
  
  # Apply ratesCalc_rowise function to obtain row-wise rates
  BFbMA_rateslists_0.2null <- ratesCalc_rowise(dataset, pcutoff = pcutoff, EScutoff = 0, BFcutoff = BFcutoff)
  
  # Save the results in their respective objects for analysis
  BFbMA_rates_null <- BFbMA_rateslists_0.2null$results_null_df %>%
    mutate(method = rep("BFbMA", n())) %>%
    select(method, everything())
  
  BFbMA_rates_0.2 <- BFbMA_rateslists_0.2null$results_true_df %>%
    mutate(method = rep("BFbMA", n())) %>%
    select(method, everything())
  
  ######
  #EUBF#
  ######
  # Specify which RDS file is being used
  dataset <- EUBF_lists_0.2null_regrouped_deltap
  
  # Apply ratesCalc_rowise function to obtain row-wise rates
  EUBF_rateslists_0.2null <- ratesCalc_rowise(dataset, pcutoff = pcutoff, EScutoff = 0, BFcutoff = BFcutoff)
  
  # Save the results in their respective objects for analysis
  EUBF_rates_null <- EUBF_rateslists_0.2null$results_null_df %>%
    mutate(method = rep("EUBF", n())) %>%
    select(method, everything())
  
  EUBF_rates_0.2 <- EUBF_rateslists_0.2null$results_true_df %>%
    mutate(method = rep("EUBF", n())) %>%
    select(method, everything())
  
  ######
  #iBF#
  ######
  # Specify which RDS file is being used
  dataset <- iBF_lists_0.2null_regrouped_deltap
  
  # Apply ratesCalc_rowise function to obtain row-wise rates
  iBF_rateslists_0.2null <- ratesCalc_rowise(dataset, pcutoff = pcutoff, EScutoff = 0, BFcutoff = BFcutoff)
  
  # Save the results in their respective objects for analysis
  iBF_rates_null <- iBF_rateslists_0.2null$results_null_df %>%
    mutate(method = rep("iBF", n())) %>%
    select(method, everything())
  
  iBF_rates_0.2 <- iBF_rateslists_0.2null$results_true_df %>%
    mutate(method = rep("iBF", n())) %>%
    select(method, everything())
  
  
  # ######
  # #iBF2# here
  # ######
  # # Specify which RDS file is being used
  # dataset <- iBF2_lists_0.2null_regrouped_deltap
  # 
  # # Apply ratesCalc_rowise function to obtain row-wise rates
  # iBF2_rateslists_0.2null <- ratesCalc_rowise(dataset, pcutoff = pcutoff, EScutoff = 0, BFcutoff = BFcutoff)
  # 
  # # Save the results in their respective objects for analysis
  # iBF2_rates_null <- iBF2_rateslists_0.2null$results_null_df %>%
  #   mutate(method = rep("iBF2", n())) %>%
  #   select(method, everything())
  # 
  # iBF2_rates_0.2 <- iBF2_rateslists_0.2null$results_true_df %>%
  #   mutate(method = rep("iBF2", n())) %>%
  #   select(method, everything())
  
  # Combine MABF methods row-wise rates results #here
  MABF_rates_0.2_rowise <- rbind(FEMABF_rates_0.2, BFbMA_rates_0.2, EUBF_rates_0.2, iBF_rates_0.2)
  MABF_rates_null_rowise <- rbind(FEMABF_rates_null, BFbMA_rates_null, EUBF_rates_null, iBF_rates_null)
  
  # Save calculated rates
  saveRDS(MABF_rates_0.2_rowise, paste0("./RDS/MABFanalyses/row-wise/MABF_rates_0.2_rowise_pcutoff", pcutoff, "_BFcutoff", BFcutoff, ".RDS"))
  saveRDS(MABF_rates_null_rowise, paste0("./RDS/MABFanalyses/row-wise/MABF_rates_null_rowise_pcutoff", pcutoff, "_BFcutoff", BFcutoff, ".RDS"))
}

# Loop over different cutoff values for anecdotal evidence
BFcutoff_values <- c(1, 3, 10, 30)
pcutoff_values <- c(0.01, 0.05) #here
for (BFcutoff in BFcutoff_values) {
  for (pcutoff in pcutoff_values) {
    process_datasets(BFcutoff, pcutoff)
  }
}
