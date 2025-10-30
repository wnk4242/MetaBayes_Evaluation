# This script will create a smaller sample size for Shiny demonstration (only 500 thresholds).
# This script will create true/false success, true/false failure, true/false positive, true/false negative
# The successes, failures, positives, negatives will be used to calculate rates for visualization.
######################################
#MABF method: iBF; dataset: 0.5null#
######################################
# Load necessary data and packages
library(tidyverse)
library(ggplot2)
library(scales) # scale_y_continuous(labels = percent) +
library(purrr)
library(furrr)

# Load needed helper functions
source("./Rfiles/helper functions_XXR.R")

options(future.globals.maxSize= 2097152000) #To avoid memory limit issues in the future package: 2000MB*1024^2

# Read in needed RDS files
iBF_lists_0.5null_regrouped_deltap <- readRDS("./RDS/MABF4ROC/4TSFS/iBF_lists_0.5null_regrouped_deltap.RDS")

# Set up parameters
num_deltai = 500
rep_times = 500
BFcutoff_null <- seq(from = 0.001, by = 0.004, length.out = 250)
BFcutoff_alt <- seq(from = 1.01, by = 0.12, length.out = 250)
BFcutoff <- c(BFcutoff_null, BFcutoff_alt) # min BF in the dataset is 0.00334158, so to make TPR=FPR=1, the lowest cutoff must be smaller than that
#BFcutoff <- c(BFcutoff, Inf) #max BF is 7.676577e+80, to make TPR=FPR=0, the highest cutoff must be larger than that
pcutoff <- c(0.01, 0.05, 0.1)
EScutoff <- c(0)

#Create an empty dataframe to contain data
iBF_TFPS4Viz_0.5null <- data.frame()

# Create a data frame of all combinations of cutoffs using expand.grid
cutoff_combinations <- expand.grid(BFcutoff = BFcutoff, pcutoff = pcutoff, EScutoff = EScutoff)
# Convert the data frame to a list of lists
cutoff <- split(cutoff_combinations, seq(nrow(cutoff_combinations)))

# Set up parallel processing plan
plan(multisession, workers = parallel::detectCores())

# Apply ratesCalc_ROC and thresholds to every sublist of MABF values
list_of_sublists <- iBF_lists_0.5null_regrouped_deltap

# Process all sublists in parallel
iBF_TFPS4Viz_0.5null <- future_map_dfr(names(list_of_sublists), function(sublist_name) {
  message("Processing sublist: ", sublist_name)
  sublist <- list_of_sublists[[sublist_name]]
  process_sublist_AiO(sublist, sublist_name)
})

# End parallel processing plan
plan(sequential)

# Reorder columns and break column names into factors
iBF_TFPS4Viz_0.5null <- iBF_TFPS4Viz_0.5null %>% 
  select(Sublist, Matrix, threshold_BF, threshold_p, threshold_ES, 
         Successes,Failures,Successes1,Failures1,Successes2,Failures2,
         TS1, TS2, FF1, FF2, TF1, TF2, FS1, FS2,
         TP, TN, FP, FN) %>%
  separate(col = "Sublist", into = c("true effect", "orig.n", "QRP level", "PB level"), sep = "_") %>% 
  separate(col = "Matrix", into = c("rep number","rep.n"), sep = "_") 

# Add MABF method name column
iBF_TFPS4Viz_0.5null <- iBF_TFPS4Viz_0.5null %>%
  mutate(method = "iBF") %>% 
  relocate(method)

# Assign group number
iBF_TFPS4Viz_0.5null <- iBF_TFPS4Viz_0.5null %>% 
  mutate(group = (row_number()-1) %/% length(BFcutoff) + 1) %>% 
  select(method, group, everything())

# Save data as RDS file
saveRDS(iBF_TFPS4Viz_0.5null, "./RDS/MABFanalyses/matrix-wise/rates4Viz/iBF_TFPS4Viz_0.5null.RDS")
