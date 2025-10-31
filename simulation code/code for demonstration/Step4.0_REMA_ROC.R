# This R script prepares data to create ROC curves for REMA outcomes.
# Note: "FEMA" was accidentally used as the object name for the random-effects meta-analysis method.
################################
#Method: REMA; dataset: 0.2null#
################################
# Load necessary data and packages
library(tidyverse)
library(purrr)
library(furrr)
source("helper functions_XXR.R")
FEMA_lists_0.2null_regrouped <- readRDS("FEMA_lists_0.2null_regrouped.RDS")
options(future.globals.maxSize= 891289600) #To avoid memory limit issues in the future package

# An argument that will be used in ratesCalc_ROC()
num_deltai = 500

# Define the thresholds to run through to calculate TPR and FPR
pcutoff_r <- seq(from = 0.00001, to = 1, by = 0.0001)

# Initialize an empty dataframe to store the results
FEMA_TPFP4ROC_0.2null <- data.frame()

# Set up parallel processing plan
plan(multisession, workers = parallel::detectCores())

# Apply ratesCalcMA_ROC and thresholds to every sublist of MABF values
list_of_sublists <- FEMA_lists_0.2null_regrouped

# Process all sublists in parallel
FEMA_TPFP4ROC_0.2null <- future_map_dfr(names(list_of_sublists), function(sublist_name) {
  message("Processing sublist: ", sublist_name)
  sublist <- list_of_sublists[[sublist_name]]
  processMA_sublist(sublist, sublist_name)
})

# End parallel processing plan
plan(sequential)

# Reorder columns and break column names into factors
FEMA_TPFP4ROC_0.2null <- FEMA_TPFP4ROC_0.2null %>% 
  select(Sublist, Matrix, TPR, FPR, TNR, FNR, TP, FN, FP, TN, threshold) %>% 
  separate(col = "Sublist", into = c("true effect", "orig.n", "QRP level", "PB level"), sep = "_") %>% 
  separate(col = "Matrix", into = c("rep number","rep.n"), sep = "_") 

# Assign group number
FEMA_TPFP4ROC_0.2null <- FEMA_TPFP4ROC_0.2null %>% 
  mutate(group = (row_number()-1) %/% length(pcutoff_r) + 1)

# Add MABF method name column
FEMA_TPFP4ROC_0.2null <- FEMA_TPFP4ROC_0.2null %>%
  mutate(method = "FEMA") %>% 
  relocate(method)

# Save data as RDS file
saveRDS(FEMA_TPFP4ROC_0.2null, "FEMA_0.2null.RDS")

