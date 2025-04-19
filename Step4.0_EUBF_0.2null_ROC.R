# This R script prepares data to create ROC curves for EUBF outcomes.
#####################################
#MABF method: EUBF; dataset: 0.2null#
#####################################
# Load necessary data and packages
library(tidyverse)
library(purrr)
library(furrr)
source("helper functions_XXR.R")
EUBF_lists_0.2null_regrouped <- readRDS("EUBF_lists_0.2null_regrouped.RDS")

# An argument that will be used in ratesCalc_ROC()
num_deltai = 500

# Define the thresholds to run through to calculate TPR and FPR
BFcutoff <- seq(from = 0.005, to = 100, by = 0.01)

# Initialize an empty dataframe to store the results
EUBF_TPFP4ROC_0.2null <- data.frame()

# Set up parallel processing plan
plan(multisession, workers = parallel::detectCores())

# Apply ratesCalc_ROC and thresholds to every sublist of MABF values
list_of_sublists <- EUBF_lists_0.2null_regrouped

# Process all sublists in parallel
EUBF_TPFP4ROC_0.2null <- future_map_dfr(names(list_of_sublists), function(sublist_name) {
  message("Processing sublist: ", sublist_name)
  sublist <- list_of_sublists[[sublist_name]]
  process_sublist(sublist, sublist_name)
})

# End parallel processing plan
plan(sequential)

# Reorder columns and break column names into factors
EUBF_TPFP4ROC_0.2null <- EUBF_TPFP4ROC_0.2null %>% 
  select(Sublist, Matrix, TPR, FPR, TNR, FNR, TP, FN, FP, TN, threshold) %>% 
  separate(col = "Sublist", into = c("true effect", "orig.n", "QRP level", "PB level"), sep = "_") %>% 
  separate(col = "Matrix", into = c("rep number","rep.n"), sep = "_") 

# Assign group number
EUBF_TPFP4ROC_0.2null <- EUBF_TPFP4ROC_0.2null %>% 
  mutate(group = (row_number()-1) %/% length(BFcutoff) + 1)

# Add MABF method name column
EUBF_TPFP4ROC_0.2null <- EUBF_TPFP4ROC_0.2null %>%
  mutate(method = "EUBF") %>% 
  relocate(method)

# Save data as RDS file
saveRDS(EUBF_TPFP4ROC_0.2null, "EUBF_0.2null.RDS")

