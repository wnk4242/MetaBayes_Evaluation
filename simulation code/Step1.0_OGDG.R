# This R script generates original study results.
# OGDG stands for original data generation.
######################################################
#Step 0: Load necessary packages and custom functions#
######################################################
source("libraries.R")
source("QRP and PB functions.R")
source("helper functions.R")
source("MABFs.R")
source("BAMA functions.R")
source("rows2keep.R")
##################################
#Step 1: Generate original studies#
##################################
# Define parameter levels for original studies
# Use parallel processing
# Determine the number of cores to use (total cores - 1)
num_cores <- parallel::detectCores()
# Set up parallel processing to leave one core free
# plan() is from the future package
plan(multisession, workers = num_cores)

set.seed(123)
params_o <- expand.grid(
  num_deltai = 500,
  delta = c(0, 0.2, 0.5),
  tau = c(0.02, 0.05, 0.125),
  fixed.n = c(20, 50, 200),
  qrpEnv = c("none", "medium","high"),
  censorFunc = c("low", "medium","high"),
  stringsAsFactors = FALSE
)
# Remove some rows that do not represent realistic scenarios
params_o <- params_o[rows_to_keep,] #use rows2keep file to create rows_to_keep
# Remove unwanted rows
row.names(params_o) <- NULL
params_o <- params_o[-c(10:36,46:72),]
row.names(params_o) <- NULL
# Apply simORIGs to parameter combinations
df_lists <- purrr::pmap(params_o, simORIGs)
all_orig_results <- do.call(rbind, df_lists)
# Convert df_lists to named lists for easier access if needed
IDs <- purrr::pmap(params_o,naming_orig)
names(df_lists) <- IDs
# Reset to sequential plan
plan(sequential)
# Save results
save.image("Step1.0_500_OGDG.RData")
