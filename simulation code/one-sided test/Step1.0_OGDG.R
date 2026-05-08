
#This trial set num_deltai =500, with multi-core parallel for all methods
#0415 Update: old name Step1_DataGenerationPart1.R changed to Step1_OGDG.R (original study data generation)
#0413 Update: separate generation of original and replication studies from feeding data to MABF methods
#0411 Update: added save.image at the end for HPRC output
#0410 Update: added source for rows2keep.R
#0409 Update: added my own prior for ES; changed params fro original studies (params_o); added rows_to_keep
#rm(list = ls())
######################################################
#Step 0: Load necessary packages and custom functions#
######################################################
source("./Rfiles/libraries.R")
source("./Rfiles/QRP and PB functions.R")
source("./Rfiles/helper functions.R")
source("./Rfiles/MABFs.R")
source("./Rfiles/BAMA functions.R")
source("./Rfiles/rows2keep.R")
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
# If you only save df_lists and params_o and import them to generate replications,
# you will not use seed and will obtain different generated replication results
# every time you rerun Step1.1 R files. This is because set.seed(123) is not 
# imported into Step1.1 R files, but it is imported into Step1.1 R files if you
# import Step1_500_OGDG.RData into Step1.1 R files.
# saveRDS(df_lists, file = "df_lists.RDS")
# saveRDS(params_o, file = "params_o.RDS")
# If you only want to load df_lists and parama_o in Step1.1 R files (not here), 
# use readRDS():
# df_lists<-readRDS(file = "./df_lists.RDS")
# params_o<-readRDS(file = "./params_o.RDS")
save.image("Step1.0_500_OGDG.RData")
