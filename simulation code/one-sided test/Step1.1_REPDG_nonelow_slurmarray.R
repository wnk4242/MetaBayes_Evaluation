# This code integrates the Slurm array command to divide the workload to nine individual arrays.
# This code only generate replications of original studies not affected by QRP or PB
######################################################
#Step 0: Load necessary packages and custom functions#
######################################################
# Load workspace of original studies (set.seed will be effective)
# Load original studies data where 500 original studies are sampled for a condition
load("./RDS/OGDG/Step1.0_500_OGDG.RData")
source("./Rfiles/libraries.R")
source("./Rfiles/QRP and PB functions.R")
source("./Rfiles/helper functions.R")
source("./Rfiles/MABFs.R")
source("./Rfiles/BAMA functions.R")
source("./Rfiles/rows2keep.R")
# Use parallel processing
# Determine the number of cores to use
num_cores <- parallel::detectCores()
# plan() is from the future package
plan(multisession, workers = num_cores)
####################################################################
#Step 1: Divide the list of original study sublists into nine parts#
####################################################################
# Read the environment variable to determine which segment to process
segment_index <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# Only choose og studies not affected by QRP or PB
# orig_nonelow_lists <- df_lists[1:9]
# Define the list segments based on the index
lists_segments <- list(df_lists[1],df_lists[2],df_lists[3],
                       df_lists[4],df_lists[5],df_lists[6],
                       df_lists[7],df_lists[8],df_lists[9])
params_o_slices <- list(1,2,3,4,5,6,7,8,9) # Corresponding indices in params_o
# Select the segment
orig_nonelow_lists <- lists_segments[[segment_index]]
#######################################
#Step 2: Generate replication studies#
#######################################
nonelow_rep_results <- list()
nonelow_IDs <- list()
nonelow_rep_params <- list()
# To show progress message:
# Define the total number of iterations at the beginning
total_iterations <- length(orig_nonelow_lists)
# Report progress every 10%
progress_threshold <- 10
for (i in seq_along(orig_nonelow_lists)) {
  delta_i_values <- orig_nonelow_lists[[i]][, 10]
  num_reps <- c(2,5,10) # number of replications
  samplesize_reps <- c(40,100,400) # sample size (for one arm, total sample size for a replication is twice as large)
  rep_times <- 500 # repetition times, not how many replication studies, which is num_reps
  params_r <- expand.grid(times = rep_times, delta_i = delta_i_values,
                          num_deltai = num_reps,
                          fixed.n = samplesize_reps)
  results_lists <- purrr::pmap(params_r, pass_params)
  IDs <- purrr::pmap(params_r, naming_rep)
  names(results_lists) <- IDs
  nonelow_IDs[[i]] <- IDs
  nonelow_rep_results[[i]] <- results_lists
  nonelow_rep_params[[i]] <- params_r
  # Progress reporting
  current_progress <- round((i / total_iterations) * 100)
  if (current_progress >= progress_threshold) {
    message(sprintf("Progress: %d%% completed", progress_threshold))
    progress_threshold <- progress_threshold + 10
  }
  # Ensure the 100% completion is captured
  if (i == total_iterations && current_progress == 100) {
    message("Progress: 100% completed")
  }
}
params_o_nonelow <- params_o[params_o_slices[[segment_index]], ]
IDs <- purrr::pmap(params_o_nonelow,naming_orig)
names(nonelow_rep_results) <- IDs
# Reset to sequential plan
plan(sequential)
# Save results
saveRDS(nonelow_rep_results, sprintf("./RDS/new/nonelow_rep_results_part%d.RDS", segment_index))
names(nonelow_IDs) <- names(orig_nonelow_lists)
saveRDS(nonelow_IDs, sprintf("./RDS/new/nonelow_IDs_part%d.RDS", segment_index))
