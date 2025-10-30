# This code integrates the Slurm array command to divide the workload to nine individual arrays.
# This code only generate replications of original studies affected by high QRP or PB
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
# Only choose og studies affected by high level QRP and PB
# orig_hihi_lists <- df_lists[19:27]
# Define the list segments based on the index
lists_segments <- list(df_lists[19], df_lists[20], df_lists[21],
                       df_lists[22], df_lists[23], df_lists[24],
                       df_lists[25], df_lists[26], df_lists[27])
params_o_slices <- list(19,20,21,22,23,24,25,26,27)  # Corresponding indices in params_o
# Select the segment
orig_hihi_lists <- lists_segments[[segment_index]]
#######################################
#Step 2: Generate replication studies#
#######################################
hihi_rep_results <- list()
hihi_IDs <- list()
hihi_rep_params <- list()
# To show progress message:
# Define the total number of iterations at the beginning
total_iterations <- length(orig_hihi_lists)
# Report progress every 10%
progress_threshold <- 10
for (i in seq_along(orig_hihi_lists)) {
  delta_i_values <- orig_hihi_lists[[i]][, 10]
  num_reps <- c(2,5,10) # number of replications
  samplesize_reps <- c(40,100,400) # sample size (for one arm, total sample size for a replication is twice as large)
  rep_times <- 500 # repetition times, not how many replication studies, which is num_reps
  params_r <- expand.grid(times = rep_times, delta_i = delta_i_values,
                          num_deltai = num_reps,
                          fixed.n = samplesize_reps)
  results_lists <- purrr::pmap(params_r, pass_params)
  IDs <- purrr::pmap(params_r, naming_rep)
  names(results_lists) <- IDs
  hihi_IDs[[i]] <- IDs
  hihi_rep_results[[i]] <- results_lists
  hihi_rep_params[[i]] <- params_r
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
params_o_hihi <- params_o[params_o_slices[[segment_index]], ]
IDs <- purrr::pmap(params_o_hihi,naming_orig)
names(hihi_rep_results) <- IDs
# Reset to sequential plan
plan(sequential)
# Save results
saveRDS(hihi_rep_results, sprintf("./RDS/new/hihi_200500_9parts/hihi_rep_results_part%d.RDS", segment_index))
names(hihi_IDs) <- names(orig_hihi_lists)
saveRDS(hihi_IDs, sprintf("./RDS/new/hihi_200500_9parts/hihi_IDs_part%d.RDS", segment_index))
