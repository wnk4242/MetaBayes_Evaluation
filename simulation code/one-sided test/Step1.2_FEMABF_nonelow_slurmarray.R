# This code use inner loop to show progress message.
# This code integrates the Slurm array command to divide the workload to nine individual arrays.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#Fixed-Effect Meta-Analysis Bayes Factor-evaluate replication data (QRP: none, PB: low)#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
######################################################
#Step 0: Load necessary packages and custom functions#
######################################################
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

#################################################
#Step 1: Load assigned parts of replication data#
#################################################
# Set the path to the directory containing RDS files
folder_path <- "./RDS/REP4EVAL/500500/nonelow_500500_9parts"
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

##############################################
#Step 2: Combine replication data into a list#
##############################################
nonelow_rep_results <- list(nonelow_rep_results_part1,nonelow_rep_results_part2,nonelow_rep_results_part3,
                            nonelow_rep_results_part4,nonelow_rep_results_part5,nonelow_rep_results_part6,
                            nonelow_rep_results_part7,nonelow_rep_results_part8,nonelow_rep_results_part9)
nonelow_IDs <- list(nonelow_IDs_part1,nonelow_IDs_part2,nonelow_IDs_part3,
                    nonelow_IDs_part4,nonelow_IDs_part5,nonelow_IDs_part6,
                    nonelow_IDs_part7,nonelow_IDs_part8,nonelow_IDs_part9)
remove(nonelow_rep_results_part1,nonelow_rep_results_part2,nonelow_rep_results_part3,
       nonelow_rep_results_part4,nonelow_rep_results_part5,nonelow_rep_results_part6,
       nonelow_rep_results_part7,nonelow_rep_results_part8,nonelow_rep_results_part9)
remove(nonelow_IDs_part1,nonelow_IDs_part2,nonelow_IDs_part3,
       nonelow_IDs_part4,nonelow_IDs_part5,nonelow_IDs_part6,
       nonelow_IDs_part7,nonelow_IDs_part8,nonelow_IDs_part9)
# Read the environment variable to determine which segment to process
segment_index <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
current_rep_results <- nonelow_rep_results[[segment_index]]
current_IDs <- nonelow_IDs[[segment_index]]

########################################################
#Step 3: Feed replication data into MABF for evaluation#
########################################################
# To show progress message:
# Because current_rep_results only contains one list, so using outer loop for progress message doesnt work.
# Find the total number of iterations for the inner loop across the single sublist works.
total_iterations <- length(current_rep_results[[1]])
# Report progress every 10%
progress_threshold <- 0.1*total_iterations
percentage <- 10
increment <- ceiling(0.1 * total_iterations) # Increment to increase the progress threshold
# Initialize an empty list to store results
FEMABF_lists_nonelow <- list()
# Outer loop, although only one iteration now
for(j in seq_along(current_rep_results)) {
  # Initialize a temporary list to store results for this iteration
  FEMABF_list <- list()
  # Variable to track progress based on k
  count_k <- 0
  for(k in seq_along(current_rep_results[[j]])) {
    FEMABF_templist <- suppressMessages(FEMABFcalc(current_rep_results[[j]][[k]]))
    FEMABF_list[[k]] <- FEMABF_templist  # Store result in the temporary list
    # Update progress count
    count_k <- count_k + 1
    # Check for progress update
    if (count_k >= progress_threshold) {
      message(sprintf("Progress: %d%% completed", percentage))
      progress_threshold <- progress_threshold + increment
      percentage <- percentage + 10
    }
  }
  FEMABF_lists_nonelow[[j]] <- FEMABF_list  # Store the temporary list in the main list
  names(FEMABF_lists_nonelow[[j]]) <- current_IDs[[j]]
}
# Reset to sequential plan
plan(sequential)
names(FEMABF_lists_nonelow) <- names(current_IDs)
# Save results
saveRDS(FEMABF_lists_nonelow, sprintf("./RDS/MABFoutcomes/500500/FEMABF/nonelow/FEMABF_lists_nonelow_part%d.RDS", segment_index))
