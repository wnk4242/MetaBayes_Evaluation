#++++++++++++++++++++++++++++++++++++++++#
#Evidence Updating Bayes Factor (medmed) #
#++++++++++++++++++++++++++++++++++++++++#
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
folder_path <- "./RDS/REP4EVAL/500500/medmed_500500_9parts" 
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
medmed_rep_results <- list(medmed_rep_results_part1,medmed_rep_results_part2,medmed_rep_results_part3,
                            medmed_rep_results_part4,medmed_rep_results_part5,medmed_rep_results_part6,
                            medmed_rep_results_part7,medmed_rep_results_part8,medmed_rep_results_part9)
medmed_IDs <- list(medmed_IDs_part1,medmed_IDs_part2,medmed_IDs_part3,
                    medmed_IDs_part4,medmed_IDs_part5,medmed_IDs_part6,
                    medmed_IDs_part7,medmed_IDs_part8,medmed_IDs_part9)
remove(medmed_rep_results_part1,medmed_rep_results_part2,medmed_rep_results_part3,
       medmed_rep_results_part4,medmed_rep_results_part5,medmed_rep_results_part6,
       medmed_rep_results_part7,medmed_rep_results_part8,medmed_rep_results_part9)
remove(medmed_IDs_part1,medmed_IDs_part2,medmed_IDs_part3,
       medmed_IDs_part4,medmed_IDs_part5,medmed_IDs_part6,
       medmed_IDs_part7,medmed_IDs_part8,medmed_IDs_part9)
# Read the environment variable to determine which segment to process
segment_index <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
current_rep_results <- medmed_rep_results[[segment_index]]
current_IDs <- medmed_IDs[[segment_index]]

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
EUBF_lists_medmed <- list()
# Nested loops to calculate and store results
for(j in seq_along(current_rep_results)) {
  # Initialize a temporary list to store results for this iteration
  EUBF_list <- list()
  # Variable to track progress based on k
  count_k <- 0
  for(k in seq_along(current_rep_results[[j]])) {
    EUBF_templist <- EUBFcalc(current_rep_results[[j]][[k]])
    EUBF_templist <- Reduce(c,EUBF_templist)
    EUBF_list[[k]] <- EUBF_templist  # Store result in the temporary list
    # Update progress count
    count_k <- count_k + 1
    # Check for progress update
    if (count_k >= progress_threshold) {
      message(sprintf("Progress: %d%% completed", percentage))
      progress_threshold <- progress_threshold + increment
      percentage <- percentage + 10
    }
  }
  EUBF_lists_medmed[[j]] <- EUBF_list  # Store the temporary list in the main list
  names(EUBF_lists_medmed[[j]]) <- current_IDs[[j]]
}
# Reset to sequential plan
plan(sequential)
names(EUBF_lists_medmed) <- names(current_IDs)
# Save results
saveRDS(EUBF_lists_medmed, sprintf("./RDS/MABFoutcomes/500500/EUBF/medmed/EUBF_lists_medmed_part%d.RDS", segment_index))
