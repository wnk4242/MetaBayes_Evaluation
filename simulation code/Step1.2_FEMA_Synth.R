# This R script asks the FEMA method to synthesize the generated replication data
# that are not affected by p-hacking or publication bias.
# High-performance computers are required to run this R script.
#+++++++++++++++++++++++++++#
#Fixed-Effect Meta-Analysis #
#+++++++++++++++++++++++++++#
######################################################
#Step 0: Load necessary packages and custom functions#
######################################################
source("libraries.R")
source("QRP and PB functions.R")
source("helper functions.R")
source("MABFs.R")
source("BAMA functions.R")
source("rows2keep.R")
# Use parallel processing
# Determine the number of cores to use
num_cores <- parallel::detectCores()
# plan() is from the future package
plan(multisession, workers = num_cores)

#################################################
#Step 1: Load assigned parts of replication data#
#################################################
# Set the path to the directory containing RDS files
# nonelow stands for no p-hacking or publication bias
folder_path <- "nonelow_500500_9parts"
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
# Initialize empty lists to store results: observed effect sizes and their p-values
MA_beta_lists_nonelow <- list()
MA_p_lists_nonelow <- list()
# Nested loops to calculate and store results
for(j in seq_along(current_rep_results)) {
  # Initialize a temporary lists to store observed effect sizes and their p-values for this iteration
  MA_beta_list <- list()
  MA_p_list <- list()
  # Variable to track progress based on k
  count_k <- 0
  for(k in seq_along(current_rep_results[[j]])) {
    MA_templist <- MAcalc_betap(current_rep_results[[j]][[k]])
    MA_beta_templist <- Reduce(c,lapply(MA_templist, function(sublist) sublist[[1]]))
    MA_p_templist <- Reduce(c,lapply(MA_templist, function(sublist) sublist[[2]]))
    MA_beta_list[[k]] <- MA_beta_templist  # Store result in the temporary list
    MA_p_list[[k]] <- MA_p_templist
    # Update progress count
    count_k <- count_k + 1
    # Check for progress update
    if (count_k >= progress_threshold) {
      message(sprintf("Progress: %d%% completed", percentage))
      progress_threshold <- progress_threshold + increment
      percentage <- percentage + 10
    }
  }
  MA_beta_lists_nonelow[[j]] <- MA_beta_list  # Store the temporary list in the main list
  MA_p_lists_nonelow[[j]] <- MA_p_list
  names(MA_beta_lists_nonelow[[j]]) <- current_IDs[[j]]
  names(MA_p_lists_nonelow[[j]]) <- current_IDs[[j]]
}
# Reset to sequential plan
plan(sequential)
names(MA_beta_lists_nonelow) <- names(current_IDs)
names(MA_p_lists_nonelow) <- names(current_IDs)
# Save results
saveRDS(MA_beta_lists_nonelow, sprintf("FEMA_beta_lists_nonelow_part%d.RDS", segment_index))
saveRDS(MA_p_lists_nonelow, sprintf("FEMA_p_lists_nonelow_part%d.RDS", segment_index))


