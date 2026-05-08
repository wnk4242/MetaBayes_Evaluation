#+++++++++++++++++++++++++++++++++++++++++++#
#Bayes Factor based on Meta-Analysis (hihi) #
#+++++++++++++++++++++++++++++++++++++++++++#
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
folder_path <- "./RDS/REP4EVAL/500500/hihi_500500_9parts"
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
hihi_rep_results <- list(hihi_rep_results_part1,hihi_rep_results_part2,hihi_rep_results_part3,
                            hihi_rep_results_part4,hihi_rep_results_part5,hihi_rep_results_part6,
                            hihi_rep_results_part7,hihi_rep_results_part8,hihi_rep_results_part9)
hihi_IDs <- list(hihi_IDs_part1,hihi_IDs_part2,hihi_IDs_part3,
                    hihi_IDs_part4,hihi_IDs_part5,hihi_IDs_part6,
                    hihi_IDs_part7,hihi_IDs_part8,hihi_IDs_part9)
remove(hihi_rep_results_part1,hihi_rep_results_part2,hihi_rep_results_part3,
       hihi_rep_results_part4,hihi_rep_results_part5,hihi_rep_results_part6,
       hihi_rep_results_part7,hihi_rep_results_part8,hihi_rep_results_part9)
remove(hihi_IDs_part1,hihi_IDs_part2,hihi_IDs_part3,
       hihi_IDs_part4,hihi_IDs_part5,hihi_IDs_part6,
       hihi_IDs_part7,hihi_IDs_part8,hihi_IDs_part9)
# Read the environment variable to determine which segment to process
segment_index <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
current_rep_results <- hihi_rep_results[[segment_index]]
current_IDs <- hihi_IDs[[segment_index]]

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
MA_lists_hihi <- list()
# Nested loops to calculate and store results
for(j in seq_along(current_rep_results)) {
  # Initialize a temporary list to store results for this iteration
  MA_list <- list()
  # Variable to track progress based on k
  count_k <- 0
  for(k in seq_along(current_rep_results[[j]])) {
    MA_templist <- MAcalc(current_rep_results[[j]][[k]])
    MA_list[[k]] <- MA_templist  # Store result in the temporary list
    # Update progress count
    count_k <- count_k + 1
    # Check for progress update
    if (count_k >= progress_threshold) {
      message(sprintf("[1/2 sim: MA]Progress: %d%% completed", percentage))
      progress_threshold <- progress_threshold + increment
      percentage <- percentage + 10
    }
  }
  MA_lists_hihi[[j]] <- MA_list  # Store the temporary list in the main list
  names(MA_lists_hihi[[j]]) <- current_IDs[[j]]
}

# Second calculate Bayes factors based on meta-analytic effect sizes
# Find the total number of iterations for the inner loop across the single sublist works.
total_iterations <- length(MA_lists_hihi[[1]])
# Report progress every 10%
progress_threshold <- 0.1*total_iterations
percentage <- 10
increment <- ceiling(0.1 * total_iterations) # Increment to increase the progress threshold
# Initialize an empty list to store results
BFbMA_lists_hihi <- list()
# Nested loops to calculate and store results
for(j in seq_along(MA_lists_hihi)) {
  # Initialize a temporary list to store results for this iteration
  BFbMA_list <- list()
  # Variable to track progress based on k
  count_k <- 0
  for(k in seq_along(MA_lists_hihi[[j]])) {
    BFbMA_templist <- BFbMAcalc(MA_lists_hihi[[j]][[k]])
    BFbMA_templist <- Reduce(c,BFbMA_templist)
    BFbMA_list[[k]] <- BFbMA_templist  # Store result in the temporary list
    # Update progress count
    count_k <- count_k + 1
    # Check for progress update
    if (count_k >= progress_threshold) {
      message(sprintf("[2/2 sim: BF]Progress: %d%% completed", percentage))
      progress_threshold <- progress_threshold + increment
      percentage <- percentage + 10
    }
  }
  BFbMA_lists_hihi[[j]] <- BFbMA_list  # Store the temporary list in the main list
  names(BFbMA_lists_hihi[[j]]) <- current_IDs[[j]]
}
# Reset to sequential plan
plan(sequential)
names(BFbMA_lists_hihi) <- names(current_IDs)
# Save results
saveRDS(BFbMA_lists_hihi, sprintf("./RDS/MABFoutcomes/500500/BFbMA/hihi/BFbMA_lists_hihi_part%d.RDS", segment_index))
