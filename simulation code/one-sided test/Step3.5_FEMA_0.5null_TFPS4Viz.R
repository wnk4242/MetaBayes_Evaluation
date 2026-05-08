# This script will create a smaller sample size for Shiny demonstration (only 1999 thresholds).
# This script will create true/false success, true/false failure, true/false positive, true/false negative
# The successes, failures, positives, negatives will be used to calculate rates for visualization.
######################################
#MABF method: FEMA; dataset: 0.5null#
######################################
# Load necessary data and packages
# library(tidyverse)
# library(ggplot2)
# library(scales) # scale_y_continuous(labels = percent) +
# library(purrr)
# library(furrr)
# 
# # Load needed helper functions
# source("./Rfiles/helper functions_XXR.R")
# 
# options(future.globals.maxSize= 2097152000) #To avoid memory limit issues in the future package: 2000MB*1024^2
# 
# # Read in needed RDS files
# FEMA_lists_0.5null_regrouped_deltap <- readRDS("./RDS/MABF4ROC/4TSFS/FEMA_lists_0.5null_regrouped_deltap.RDS")
# 
# # Set up parameters
# num_deltai = 500
# rep_times = 500
# pcutoff_r <- seq(from = 0, by = 0.0005, length.out = 1999) #seq(from = 0.0005, by = 0.0005, length.out = 2000) 
# pcutoff_o <- c(0.01, 0.05, 0.1) 
# EScutoff_o <- c(0)
# 
# #Create an empty dataframe to contain data
# FEMA_TFPS4Viz_0.5null <- data.frame()
# 
# # Create a data frame of all combinations of cutoffs using expand.grid
# cutoff_combinations <- expand.grid(pcutoff_r = pcutoff_r, pcutoff_o = pcutoff_o, EScutoff_o = EScutoff_o)
# # Convert the data frame to a list of lists
# cutoffs <- split(cutoff_combinations, seq(nrow(cutoff_combinations)))
# 
# # Set up parallel processing plan
# plan(multisession, workers = parallel::detectCores())
# 
# # Apply ratesCalc_ROC and thresholds to every sublist of MABF values
# list_of_sublists <- FEMA_lists_0.5null_regrouped_deltap
# 
# # Process all sublists in parallel
# FEMA_TFPS4Viz_0.5null <- future_map_dfr(names(list_of_sublists), function(sublist_name) {
#   message("Processing sublist: ", sublist_name)
#   sublist <- list_of_sublists[[sublist_name]]
#   processMA_sublist_AiO(sublist, sublist_name)
# })
# 
# # End parallel processing plan
# plan(sequential)
# 
# # Reorder columns and break column names into factors
# # Remember to change threshold_p_r to threshold_BF when combining it with other MABF methods outcomes
# FEMA_TFPS4Viz_0.5null <- FEMA_TFPS4Viz_0.5null %>% 
#   select(Sublist, Matrix, threshold_p_r, threshold_p, threshold_ES,
#          Successes,Failures,Successes1,Failures1,Successes2,Failures2,
#          TS1, TS2, FF1, FF2, TF1, TF2, FS1, FS2,
#          TP, TN, FP, FN) %>% 
#   separate(col = "Sublist", into = c("true effect", "orig.n", "QRP level", "PB level"), sep = "_") %>% 
#   separate(col = "Matrix", into = c("rep number","rep.n"), sep = "_") 
# 
# 
# # Add MABF method name column
# FEMA_TFPS4Viz_0.5null <- FEMA_TFPS4Viz_0.5null %>%
#   mutate(method = "FEMA") %>% 
#   relocate(method)
# 
# # Assign group number
# FEMA_TFPS4Viz_0.5null <- FEMA_TFPS4Viz_0.5null %>% 
#   mutate(group = (row_number()-1) %/% length(pcutoff_r) + 1) %>% 
#   select(method, group, everything())
# 
# # Save data as RDS file
# saveRDS(FEMA_TFPS4Viz_0.5null, "./RDS/MABFanalyses/matrix-wise/rates4Viz/FEMA_TFPS4Viz_0.5null.RDS")

#####4/7/2026 Update version: a local faster version + proper progress message
############################################################
# Method: FEMA
# Dataset: 0.5null
############################################################

## =========================
## 1. Load packages
## =========================
library(tidyverse)
library(future)
library(parallelly)

## =========================
## 2. Load helper functions
## =========================
source("./helper functions_XXR.R")

## =========================
## 3. Future / memory settings
## =========================
options(future.globals.maxSize = 2097152000)  # about 2 GB

## =========================
## 4. Read input data
## =========================
FEMA_lists_0.5null_regrouped_deltap <- readRDS(
  "./MABF4ROC/4TSFS/FEMA_lists_0.5null_regrouped_deltap.RDS"
)

## =========================
## 5. Set parameters
## =========================
num_deltai <- 500
rep_times  <- 500

pcutoff_r  <- seq(from = 0, by = 0.0005, length.out = 1999)
pcutoff_o  <- c(0.001, 0.01, 0.05)
EScutoff_o <- 0

## =========================
## 6. Build cutoff combinations
##    Keep as rowwise list because your helper function
##    expects cutoff$pcutoff_r style access
## =========================
cutoff_combinations <- tidyr::expand_grid(
  pcutoff_r  = pcutoff_r,
  pcutoff_o  = pcutoff_o,
  EScutoff_o = EScutoff_o
)

cutoffs <- split(cutoff_combinations, seq_len(nrow(cutoff_combinations)))

## =========================
## 7. Prepare parallel processing
## =========================
workers_to_use <- max(1, availableCores() - 1)
plan(multisession, workers = workers_to_use)

## =========================
## 8. Set up sublists
## =========================
list_of_sublists <- FEMA_lists_0.5null_regrouped_deltap
sublist_names <- names(list_of_sublists)

process_one_sublist <- function(sublist_name) {
  sublist <- list_of_sublists[[sublist_name]]
  processMA_sublist_AiO(sublist, sublist_name)
}

## =========================
## 9. Launch all futures first, with visible submission log
## =========================
cat("Submitting jobs...\n")
flush.console()

start_time <- Sys.time()

futures <- vector("list", length(sublist_names))
names(futures) <- sublist_names

for (i in seq_along(sublist_names)) {
  sublist_name <- sublist_names[i]
  
  cat(sprintf("[%d/%d] Submitting: %s\n", i, length(sublist_names), sublist_name))
  flush.console()
  
  futures[[sublist_name]] <- future::future({
    t0 <- Sys.time()
    result <- process_one_sublist(sublist_name)
    t1 <- Sys.time()
    
    list(
      sublist_name = sublist_name,
      result = result,
      elapsed_sec = as.numeric(difftime(t1, t0, units = "secs"))
    )
  })
}

cat("\nAll jobs submitted.\n\n")
flush.console()

## =========================
## 10. Monitor completion and print permanent log lines
## =========================
results_list <- vector("list", length(futures))
names(results_list) <- names(futures)

completed <- rep(FALSE, length(futures))
names(completed) <- names(futures)

n_total <- length(futures)
n_done <- 0

cat("Waiting for completed datasets:\n")
flush.console()

while (n_done < n_total) {
  for (nm in names(futures)) {
    if (!completed[[nm]] && future::resolved(futures[[nm]])) {
      out <- future::value(futures[[nm]])
      results_list[[nm]] <- out$result
      completed[[nm]] <- TRUE
      n_done <- n_done + 1
      
      current_time <- format(Sys.time(), "%H:%M:%S")
      
      cat(sprintf(
        "[%d/%d] Completed: %s | elapsed: %.2f sec | time: %s\n",
        n_done, n_total, out$sublist_name, out$elapsed_sec, current_time
      ))
      flush.console()
    }
  }
  
  Sys.sleep(0.5)
}
## =========================
## 11. End parallel processing
## =========================
plan(sequential)

total_elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat(sprintf("\nAll datasets completed. Total elapsed time: %.2f sec\n", total_elapsed))
flush.console()

## =========================
## 12. Bind results once at the end
## =========================
FEMA_TFPS4Viz_0.5null <- dplyr::bind_rows(results_list)

## =========================
## 13. Reorder columns and split labels
## =========================
FEMA_TFPS4Viz_0.5null <- FEMA_TFPS4Viz_0.5null %>%
  dplyr::select(
    Sublist, Matrix, threshold_p_r, threshold_p, threshold_ES,
    Successes, Failures, Successes1, Failures1, Successes2, Failures2,
    TS1, TS2, FF1, FF2, TF1, TF2, FS1, FS2,
    TP, TN, FP, FN
  ) %>%
  tidyr::separate(
    col  = "Sublist",
    into = c("true effect", "orig.n", "QRP level", "PB level"),
    sep  = "_"
  ) %>%
  tidyr::separate(
    col  = "Matrix",
    into = c("rep number", "rep.n"),
    sep  = "_"
  )

## =========================
## 14. Add method column
## =========================
FEMA_TFPS4Viz_0.5null <- FEMA_TFPS4Viz_0.5null %>%
  dplyr::mutate(method = "FEMA") %>%
  dplyr::relocate(method)

## =========================
## 15. Assign group number
## =========================
FEMA_TFPS4Viz_0.5null <- FEMA_TFPS4Viz_0.5null %>%
  dplyr::mutate(group = (dplyr::row_number() - 1) %/% length(pcutoff_r) + 1) %>%
  dplyr::select(method, group, dplyr::everything())

## =========================
## 16. Save output
## =========================
saveRDS(
  FEMA_TFPS4Viz_0.5null,
  "./MABFanalyses/matrix-wise/rates4Viz/FEMA_TFPS4Viz_0.5null.RDS"
)

cat("\nDone. Output saved successfully.\n")
flush.console()
