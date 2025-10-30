#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Inclusion Bayes Factor of Bayesian Averaged Meta-Analysis (nonelow)
# Chunked parallel version to avoid OOM errors on SLURM
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

######################################################
# Step 0: Load necessary packages and custom functions
######################################################
source("./Rfiles/libraries.R")
source("./Rfiles/QRP and PB functions.R")
source("./Rfiles/helper functions.R")
source("./Rfiles/MABFs.R")
source("./Rfiles/BAMA functions.R")
source("./Rfiles/rows2keep.R")

library(future)
library(future.apply)
library(progressr)

options(future.stdout = TRUE)  # Allow messages from workers
plan(multisession, workers = 6)  # Safe number of workers
options(future.globals.maxSize = 10 * 1024^3)

handlers(global = TRUE)
handlers("txtprogressbar")

#################################################
# Step 1: Load assigned part of replication data
#################################################
folder_path <- "./RDS/REP4EVAL/500500/nonelow_500500_9parts"
segment_index <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

current_rep_results <- readRDS(sprintf("%s/nonelow_rep_results_part%d.RDS", folder_path, segment_index))
current_IDs <- readRDS(sprintf("%s/nonelow_IDs_part%d.RDS", folder_path, segment_index))

########################################################
# Step 2: Feed replication data into MABF for evaluation
########################################################

# prior for true effect size (default, Gronau et al., 2017)
priorEStesting <- prior(
  family = "custom",
  function(x) dcauchy(x, location = 0, scale = 0.707),
  lower = -Inf,
  upper = Inf
)

# Output container
iBF_lists_nonelow <- list()

# Outer loop over simulation batches
for (j in seq_along(current_rep_results)) {
  current_batch <- current_rep_results[[j]]
  current_ids <- current_IDs[[j]]
  
  # Clean up large globals
  rm(current_rep_results, current_IDs)
  gc()
  
  # Set up chunked batching
  num_chunks <- 6  # Match to number of workers
  chunk_size <- ceiling(length(current_batch) / num_chunks)
  chunk_indices <- split(seq_along(current_batch), ceiling(seq_along(current_batch) / chunk_size))
  
  # Parallel over chunks
  with_progress({
    chunked_results <- future_lapply(chunk_indices, function(chunk_k) {
      message(sprintf("Job %d: Starting chunk with %d replications", segment_index, length(chunk_k)))
      chunk_results <- vector("list", length(chunk_k))
      total <- length(chunk_k)
      progress_step <- ceiling(total / 10)  # update every 10%
      
      for (i in seq_along(chunk_k)) {
        idx <- chunk_k[i]
        chunk_results[[i]] <- Reduce(c, iBFcalc(current_batch[[idx]]))
        
        if (i %% progress_step == 0 || i == total) {
          pct <- round(100 * i / total)
          message(sprintf("Job %d: Chunk progress: %d%% (%d of %d)", segment_index, pct, i, total))
        }
      }
      
      chunk_results
    })
  })
  
  # Flatten and name results
  iBF_list <- unlist(chunked_results, recursive = FALSE)
  names(iBF_list) <- current_ids
  iBF_lists_nonelow[[j]] <- iBF_list
  
  # Reload for next batch
  current_rep_results <- readRDS(sprintf("%s/nonelow_rep_results_part%d.RDS", folder_path, segment_index))
  current_IDs <- readRDS(sprintf("%s/nonelow_IDs_part%d.RDS", folder_path, segment_index))
  
  rm(current_batch, current_ids, iBF_list, chunked_results)
  gc()
}

plan(sequential)
names(iBF_lists_nonelow) <- names(current_IDs)

# Save result
saveRDS(iBF_lists_nonelow,
        sprintf("./RDS/MABFoutcomes/500500/iBF/nonelow/iBF_lists_nonelow_part%d.RDS", segment_index))
