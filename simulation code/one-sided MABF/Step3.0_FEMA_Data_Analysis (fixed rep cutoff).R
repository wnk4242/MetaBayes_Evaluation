# In this data analysis file, we calculate all kinds of rates for FEMA outcomes
# This script delivers the same results as the old Step 3.0 R script, but this one is better thanks to for loops.
# Step 0: Parameter prep
rm(list = ls()) # housekeeping
source("./Step2.0_FEMA_Data_Combination.R")
num_reps <- c(2, 5, 10) 
samplesize_reps <- c(40, 100, 400)
rep_times <- 500  

directories <- list(
  "fixed rep cutoff/BFcutoff=1, pcutoff_r=0.1, EScutoff_r=0/",
  "fixed rep cutoff/BFcutoff=3, pcutoff_r=0.05, EScutoff_r=0/",
  "fixed rep cutoff/BFcutoff=10, pcutoff_r=0.01, EScutoff_r=0/",
  "fixed rep cutoff/BFcutoff=30, pcutoff_r=0.001, EScutoff_r=0/"
)

pcutoff_r_values <- c(0.1, 0.05, 0.01, 0.001)

# Define criteria parameters
criteria_list <- list(
  list(pcutoff_o = 0.1, EScutoff_o = 0, EScutoff_r = 0, num_deltai = 500, folder = "c0"),
  list(pcutoff_o = 0.05, EScutoff_o = 0, EScutoff_r = 0, num_deltai = 500, folder = "c1"),
  list(pcutoff_o = 0.01, EScutoff_o = 0, EScutoff_r = 0, num_deltai = 500, folder = "c2"),
  list(pcutoff_o = 0.001, EScutoff_o = 0, EScutoff_r = 0, num_deltai = 500, folder = "c3")
)

total_directories <- length(directories)

for (dir_index in 1:total_directories) {
  directory <- directories[[dir_index]]
  pcutoff_r_value <- pcutoff_r_values[dir_index]
  
  # Update the criteria list with the current pcutoff_r_value
  updated_criteria_list <- lapply(criteria_list, function(criteria) {
    criteria$pcutoff_r <- pcutoff_r_value
    return(criteria)
  })
  
  for (criteria in updated_criteria_list) {
    pcutoff_o <- criteria$pcutoff_o
    EScutoff_o <- criteria$EScutoff_o
    pcutoff_r <- criteria$pcutoff_r
    EScutoff_r <- criteria$EScutoff_r
    num_deltai <- criteria$num_deltai
    folder <- criteria$folder
    
    # Step 1.1: Calculate diagnostic rates: TPR, FPR, TNR, FNR (delta=null and 0.2)
    FEMA_TFPNrates_0.2null <- do.call(rbind, lapply(FEMA_lists_0.2null_regrouped, function(sublist) {
      do.call(rbind, lapply(sublist, function(x) {
        ratesCalcMA(x, pcutoff_r, EScutoff_r, num_deltai)
      }))
    }))
    
    # Step 1.2: Calculate diagnostic rates: SR, FSR, TFR, FFR (delta=null and 0.2)
    orig_p_matrix_0.2null <- matrix(sapply(df_lists_0.2null, function(x) x$p), nrow = 2 * num_deltai, byrow = FALSE) 
    obs_d_matrix_0.2null <- matrix(sapply(df_lists_0.2null, function(x) x$d), nrow = 2 * num_deltai, byrow = FALSE)
    num_outerlists <- nrow(params_FEMAplists_002)
    num_innermatrixs <- length(num_reps) * length(samplesize_reps)
    FEMA_lists_0.2null_regrouped_deltap <- vector("list", num_outerlists)
    
    for (outerlist_index in 1:num_outerlists) {
      FEMA_lists_0.2null_regrouped_deltap[[outerlist_index]] <- list()
      for (innermatrix_index in 1:num_innermatrixs) {
        matrix_BFs_deltap <- cbind(rep(c(0, 0.2), each = num_deltai), orig_p_matrix_0.2null[, outerlist_index], obs_d_matrix_0.2null[, outerlist_index], FEMA_lists_0.2null_regrouped[[outerlist_index]][[innermatrix_index]])
        FEMA_lists_0.2null_regrouped_deltap[[outerlist_index]][[innermatrix_index]] <- matrix_BFs_deltap
      }
    }
    
    FEMA_lists_0.2null_regrouped_deltap <- as.list(FEMA_lists_0.2null_regrouped_deltap)
    names(FEMA_lists_0.2null_regrouped_deltap) <- IDs_FEMAplists_002
    FEMA_lists_0.2null_regrouped_deltap <- lapply(FEMA_lists_0.2null_regrouped_deltap, function(sublist) {
      names(sublist) <- IDs_FEMApsubls
      return(sublist)
    })
    
    for (outerlist_index in 1:num_outerlists) {
      for (innermatrix_index in 1:num_innermatrixs) {
        colnames(FEMA_lists_0.2null_regrouped_deltap[[outerlist_index]][[innermatrix_index]]) <- c("delta", "original_p", "observed_es", 1:c(rep_times * 2))
      }
    }
    
    FEMA_TFSFrates_0.2null <- do.call(rbind, lapply(FEMA_lists_0.2null_regrouped_deltap, function(sublist) {
      do.call(rbind, lapply(sublist, function(x) {
        ratesCalcMA2(x, pcutoff_o, EScutoff_o, pcutoff_r, EScutoff_r, num_deltai)
      }))
    }))
    
    FEMA_rates_all_0.2null <- cbind(FEMA_TFPNrates_0.2null, FEMA_TFSFrates_0.2null)
    
    params_FEMA_0.2null <- expand.grid(
      delta = "0,0.2",
      rep.number = c(2, 5, 10),
      rep.n = c(40, 100, 400),
      orig.n = c(20, 50, 200),
      qrpEnv = c("none", "medium", "high"),
      censorFunc = c("low", "medium", "high"),
      stringsAsFactors = TRUE
    )
    
    params_FEMA_0.2null <- params_FEMA_0.2null[-c(28:108, 136:216), ]
    row.names(params_FEMA_0.2null) <- NULL
    rates_FEMA_0.2null <- cbind(params_FEMA_0.2null, FEMA_rates_all_0.2null)
    
    # Step 2.1: Calculate diagnostic rates: TPR, FPR, TNR, FNR (delta=null and 0.5)
    FEMA_TFPNrates_0.5null <- do.call(rbind, lapply(FEMA_lists_0.5null_regrouped, function(sublist) {
      do.call(rbind, lapply(sublist, function(x) {
        ratesCalcMA(x, pcutoff_r, EScutoff_r, num_deltai)
      }))
    }))
    
    # Step 2.2: Calculate diagnostic rates: TSR, FSR, TFR, FFR (delta=null and 0.5)
    orig_p_matrix_0.5null <- matrix(sapply(df_lists_0.5null, function(x) x$p), nrow = 2 * num_deltai, byrow = FALSE)
    obs_d_matrix_0.5null <- matrix(sapply(df_lists_0.5null, function(x) x$d), nrow = 2 * num_deltai, byrow = FALSE)
    num_outerlists <- nrow(params_FEMAplists_005)
    num_innermatrixs <- length(num_reps) * length(samplesize_reps)
    FEMA_lists_0.5null_regrouped_deltap <- vector("list", num_outerlists)
    
    for (outerlist_index in 1:num_outerlists) {
      FEMA_lists_0.5null_regrouped_deltap[[outerlist_index]] <- list()
      for (innermatrix_index in 1:num_innermatrixs) {
        matrix_BFs_deltap <- cbind(rep(c(0, 0.5), each = num_deltai), orig_p_matrix_0.5null[, outerlist_index], obs_d_matrix_0.5null[, outerlist_index], FEMA_lists_0.5null_regrouped[[outerlist_index]][[innermatrix_index]])
        FEMA_lists_0.5null_regrouped_deltap[[outerlist_index]][[innermatrix_index]] <- matrix_BFs_deltap
      }
    }
    
    FEMA_lists_0.5null_regrouped_deltap <- as.list(FEMA_lists_0.5null_regrouped_deltap)
    names(FEMA_lists_0.5null_regrouped_deltap) <- IDs_FEMAplists_005
    FEMA_lists_0.5null_regrouped_deltap <- lapply(FEMA_lists_0.5null_regrouped_deltap, function(sublist) {
      names(sublist) <- IDs_FEMApsubls
      return(sublist)
    })
    
    for (outerlist_index in 1:num_outerlists) {
      for (innermatrix_index in 1:num_innermatrixs) {
        colnames(FEMA_lists_0.5null_regrouped_deltap[[outerlist_index]][[innermatrix_index]]) <- c("delta", "original_p", "observed_es", 1:c(rep_times * 2))
      }
    }
    
    FEMA_TFSFrates_0.5null <- do.call(rbind, lapply(FEMA_lists_0.5null_regrouped_deltap, function(sublist) {
      do.call(rbind, lapply(sublist, function(x) {
        ratesCalcMA2(x, pcutoff_o, EScutoff_o, pcutoff_r, EScutoff_r, num_deltai)
      }))
    }))
    
    FEMA_rates_all_0.5null <- cbind(FEMA_TFPNrates_0.5null, FEMA_TFSFrates_0.5null)
    
    params_FEMA_0.5null <- expand.grid(
      delta = "0,0.5",
      rep.number = c(2, 5, 10),
      rep.n = c(40, 100, 400),
      orig.n = c(20, 50, 200),
      qrpEnv = c("none", "medium", "high"),
      censorFunc = c("low", "medium", "high"),
      stringsAsFactors = TRUE
    )
    
    params_FEMA_0.5null <- params_FEMA_0.5null[-c(28:108, 136:216), ]
    row.names(params_FEMA_0.5null) <- NULL
    rates_FEMA_0.5null <- cbind(params_FEMA_0.5null, FEMA_rates_all_0.5null)
    
    # Step 3: Save analysis outcomes
    saveRDS(rates_FEMA_0.2null, paste0("./MABFanalyses/matrix-wise/rates4Plot/", directory, folder, "/rates_FEMA_0.2null_", folder, ".RDS"))
    saveRDS(rates_FEMA_0.5null, paste0("./MABFanalyses/matrix-wise/rates4Plot/", directory, folder, "/rates_FEMA_0.5null_", folder, ".RDS"))
  }
  
  # Print progress update
  cat(paste(dir_index, "of", total_directories, "directories complete.\n"))
}


saveRDS(FEMA_lists_0.2null_regrouped_deltap, paste0("./MABF4ROC/4TSFS/FEMA_lists_0.2null_regrouped_deltap.RDS"))
saveRDS(FEMA_lists_0.5null_regrouped_deltap, paste0("./MABF4ROC/4TSFS/FEMA_lists_0.5null_regrouped_deltap.RDS"))