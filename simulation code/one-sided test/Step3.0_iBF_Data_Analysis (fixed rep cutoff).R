# In this data analysis file, we calculate all kinds of rates for MABF outcomes
# This script delivers the same results as the old Step 3.0 R script, but this one is better thanks to for loops.
# Step 0: Parameter prep
rm(list = ls()) #housekeeping
source("./Step2.0_iBF_Data_Combination.R")
num_reps <- c(2,5,10) 
samplesize_reps <- c(40,100,400)
rep_times <- 500  

directories <- list(
  "fixed rep cutoff/BFcutoff=1, pcutoff_r=0.1, EScutoff_r=0/",
  "fixed rep cutoff/BFcutoff=3, pcutoff_r=0.05, EScutoff_r=0/",
  "fixed rep cutoff/BFcutoff=10, pcutoff_r=0.01, EScutoff_r=0/",
  "fixed rep cutoff/BFcutoff=30, pcutoff_r=0.001, EScutoff_r=0/"
)

BFcutoff_values <- c(1, 3, 10, 30)

# Define criteria parameters
criteria_list <- list(
  list(pcutoff = 0.1, EScutoff = 0, num_deltai = 500, folder = "c0"),
  list(pcutoff = 0.05, EScutoff = 0, num_deltai = 500, folder = "c1"),
  list(pcutoff = 0.01, EScutoff = 0, num_deltai = 500, folder = "c2"),
  list(pcutoff = 0.001, EScutoff = 0, num_deltai = 500, folder = "c3")
)

total_directories <- length(directories)

for (dir_index in 1:total_directories) {
  directory <- directories[[dir_index]]
  BFcutoff_value <- BFcutoff_values[dir_index]
  
  # Update the criteria list with the current BFcutoff_value
  updated_criteria_list <- lapply(criteria_list, function(criteria) {
    criteria$BFcutoff <- BFcutoff_value
    return(criteria)
  })
  
  for (criteria in updated_criteria_list) {
    pcutoff <- criteria$pcutoff
    EScutoff <- criteria$EScutoff
    BFcutoff <- criteria$BFcutoff
    num_deltai <- criteria$num_deltai
    folder <- criteria$folder
    
    # Step 1.1: Calculate diagnostic rates: TPR, FPR, TNR, FNR (delta=null and 0.2)
    iBF_TFPNrates_0.2null <- do.call(rbind, lapply(iBF_lists_0.2null_regrouped, function(sublist) {
      do.call(rbind, lapply(sublist, function(x) {
        ratesCalc(x, BFcutoff, num_deltai)
      }))
    }))
    
    # Step 1.2: Calculate diagnostic rates: SR,FSR,TFR,FFR (delta=null and 0.2)
    orig_p_matrix_0.2null <- matrix(sapply(df_lists_0.2null, function(x) x$p), nrow = 2*num_deltai, byrow = FALSE) 
    obs_d_matrix_0.2null <- matrix(sapply(df_lists_0.2null, function(x) x$d), nrow = 2*num_deltai, byrow = FALSE)
    num_outerlists <- nrow(params_MABFlists_002)
    num_innermatrixs <- length(num_reps)*length(samplesize_reps)
    iBF_lists_0.2null_regrouped_deltap <- vector("list", num_outerlists)
    
    for (outerlist_index in 1:num_outerlists) {
      iBF_lists_0.2null_regrouped_deltap[[outerlist_index]] <- list()
      for (innermatrix_index in 1:num_innermatrixs) {
        matrix_BFs_deltap <- cbind(rep(c(0,0.2), each = num_deltai), orig_p_matrix_0.2null[, outerlist_index], obs_d_matrix_0.2null[, outerlist_index], iBF_lists_0.2null_regrouped[[outerlist_index]][[innermatrix_index]])
        iBF_lists_0.2null_regrouped_deltap[[outerlist_index]][[innermatrix_index]] <- matrix_BFs_deltap
      }
    }
    
    iBF_lists_0.2null_regrouped_deltap <- as.list(iBF_lists_0.2null_regrouped_deltap)
    names(iBF_lists_0.2null_regrouped_deltap) <- IDs_MABFlists_002
    iBF_lists_0.2null_regrouped_deltap <- lapply(iBF_lists_0.2null_regrouped_deltap, function(sublist) {
      names(sublist) <- IDs_MABFsubls
      return(sublist)
    })
    
    for (outerlist_index in 1:num_outerlists) {
      for (innermatrix_index in 1:num_innermatrixs) {
        colnames(iBF_lists_0.2null_regrouped_deltap[[outerlist_index]][[innermatrix_index]]) <- c("delta","original_p","observed_es",1:rep_times)
      }
    }
    
    iBF_TFSFrates_0.2null <- do.call(rbind, lapply(iBF_lists_0.2null_regrouped_deltap, function(sublist) {
      do.call(rbind, lapply(sublist, function(x) {
        ratesCalc2(x, pcutoff, EScutoff, BFcutoff, num_deltai)
      }))
    }))
    
    iBF_rates_all_0.2null <- cbind(iBF_TFPNrates_0.2null, iBF_TFSFrates_0.2null)
    
    params_iBF_0.2null <- expand.grid(
      delta = "0,0.2",
      rep.number = c(2, 5, 10),
      rep.n = c(40, 100, 400),
      orig.n = c(20, 50, 200),
      qrpEnv = c("none", "medium","high"),
      censorFunc = c("low", "medium","high"),
      stringsAsFactors = TRUE
    )
    
    params_iBF_0.2null <- params_iBF_0.2null[-c(28:108,136:216),]
    row.names(params_iBF_0.2null) <- NULL
    rates_iBF_0.2null <- cbind(params_iBF_0.2null, iBF_rates_all_0.2null)
    
    # Step 2.1: Calculate diagnostic rates: TPR, FPR, TNR, FNR (delta=null and 0.5)
    iBF_TFPNrates_0.5null <- do.call(rbind, lapply(iBF_lists_0.5null_regrouped, function(sublist) {
      do.call(rbind, lapply(sublist, function(x) {
        ratesCalc(x, BFcutoff, num_deltai)
      }))
    }))
    
    # Step 2.2: Calculate diagnostic rates: TSR,FSR,TFR,FFR (delta=null and 0.5)
    orig_p_matrix_0.5null <- matrix(sapply(df_lists_0.5null, function(x) x$p), nrow = 2*num_deltai, byrow = FALSE)
    obs_d_matrix_0.5null <- matrix(sapply(df_lists_0.5null, function(x) x$d), nrow = 2*num_deltai, byrow = FALSE)
    num_outerlists <- nrow(params_MABFlists_005)
    num_innermatrixs <- length(num_reps)*length(samplesize_reps)
    iBF_lists_0.5null_regrouped_deltap <- vector("list", num_outerlists)
    
    for (outerlist_index in 1:num_outerlists) {
      iBF_lists_0.5null_regrouped_deltap[[outerlist_index]] <- list()
      for (innermatrix_index in 1:num_innermatrixs) {
        matrix_BFs_deltap <- cbind(rep(c(0,0.5), each = num_deltai), orig_p_matrix_0.5null[, outerlist_index], obs_d_matrix_0.5null[, outerlist_index], iBF_lists_0.5null_regrouped[[outerlist_index]][[innermatrix_index]])
        iBF_lists_0.5null_regrouped_deltap[[outerlist_index]][[innermatrix_index]] <- matrix_BFs_deltap
      }
    }
    
    iBF_lists_0.5null_regrouped_deltap <- as.list(iBF_lists_0.5null_regrouped_deltap)
    names(iBF_lists_0.5null_regrouped_deltap) <- IDs_MABFlists_005
    iBF_lists_0.5null_regrouped_deltap <- lapply(iBF_lists_0.5null_regrouped_deltap, function(sublist) {
      names(sublist) <- IDs_MABFsubls
      return(sublist)
    })
    
    for (outerlist_index in 1:num_outerlists) {
      for (innermatrix_index in 1:num_innermatrixs) {
        colnames(iBF_lists_0.5null_regrouped_deltap[[outerlist_index]][[innermatrix_index]]) <- c("delta","original_p", "observed_es", 1:rep_times)
      }
    }
    
    iBF_TFSFrates_0.5null <- do.call(rbind, lapply(iBF_lists_0.5null_regrouped_deltap, function(sublist) {
      do.call(rbind, lapply(sublist, function(x) {
        ratesCalc2(x, pcutoff, EScutoff, BFcutoff, num_deltai)
      }))
    }))
    
    iBF_rates_all_0.5null <- cbind(iBF_TFPNrates_0.5null, iBF_TFSFrates_0.5null)
    
    params_iBF_0.5null <- expand.grid(
      delta = "0,0.5",
      rep.number = c(2, 5, 10),
      rep.n = c(40, 100, 400),
      orig.n = c(20, 50, 200),
      qrpEnv = c("none", "medium","high"),
      censorFunc = c("low", "medium","high"),
      stringsAsFactors = TRUE
    )
    
    params_iBF_0.5null <- params_iBF_0.5null[-c(28:108,136:216),]
    row.names(params_iBF_0.5null) <- NULL
    rates_iBF_0.5null <- cbind(params_iBF_0.5null, iBF_rates_all_0.5null)
    
    # Step 3: Save analysis outcomes
    saveRDS(rates_iBF_0.2null, paste0("./MABFanalyses/matrix-wise/rates4Plot/", directory, folder, "/rates_iBF_0.2null_", folder, ".RDS"))
    saveRDS(rates_iBF_0.5null, paste0("./MABFanalyses/matrix-wise/rates4Plot/", directory, folder, "/rates_iBF_0.5null_", folder, ".RDS"))
  }
  
  # Print progress update
  cat(paste(dir_index, "of", total_directories, "directories complete.\n"))
}

#saveRDS(iBF_lists_0.2null_regrouped_deltap, file = "./MABF4ROC/4TSFS/iBF_lists_0.2null_regrouped_deltap.RDS")
#saveRDS(iBF_lists_0.5null_regrouped_deltap, file = "./MABF4ROC/4TSFS/iBF_lists_0.5null_regrouped_deltap.RDS")
