#In this data analysis file, we calculate all kinds of rates for FEMA outcomes
########################
#Step 0: Parameter prep#
########################
# To analyze data, first run Step2.0_FEMA_Data_Combination.R
source("./Step2.0_FEMA_Data_Combination.R")
# These parameters will be used to create FEMA analysis tables for both the 
# delta=null and 0.2 dataset and the delta=null and 0.5 dataset
num_reps <- c(2,5,10) # Must be same as the values in Step 1.1
samplesize_reps <- c(40,100,400) # Must be same as the values in Step 1.1
rep_times <- 500  #  Must be same as the values in Step 1.1

# Define the criteria and pcutoff values
criteria <- list(
    list(EScutoff_o = 0, pcutoff_r = 0.1, EScutoff_r = 0, num_deltai = 500), #c0
    list(EScutoff_o = 0, pcutoff_r = 0.05, EScutoff_r = 0, num_deltai = 500), #c1
    list(EScutoff_o = 0, pcutoff_r = 0.01, EScutoff_r = 0, num_deltai = 500), #c2
    list(EScutoff_o = 0, pcutoff_r = 0.001, EScutoff_r = 0, num_deltai = 500) #c3
)

# This is the variable you want it to be fixed at a value
# pcutoff_o_values is the fixed original p values
pcutoff_o_values <- c(0.1, 0.05, 0.01, 0.001)

# Loop over each pcutoff and each criterion
for (pcutoff_o in pcutoff_o_values) {
    for (criterion in seq_along(criteria)) {
        EScutoff_o <- criteria[[criterion]]$EScutoff_o
        pcutoff_r <- criteria[[criterion]]$pcutoff_r
        EScutoff_r <- criteria[[criterion]]$EScutoff_r
        num_deltai <- criteria[[criterion]]$num_deltai
        
        cat("Processing pcutoff:", pcutoff_o, "Criterion:", criterion, "\n")
        
        ###############################################################################
        #Step 1.1: Calculate diagnostic rates: TPR, FPR, TNR, FNR (delta=null and 0.2)#
        ###############################################################################
        # TPR, FPR, TNR, FNR
        # Apply ratesCalcMA to each sublist, combine the results with rbind, and then combine everything into a dataframe
        FEMA_TFPNrates_0.2null <- do.call(rbind, lapply(FEMA_lists_0.2null_regrouped, function(sublist) {
            do.call(rbind, lapply(sublist, function(x) {
                ratesCalcMA(x, pcutoff_r, EScutoff_r, num_deltai)
            }))
        }))
        
        ###########################################################################
        #Step 1.2: Calculate diagnostic rates: TSR,FSR,TFR,FFR (delta=null and 0.2)#
        ###########################################################################
        # Extracting original p values and converting to a matrix
        orig_p_matrix_0.2null <- matrix(sapply(df_lists_0.2null, function(x) x$p), nrow = 2*num_deltai, byrow = FALSE) 
        obs_d_matrix_0.2null <- matrix(sapply(df_lists_0.2null, function(x) x$d), nrow = 2*num_deltai, byrow = FALSE)
        # Initialize a numerically indexed list for later calculations of rates
        num_outerlists <- nrow(params_FEMAplists_002)
        num_innermatrixs <- length(num_reps)*length(samplesize_reps)
        FEMA_lists_0.2null_regrouped_deltap <- vector("list", num_outerlists)
        # Fill the list of deltas, p values, and MABFs
        for(i in 1:num_outerlists) {
            FEMA_lists_0.2null_regrouped_deltap[[i]] <- list() # Initialize each trial as an empty list
            for(j in 1:num_innermatrixs) {
                # Dynamically create matrix values
                matrix_BFs_deltap <- cbind(rep(c(0,0.2), each = num_deltai), orig_p_matrix_0.2null[, i], obs_d_matrix_0.2null[, i], FEMA_lists_0.2null_regrouped[[i]][[j]])
                # Assign the created matrix to the correct position in FEMA_lists_0.2null_regrouped_deltap
                FEMA_lists_0.2null_regrouped_deltap[[i]][[j]] <- matrix_BFs_deltap
            }
        }
        # Name the lists and matricies 
        FEMA_lists_0.2null_regrouped_deltap <- as.list(FEMA_lists_0.2null_regrouped_deltap)
        # Naming the main sublists
        names(FEMA_lists_0.2null_regrouped_deltap) <- IDs_FEMAplists_002
        # Naming contained matricies
        FEMA_lists_0.2null_regrouped_deltap <- lapply(FEMA_lists_0.2null_regrouped_deltap, function(sublist) {
            names(sublist) <- IDs_FEMApsubls
            return(sublist)
        })
        # Name the columns of the matricies
        for(i in 1:num_outerlists) {
            for(j in 1:num_innermatrixs) {
                colnames(FEMA_lists_0.2null_regrouped_deltap[[i]][[j]]) <- c("delta","original_p","observed_es",1:c(rep_times * 2))
            }
        }
        # Apply ratesCalcMA2 to each sublist, combine the results with rbind, and then combine everything into a dataframe
        FEMA_TFSFrates_0.2null <- do.call(rbind, lapply(FEMA_lists_0.2null_regrouped_deltap, function(sublist) {
            do.call(rbind, lapply(sublist, function(x) {
                ratesCalcMA2(x, pcutoff_o, EScutoff_o, pcutoff_r, EScutoff_r, num_deltai)
            }))
        }))
        # Combine all rates to a table
        FEMA_rates_all_0.2null <- cbind(FEMA_TFPNrates_0.2null,FEMA_TFSFrates_0.2null)
        # Compile reportable table
        params_FEMA_0.2null <- expand.grid(
            delta = "0,0.2",
            rep.number = c(2, 5, 10),
            rep.n = c(40, 100, 400),
            orig.n = c(20, 50, 200),
            qrpEnv = c("none", "medium","high"),
            censorFunc = c("low", "medium","high"),
            stringsAsFactors = TRUE
        )
        # Remove unused conditions
        params_FEMA_0.2null <- params_FEMA_0.2null[-c(28:108,136:216),]
        row.names(params_FEMA_0.2null) <- NULL
        # analysis_matrix
        rates_FEMA_0.2null <- cbind(params_FEMA_0.2null, FEMA_rates_all_0.2null)
        
        ###############################################################################
        #Step 2.1: Calculate diagnostic rates: TPR, FPR, TNR, FNR (delta=null and 0.5)#
        ###############################################################################
        # TPR, FPR, TNR, FNR
        # Apply ratesCalcMA to each sublist, combine the results with rbind, and then combine everything into a dataframe
        FEMA_TFPNrates_0.5null <- do.call(rbind, lapply(FEMA_lists_0.5null_regrouped, function(sublist) {
            do.call(rbind, lapply(sublist, function(x) {
                ratesCalcMA(x, pcutoff_r, EScutoff_r, num_deltai)
            }))
        }))
        
        ############################################################################
        #Step 2.2: Calculate diagnostic rates: TSR,FSR,TFR,FFR (delta=null and 0.5)#
        ############################################################################
        # Extracting original p values and converting to a matrix
        orig_p_matrix_0.5null <- matrix(sapply(df_lists_0.5null, function(x) x$p), nrow = 2*num_deltai, byrow = FALSE)
        obs_d_matrix_0.5null <- matrix(sapply(df_lists_0.5null, function(x) x$d), nrow = 2*num_deltai, byrow = FALSE)
        # Initialize a numerically indexed list for later calculations of rates
        num_outerlists <- nrow(params_FEMAplists_005)
        num_innermatrixs <- length(num_reps)*length(samplesize_reps)
        FEMA_lists_0.5null_regrouped_deltap <- vector("list", num_outerlists)
        # Fill the list of deltas, p values, and MABFs
        for(i in 1:num_outerlists) {
            FEMA_lists_0.5null_regrouped_deltap[[i]] <- list() # Initialize each trial as an empty list
            for(j in 1:num_innermatrixs) {
                # Dynamically create matrix values
                matrix_BFs_deltap <- cbind(rep(c(0,0.5), each = num_deltai), orig_p_matrix_0.5null[, i], obs_d_matrix_0.5null[, i], FEMA_lists_0.5null_regrouped[[i]][[j]])
                # Assign the created matrix to the correct position in FEMA_lists_0.5null_regrouped_deltap
                FEMA_lists_0.5null_regrouped_deltap[[i]][[j]] <- matrix_BFs_deltap
            }
        }
        # Name the lists and matricies 
        FEMA_lists_0.5null_regrouped_deltap <- as.list(FEMA_lists_0.5null_regrouped_deltap)
        # Naming the main sublists
        names(FEMA_lists_0.5null_regrouped_deltap) <- IDs_FEMAplists_005
        # Naming contained matricies
        FEMA_lists_0.5null_regrouped_deltap <- lapply(FEMA_lists_0.5null_regrouped_deltap, function(sublist) {
            names(sublist) <- IDs_FEMApsubls
            return(sublist)
        })
        # Name the columns of the matricies
        for(i in 1:num_outerlists) {
            for(j in 1:num_innermatrixs) {
                colnames(FEMA_lists_0.5null_regrouped_deltap[[i]][[j]]) <- c("delta","original_p", "observed_es", 1:c(rep_times * 2))
            }
        }
        # Apply ratesCalcMA2 to each sublist, combine the results with rbind, and then combine everything into a dataframe
        FEMA_TFSFrates_0.5null <- do.call(rbind, lapply(FEMA_lists_0.5null_regrouped_deltap, function(sublist) {
            do.call(rbind, lapply(sublist, function(x) {
                ratesCalcMA2(x, pcutoff_o, EScutoff_o, pcutoff_r, EScutoff_r, num_deltai)
            }))
        }))
        # Combine all rates to a table
        FEMA_rates_all_0.5null <- cbind(FEMA_TFPNrates_0.5null,FEMA_TFSFrates_0.5null)
        # Compile reportable table
        params_FEMA_0.5null <- expand.grid(
            delta = "0,0.5",
            rep.number = c(2, 5, 10),
            rep.n = c(40, 100, 400),
            orig.n = c(20, 50, 200),
            qrpEnv = c("none", "medium","high"),
            censorFunc = c("low", "medium","high"),
            stringsAsFactors = TRUE
        )
        # Remove unused conditions
        params_FEMA_0.5null <- params_FEMA_0.5null[-c(28:108,136:216),]
        row.names(params_FEMA_0.5null) <- NULL
        # analysis_matrix
        rates_FEMA_0.5null <- cbind(params_FEMA_0.5null, FEMA_rates_all_0.5null)
        
        ################################
        #Step 3: Save analysis outcomes#
        ################################
        # Create directory if it does not exist
        dir_path <- paste0("./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=", pcutoff_o, ", EScutoff_o=0", "/c", criterion - 1)
        if (!dir.exists(dir_path)) {
            dir.create(dir_path, recursive = TRUE)
        }
        
        # Save results
        saveRDS(rates_FEMA_0.2null, paste0(dir_path, "/rates_FEMA_0.2null_c", criterion - 1, ".RDS"))
        saveRDS(rates_FEMA_0.5null, paste0(dir_path, "/rates_FEMA_0.5null_c", criterion - 1, ".RDS"))
    }
}

#####
# # A different way to obtain the same results:
# # In this data analysis file, we calculate all kinds of rates for FEMA outcomes
# # This script delivers the same results as the old Step 3.0 R script, but this one is better thanks to for loops.
# # Step 0: Parameter prep
# rm(list = ls()) # housekeeping
# source("./Step2.0_FEMA_Data_Combination.R")
# num_reps <- c(2, 5, 10) 
# samplesize_reps <- c(40, 100, 400)
# rep_times <- 500  
# 
# directories <- list(
#   "fixed original cutoff/pcutoff_o=0.1, EScutoff_o=0/",
#   "fixed original cutoff/pcutoff_o=0.05, EScutoff_o=0/",
#   "fixed original cutoff/pcutoff_o=0.01, EScutoff_o=0/",
#   "fixed original cutoff/pcutoff_o=0.001, EScutoff_o=0/"
# )
# 
# # This is the variable you want it to be fixed at a value
# pcutoff_o_values <- c(0.1, 0.05, 0.01, 0.001)
# 
# # Define criteria parameters
# criteria_list <- list(
#   list(EScutoff_o = 0, pcutoff_r = 0.1, EScutoff_r = 0, num_deltai = 500, folder = "c0"),
#   list(EScutoff_o = 0, pcutoff_r = 0.05, EScutoff_r = 0, num_deltai = 500, folder = "c1"),
#   list(EScutoff_o = 0, pcutoff_r = 0.01, EScutoff_r = 0, num_deltai = 500, folder = "c2"),
#   list(EScutoff_o = 0, pcutoff_r = 0.001, EScutoff_r = 0, num_deltai = 500, folder = "c3")
# )
# 
# total_directories <- length(directories)
# 
# for (dir_index in 1:total_directories) {
#   directory <- directories[[dir_index]]
#   pcutoff_o_value <- pcutoff_o_values[dir_index]
#   
#   # Update the criteria list with the current pcutoff_o_value
#   updated_criteria_list <- lapply(criteria_list, function(criteria) {
#     criteria$pcutoff_o <- pcutoff_o_value
#     return(criteria)
#   })
#   
#   for (criteria in updated_criteria_list) {
#     pcutoff_o <- criteria$pcutoff_o
#     EScutoff_o <- criteria$EScutoff_o
#     pcutoff_r <- criteria$pcutoff_r
#     EScutoff_r <- criteria$EScutoff_r
#     num_deltai <- criteria$num_deltai
#     folder <- criteria$folder
#     
#     # Step 1.1: Calculate diagnostic rates: TPR, FPR, TNR, FNR (delta=null and 0.2)
#     FEMA_TFPNrates_0.2null <- do.call(rbind, lapply(FEMA_lists_0.2null_regrouped, function(sublist) {
#       do.call(rbind, lapply(sublist, function(x) {
#         ratesCalcMA(x, pcutoff_r, EScutoff_r, num_deltai)
#       }))
#     }))
#     
#     # Step 1.2: Calculate diagnostic rates: SR, FSR, TFR, FFR (delta=null and 0.2)
#     orig_p_matrix_0.2null <- matrix(sapply(df_lists_0.2null, function(x) x$p), nrow = 2 * num_deltai, byrow = FALSE) 
#     obs_d_matrix_0.2null <- matrix(sapply(df_lists_0.2null, function(x) x$d), nrow = 2 * num_deltai, byrow = FALSE)
#     num_outerlists <- nrow(params_FEMAplists_002)
#     num_innermatrixs <- length(num_reps) * length(samplesize_reps)
#     FEMA_lists_0.2null_regrouped_deltap <- vector("list", num_outerlists)
#     
#     for (outerlist_index in 1:num_outerlists) {
#       FEMA_lists_0.2null_regrouped_deltap[[outerlist_index]] <- list()
#       for (innermatrix_index in 1:num_innermatrixs) {
#         matrix_BFs_deltap <- cbind(rep(c(0, 0.2), each = num_deltai), orig_p_matrix_0.2null[, outerlist_index], obs_d_matrix_0.2null[, outerlist_index], FEMA_lists_0.2null_regrouped[[outerlist_index]][[innermatrix_index]])
#         FEMA_lists_0.2null_regrouped_deltap[[outerlist_index]][[innermatrix_index]] <- matrix_BFs_deltap
#       }
#     }
#     
#     FEMA_lists_0.2null_regrouped_deltap <- as.list(FEMA_lists_0.2null_regrouped_deltap)
#     names(FEMA_lists_0.2null_regrouped_deltap) <- IDs_FEMAplists_002
#     FEMA_lists_0.2null_regrouped_deltap <- lapply(FEMA_lists_0.2null_regrouped_deltap, function(sublist) {
#       names(sublist) <- IDs_FEMApsubls
#       return(sublist)
#     })
#     
#     for (outerlist_index in 1:num_outerlists) {
#       for (innermatrix_index in 1:num_innermatrixs) {
#         colnames(FEMA_lists_0.2null_regrouped_deltap[[outerlist_index]][[innermatrix_index]]) <- c("delta", "original_p", "observed_es", 1:c(rep_times * 2))
#       }
#     }
#     
#     FEMA_TFSFrates_0.2null <- do.call(rbind, lapply(FEMA_lists_0.2null_regrouped_deltap, function(sublist) {
#       do.call(rbind, lapply(sublist, function(x) {
#         ratesCalcMA2(x, pcutoff_o, EScutoff_o, pcutoff_r, EScutoff_r, num_deltai)
#       }))
#     }))
#     
#     FEMA_rates_all_0.2null <- cbind(FEMA_TFPNrates_0.2null, FEMA_TFSFrates_0.2null)
#     
#     params_FEMA_0.2null <- expand.grid(
#       delta = "0,0.2",
#       rep.number = c(2, 5, 10),
#       rep.n = c(40, 100, 400),
#       orig.n = c(20, 50, 200),
#       qrpEnv = c("none", "medium", "high"),
#       censorFunc = c("low", "medium", "high"),
#       stringsAsFactors = TRUE
#     )
#     
#     params_FEMA_0.2null <- params_FEMA_0.2null[-c(28:108, 136:216), ]
#     row.names(params_FEMA_0.2null) <- NULL
#     rates_FEMA_0.2null <- cbind(params_FEMA_0.2null, FEMA_rates_all_0.2null)
#     
#     # Step 2.1: Calculate diagnostic rates: TPR, FPR, TNR, FNR (delta=null and 0.5)
#     FEMA_TFPNrates_0.5null <- do.call(rbind, lapply(FEMA_lists_0.5null_regrouped, function(sublist) {
#       do.call(rbind, lapply(sublist, function(x) {
#         ratesCalcMA(x, pcutoff_r, EScutoff_r, num_deltai)
#       }))
#     }))
#     
#     # Step 2.2: Calculate diagnostic rates: TSR, FSR, TFR, FFR (delta=null and 0.5)
#     orig_p_matrix_0.5null <- matrix(sapply(df_lists_0.5null, function(x) x$p), nrow = 2 * num_deltai, byrow = FALSE)
#     obs_d_matrix_0.5null <- matrix(sapply(df_lists_0.5null, function(x) x$d), nrow = 2 * num_deltai, byrow = FALSE)
#     num_outerlists <- nrow(params_FEMAplists_005)
#     num_innermatrixs <- length(num_reps) * length(samplesize_reps)
#     FEMA_lists_0.5null_regrouped_deltap <- vector("list", num_outerlists)
#     
#     for (outerlist_index in 1:num_outerlists) {
#       FEMA_lists_0.5null_regrouped_deltap[[outerlist_index]] <- list()
#       for (innermatrix_index in 1:num_innermatrixs) {
#         matrix_BFs_deltap <- cbind(rep(c(0, 0.5), each = num_deltai), orig_p_matrix_0.5null[, outerlist_index], obs_d_matrix_0.5null[, outerlist_index], FEMA_lists_0.5null_regrouped[[outerlist_index]][[innermatrix_index]])
#         FEMA_lists_0.5null_regrouped_deltap[[outerlist_index]][[innermatrix_index]] <- matrix_BFs_deltap
#       }
#     }
#     
#     FEMA_lists_0.5null_regrouped_deltap <- as.list(FEMA_lists_0.5null_regrouped_deltap)
#     names(FEMA_lists_0.5null_regrouped_deltap) <- IDs_FEMAplists_005
#     FEMA_lists_0.5null_regrouped_deltap <- lapply(FEMA_lists_0.5null_regrouped_deltap, function(sublist) {
#       names(sublist) <- IDs_FEMApsubls
#       return(sublist)
#     })
#     
#     for (outerlist_index in 1:num_outerlists) {
#       for (innermatrix_index in 1:num_innermatrixs) {
#         colnames(FEMA_lists_0.5null_regrouped_deltap[[outerlist_index]][[innermatrix_index]]) <- c("delta", "original_p", "observed_es", 1:c(rep_times * 2))
#       }
#     }
#     
#     FEMA_TFSFrates_0.5null <- do.call(rbind, lapply(FEMA_lists_0.5null_regrouped_deltap, function(sublist) {
#       do.call(rbind, lapply(sublist, function(x) {
#         ratesCalcMA2(x, pcutoff_o, EScutoff_o, pcutoff_r, EScutoff_r, num_deltai)
#       }))
#     }))
#     
#     FEMA_rates_all_0.5null <- cbind(FEMA_TFPNrates_0.5null, FEMA_TFSFrates_0.5null)
#     
#     params_FEMA_0.5null <- expand.grid(
#       delta = "0,0.5",
#       rep.number = c(2, 5, 10),
#       rep.n = c(40, 100, 400),
#       orig.n = c(20, 50, 200),
#       qrpEnv = c("none", "medium", "high"),
#       censorFunc = c("low", "medium", "high"),
#       stringsAsFactors = TRUE
#     )
#     
#     params_FEMA_0.5null <- params_FEMA_0.5null[-c(28:108, 136:216), ]
#     row.names(params_FEMA_0.5null) <- NULL
#     rates_FEMA_0.5null <- cbind(params_FEMA_0.5null, FEMA_rates_all_0.5null)
#     
#     # Step 3: Save analysis outcomes
#     saveRDS(rates_FEMA_0.2null, paste0("./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/", directory, folder, "/rates_FEMA_0.2null_", folder, ".RDS"))
#     saveRDS(rates_FEMA_0.5null, paste0("./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/", directory, folder, "/rates_FEMA_0.5null_", folder, ".RDS"))
#   }
#   # Print progress update
#   cat(paste(dir_index, "of", total_directories, "directories complete.", "\n"))
# }
# 
# 
