##For TPR FPR

library(dplyr)

# read the dataset
shinyROC_0.2null <- readRDS("./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.2null.RDS")

# add MCSE and 95% CI columns
shinyROC_0.2null <- shinyROC_0.2null %>%
  mutate(
    n_tpr = TP + FN,
    n_fpr = FP + TN,
    
    mcse_tpr = sqrt(TPR * (1 - TPR) / n_tpr),
    mcse_fpr = sqrt(FPR * (1 - FPR) / n_fpr),
    
    tpr_ci_lo = pmax(0, TPR - 1.96 * mcse_tpr),
    tpr_ci_hi = pmin(1, TPR + 1.96 * mcse_tpr),
    
    fpr_ci_lo = pmax(0, FPR - 1.96 * mcse_fpr),
    fpr_ci_hi = pmin(1, FPR + 1.96 * mcse_fpr)
  )

# save it back
saveRDS(shinyROC_0.2null, "./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.2null.RDS")


shinyROC_0.5null <- readRDS("./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.5null.RDS")

# add MCSE and 95% CI columns
shinyROC_0.5null <- shinyROC_0.5null %>%
  mutate(
    n_tpr = TP + FN,
    n_fpr = FP + TN,
    
    mcse_tpr = sqrt(TPR * (1 - TPR) / n_tpr),
    mcse_fpr = sqrt(FPR * (1 - FPR) / n_fpr),
    
    tpr_ci_lo = pmax(0, TPR - 1.96 * mcse_tpr),
    tpr_ci_hi = pmin(1, TPR + 1.96 * mcse_tpr),
    
    fpr_ci_lo = pmax(0, FPR - 1.96 * mcse_fpr),
    fpr_ci_hi = pmin(1, FPR + 1.96 * mcse_fpr)
  )





add_mcse_ci <- function(df) {
  df %>%
    mutate(
      n_tpr = TP + FN,
      n_fpr = FP + TN,
      
      mcse_tpr = sqrt(TPR * (1 - TPR) / n_tpr),
      mcse_fpr = sqrt(FPR * (1 - FPR) / n_fpr),
      
      tpr_ci_lo = pmax(0, TPR - 1.96 * mcse_tpr),
      tpr_ci_hi = pmin(1, TPR + 1.96 * mcse_tpr),
      
      fpr_ci_lo = pmax(0, FPR - 1.96 * mcse_fpr),
      fpr_ci_hi = pmin(1, FPR + 1.96 * mcse_fpr)
    )
}

shinyROC_0.2null <- readRDS("./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.2null.RDS")
shinyROC_0.5null <- readRDS("./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.5null.RDS")

shinyROC_0.2null <- add_mcse_ci(shinyROC_0.2null)
shinyROC_0.5null <- add_mcse_ci(shinyROC_0.5null)

saveRDS(shinyROC_0.2null, "./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.2null.RDS")
saveRDS(shinyROC_0.5null, "./MABFanalyses/matrix-wise/rates4ROC/shinyROC_0.5null.RDS")






################
#For TSR FSR
# read the dataset
rates_FEMA_0.2null_c1 <- readRDS("./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.01, EScutoff_o=0/c1/rates_FEMA2_0.2null_c1.RDS")

# add MCSE and 95% CI columns
a <- rates_FEMA_0.2null_c1 %>%
  mutate(
    n_tsr1  = TS1 + FS1,
    n_fsr2 = TS2 + FS2,
    
    mcse_tsr1 = sqrt(TSR1 * (1 - TSR1) / n_tsr1),
    mcse_fsr2 = sqrt(FSR2 * (1 - FSR2) / n_fsr2),
    
    tsr1_ci_lo = pmax(0, TSR1 - 1.96 * mcse_tsr1),
    tsr1_ci_hi = pmin(1, TSR1 + 1.96 * mcse_tsr1),
    
    fsr2_ci_lo = pmax(0, FSR2 - 1.96 * mcse_fsr2),
    fsr2_ci_hi = pmin(1, FSR2 + 1.96 * mcse_fsr2)
  )



library(dplyr)

# folder containing the RDS files
folder_path <- "your/folder/path"

# get all RDS files
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)

# files for each group
files_02 <- rds_files[grepl("0\\.2null", basename(rds_files))]
files_05 <- rds_files[grepl("0\\.5null", basename(rds_files))]

# combine by matching column names
combined_02null <- bind_rows(lapply(files_02, readRDS), .id = "file_id")
combined_05null <- bind_rows(lapply(files_05, readRDS), .id = "file_id")

# optional: add actual file names instead of numeric file_id
combined_02null$file_name <- basename(files_02)[as.integer(combined_02null$file_id)]
combined_05null$file_name <- basename(files_05)[as.integer(combined_05null$file_id)]

# optional: remove file_id if you do not need it
combined_02null$file_id <- NULL
combined_05null$file_id <- NULL