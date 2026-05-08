###For FPR TPR
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


#####

#For TSR FSR

folder_001 <- "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.01, EScutoff_o=0/c1"
folder_005 <- "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.05, EScutoff_o=0/c1"

# function to read, filter, and tag files from one folder
read_rates_files <- function(folder_path, pcutoff_label) {
  rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
  
  # exclude iBF2 and BFbMA
  rds_files <- rds_files[!grepl("iBF2|BFbMA", basename(rds_files))]
  
  df_list <- lapply(rds_files, function(f) {
    dat <- readRDS(f)
    dat$file_name <- basename(f)
    dat$pcutoff_o <- pcutoff_label
    dat
  })
  
  bind_rows(df_list)
}

# read both folders
dat_001 <- read_rates_files(folder_001, "0.01")
dat_005 <- read_rates_files(folder_005, "0.05")

# combine both folders first
all_dat <- bind_rows(dat_001, dat_005)

# then split by 0.2null and 0.5null in file name
combined_02null <- all_dat %>%
  filter(grepl("0\\.2null", file_name))

combined_05null <- all_dat %>%
  filter(grepl("0\\.5null", file_name))



# add MCSE and 95% CI columns
combined_02null <- combined_02null %>%
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

combined_05null <- combined_05null %>%
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
