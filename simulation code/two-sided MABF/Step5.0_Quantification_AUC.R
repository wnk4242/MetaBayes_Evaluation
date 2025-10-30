# This script obtain quantitative values in the visualizations in Step 4.0 so we can compare the methods formally
# I used two approaches: matrix-wise approach (obtain an average AUC, but no SD statistics); row-wise approach (can obtain 500 individual AUC for each scenario and thier SDs. We will report the second approach results)
# Matrix-wise approach use trapz() to calculate AUC, the row-wise approach use pROC package to calculate AUC. AUC results are highly similar.
# This script converts character variables in rates_0.2null to factors

############# Matrix-wise approach ##############
# Load necessary packages
library(dplyr)
library(tidyr)
library(caTools)#AUC for MABFs
library(purrr)
rm(list = ls()) # housekeeping

# Summarize simulation scenarios showing group number and associated variables
scenario_df <- data.frame(
  scenario = c(1:81),
  true_es = 0.2,
  orig.n = rep(c(20, 50, 200), each=9, times=3),
  rep.number = rep(c(2, 5, 10), times = 27),
  rep.n = rep(c(40,100,400), each = 3, times = 9),
  bias = rep(c('none','medium','high'), each=27)
)
scenario_df

# Set the path to the directory containing RDS files
folder_path <- "./MABFanalyses/matrix-wise/rates4ROC"

# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

# Combine rate outcomes of different methods in order to compare them using AUC
rates_0.2null <- rbind(EUBF_TPFP4ROC_0.2null, FEMABF_TPFP4ROC_0.2null, iBF_TPFP4ROC_0.2null, BFbMA_TPFP4ROC_0.2null, FEMA_TPFP4ROC_0.2null) 

# Convert character to factor
rates_0.2null <- rates_0.2null %>%
  mutate(true_es = 0.2) %>% 
  mutate(across(c(method, true_es, orig.n, `QRP level`, `PB level`, `rep number`,rep.n), as.factor)) %>% 
  relocate(true_es, .after = `true effect`) %>% 
  relocate(group, .before = `method`) %>% 
  rename(scenario = group)



# Calculate AUC of each method across levels of a variable (e.g., orig.n, rep number)
auc_comparison <- rates_0.2null %>%
  filter(`rep number`== 10) %>% 
  filter(rep.n == 40) %>% 
  group_by(method) %>%
  arrange(FPR, .by_group = TRUE) %>%
  summarise(auc = trapz(FPR, TPR), .groups = "drop") 

# View results
print(auc_comparison)

# Obtain AUC for all methods for all scenarios
auc_comparison_all <- map_dfr(1:81, function(scen) {
  rates_0.2null %>%
    filter(scenario == scen) %>%
    group_by(method) %>%
    arrange(FPR, .by_group = TRUE) %>%
    summarise(
      auc = trapz(FPR, TPR),
      .groups = "drop"
    ) %>%
    mutate(scenario = scen)  # add scenario as a column
})

# long format
scenario_with_auc <- scenario_df %>%
  left_join(auc_comparison_all, by = "scenario")

# wide format
auc_wide <- auc_comparison_all %>%
  pivot_wider(
    names_from = method,
    values_from = auc,
    names_prefix = "auc_"
  )

scenario_with_auc <- scenario_df %>%
  left_join(auc_wide, by = "scenario")
print(scenario_with_auc)
#If you look at TPR/FPR AUC, you will find the rep.n do not make any impact on AUC. The original and rep sample size make a lot of impact on AUC
#If you use trapz() to calculate AUC for rep curve (TSR1vsFSR2), it will inflate AUC for bad methods (e.g. BFbMA). The reason the AUC is inflated is because trapz() connects the point before the y axis remain horizontal (which is around (x=0.03, y=0.26) and the point (1,1), then it calculate the AUC

############# Row-wise approach (use this results for final analysis)##############
# install.packages("pROC")
library(pROC)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)

# Set the path to the directory containing RDS files
folder_path <- "./MABF4ROC/4TPFP"

# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}


# ---- AUC per replication -------------------------------------
auc_per_column <- function(M,
                           n_null = 500, n_true = 500,
                           use_log = TRUE, log_base = 10,
                           direction = "auto", quiet = TRUE) {
  M <- as.matrix(M)
  stopifnot(nrow(M) == (n_null + n_true))
  y_true <- c(rep(0L, n_null), rep(1L, n_true))
  
  map_dbl(seq_len(ncol(M)), function(j) {
    scores <- M[, j]
    
    # optional log transform
    if (use_log) scores <- log(scores, base = log_base)
    
    # handle problematic values
    scores[!is.finite(scores)] <- NA
    
    # skip column if all scores invalid
    if (all(is.na(scores))) return(NA_real_)
    
    # safely compute ROC AUC
    tryCatch({
      roc_obj <- roc(response = y_true,
                     predictor = scores,
                     direction = if (direction == "auto") NULL else direction,
                     quiet = quiet, na.rm = TRUE)
      as.numeric(auc(roc_obj))
    }, error = function(e) NA_real_)
  })
}

# ---- Summarize a numeric AUC vector ---------------------------------------
summarize_auc <- function(aucs) {
  aucs <- aucs[is.finite(aucs)]
  n <- length(aucs)
  mean_AUC <- mean(aucs)
  sd_AUC   <- sd(aucs)
  mc_se    <- sd_AUC / sqrt(n)
  tibble(
    mean_AUC = mean_AUC,
    sd_AUC   = sd_AUC,
    mc_se    = mc_se,
    mc_ci_lo = mean_AUC - 1.96 * mc_se,
    mc_ci_hi = mean_AUC + 1.96 * mc_se
  )
}

# ---- run across ALL scenarios & sublists ------------------------------
# list_obj: top-level list (e.g., iBF_lists_0.2null_regrouped)
# return = "summary" (default) or "aucs" (per-replication AUCs)
summarize_auc_all <- function(list_obj,
                              use_log = TRUE, log_base = 10,
                              direction = "auto") {
  
  # Fixed sample sizes
  n_null <- 500
  n_true <- 500
  
  out <- imap(list_obj, function(scen_list, scenario_name) {
    imap_dfr(scen_list, function(M, cell_name) {
      aucs <- auc_per_column(M, n_null, n_true, use_log, log_base, direction)
      tibble(scenario  = scenario_name,
             cell      = cell_name,
             replicate = seq_along(aucs),
             AUC       = aucs)
    })
  }) %>% bind_rows()
  
  # --- collapse across scenarios -------------
  out %>%
    tidyr::separate(cell, into = c("rep.n","rep.number"), sep = "_", convert = TRUE) %>%
    group_by(rep.n, rep.number) %>%
    summarise(summarize_auc(AUC), .groups = "drop") %>%
    arrange(rep.n, rep.number)
}

# Calculate AUC statistics for every scenario
# 0.2null
# FEMA
FEMA_0.2null_AUC <- summarize_auc_all(FEMA_lists_0.2null_regrouped)%>% 
  mutate(method = 'FEMA') %>% 
  relocate(method,.before = 1)

# BFbMA
BFbMA_0.2null_AUC <- summarize_auc_all(BFbMA_lists_0.2null_regrouped)%>% 
  mutate(method = 'BFbMA') %>% 
  relocate(method,.before = 1)

# FEMABF
FEMABF_0.2null_AUC <- summarize_auc_all(FEMABF_lists_0.2null_regrouped)%>% 
  mutate(method = 'FEMABF') %>% 
  relocate(method,.before = 1)

# EUBF
EUBF_0.2null_AUC <- summarize_auc_all(EUBF_lists_0.2null_regrouped)%>% 
  mutate(method = 'EUBF') %>% 
  relocate(method,.before = 1)

# iBF
iBF_0.2null_AUC <- summarize_auc_all(iBF_lists_0.2null_regrouped)%>% 
  mutate(method = 'iBF') %>% 
  relocate(method,.before = 1)

save(BFbMA_0.2null_AUC, EUBF_0.2null_AUC, FEMA_0.2null_AUC, FEMABF_0.2null_AUC, iBF_0.2null_AUC,
     file = "./MABFanalyses/matrix-wise/AUC/MABF_0.2null_AUC.RData")


# 0.5null
# FEMA
FEMA_0.5null_AUC <- summarize_auc_all(FEMA_lists_0.5null_regrouped)%>% 
  mutate(method = 'FEMA') %>% 
  relocate(method,.before = 1)

# BFbMA
BFbMA_0.5null_AUC <- summarize_auc_all(BFbMA_lists_0.5null_regrouped)%>% 
  mutate(method = 'BFbMA') %>% 
  relocate(method,.before = 1)

# FEMABF
FEMABF_0.5null_AUC <- summarize_auc_all(FEMABF_lists_0.5null_regrouped)%>% 
  mutate(method = 'FEMABF') %>% 
  relocate(method,.before = 1)

# EUBF
EUBF_0.5null_AUC <- summarize_auc_all(EUBF_lists_0.5null_regrouped)%>% 
  mutate(method = 'EUBF') %>% 
  relocate(method,.before = 1)

# iBF
iBF_0.5null_AUC <- summarize_auc_all(iBF_lists_0.5null_regrouped)%>% 
  mutate(method = 'iBF') %>% 
  relocate(method,.before = 1)

save(BFbMA_0.5null_AUC, EUBF_0.5null_AUC, FEMA_0.5null_AUC, FEMABF_0.5null_AUC, iBF_0.5null_AUC,
     file = "./MABFanalyses/matrix-wise/AUC/MABF_0.5null_AUC.RData")
