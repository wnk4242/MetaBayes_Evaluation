# Here, we calculate the classification rates (TPR, FPR, TNR, FNR) for the original study, 
# namely we calculate TFPN rates for the frequentist hypothesis testing (p value)
# We also calculate the mean and median of original study observed ES based on true ES = 0, 0.2, 0.5

################################################
#Step 0.1: Load original study dataset df_lists#
################################################
library(ggplot2)
library(dplyr)
library(tidyr)
# Load necessary helper functions: ratesCalc_TNFP and ratesCalc_TPFN 
source("./helper functions.R")
source("./helper functions_XXR.R")
# Set the path to the directory containing RDS files
folder_path <- "./OGDG/10000" 
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}
df_lists <- df_lists_10000 #for convenience purpose
remove(df_lists_10000)
############################
#Step 2.0: Regroup outcomes#
############################
# df_lists with only null effects
name_with_0 <- startsWith(names(df_lists), c("0_"))
df_lists_0 <- df_lists[name_with_0]
# df_lists with only underlying true effect = 0.2
name_with_0.2 <- startsWith(names(df_lists), c("0.2_"))
df_lists_0.2 <- df_lists[name_with_0.2]
# df_lists with only underlying true effect = 0.5
name_with_0.5 <- startsWith(names(df_lists), c("0.5_"))
df_lists_0.5 <- df_lists[name_with_0.5]
################################################################
#Step 3.0: Calculate TPR, FPR, TNR, FNR when true effect is 0.2#
################################################################
# Apply ratesCalc functions to calculate TPR, FPR, TNR, FNR
# Calculate TNR and FPR
rates_orig_null <- lapply(df_lists_0, ratesCalc_TNFP)
# Combine the results into a single dataframe
rates_orig_null <- do.call(rbind, rates_orig_null)
# Calculate TPR and FNR
rates_orig_0.2 <- lapply(df_lists_0.2, ratesCalc_TPFN)
# Combine the results into a single dataframe
rates_orig_0.2 <- do.call(rbind, rates_orig_0.2)
# Combine the dataframes of rates 
rates_orig_0.2null <- cbind(rates_orig_null, rates_orig_0.2)
# Compile reportable table
params_orig_0.2null <- expand.grid(
  method = "Significance",
  delta = "0,0.2",
  orig.n = c(20, 50, 200),
  qrpEnv = c("none", "medium","high"),
  censorFunc = c("low", "medium","high"),
  stringsAsFactors = TRUE
)
# Remove unused conditions
params_orig_0.2null <- params_orig_0.2null[-c(4:12,16:24),]
row.names(params_orig_0.2null) <- NULL
# Combine parameters and rates
rates_orig_0.2null <- cbind(params_orig_0.2null, rates_orig_0.2null)
row.names(rates_orig_0.2null) <- NULL
# Convert orig.n (original sample size) from numeric to categorical
rates_orig_0.2null <- rates_orig_0.2null %>% 
  mutate(orig.n.level = factor(orig.n,
                               levels = c(20, 50, 200),
                               labels = c("small", "medium", "large")))
# Relocate the new column directly after
rates_orig_0.2null <- rates_orig_0.2null %>% 
  select(-orig.n.level, everything()) %>% 
  relocate(orig.n.level, .after = orig.n)

################################################################
#Step 3.0: Calculate TPR, FPR, TNR, FNR when true effect is 0.5#
################################################################
# Apply ratesCalc functions to calculate TPR, FPR, TNR, FNR
# Calculate TNR and FPR
rates_orig_null <- lapply(df_lists_0, ratesCalc_TNFP)
# Combine the results into a single dataframe
rates_orig_null <- do.call(rbind, rates_orig_null)
# Calculate TPR and FNR
rates_orig_0.5 <- lapply(df_lists_0.5, ratesCalc_TPFN)
# Combine the results into a single dataframe
rates_orig_0.5 <- do.call(rbind, rates_orig_0.5)
# Combine the dataframes of rates 
rates_orig_0.5null <- cbind(rates_orig_null, rates_orig_0.5)
# Compile reportable table
params_orig_0.5null <- expand.grid(
  method = "Significance",
  delta = "0,0.5",
  orig.n = c(20, 50, 200),
  qrpEnv = c("none", "medium","high"),
  censorFunc = c("low", "medium","high"),
  stringsAsFactors = TRUE
)
# Remove unused conditions
params_orig_0.5null <- params_orig_0.5null[-c(4:12,16:24),]
row.names(params_orig_0.5null) <- NULL
# Combine parameters and rates
rates_orig_0.5null <- cbind(params_orig_0.5null, rates_orig_0.5null)
row.names(rates_orig_0.5null) <- NULL
# Convert orig.n (original sample size) from numeric to categorical
rates_orig_0.5null <- rates_orig_0.5null %>% 
  mutate(orig.n.level = factor(orig.n,
                               levels = c(20, 50, 200),
                               labels = c("small", "medium", "large")))
# Relocate the new column directly after
rates_orig_0.5null <- rates_orig_0.5null %>% 
  select(-orig.n.level, everything()) %>% 
  relocate(orig.n.level, .after = orig.n)


####################################################################################################
#In the following script, we extract and calculate mean and median effect sizes of original studies#
####################################################################################################
################################################
#Step 0.1: Load original study dataset df_lists#
################################################
library(ggplot2)
library(dplyr)
library(tidyr)
# Set the path to the directory containing RDS files
folder_path <- "./OGDG/10000" 
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}
df_lists <- df_lists_10000 #for convenience purpose
remove(df_lists_10000)
############################
#Step 2.0: Regroup outcomes#
############################
# df_lists with only null effects
name_with_0 <- startsWith(names(df_lists), c("0_"))
df_lists_0 <- df_lists[name_with_0]
# df_lists with only underlying true effect = 0.2
name_with_0.2 <- startsWith(names(df_lists), c("0.2_"))
df_lists_0.2 <- df_lists[name_with_0.2]
# df_lists with only underlying true effect = 0.5
name_with_0.5 <- startsWith(names(df_lists), c("0.5_"))
df_lists_0.5 <- df_lists[name_with_0.5]
#########################################################################
#Step 3.0: Calculate mean and median ES when true effect is 0.2 and null#
#########################################################################
# Apply ratesCalc functions to calculate TPR, FPR, TNR, FNR
# Calculate TNR and FPR
es_orig_null <- lapply(df_lists_0, esCalc_null)
# Combine the results into a single dataframe
es_orig_null <- do.call(rbind, es_orig_null)
# Calculate TPR and FNR
es_orig_0.2 <- lapply(df_lists_0.2, esCalc_true)
# Combine the results into a single dataframe
es_orig_0.2 <- do.call(rbind, es_orig_0.2)
# Combine the dataframes of rates 
es_orig_0.2null <- cbind(es_orig_null, es_orig_0.2)
# Compile reportable table
params_orig_0.2null <- expand.grid(
  method = "Significance",
  delta = "0,0.2",
  orig.n = c(20, 50, 200),
  qrpEnv = c("none", "medium","high"),
  censorFunc = c("low", "medium","high"),
  stringsAsFactors = TRUE
)
# Remove unused conditions
params_orig_0.2null <- params_orig_0.2null[-c(4:12,16:24),]
row.names(params_orig_0.2null) <- NULL
# Combine parameters and rates
es_orig_0.2null <- cbind(params_orig_0.2null, es_orig_0.2null)
row.names(es_orig_0.2null) <- NULL


#########################################################################
#Step 3.0: Calculate mean and median ES when true effect is 0.5 and null#
#########################################################################
# Apply ratesCalc functions to calculate TPR, FPR, TNR, FNR
# Calculate TNR and FPR
es_orig_null <- lapply(df_lists_0, esCalc_null)
# Combine the results into a single dataframe
es_orig_null <- do.call(rbind, es_orig_null)
# Calculate TPR and FNR
es_orig_0.5 <- lapply(df_lists_0.5, esCalc_true)
# Combine the results into a single dataframe
es_orig_0.5 <- do.call(rbind, es_orig_0.5)
# Combine the dataframes of rates 
es_orig_0.5null <- cbind(es_orig_null, es_orig_0.5)
# Compile reportable table
params_orig_0.5null <- expand.grid(
  method = "Significance",
  delta = "0,0.5",
  orig.n = c(20, 50, 200),
  qrpEnv = c("none", "medium","high"),
  censorFunc = c("low", "medium","high"),
  stringsAsFactors = TRUE
)
# Remove unused conditions
params_orig_0.5null <- params_orig_0.5null[-c(4:12,16:24),]
row.names(params_orig_0.5null) <- NULL
# Combine parameters and rates
es_orig_0.5null <- cbind(params_orig_0.5null, es_orig_0.5null)
row.names(es_orig_0.5null) <- NULL

