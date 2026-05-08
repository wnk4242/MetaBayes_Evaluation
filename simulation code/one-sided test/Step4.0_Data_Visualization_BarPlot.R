# This is matrix-wise rates, not row-wise rates
################################################
#Step 0: Load necessary libraries and functions#
################################################
library(tidyverse)
library(scales) #for creating plots
setwd("D:/wnk/PhD/PhD research/dissertation/code/hprc")
source("./helper functions.R") #barplot_XXRbyMABF_overall and barplot_XXRbyMABF_cond
#############################################################################
#Step 1: Load R data files containing classification accuracy rates of MABFs#
#############################################################################
# Set the path to the directory containing RDS files
folder_path <- "./MABFanalyses/matrix-wise/rates4Plot/fixed rep cutoff/BFcutoff=3, pcutoff_r=0.05, EScutoff_r=0/c3" 
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}
#######################################
#Step 2: Combine datasets for analyses#
#######################################
# combine data sets of all four methods
# 0.2null
rates_0.2null <- rbind(rates_FEMABF_0.2null_c3[,-c(11,12,13,18,19,20)],rates_BFbMA_0.2null_c3[,-c(11,12,13,18,19,20)],rates_EUBF_0.2null_c3[,-c(11,12,13,18,19,20)],rates_iBF_0.2null_c3[,-c(11,12,13,18,19,20)],rates_FEMA_0.2null_c3)
rates_0.2null$method <- rep(c("FEMABF","BFbMA","EUBF","iBF","FEMA"), each = nrow(rates_iBF_0.2null_c3))
rates_0.2null <- rates_0.2null %>%
  select(method, everything())
row.names(rates_0.2null) <- NULL
# Convert rep.n (replication sample size) from numeric to categorical
rates_0.2null <- rates_0.2null %>% 
  mutate(rep.n.level = factor(rep.n,
                            levels = c(40, 100, 400),
                            labels = c("small", "medium", "large"))) %>% 
  mutate(orig.n.level = factor(orig.n,
                             levels = c(20, 50, 200),
                             labels = c("small", "medium", "large"))) %>% 
  mutate(rep.number.level = factor(rep.number,
                                 levels = c(2, 5, 10),
                                 labels = c("small", "medium", "large")))
# Relocate the new column directly after rep.n
rates_0.2null <- rates_0.2null %>% 
  select(-rep.n.level, everything()) %>% 
  relocate(rep.n.level, .after = rep.n) %>% 
  select(-orig.n.level, everything()) %>% 
  relocate(orig.n.level, .after = orig.n) %>% 
  select(-rep.number.level, everything()) %>% 
  relocate(rep.number.level, .after = rep.number) 
# 0.5null
rates_0.5null <- rbind(rates_FEMABF_0.5null_c3[,-c(11,12,13,18,19,20)],rates_BFbMA_0.5null_c3[,-c(11,12,13,18,19,20)],rates_EUBF_0.5null_c3[,-c(11,12,13,18,19,20)],rates_iBF_0.5null_c3[,-c(11,12,13,18,19,20)],rates_FEMA_0.5null_c3)
rates_0.5null$method <- rep(c("FEMABF","BFbMA","EUBF","iBF","FEMA"), each = nrow(rates_iBF_0.5null_c3))
rates_0.5null <- rates_0.5null %>%
  select(method, everything())
row.names(rates_0.5null) <- NULL
# Convert rep.n (replication sample size) from numeric to categorical
rates_0.5null <- rates_0.5null %>% 
  mutate(rep.n.level = factor(rep.n,
                              levels = c(40, 100, 400),
                              labels = c("small", "medium", "large"))) %>% 
  mutate(orig.n.level = factor(orig.n,
                               levels = c(20, 50, 200),
                               labels = c("small", "medium", "large"))) %>% 
  mutate(rep.number.level = factor(rep.number,
                                   levels = c(2, 5, 10),
                                   labels = c("small", "medium", "large")))
# Relocate the new column directly after rep.n
rates_0.5null <- rates_0.5null %>% 
  select(-rep.n.level, everything()) %>% 
  relocate(rep.n.level, .after = rep.n) %>% 
  select(-orig.n.level, everything()) %>% 
  relocate(orig.n.level, .after = orig.n) %>% 
  select(-rep.number.level, everything()) %>% 
  relocate(rep.number.level, .after = rep.number) 

#########################################################
#Step 3: Create bar plots for all types of overall rates#
#########################################################
# For overall rates, we need to set three parameters:
# TrueES (0.2 or 0.5), 
# data2plot (rates_0.5null), and 
# ratetype2calc (XXR)

##########################################################################################################
#Overall rate (TPR and FPR) of the significance method used by the original study when true effect is 0.2#
##########################################################################################################
# TPR of the significance method at different QRP levels
barplot_XXRorig_overall(trueES = 0.2, ratetype = "TPR", param = "method", dataset = rates_0.2null)

# FPR of the significance method at different QRP levels
barplot_XXRorig_overall(trueES = 0.2, ratetype = "FPR", param = "method", dataset = rates_0.2null)
##########################################################################################################
#Overall rate (TPR and FPR) of the significance method used by the original study when true effect is 0.5#
##########################################################################################################
# TPR of the significance method at different QRP levels
barplot_XXRorig_overall(trueES = 0.5, ratetype = "TPR", param = "method", dataset = rates_0.5null)

# FPR of the significance method at different QRP levels
# Which is the same as when true effect size is 0.2 because we use the same null data set
barplot_XXRorig_overall(trueES = 0.5, ratetype = "FPR", param = "method", dataset = rates_0.5null)

#####################################################
#Overall rate by MABF method when true effect is 0.2#
#####################################################
# Set up parameters
TrueES <- 0.2
data2plot <- rates_0.2null
ratetype2calc <- c("ADR", "ADRTE", "ADRNE", "SAR", "FAR", "FCR", "SCR", "SAR1", "SAR2", "FAR1", "FAR2", "FCR1", "FCR2", "SCR1", "SCR2", "TPR", "FPR", "TSR", "FSR", "TFR", "FFR", "TSR1", "TSR2", "FSR1", "FSR2", "TFR1", "TFR2", "FFR1", "FFR2")
# Use the lambda syntax in purrr::walk to save plots
# .x replaces ratetype2calc
walk(ratetype2calc, ~ {
  plot <- barplot_XXRbyMABF_overall(TrueES, .x, data2plot)
  filename <- paste0("plots/bar plot/by MABF methods/overall rate/", "Overall_", .x, "_", TrueES, ".jpg")
  ggsave(filename = filename, plot = plot, width = 8, height = 6)
})

#####################################################
#Overall rate by MABF method when true effect is 0.5#
#####################################################
# Set up parameters
TrueES <- 0.5
data2plot <- rates_0.5null
ratetype2calc <- c("ADR", "ADRTE", "ADRNE", "SAR", "FAR", "FCR", "SCR", "SAR1", "SAR2", "FAR1", "FAR2", "FCR1", "FCR2", "SCR1", "SCR2", "TPR", "FPR", "TSR", "FSR", "TFR", "FFR", "TSR1", "TSR2", "FSR1", "FSR2", "TFR1", "TFR2", "FFR1", "FFR2")
# Use the lambda syntax in purrr::walk to save plots
walk(ratetype2calc, ~ {
  plot <- barplot_XXRbyMABF_overall(TrueES, .x, data2plot)
  filename <- paste0("plots/bar plot/by MABF methods/overall rate/", "Overall_", .x, "_", TrueES, ".jpg")
  ggsave(filename = filename, plot = plot, width = 8, height = 6)
})

########################################################################################
#Step 4: Create bar plots for all types of conditional rates for parameter combinations#
########################################################################################
# For conditional rates, we need to set three parameters:
# TrueES (0.2 or 0.5), 
# data2plot (rates_0.5null), 
# ratetype2calc (XXR), and
# dir2save (0.2null or 0.5null)
#####################################################################################
#Rates conditioned on one parameter when true effect is 0.2: replication sample size#
#####################################################################################
# Set up parameters
paramcomb2plot <- expand.grid(repsampsize = c("small", "medium", "large")) #paramters to plot
TrueES <- 0.2
data2plot <- rates_0.2null
ratetype2calc <- c("ADR", "ADRTE", "ADRNE", "SAR", "FAR", "FCR", "SCR", "SAR1", "SAR2", "FAR1", "FAR2", "FCR1", "FCR2", "SCR1", "SCR2", "TPR", "FPR", "TSR", "FSR", "TFR", "FFR", "TSR1", "TSR2", "FSR1", "FSR2", "TFR1", "TFR2", "FFR1", "FFR2")
dir2save <- c("0.2null")
# Generate directory name based on parameters
filename_parts <- lapply(names(paramcomb2plot), function(name) {
  if (!is.null(paramcomb2plot[[name]])) {
    paste0(name)
  }
})
# Combine factor names generated by filename_parts into a single directory name
dir.name <- paste0("D:/wnk/PhD/PhD research/dissertation/code/hprc/plots/bar plot/by MABF methods/conditional rate/", length(paramcomb2plot), " factor/", paste0(filename_parts, collapse = "_"))
# Create the directory
dir.create(dir.name, recursive = TRUE)  # Set recursive = TRUE to create intermediate directories if they don't exist
# Set the working directory
setwd(dir.name)
# Generate and save plots for each combination
# Iterate over each rate type
for (ratetype2calc in ratetype2calc) {
  pwalk(paramcomb2plot, savebarplot_XXRbyMABF_cond)
}
# Change back to original working directory
setwd("D:/wnk/PhD/PhD research/dissertation/code/hprc")
################################################################################
#Rates conditioned on one parameter when true effect is 0.2: replication number#
################################################################################
# Set up parameters
paramcomb2plot <- expand.grid(repnumber = c("small", "medium", "large"))
TrueES <- 0.2
data2plot <- rates_0.2null
ratetype2calc <- c("ADR", "ADRTE", "ADRNE", "SAR", "FAR", "FCR", "SCR", "SAR1", "SAR2", "FAR1", "FAR2", "FCR1", "FCR2", "SCR1", "SCR2", "TPR", "FPR", "TSR", "FSR", "TFR", "FFR", "TSR1", "TSR2", "FSR1", "FSR2", "TFR1", "TFR2", "FFR1", "FFR2")
dir2save <- c("0.2null")
# Generate directory name based on parameters
filename_parts <- lapply(names(paramcomb2plot), function(name) {
  if (!is.null(paramcomb2plot[[name]])) {
    paste0(name)
  }
})
# Combine factor names generated by filename_parts into a single directory name
dir.name <- paste0("D:/wnk/PhD/PhD research/dissertation/code/hprc/plots/bar plot/by MABF methods/conditional rate/", length(paramcomb2plot), " factor/", paste0(filename_parts, collapse = "_"))
# Create the directory
dir.create(dir.name, recursive = TRUE)  # Set recursive = TRUE to create intermediate directories if they don't exist
# Set the working directory
setwd(dir.name)
# Generate and save plots for each combination
# Iterate over each rate type
for (ratetype2calc in ratetype2calc) {
  pwalk(paramcomb2plot, savebarplot_XXRbyMABF_cond)
}
# Change back to original working directory
setwd("D:/wnk/PhD/PhD research/dissertation/code/hprc")
##################################################################################
#Rates conditioned on one parameter when true effect is 0.2: original sample size#
##################################################################################
# Set up parameters
paramcomb2plot <- expand.grid(origsampsize = c("small", "medium", "large"))
TrueES <- 0.2
data2plot <- rates_0.2null
ratetype2calc <- c("ADR", "ADRTE", "ADRNE", "SAR", "FAR", "FCR", "SCR", "SAR1", "SAR2", "FAR1", "FAR2", "FCR1", "FCR2", "SCR1", "SCR2", "TPR", "FPR", "TSR", "FSR", "TFR", "FFR", "TSR1", "TSR2", "FSR1", "FSR2", "TFR1", "TFR2", "FFR1", "FFR2")
dir2save <- c("0.2null")
# Generate directory name based on parameters
filename_parts <- lapply(names(paramcomb2plot), function(name) {
  if (!is.null(paramcomb2plot[[name]])) {
    paste0(name)
  }
})
# Combine factor names generated by filename_parts into a single directory name
dir.name <- paste0("D:/wnk/PhD/PhD research/dissertation/code/hprc/plots/bar plot/by MABF methods/conditional rate/", length(paramcomb2plot), " factor/", paste0(filename_parts, collapse = "_"))
# Create the directory
dir.create(dir.name, recursive = TRUE)  # Set recursive = TRUE to create intermediate directories if they don't exist
# Set the working directory
setwd(dir.name)
# Generate and save plots for each combination
# Iterate over each rate type
for (ratetype2calc in ratetype2calc) {
  pwalk(paramcomb2plot, savebarplot_XXRbyMABF_cond)
}
# Change back to original working directory
setwd("D:/wnk/PhD/PhD research/dissertation/code/hprc")
#############################################################################
#Rates conditioned on one parameter when true effect is 0.2: QRP environment#
#############################################################################
# Set up parameters
paramcomb2plot <- expand.grid(qrplevel= c("none", "medium", "high"))
TrueES <- 0.2
data2plot <- rates_0.2null
ratetype2calc <- c("ADR", "ADRTE", "ADRNE", "SAR", "FAR", "FCR", "SCR", "SAR1", "SAR2", "FAR1", "FAR2", "FCR1", "FCR2", "SCR1", "SCR2", "TPR", "FPR", "TSR", "FSR", "TFR", "FFR", "TSR1", "TSR2", "FSR1", "FSR2", "TFR1", "TFR2", "FFR1", "FFR2")
dir2save <- c("0.2null")
# Generate directory name based on parameters
filename_parts <- lapply(names(paramcomb2plot), function(name) {
  if (!is.null(paramcomb2plot[[name]])) {
    paste0(name)
  }
})
# Combine factor names generated by filename_parts into a single directory name
dir.name <- paste0("D:/wnk/PhD/PhD research/dissertation/code/hprc/plots/bar plot/by MABF methods/conditional rate/", length(paramcomb2plot), " factor/", paste0(filename_parts, collapse = "_"))
# Create the directory
dir.create(dir.name, recursive = TRUE)  # Set recursive = TRUE to create intermediate directories if they don't exist
# Set the working directory
setwd(dir.name)
# Generate and save plots for each combination
# Iterate over each rate type
for (ratetype2calc in ratetype2calc) {
  pwalk(paramcomb2plot, savebarplot_XXRbyMABF_cond)
}
# Change back to original working directory
setwd("D:/wnk/PhD/PhD research/dissertation/code/hprc")
#####################################################################################
#Rates conditioned on one parameter when true effect is 0.5: replication sample size#
#####################################################################################
# Set up parameters
paramcomb2plot <- expand.grid(repsampsize = c("small", "medium", "large"))
TrueES <- 0.5
data2plot <- rates_0.5null
ratetype2calc <- c("ADR", "ADRTE", "ADRNE", "SAR", "FAR", "FCR", "SCR", "SAR1", "SAR2", "FAR1", "FAR2", "FCR1", "FCR2", "SCR1", "SCR2", "TPR", "FPR", "TSR", "FSR", "TFR", "FFR", "TSR1", "TSR2", "FSR1", "FSR2", "TFR1", "TFR2", "FFR1", "FFR2")
dir2save <- c("0.5null")
# Generate directory name based on parameters
filename_parts <- lapply(names(paramcomb2plot), function(name) {
  if (!is.null(paramcomb2plot[[name]])) {
    paste0(name)
  }
})
# Combine factor names generated by filename_parts into a single directory name
dir.name <- paste0("D:/wnk/PhD/PhD research/dissertation/code/hprc/plots/bar plot/by MABF methods/conditional rate/", length(paramcomb2plot), " factor/", paste0(filename_parts, collapse = "_"))
# Create the directory
dir.create(dir.name, recursive = TRUE)  # Set recursive = TRUE to create intermediate directories if they don't exist
# Set the working directory
setwd(dir.name)
# Generate and save plots for each combination
# Iterate over each rate type
for (ratetype2calc in ratetype2calc) {
  pwalk(paramcomb2plot, savebarplot_XXRbyMABF_cond)
}
# Change back to original working directory
setwd("D:/wnk/PhD/PhD research/dissertation/code/hprc")
################################################################################
#Rates conditioned on one parameter when true effect is 0.5: replication number#
################################################################################
# Set up parameters
paramcomb2plot <- expand.grid(repnumber = c("small", "medium", "large"))
TrueES <- 0.5
data2plot <- rates_0.5null
ratetype2calc <- c("ADR", "ADRTE", "ADRNE", "SAR", "FAR", "FCR", "SCR", "SAR1", "SAR2", "FAR1", "FAR2", "FCR1", "FCR2", "SCR1", "SCR2", "TPR", "FPR", "TSR", "FSR", "TFR", "FFR", "TSR1", "TSR2", "FSR1", "FSR2", "TFR1", "TFR2", "FFR1", "FFR2")
dir2save <- c("0.5null")
# Generate directory name based on parameters
filename_parts <- lapply(names(paramcomb2plot), function(name) {
  if (!is.null(paramcomb2plot[[name]])) {
    paste0(name)
  }
})
# Combine factor names generated by filename_parts into a single directory name
dir.name <- paste0("D:/wnk/PhD/PhD research/dissertation/code/hprc/plots/bar plot/by MABF methods/conditional rate/", length(paramcomb2plot), " factor/", paste0(filename_parts, collapse = "_"))
# Create the directory
dir.create(dir.name, recursive = TRUE)  # Set recursive = TRUE to create intermediate directories if they don't exist
# Set the working directory
setwd(dir.name)
# Generate and save plots for each combination
# Iterate over each rate type
for (ratetype2calc in ratetype2calc) {
  pwalk(paramcomb2plot, savebarplot_XXRbyMABF_cond)
}
# Change back to original working directory
setwd("D:/wnk/PhD/PhD research/dissertation/code/hprc")
##################################################################################
#Rates conditioned on one parameter when true effect is 0.5: original sample size#
##################################################################################
# Set up parameters
paramcomb2plot <- expand.grid(origsampsize = c("small", "medium", "large"))
TrueES <- 0.5
data2plot <- rates_0.5null
ratetype2calc <- c("ADR", "ADRTE", "ADRNE", "SAR", "FAR", "FCR", "SCR", "SAR1", "SAR2", "FAR1", "FAR2", "FCR1", "FCR2", "SCR1", "SCR2", "TPR", "FPR", "TSR", "FSR", "TFR", "FFR", "TSR1", "TSR2", "FSR1", "FSR2", "TFR1", "TFR2", "FFR1", "FFR2")
dir2save <- c("0.5null")
# Generate directory name based on parameters
filename_parts <- lapply(names(paramcomb2plot), function(name) {
  if (!is.null(paramcomb2plot[[name]])) {
    paste0(name)
  }
})
# Combine factor names generated by filename_parts into a single directory name
dir.name <- paste0("D:/wnk/PhD/PhD research/dissertation/code/hprc/plots/bar plot/by MABF methods/conditional rate/", length(paramcomb2plot), " factor/", paste0(filename_parts, collapse = "_"))
# Create the directory
dir.create(dir.name, recursive = TRUE)  # Set recursive = TRUE to create intermediate directories if they don't exist
# Set the working directory
setwd(dir.name)
# Generate and save plots for each combination
# Iterate over each rate type
for (ratetype2calc in ratetype2calc) {
  pwalk(paramcomb2plot, savebarplot_XXRbyMABF_cond)
}
# Change back to original working directory
setwd("D:/wnk/PhD/PhD research/dissertation/code/hprc")
#############################################################################
#Rates conditioned on one parameter when true effect is 0.5: QRP environment#
#############################################################################
# Set up parameters
# You can add up to 4 parameters to plot, for example:
# paramcomb2plot <- expand.grid(qrplevel= c("none", "medium", "high"),
#                               origsampsize = c("small", "medium", "large"))
paramcomb2plot <- expand.grid(qrplevel= c("none", "medium", "high"))
TrueES <- 0.5
data2plot <- rates_0.5null
ratetype2calc <- c("ADR", "ADRTE", "ADRNE", "SAR", "FAR", "FCR", "SCR", "SAR1", "SAR2", "FAR1", "FAR2", "FCR1", "FCR2", "SCR1", "SCR2", "TPR", "FPR", "TSR", "FSR", "TFR", "FFR", "TSR1", "TSR2", "FSR1", "FSR2", "TFR1", "TFR2", "FFR1", "FFR2")
dir2save <- c("0.5null")
# Generate directory name based on parameters
filename_parts <- lapply(names(paramcomb2plot), function(name) {
  if (!is.null(paramcomb2plot[[name]])) {
    paste0(name)
  }
})
# Combine factor names generated by filename_parts into a single directory name
dir.name <- paste0("D:/wnk/PhD/PhD research/dissertation/code/hprc/plots/bar plot/by MABF methods/conditional rate/", length(paramcomb2plot), " factor/", paste0(filename_parts, collapse = "_"))
# Create the directory
dir.create(dir.name, recursive = TRUE)  # Set recursive = TRUE to create intermediate directories if they don't exist
# Set the working directory
setwd(dir.name)
# Generate and save plots for each combination
# Iterate over each rate type
for (ratetype2calc in ratetype2calc) {
  pwalk(paramcomb2plot, savebarplot_XXRbyMABF_cond)
}
# Change back to original working directory
setwd("D:/wnk/PhD/PhD research/dissertation/code/hprc")
