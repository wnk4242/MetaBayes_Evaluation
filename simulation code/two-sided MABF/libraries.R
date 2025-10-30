# 08042025 Update: Installed rstan, metaBMA, etc in the scratch directory (see installing libraries on servers.txt). Set up a new package directory here to load these packages in addtion to those installed in the home directory.
# 08042025 Update: I used to have code in here to install packages but now followed best practice to remove them. One should install packages before SLURM jobs, not during.
# 04172024 Update: Added foreach and doParallel packages
# Define the packages
packages <- c("metafor", "MASS", "pwr", "BayesFactor", "metaBMA", 
              "openxlsx", "tidyverse", "ggplot2", "LaplacesDemon", 
              "future", "future.apply", "progressr", "foreach", "doParallel")

# Tell R where to find packages installed in /scratch
.libPaths(c("/scratch/user/wangnaike/wnk4242finalsim/Rlibs", .libPaths()))

# Load packages
for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    stop(paste("Package", package, "is not installed in any accessible library paths."))
  }
}

# 
# #for conducting meta-analysis
# library(metafor)
# 
# #for general stats functions
# library(MASS)
# library(pwr) #used for simMA
# 
# #for Bayes factor functions
# library(BayesFactor)
# library(metaBMA)
# 
# #for exporting data
# library(openxlsx)#export dataframe to excel
# 
# #for data manipulation
# library(tidyverse)
# 
# #for visualization
# library(ggplot2)
# 
# #for extra distributions
# library(LaplacesDemon) #used for Bayesian averaged meta-analysis
# 
# #for parallel process
# library(future)
# library(future.apply)
# library(progressr)
# #if parallel package is used, it may mask certain functions in purrr
# #therefore, use parallel::function and purrr::function, instead of loading these packages
# #library(parallel) 
# library(foreach)
# library(doParallel)