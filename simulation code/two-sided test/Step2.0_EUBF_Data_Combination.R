#########################################
#Step 0: Load necessary custom functions#
#########################################
source("./helper functions.R")
source("./helper functions_XXR.R")
################################################
#Step 0.1: Load original study dataset df_lists#
################################################
# Set the path to the directory containing RDS files
folder_path <- "./OGDG/500" 
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}
df_lists <- df_lists_500 #for convenience purpose
remove(df_lists_500)
###################################
#Step 0.2: Load MABF outcome lists#
###################################
# Set the path to the directory containing RDS files
folder_path <- "./MABFoutcomes/500500/EUBF" 
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}
#####################################
#Step 1: Combine splitted MABF lists#
#####################################
# Combine 9 individual lists into an overall list
EUBF_lists <- c(EUBF_lists_nonelow_part1,EUBF_lists_nonelow_part2,EUBF_lists_nonelow_part3,
                  EUBF_lists_nonelow_part4,EUBF_lists_nonelow_part5,EUBF_lists_nonelow_part6,
                  EUBF_lists_nonelow_part7,EUBF_lists_nonelow_part8,EUBF_lists_nonelow_part9,
                  EUBF_lists_medmed_part1,EUBF_lists_medmed_part2,EUBF_lists_medmed_part3,
                  EUBF_lists_medmed_part4,EUBF_lists_medmed_part5,EUBF_lists_medmed_part6,
                  EUBF_lists_medmed_part7,EUBF_lists_medmed_part8,EUBF_lists_medmed_part9,
                  EUBF_lists_hihi_part1,EUBF_lists_hihi_part2,EUBF_lists_hihi_part3,
                  EUBF_lists_hihi_part4,EUBF_lists_hihi_part5,EUBF_lists_hihi_part6,
                  EUBF_lists_hihi_part7,EUBF_lists_hihi_part8,EUBF_lists_hihi_part9)
# Remove wont-be-used lists
remove(EUBF_lists_nonelow_part1,EUBF_lists_nonelow_part2,EUBF_lists_nonelow_part3,
       EUBF_lists_nonelow_part4,EUBF_lists_nonelow_part5,EUBF_lists_nonelow_part6,
       EUBF_lists_nonelow_part7,EUBF_lists_nonelow_part8,EUBF_lists_nonelow_part9,
       EUBF_lists_medmed_part1,EUBF_lists_medmed_part2,EUBF_lists_medmed_part3,
       EUBF_lists_medmed_part4,EUBF_lists_medmed_part5,EUBF_lists_medmed_part6,
       EUBF_lists_medmed_part7,EUBF_lists_medmed_part8,EUBF_lists_medmed_part9,
       EUBF_lists_hihi_part1,EUBF_lists_hihi_part2,EUBF_lists_hihi_part3,
       EUBF_lists_hihi_part4,EUBF_lists_hihi_part5,EUBF_lists_hihi_part6,
       EUBF_lists_hihi_part7,EUBF_lists_hihi_part8,EUBF_lists_hihi_part9)

################################################
#Step 2.0: Regroup and restucture MABF outcomes#
################################################
# The complete MABF list contains sublists of MABFs based on delta 0, 0.2, 0.5, separately
# the following code will divide these sublists into two groups: one with 0 and 0.2, one with 0 and 0.5
# By doing so, we can combine sublists of 0 and 0.2 together and combine sublists of 0 and 0.5 together
name_witht_0.2 <- !startsWith(names(EUBF_lists), "0.2")
EUBF_lists_0.5null <- EUBF_lists[name_witht_0.2]
name_witht_0.5 <- !startsWith(names(EUBF_lists), "0.5")
EUBF_lists_0.2null <- EUBF_lists[name_witht_0.5]
remove(EUBF_lists)
# We also need to split up df_lists, otherwise p values will be incorrectly assigned
df_lists_0.5null <- df_lists[name_witht_0.2]
df_lists_0.2null <- df_lists[name_witht_0.5]
remove(df_lists)
############################################################################
#Step 2.1: Create IDs for sublists and contained matricies in MABF outcomes#
############################################################################
num_deltai = 500 # This has to match the parameter setting in Step 1 Data Genearation
# Create IDs used to name the main sublists for delta=0 and 0.2
params_MABFlists_002 <- expand.grid(
  delta = "0,0.2",
  orig.n = c(20,50,200),
  qrpEnv = c("none","medium","high"),
  censorFunc = c("low","medium","high"),
  stringsAsFactors = FALSE
)
# Remove unused conditions
params_MABFlists_002 <- params_MABFlists_002[-c(4:12,16:24),]
row.names(params_MABFlists_002) <- NULL
IDs_MABFlists_002 <- purrr::pmap(params_MABFlists_002,naming_MABFlists)
# Create IDs used to name the main sublists for delta=0 and 0.5
params_MABFlists_005 <- expand.grid(
  delta = "0,0.5",
  orig.n = c(20,50,200),
  qrpEnv = c("none","medium","high"),
  censorFunc = c("low","medium","high"),
  stringsAsFactors = FALSE
)
# Remove unused conditions
params_MABFlists_005 <- params_MABFlists_005[-c(4:12,16:24),]
row.names(params_MABFlists_005) <- NULL
IDs_MABFlists_005 <- purrr::pmap(params_MABFlists_005,naming_MABFlists)
# Create IDs used to name contained matricies
params_MABFsubls <- expand.grid(
  rep.number = c(2, 5, 10),
  rep.n = c(40, 100, 400),
  stringsAsFactors = FALSE
)
IDs_MABFsubls <- purrr::pmap(params_MABFsubls,naming_MABFsubls)

########################################################
#Step 2.2: Combine null and true effects of delta = 0.2#
########################################################
# Initialize an empty list to store the final groups
# EUBF_lists_0.2null_regrouped combines sublists with delta of 0 and 0.2
EUBF_lists_0.2null_regrouped = list()
# Iterate through the list of sublists two at a time
for(i in seq(1, length(EUBF_lists_0.2null), by = 2)) {
  # Combine vectors from the current pair of sublists
  combined = combine_vectors(EUBF_lists_0.2null[[i]], EUBF_lists_0.2null[[i+1]], chunk_size=num_deltai)
  # Add the combined result to the final groups list
  EUBF_lists_0.2null_regrouped[[length(EUBF_lists_0.2null_regrouped) + 1]] = combined
}
# Naming EUBF sublists and matricies within
# Naming the main sublists
names(EUBF_lists_0.2null_regrouped) <- IDs_MABFlists_002
# Naming contained matricies
EUBF_lists_0.2null_regrouped <- lapply(EUBF_lists_0.2null_regrouped, function(sublist) {
  names(sublist) <- IDs_MABFsubls
  return(sublist)
})

########################################################
#Step 2.3: Combine null and true effects of delta = 0.5#
########################################################
# Initialize an empty list to store the final groups
# EUBF_lists_0.5null_regrouped combines sublists with delta of 0 and 0.5
EUBF_lists_0.5null_regrouped = list()
# Iterate through the list of sublists two at a time
for(i in seq(1, length(EUBF_lists_0.5null), by = 2)) {
  # Combine vectors from the current pair of sublists
  combined = combine_vectors(EUBF_lists_0.5null[[i]], EUBF_lists_0.5null[[i+1]], chunk_size=num_deltai)
  # Add the combined result to the final groups list
  EUBF_lists_0.5null_regrouped[[length(EUBF_lists_0.5null_regrouped) + 1]] = combined
}
# Naming EUBF sublists and matricies within
# Naming the main sublists
names(EUBF_lists_0.5null_regrouped) <- IDs_MABFlists_005
# Naming contained matricies
EUBF_lists_0.5null_regrouped <- lapply(EUBF_lists_0.5null_regrouped, function(sublist) {
  names(sublist) <- IDs_MABFsubls
  return(sublist)
})


#saveRDS(EUBF_lists_0.2null_regrouped, file = "./MABF4ROC/4TPFP/EUBF_lists_0.2null_regrouped.RDS")
#saveRDS(EUBF_lists_0.5null_regrouped, file = "./MABF4ROC/4TPFP/EUBF_lists_0.5null_regrouped.RDS")
