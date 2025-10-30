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
#Step 0.2: Load FEMA outcome lists#
###################################
# Set the path to the directory containing RDS files
folder_path <- "./MABFoutcomes/500500/MA/FEMA" 
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}

# Combine beta matricies
#####################################
#Step 1: Combine splitted FEMA lists#
#####################################
# Combine 9 individual lists into an overall list
FEMA_beta_lists <- c(FEMA_beta_lists_nonelow_part1,FEMA_beta_lists_nonelow_part2,FEMA_beta_lists_nonelow_part3,
                     FEMA_beta_lists_nonelow_part4,FEMA_beta_lists_nonelow_part5,FEMA_beta_lists_nonelow_part6,
                     FEMA_beta_lists_nonelow_part7,FEMA_beta_lists_nonelow_part8,FEMA_beta_lists_nonelow_part9,
                     FEMA_beta_lists_medmed_part1,FEMA_beta_lists_medmed_part2,FEMA_beta_lists_medmed_part3,
                     FEMA_beta_lists_medmed_part4,FEMA_beta_lists_medmed_part5,FEMA_beta_lists_medmed_part6,
                     FEMA_beta_lists_medmed_part7,FEMA_beta_lists_medmed_part8,FEMA_beta_lists_medmed_part9,
                     FEMA_beta_lists_hihi_part1,FEMA_beta_lists_hihi_part2,FEMA_beta_lists_hihi_part3,
                     FEMA_beta_lists_hihi_part4,FEMA_beta_lists_hihi_part5,FEMA_beta_lists_hihi_part6,
                     FEMA_beta_lists_hihi_part7,FEMA_beta_lists_hihi_part8,FEMA_beta_lists_hihi_part9)
# Remove wont-be-used lists
remove(FEMA_beta_lists_nonelow_part1,FEMA_beta_lists_nonelow_part2,FEMA_beta_lists_nonelow_part3,
       FEMA_beta_lists_nonelow_part4,FEMA_beta_lists_nonelow_part5,FEMA_beta_lists_nonelow_part6,
       FEMA_beta_lists_nonelow_part7,FEMA_beta_lists_nonelow_part8,FEMA_beta_lists_nonelow_part9,
       FEMA_beta_lists_medmed_part1,FEMA_beta_lists_medmed_part2,FEMA_beta_lists_medmed_part3,
       FEMA_beta_lists_medmed_part4,FEMA_beta_lists_medmed_part5,FEMA_beta_lists_medmed_part6,
       FEMA_beta_lists_medmed_part7,FEMA_beta_lists_medmed_part8,FEMA_beta_lists_medmed_part9,
       FEMA_beta_lists_hihi_part1,FEMA_beta_lists_hihi_part2,FEMA_beta_lists_hihi_part3,
       FEMA_beta_lists_hihi_part4,FEMA_beta_lists_hihi_part5,FEMA_beta_lists_hihi_part6,
       FEMA_beta_lists_hihi_part7,FEMA_beta_lists_hihi_part8,FEMA_beta_lists_hihi_part9)

################################################
#Step 2.0: Regroup and restucture Fixed-Effect Meta-Analysis' beta outcomes#
################################################
# The complete FEMA beta list contains sublists of FEMA betas based on delta 0, 0.2, 0.5, separately
# the following code will divide these sublists into two groups: one with 0 and 0.2, one with 0 and 0.5
# By doing so, we can combine sublists of 0 and 0.2 together and combine sublists of 0 and 0.5 together
name_witht_0.2 <- !startsWith(names(FEMA_beta_lists), "0.2")
FEMA_beta_lists_0.5null <- FEMA_beta_lists[name_witht_0.2]
name_witht_0.5 <- !startsWith(names(FEMA_beta_lists), "0.5")
FEMA_beta_lists_0.2null <- FEMA_beta_lists[name_witht_0.5]
remove(FEMA_beta_lists)
# We also need to split up df_lists, otherwise p values will be incorrectly assigned
df_lists_0.5null <- df_lists[name_witht_0.2]
df_lists_0.2null <- df_lists[name_witht_0.5]
############################################################################
#Step 2.1: Create IDs for sublists and contained matricies in FEMA beta outcomes#
############################################################################
num_deltai = 500 # This has to match the parameter setting in Step 1 Data Genearation
# Create IDs used to name the main sublists for delta=0 and 0.2
params_FEMAbetalists_002 <- expand.grid(
  delta = "0,0.2",
  orig.n = c(20,50,200),
  qrpEnv = c("none","medium","high"),
  censorFunc = c("low","medium","high"),
  stringsAsFactors = FALSE
)
# Remove unused conditions
params_FEMAbetalists_002 <- params_FEMAbetalists_002[-c(4:12,16:24),]
row.names(params_FEMAbetalists_002) <- NULL
IDs_FEMAbetalists_002 <- purrr::pmap(params_FEMAbetalists_002,naming_MABFlists)
# Create IDs used to name the main sublists for delta=0 and 0.5
params_FEMAbetalists_005 <- expand.grid(
  delta = "0,0.5",
  orig.n = c(20,50,200),
  qrpEnv = c("none","medium","high"),
  censorFunc = c("low","medium","high"),
  stringsAsFactors = FALSE
)
# Remove unused conditions
params_FEMAbetalists_005 <- params_FEMAbetalists_005[-c(4:12,16:24),]
row.names(params_FEMAbetalists_005) <- NULL
IDs_FEMAbetalists_005 <- purrr::pmap(params_FEMAbetalists_005,naming_MABFlists)
# Create IDs used to name contained matricies
params_FEMAbetasubls <- expand.grid(
  rep.number = c(2, 5, 10),
  rep.n = c(40, 100, 400),
  stringsAsFactors = FALSE
)
IDs_FEMAbetasubls <- purrr::pmap(params_FEMAbetasubls,naming_MABFsubls)

########################################################
#Step 2.2: Combine null and true effects of delta = 0.2#
########################################################
# Initialize an empty list to store the final groups
# FEMA_beta_lists_0.2null_regrouped combines sublists with delta of 0 and 0.2
FEMA_beta_lists_0.2null_regrouped = list()
# Iterate through the list of sublists two at a time
for(i in seq(1, length(FEMA_beta_lists_0.2null), by = 2)) {
  # Combine vectors from the current pair of sublists
  combined = combine_vectors(FEMA_beta_lists_0.2null[[i]], FEMA_beta_lists_0.2null[[i+1]], chunk_size=num_deltai)
  # Add the combined result to the final groups list
  FEMA_beta_lists_0.2null_regrouped[[length(FEMA_beta_lists_0.2null_regrouped) + 1]] = combined
}
# Naming FEMA_beta sublists and matricies within
# Naming the main sublists
names(FEMA_beta_lists_0.2null_regrouped) <- IDs_FEMAbetalists_002
# Naming contained matricies
FEMA_beta_lists_0.2null_regrouped <- lapply(FEMA_beta_lists_0.2null_regrouped, function(sublist) {
  names(sublist) <- IDs_FEMAbetasubls
  return(sublist)
})

########################################################
#Step 2.3: Combine null and true effects of delta = 0.5#
########################################################
# Initialize an empty list to store the final groups
# FEMA_beta_lists_0.5null_regrouped combines sublists with delta of 0 and 0.5
FEMA_beta_lists_0.5null_regrouped = list()
# Iterate through the list of sublists two at a time
for(i in seq(1, length(FEMA_beta_lists_0.5null), by = 2)) {
  # Combine vectors from the current pair of sublists
  combined = combine_vectors(FEMA_beta_lists_0.5null[[i]], FEMA_beta_lists_0.5null[[i+1]], chunk_size=num_deltai)
  # Add the combined result to the final groups list
  FEMA_beta_lists_0.5null_regrouped[[length(FEMA_beta_lists_0.5null_regrouped) + 1]] = combined
}
# Naming FEMA_beta sublists and matricies within
# Naming the main sublists
names(FEMA_beta_lists_0.5null_regrouped) <- IDs_FEMAbetalists_005
# Naming contained matricies
FEMA_beta_lists_0.5null_regrouped <- lapply(FEMA_beta_lists_0.5null_regrouped, function(sublist) {
  names(sublist) <- IDs_FEMAbetasubls
  return(sublist)
})

# Combine p matricies
#####################################
#Step 1: Combine splitted FEMA lists#
#####################################
# Combine 9 individual lists into an overall list
FEMA_p_lists <- c(FEMA_p_lists_nonelow_part1,FEMA_p_lists_nonelow_part2,FEMA_p_lists_nonelow_part3,
                  FEMA_p_lists_nonelow_part4,FEMA_p_lists_nonelow_part5,FEMA_p_lists_nonelow_part6,
                  FEMA_p_lists_nonelow_part7,FEMA_p_lists_nonelow_part8,FEMA_p_lists_nonelow_part9,
                  FEMA_p_lists_medmed_part1,FEMA_p_lists_medmed_part2,FEMA_p_lists_medmed_part3,
                  FEMA_p_lists_medmed_part4,FEMA_p_lists_medmed_part5,FEMA_p_lists_medmed_part6,
                  FEMA_p_lists_medmed_part7,FEMA_p_lists_medmed_part8,FEMA_p_lists_medmed_part9,
                  FEMA_p_lists_hihi_part1,FEMA_p_lists_hihi_part2,FEMA_p_lists_hihi_part3,
                  FEMA_p_lists_hihi_part4,FEMA_p_lists_hihi_part5,FEMA_p_lists_hihi_part6,
                  FEMA_p_lists_hihi_part7,FEMA_p_lists_hihi_part8,FEMA_p_lists_hihi_part9)
# Remove wont-be-used lists
remove(FEMA_p_lists_nonelow_part1,FEMA_p_lists_nonelow_part2,FEMA_p_lists_nonelow_part3,
       FEMA_p_lists_nonelow_part4,FEMA_p_lists_nonelow_part5,FEMA_p_lists_nonelow_part6,
       FEMA_p_lists_nonelow_part7,FEMA_p_lists_nonelow_part8,FEMA_p_lists_nonelow_part9,
       FEMA_p_lists_medmed_part1,FEMA_p_lists_medmed_part2,FEMA_p_lists_medmed_part3,
       FEMA_p_lists_medmed_part4,FEMA_p_lists_medmed_part5,FEMA_p_lists_medmed_part6,
       FEMA_p_lists_medmed_part7,FEMA_p_lists_medmed_part8,FEMA_p_lists_medmed_part9,
       FEMA_p_lists_hihi_part1,FEMA_p_lists_hihi_part2,FEMA_p_lists_hihi_part3,
       FEMA_p_lists_hihi_part4,FEMA_p_lists_hihi_part5,FEMA_p_lists_hihi_part6,
       FEMA_p_lists_hihi_part7,FEMA_p_lists_hihi_part8,FEMA_p_lists_hihi_part9)

################################################
#Step 2.0: Regroup and restucture Fixed-Effect Meta-Analysis' beta outcomes#
################################################
# The complete FEMA beta list contains sublists of FEMA betas based on delta 0, 0.2, 0.5, separately
# the following code will divide these sublists into two groups: one with 0 and 0.2, one with 0 and 0.5
# By doing so, we can combine sublists of 0 and 0.2 together and combine sublists of 0 and 0.5 together
name_witht_0.2 <- !startsWith(names(FEMA_p_lists), "0.2")
FEMA_p_lists_0.5null <- FEMA_p_lists[name_witht_0.2]
name_witht_0.5 <- !startsWith(names(FEMA_p_lists), "0.5")
FEMA_p_lists_0.2null <- FEMA_p_lists[name_witht_0.5]
remove(FEMA_p_lists)
# We also need to split up df_lists, otherwise p values will be incorrectly assigned
df_lists_0.5null <- df_lists[name_witht_0.2]
df_lists_0.2null <- df_lists[name_witht_0.5]
remove(df_lists)
############################################################################
#Step 2.1: Create IDs for sublists and contained matricies in FEMA beta outcomes#
############################################################################
num_deltai = 500 # This has to match the parameter setting in Step 1 Data Genearation
# Create IDs used to name the main sublists for delta=0 and 0.2
params_FEMAplists_002 <- expand.grid(
  delta = "0,0.2",
  orig.n = c(20,50,200),
  qrpEnv = c("none","medium","high"),
  censorFunc = c("low","medium","high"),
  stringsAsFactors = FALSE
)
# Remove unused conditions
params_FEMAplists_002 <- params_FEMAplists_002[-c(4:12,16:24),]
row.names(params_FEMAplists_002) <- NULL
IDs_FEMAplists_002 <- purrr::pmap(params_FEMAplists_002,naming_MABFlists)
# Create IDs used to name the main sublists for delta=0 and 0.5
params_FEMAplists_005 <- expand.grid(
  delta = "0,0.5",
  orig.n = c(20,50,200),
  qrpEnv = c("none","medium","high"),
  censorFunc = c("low","medium","high"),
  stringsAsFactors = FALSE
)
# Remove unused conditions
params_FEMAplists_005 <- params_FEMAplists_005[-c(4:12,16:24),]
row.names(params_FEMAplists_005) <- NULL
IDs_FEMAplists_005 <- purrr::pmap(params_FEMAplists_005,naming_MABFlists)
# Create IDs used to name contained matricies
params_FEMApsubls <- expand.grid(
  rep.number = c(2, 5, 10),
  rep.n = c(40, 100, 400),
  stringsAsFactors = FALSE
)
IDs_FEMApsubls <- purrr::pmap(params_FEMApsubls,naming_MABFsubls)

########################################################
#Step 2.2: Combine null and true effects of delta = 0.2#
########################################################
# Initialize an empty list to store the final groups
# FEMA_p_lists_0.2null_regrouped combines sublists with delta of 0 and 0.2
FEMA_p_lists_0.2null_regrouped = list()
# Iterate through the list of sublists two at a time
for(i in seq(1, length(FEMA_p_lists_0.2null), by = 2)) {
  # Combine vectors from the current pair of sublists
  combined = combine_vectors(FEMA_p_lists_0.2null[[i]], FEMA_p_lists_0.2null[[i+1]], chunk_size=num_deltai)
  # Add the combined result to the final groups list
  FEMA_p_lists_0.2null_regrouped[[length(FEMA_p_lists_0.2null_regrouped) + 1]] = combined
}
# Naming FEMA_p sublists and matricies within
# Naming the main sublists
names(FEMA_p_lists_0.2null_regrouped) <- IDs_FEMAplists_002
# Naming contained matricies
FEMA_p_lists_0.2null_regrouped <- lapply(FEMA_p_lists_0.2null_regrouped, function(sublist) {
  names(sublist) <- IDs_FEMApsubls
  return(sublist)
})

########################################################
#Step 2.3: Combine null and true effects of delta = 0.5#
########################################################
# Initialize an empty list to store the final groups
# FEMA_p_lists_0.5null_regrouped combines sublists with delta of 0 and 0.5
FEMA_p_lists_0.5null_regrouped = list()
# Iterate through the list of sublists two at a time
for(i in seq(1, length(FEMA_p_lists_0.5null), by = 2)) {
  # Combine vectors from the current pair of sublists
  combined = combine_vectors(FEMA_p_lists_0.5null[[i]], FEMA_p_lists_0.5null[[i+1]], chunk_size=num_deltai)
  # Add the combined result to the final groups list
  FEMA_p_lists_0.5null_regrouped[[length(FEMA_p_lists_0.5null_regrouped) + 1]] = combined
}
# Naming FEMA_p sublists and matricies within
# Naming the main sublists
names(FEMA_p_lists_0.5null_regrouped) <- IDs_FEMAplists_005
# Naming contained matricies
FEMA_p_lists_0.5null_regrouped <- lapply(FEMA_p_lists_0.5null_regrouped, function(sublist) {
  names(sublist) <- IDs_FEMApsubls
  return(sublist)
})


########################################################
#Last Step: Combine beta and p values into a single matrix#
########################################################
FEMA_lists_0.2null_regrouped <- mapply(combine_matrices, FEMA_beta_lists_0.2null_regrouped, FEMA_p_lists_0.2null_regrouped, SIMPLIFY = FALSE)
FEMA_lists_0.5null_regrouped <- mapply(combine_matrices, FEMA_beta_lists_0.5null_regrouped, FEMA_p_lists_0.5null_regrouped, SIMPLIFY = FALSE)

# Save data as RDS file
#saveRDS(FEMA_lists_0.2null_regrouped, "./MABF4ROC/FEMA_lists_0.2null_regrouped.RDS")
#saveRDS(FEMA_lists_0.5null_regrouped, "./MABF4ROC/FEMA_lists_0.5null_regrouped.RDS")
