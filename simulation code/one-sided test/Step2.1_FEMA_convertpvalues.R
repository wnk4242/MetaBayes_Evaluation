# Convert two-side FEMA p values into one-sided p values
# Set the folder containing the RDS files
data_dir <- "MABF4ROC/4TPFP"

# Import the two-sided FEMA meta-analysis datasets
FEMA_lists_0.2null_regrouped <- readRDS(file.path(data_dir, "FEMA_lists_0.2null_regrouped.RDS"))
FEMA_lists_0.5null_regrouped <- readRDS(file.path(data_dir, "FEMA_lists_0.5null_regrouped.RDS"))

# Convert one 1000 x 1000 matrix from two-sided to one-sided p-values
convert_two_to_one_sided <- function(mat) {
  # Columns 1:500 contain summary effect sizes
  es_mat <- mat[, 1:500]
  
  # Columns 501:1000 contain the corresponding two-sided p-values
  p2_mat <- mat[, 501:1000]
  
  # Convert to one-sided p-values for H0: theta = 0 vs H1: theta > 0
  p1_mat <- ifelse(
    es_mat > 0,
    p2_mat / 2,
    ifelse(es_mat < 0, 1 - p2_mat / 2, 0.5)
  )
  
  # Return the same structure: first 500 columns = effect sizes,
  # last 500 columns = one-sided p-values
  cbind(es_mat, p1_mat)
}

# Apply the conversion to the full nested list
convert_nested_list <- function(lst) {
  lapply(lst, function(condition_list) {
    lapply(condition_list, convert_two_to_one_sided)
  })
}

# Create one-sided versions
FEMA_lists_0.2null_regrouped_onesided <- convert_nested_list(FEMA_lists_0.2null_regrouped)
FEMA_lists_0.5null_regrouped_onesided <- convert_nested_list(FEMA_lists_0.5null_regrouped)


# Save the converted datasets
saveRDS(
  FEMA_lists_0.2null_regrouped_onesided,
  file = file.path(data_dir, "FEMA_lists_0.2null_regrouped_onesided.RDS")
)

saveRDS(
  FEMA_lists_0.5null_regrouped_onesided,
  file = file.path(data_dir, "FEMA_lists_0.5null_regrouped_onesided.RDS")
)


##########################
# NOTE: FEMA is random effects, FEMA2 is fixed effect
# Convert two-sided FEMA p values into one-sided 
# Folder containing the FEMA beta and p-value RDS files
data_dir <- "MABFoutcomes/500500/MA/FEMA/"


# Convert two-sided p-values to one-sided p-values for H0: theta = 0 vs H1: theta > 0
convert_pvec_to_onesided <- function(beta_vec, p2_vec) {
  if (length(beta_vec) != length(p2_vec)) {
    stop("beta_vec and p2_vec must have the same length.")
  }
  
  ifelse(
    beta_vec > 0,
    p2_vec / 2,
    ifelse(beta_vec < 0, 1 - p2_vec / 2, 0.5)
  )
}


# Convert one paired beta-list object and p-list object
convert_nested_p_list <- function(beta_obj, p_obj) {
  if (!identical(names(beta_obj), names(p_obj))) {
    stop("Top-level names do not match between beta and p objects.")
  }
  
  out <- vector("list", length(beta_obj))
  names(out) <- names(beta_obj)
  
  for (cond_name in names(beta_obj)) {
    beta_cond <- beta_obj[[cond_name]]
    p_cond    <- p_obj[[cond_name]]
    
    if (!identical(names(beta_cond), names(p_cond))) {
      stop(sprintf("Inner names do not match for condition: %s", cond_name))
    }
    
    out[[cond_name]] <- lapply(names(beta_cond), function(x) {
      convert_pvec_to_onesided(beta_cond[[x]], p_cond[[x]])
    })
    names(out[[cond_name]]) <- names(beta_cond)
  }
  
  out
}


# Find all beta and two-sided p-value RDS files
beta_files <- list.files(
  data_dir,
  pattern = "^FEMA_beta_lists_.*\\.RDS$",
  full.names = TRUE
)

p_files <- list.files(
  data_dir,
  pattern = "^FEMA_p_lists_.*\\.RDS$",
  full.names = TRUE
)


# Match files by shared suffix
beta_keys <- sub("^FEMA_beta_lists_", "", basename(beta_files))
p_keys    <- sub("^FEMA_p_lists_", "", basename(p_files))

common_keys <- intersect(beta_keys, p_keys)


# Generate one-sided p-value RDS files
for (key in common_keys) {
  beta_path <- file.path(data_dir, paste0("FEMA_beta_lists_", key))
  p_path    <- file.path(data_dir, paste0("FEMA_p_lists_", key))
  
  beta_obj <- readRDS(beta_path)
  p_obj    <- readRDS(p_path)
  
  p1_obj <- convert_nested_p_list(beta_obj, p_obj)
  
  out_path <- file.path(data_dir, paste0("FEMA_p1_lists_", key))
  saveRDS(p1_obj, out_path)
  
  message("Saved: ", basename(out_path))
}


# After all p1 files are created, delete original two-sided p files
for (key in common_keys) {
  p_path <- file.path(data_dir, paste0("FEMA_p_lists_", key))
  
  if (file.exists(p_path)) {
    file.remove(p_path)
    message("Deleted original: ", basename(p_path))
  }
}


# Rename p1 files back to FEMA_p_lists_*.RDS
for (key in common_keys) {
  old_path <- file.path(data_dir, paste0("FEMA_p1_lists_", key))
  new_path <- file.path(data_dir, paste0("FEMA_p_lists_", key))
  
  if (file.exists(old_path)) {
    file.rename(from = old_path, to = new_path)
    message("Renamed: ", basename(old_path), " -> ", basename(new_path))
  }
}


################################
# Convert two-sided FEMA2 p values into one-sided (this is the code for fixed effect)
# Folder containing the FEMA beta and p-value RDS files
data_dir <- "MABFoutcomes/500500/MA/FEMA2/"


# Convert two-sided p-values to one-sided p-values for H0: theta = 0 vs H1: theta > 0
convert_pvec_to_onesided <- function(beta_vec, p2_vec) {
  if (length(beta_vec) != length(p2_vec)) {
    stop("beta_vec and p2_vec must have the same length.")
  }
  
  ifelse(
    beta_vec > 0,
    p2_vec / 2,
    ifelse(beta_vec < 0, 1 - p2_vec / 2, 0.5)
  )
}


# Convert one paired beta-list object and p-list object
convert_nested_p_list <- function(beta_obj, p_obj) {
  if (!identical(names(beta_obj), names(p_obj))) {
    stop("Top-level names do not match between beta and p objects.")
  }
  
  out <- vector("list", length(beta_obj))
  names(out) <- names(beta_obj)
  
  for (cond_name in names(beta_obj)) {
    beta_cond <- beta_obj[[cond_name]]
    p_cond    <- p_obj[[cond_name]]
    
    if (!identical(names(beta_cond), names(p_cond))) {
      stop(sprintf("Inner names do not match for condition: %s", cond_name))
    }
    
    out[[cond_name]] <- lapply(names(beta_cond), function(x) {
      convert_pvec_to_onesided(beta_cond[[x]], p_cond[[x]])
    })
    names(out[[cond_name]]) <- names(beta_cond)
  }
  
  out
}


# Find all beta and two-sided p-value RDS files
beta_files <- list.files(
  data_dir,
  pattern = "^FEMA2_beta_lists_.*\\.RDS$",
  full.names = TRUE
)

p_files <- list.files(
  data_dir,
  pattern = "^FEMA2_p_lists_.*\\.RDS$",
  full.names = TRUE
)


# Match files by shared suffix
beta_keys <- sub("^FEMA2_beta_lists_", "", basename(beta_files))
p_keys    <- sub("^FEMA2_p_lists_", "", basename(p_files))

common_keys <- intersect(beta_keys, p_keys)


# Generate one-sided p-value RDS files
for (key in common_keys) {
  beta_path <- file.path(data_dir, paste0("FEMA2_beta_lists_", key))
  p_path    <- file.path(data_dir, paste0("FEMA2_p_lists_", key))
  
  beta_obj <- readRDS(beta_path)
  p_obj    <- readRDS(p_path)
  
  p1_obj <- convert_nested_p_list(beta_obj, p_obj)
  
  out_path <- file.path(data_dir, paste0("FEMA2_p1_lists_", key))
  saveRDS(p1_obj, out_path)
  
  message("Saved: ", basename(out_path))
}


# After all p1 files are created, delete original two-sided p files
for (key in common_keys) {
  p_path <- file.path(data_dir, paste0("FEMA2_p_lists_", key))
  
  if (file.exists(p_path)) {
    file.remove(p_path)
    message("Deleted original: ", basename(p_path))
  }
}


# Rename p1 files back to FEMA_p_lists_*.RDS
for (key in common_keys) {
  old_path <- file.path(data_dir, paste0("FEMA2_p1_lists_", key))
  new_path <- file.path(data_dir, paste0("FEMA2_p_lists_", key))
  
  if (file.exists(old_path)) {
    file.rename(from = old_path, to = new_path)
    message("Renamed: ", basename(old_path), " -> ", basename(new_path))
  }
}
