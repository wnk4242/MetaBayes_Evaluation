library(tidyverse)
library(httr)
library(jsonlite)
library(glmmTMB)
library(lmtest)
library(emmeans)
library(car)
library(ggplot2)
library(ggeffects)
#####
# # Step 0 Add new columns: orig.alpha, bias.level (done, no need to run again)
# folder_path <- "./MABFanalyses/row-wise/orig.alpha=0.01/for step 0 (raw)"
# output_path <- file.path(folder_path, "updated")
# # Create output folder if it doesn't exist
# if (!dir.exists(output_path)) {
#   dir.create(output_path, recursive = TRUE)
# }
# rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# # Add and relocate new columns: orig.alpha, bias.level
# for (file_path in rds_files) {
#   file_name <- tools::file_path_sans_ext(basename(file_path))
#   # Read and modify
#   df <- readRDS(file_path) %>%
#     mutate(
#       orig.alpha = 0.01,
#       bias.level = PB.level
#     ) %>%
#     relocate(orig.alpha, .after = orig.n) %>%
#     relocate(bias.level, .before = QRP.level)
#   # Assign to global environment
#   assign(file_name, df, envir = .GlobalEnv)
#   # Save to updated folder
#   saveRDS(df, file.path(output_path, paste0(file_name, ".RDS")))
# }
# # Step 0 Add new columns: orig.alpha, bias.level (done, no need to run again)
# folder_path <- "./MABFanalyses/row-wise/orig.alpha=0.05/for step 0 (raw)"
# output_path <- file.path(folder_path, "updated")
# # Create output folder if it doesn't exist
# if (!dir.exists(output_path)) {
#   dir.create(output_path, recursive = TRUE)
# }
# rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# # Add and relocate new columns: orig.alpha, bias.level
# for (file_path in rds_files) {
#   file_name <- tools::file_path_sans_ext(basename(file_path))
#   # Read and modify
#   df <- readRDS(file_path) %>%
#     mutate(
#       orig.alpha = 0.05,
#       bias.level = PB.level
#     ) %>%
#     relocate(orig.alpha, .after = orig.n) %>%
#     relocate(bias.level, .before = QRP.level)
#   # Assign to global environment
#   assign(file_name, df, envir = .GlobalEnv)
#   # Save to updated folder
#   saveRDS(df, file.path(output_path, paste0(file_name, ".RDS")))
# }
# # Step 00: Combine corresponding datasets with two alpha levels
# # Set the path to the directory containing RDS files
# library(dplyr)
# library(purrr)
# library(stringr)
# folder_path <- "./MABFanalyses/row-wise/ready to combine"
# save_path <- file.path(folder_path, "combined")  # New folder to save combined files
# # Create new folder if it doesn't exist
# if (!dir.exists(save_path)) dir.create(save_path)
# # List all RDS files in the directory
# rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# # Extract filename only (for parsing) and metadata
# meta_info <- tibble(
#   file = rds_files,
#   filename = basename(rds_files)
# ) %>%
#   mutate(
#     effect = str_extract(filename, "(?<=rates_)[^_]+"),
#     BFcutoff = str_extract(filename, "BFcutoff[0-9]+")
#   )
# # Group by effect and BFcutoff and combine rows of pcutoff 0.01 and 0.05
# combined_named <- meta_info %>%
#   group_by(effect, BFcutoff) %>%
#   group_map(~ {
#     dfs <- map(.x$file, readRDS)
#     bind_rows(dfs)
#   }) %>%
#   set_names(
#     meta_info %>%
#       distinct(effect, BFcutoff) %>%
#       transmute(name = paste0("MABF_rates_", effect, "_rowise_", BFcutoff)) %>%
#       pull(name)
#   )
# # Save each combined dataset
# walk2(combined_named, names(combined_named), ~ saveRDS(.x, file.path(save_path, paste0(.y, ".RDS"))))
#####
# # Check levels in independent variables
MABF_rates_0.2_rowise_BFcutoff1 %>%
  summarise(across(c(1,3,4,5,6,9,10), ~ list(unique(.x)))) %>%
  pivot_longer(everything(), names_to = "column", values_to = "levels") %>%
  mutate(levels = map_chr(levels, ~ paste(.x, collapse = ", ")))

# # Check unique combinations across independent variables
df <-  MABF_rates_0.2_rowise_BFcutoff1 %>%
  distinct(method, true.effect, orig.n, orig.alpha, bias.level, rep.number, rep.n) %>%
  filter(method == "FEMABF") %>%
  print()

#####
# Step 1 Load datasets
# Set the path to the directory containing RDS files
folder_path <- "./MABFanalyses/row-wise" #These dataset contain original study p and ES for calculating TSR, FSR, etc.
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}




#####
# Step 2 Modeling
#####
# For TPR and FPR, rule out useless predictors: orig.n, orig.alpha, bias level
# To better quantify practical significance, compute effect sizes or marginal means:
# library(emmeans)
# model_tpr_full     <- glmmTMB(TPR ~ rep.n+rep.number+bias.level, data = data_combined, family = beta_family())
# emmeans(model_tpr_full, ~ bias.level) #increasing bias.level yield very small change in TPR although significant
# model_tpr_full     <- glmmTMB(TPR ~ rep.n+rep.number+orig.n, data = data_combined, family = beta_family())
# emmeans(model_tpr_full, ~ orig.n) #increasing orig.n yield very small change in TPR although significant
# model_tpr_full     <- glmmTMB(TPR ~ rep.n+rep.number+orig.alpha, data = data_combined, family = beta_family())
# emmeans(model_tpr_full, ~ orig.alpha) #increasing orig.alpha yield very small change in TPR although significant

#####
# Uses beta logistic regression in betareg library, this is slow
# analyze_MABF_effects <- function(.method, BFcutoff, true.effect) {
#   # Construct dataset names
#   tpr_rds <- sprintf("MABF_rates_%.1f_rowise_pcutoff0.05_BFcutoff%d", true.effect, BFcutoff)
#   fpr_rds <- sprintf("MABF_rates_null_rowise_pcutoff0.05_BFcutoff%d", BFcutoff)
#   # Load datasets
#   data_tpr <- get(tpr_rds, envir = .GlobalEnv)
#   data_fpr <- get(fpr_rds, envir = .GlobalEnv)
#   # Select columns
#   data_tpr <- data_tpr %>% select(1, 3, 4, 5, 6, 7, 8, TPR = 21)
#   data_fpr <- data_fpr %>% select(1, 3, 4, 5, 6, 7, 8, FPR = 22)
#   # Factorize
#   data_tpr <- data_tpr %>% mutate(across(c(method, true.effect, orig.n, QRP.level, PB.level, rep.number, rep.n), as.factor))
#   data_fpr <- data_fpr %>% mutate(across(c(method, true.effect, orig.n, QRP.level, PB.level, rep.number, rep.n), as.factor))
#   # Filter
#   data_tpr <- data_tpr %>% filter(method == .method, orig.n == "20", QRP.level == "none")
#   data_fpr <- data_fpr %>% filter(method == .method, orig.n == "20", QRP.level == "none")
#   # Add scenario_id
#   data_tpr <- data_tpr %>% mutate(scenario_id = rep(1:9, each = 500)) %>% relocate(scenario_id, .before = 1)
#   data_fpr <- data_fpr %>% mutate(scenario_id = rep(1:9, each = 500)) %>% relocate(scenario_id, .before = 1)
#   # Set factor levels
#   data_tpr$rep.n <- factor(data_fpr$rep.n, levels = c("40", "100", "400"))
#   data_tpr$rep.number <- factor(data_fpr$rep.number, levels = c("2", "5", "10"))
#   data_fpr$rep.n <- factor(data_fpr$rep.n, levels = c("40", "100", "400"))
#   data_fpr$rep.number <- factor(data_fpr$rep.number, levels = c("2", "5", "10"))
#   # Merge
#   data_combined <- data_tpr %>% mutate(FPR = data_fpr$FPR)
#   # Beta regression
#   model_tpr_full <- betareg(TPR ~ rep.n * rep.number, data = data_combined)
#   model_fpr_full <- betareg(FPR ~ rep.n * rep.number, data = data_combined)
#   # Reduced models for TPR
#   model_tpr_dropn <- betareg(TPR ~ rep.number, data = data_combined)
#   model_tpr_dropnum <- betareg(TPR ~ rep.n, data = data_combined)
#   # Reduced models for FPR
#   model_fpr_dropn <- betareg(FPR ~ rep.number, data = data_combined)
#   model_fpr_dropnum <- betareg(FPR ~ rep.n, data = data_combined)
#   # Likelihood ratio tests
#   lr_tpr_dropn <- lrtest(model_tpr_full, model_tpr_dropn)
#   lr_tpr_dropnum <- lrtest(model_tpr_full, model_tpr_dropnum)
#   lr_fpr_dropn <- lrtest(model_fpr_full, model_fpr_dropn)
#   lr_fpr_dropnum <- lrtest(model_fpr_full, model_fpr_dropnum)
#   # Return results
#   list(
#     Beta_TPR = summary(model_tpr_full),
#     Beta_FPR = summary(model_fpr_full),
#     LR_Tests = list(
#       TPR = list(drop_rep.n = lr_tpr_dropn, drop_rep.number = lr_tpr_dropnum),
#       FPR = list(drop_rep.n = lr_fpr_dropn, drop_rep.number = lr_fpr_dropnum)
#     )
#   )
# }

#####
# Use glmmTMB, faster
analyze_MABF_effects <- function(.method, BFcutoff, true.effect, dv, ivs, interactions = NULL) {
  
  
  # Dataset names
  tpr_rds <- sprintf("MABF_rates_%.1f_rowise_BFcutoff%d", true.effect, BFcutoff)
  fpr_rds <- sprintf("MABF_rates_null_rowise_BFcutoff%d", BFcutoff)
  
  # Load
  data_tpr <- get(tpr_rds, envir = .GlobalEnv)
  data_fpr <- get(fpr_rds, envir = .GlobalEnv)
  
  # Select columns
  data_tpr <- data_tpr %>% select(method = 1, true.effect = 3, orig.n = 4, orig.alpha = 5,
                                  bias.level = 6, rep.number = 9, rep.n = 10,
                                  ADRt = 11, TPR = 23, TSR1 = 29)
  data_fpr <- data_fpr %>% select(method = 1, true.effect = 3, orig.n = 4, orig.alpha = 5,
                                  bias.level = 6, rep.number = 9, rep.n = 10,
                                  ADRn = 11, FPR = 24, TSR2 = 29, FSR2 = 30)
  
  # Factorize
  to_factor <- c("method", "true.effect", "orig.n", "orig.alpha", "bias.level", "rep.number", "rep.n")
  data_tpr <- data_tpr %>% mutate(across(all_of(to_factor), as.factor))
  data_fpr <- data_fpr %>% mutate(across(all_of(to_factor), as.factor))
  
  # Filter
  data_tpr <- data_tpr %>% filter(method == .method)
  data_fpr <- data_fpr %>% filter(method == .method)
  
  # Scenario ID
  data_tpr <- data_tpr %>% mutate(scenario_id = rep(1:162, each = 500)) %>% relocate(scenario_id)
  data_fpr <- data_fpr %>% mutate(scenario_id = rep(1:162, each = 500)) %>% relocate(scenario_id)
  
  # Factor levels
  factor_levels <- list(rep.n = c("40", "100", "400"),
                        rep.number = c("2", "5", "10"),
                        orig.n = c("20", "50", "200"),
                        orig.alpha = c("0.01", "0.05"),
                        bias.level = c("low", "medium", "high"))
  for (var in names(factor_levels)) {
    data_tpr[[var]] <- factor(data_tpr[[var]], levels = factor_levels[[var]])
    data_fpr[[var]] <- factor(data_fpr[[var]], levels = factor_levels[[var]])
  }
  
  # Merge
  data_combined <- data_tpr %>%
    mutate(FPR = data_fpr$FPR,
           ADRn = data_fpr$ADRn,
           FSR2 = data_fpr$FSR2)
  
  # Shrink bounded values
  epsilon <- 1e-6
  data_combined <- data_combined %>%
    mutate(across(c(TPR, FPR, ADRt, ADRn, TSR1, TSR2, FSR2),
                  ~ ifelse(. <= 0, epsilon, ifelse(. >= 1, 1 - epsilon, .))))
  
  # Build formula
  main_effects_formula <- paste(ivs, collapse = " + ")
  if (!is.null(interactions) && length(interactions) > 0) {
    interaction_formula <- paste(interactions, collapse = " + ")
    full_rhs <- paste(main_effects_formula, interaction_formula, sep = " + ")
  } else {
    full_rhs <- main_effects_formula
  }
  full_formula <- reformulate(full_rhs, response = dv)
  
  # Fit full model
  full_model <- glmmTMB(full_formula, data = data_combined, family = beta_family(link = "logit"))
  
  # drop1 table
  drop1_table <- drop1(full_model, test = "Chisq")
  
  # emmeans results
  all_terms <- c(ivs, interactions)
  emmeans_results <- list()
  for (term in all_terms) {
    # Skip invalid interactions (e.g., a typo)
    tryCatch({
      emmeans_results[[term]] <- emmeans(full_model, specs = as.formula(paste("~", term)))
    }, error = function(e) {
      emmeans_results[[term]] <- paste("Error in emmeans for term:", term)
    })
  }
  
  # Return
  list(
    Full_Model = summary(full_model),
    Drop1_Table = drop1_table,
    EMMeans = emmeans_results,
    Wald_type3 = Anova(full_model, type = 3)
    # , data_used = data_combined
  )
}


#####
# Step 3 Reporting
# Save model results
results <- analyze_MABF_effects(
  .method = "FEMABF",
  BFcutoff = 1,
  true.effect = 0.2,
  dv = "TPR",
  ivs = c("rep.n", "rep.number", "bias.level"),
  interactions = c("rep.n:bias.level", "rep.number:bias.level")
)

# Convert model outputs to a clean text block
output_text <- paste(
  "Wald teset (Type 3):\n",
  paste(capture.output(print(results$Wald_Type3)), collapse = "\n"), "\n\n",
  "===============================\n\n",
  "Beta Regression (Full Model):\n",
  paste(capture.output(print(results$Full_Model)), collapse = "\n"), "\n\n",
  "===============================\n\n",
  "Drop1 Likelihood Ratio Tests:\n",
  paste(capture.output(print(results$Drop1_Table)), collapse = "\n"), "\n\n",
  "===============================\n\n",
  "Estimated Marginal Means:\n",
  paste(
    lapply(names(results$EMMeans), function(term) {
      paste0("EMMeans for ", term, ":\n",
             paste(capture.output(print(results$EMMeans[[term]])), collapse = "\n"))
    }),
    collapse = "\n\n"
  )
)

# Print the final output
cat(output_text)




#####
# Step 4 ChatGPT interpretation
# Prepare API request
# gpt-4 or gpt-3.5-turbo both work
api_key <- readLines("openai_key.txt")[1]
response <- POST(
  url = "https://api.openai.com/v1/chat/completions",
  add_headers(
    Authorization = paste("Bearer", api_key),
    `Content-Type` = "application/json"
  ),
  body = toJSON(list(
    model = "gpt-4",
    messages = list(
      list(role = "user", content = output_text),
      list(role = "user", content = 
            "Please interpret these statistical results in APA style.
           Convert the logit to proportion in Estimated Marginal Means.")
    )
  ), auto_unbox = TRUE)
)

# Save response
raw_response <- content(response, as = "text", encoding = "UTF-8")

# Try parsing if valid
parsed <- fromJSON(raw_response, simplifyVector = FALSE)

# Print ChatGPT interpretation nicely
if (!is.null(parsed$choices)) {
  cat(parsed$choices[[1]]$message$content)
}










######
#Manual test
# Set the path to the directory containing RDS files
folder_path <- "./MABFanalyses/row-wise" #These dataset contain original study p and ES for calculating TSR, FSR, etc.
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}


# Construct dataset names
tpr_rds <- MABF_rates_0.2_rowise_BFcutoff1
fpr_rds <- MABF_rates_null_rowise_BFcutoff1
# Load datasets
data_tpr <- tpr_rds
data_fpr <- fpr_rds 
# Select relevant columns #HEREHERE #Change
data_tpr <- data_tpr %>% select(method = 1, true.effect = 3, orig.n = 4, orig.alpha = 5, bias.level = 6, rep.number = 9, rep.n = 10, ADRt = 11, TPR = 23, TSR1 = 29) #here ADRt = ADR_true
data_fpr <- data_fpr %>% select(method = 1, true.effect = 3, orig.n = 4, orig.alpha = 5, bias.level = 6, rep.number = 9, rep.n = 10, ADRn = 11, FPR = 24, TSR2 = 29, FSR2 = 30) #here ADRn = ADR_null
# Factorize variables #HEREHERE #Change
data_tpr <- data_tpr %>%
  mutate(across(c(method, true.effect, orig.n, orig.alpha, bias.level, rep.number, rep.n), as.factor))
data_fpr <- data_fpr %>%
  mutate(across(c(method, true.effect, orig.n, orig.alpha, bias.level, rep.number, rep.n), as.factor))
# Filter to target method/scenario #HEREHERE #Change
data_tpr <- data_tpr %>% filter(method == "FEMABF")
data_fpr <- data_fpr %>% filter(method == "FEMABF")
# Add scenario_id #Change
data_tpr <- data_tpr %>% mutate(scenario_id = rep(1:162, each = 500)) %>% relocate(scenario_id, .before = 1)
data_fpr <- data_fpr %>% mutate(scenario_id = rep(1:162, each = 500)) %>% relocate(scenario_id, .before = 1)
# Set factor levels explicitly #Change
data_tpr$rep.n <- factor(data_fpr$rep.n, levels = c("40", "100", "400"))
data_tpr$rep.number <- factor(data_fpr$rep.number, levels = c("2", "5", "10"))
data_tpr$orig.n <- factor(data_tpr$orig.n, levels = c("20", "50", "200"))
data_tpr$orig.alpha <- factor(data_tpr$orig.alpha, levels = c("0.01", "0.05"))
data_tpr$bias.level <- factor(data_tpr$bias.level, levels = c("low", "medium", "high"))

data_fpr$rep.n <- factor(data_fpr$rep.n, levels = c("40", "100", "400"))
data_fpr$rep.number <- factor(data_fpr$rep.number, levels = c("2", "5", "10"))
data_fpr$orig.n <- factor(data_fpr$orig.n, levels = c("20", "50", "200"))
data_fpr$orig.alpha <- factor(data_fpr$orig.alpha, levels = c("0.01", "0.05"))
data_fpr$bias.level <- factor(data_fpr$bias.level, levels = c("low", "medium", "high"))
# Merge TPR and FPR #Change
data_combined <- data_tpr %>% mutate(FPR = data_fpr$FPR, ADRn = data_fpr$ADRn, TSR2 = data_fpr$TSR2, FSR2 = data_fpr$FSR2) #here
# Apply minimal epsilon shrinkage to avoid exact 0 or 1
epsilon <- 1e-6
data_combined <- data_combined %>%
  mutate(
    TPR = ifelse(TPR <= 0, epsilon, ifelse(TPR >= 1, 1 - epsilon, TPR)),
    FPR = ifelse(FPR <= 0, epsilon, ifelse(FPR >= 1, 1 - epsilon, FPR)),
    ADRt = ifelse(ADRt<= 0, epsilon, ifelse(ADRt >= 1, 1 - epsilon, ADRt)), #here
    ADRn = ifelse(ADRn <= 0, epsilon, ifelse(ADRn >= 1, 1 - epsilon, ADRn)), #here
    TSR1 = ifelse(TSR1 <= 0, epsilon, ifelse(TSR1 >= 1, 1 - epsilon, TSR1)), #Change
    TSR2 = ifelse(TSR2 <= 0, epsilon, ifelse(TSR2 >= 1, 1 - epsilon, TSR2)),
    FSR2 = ifelse(FSR2 <= 0, epsilon, ifelse(FSR2 >= 1, 1 - epsilon, FSR2)) #Change
  )


# Fit models
model_tpr_full     <- glmmTMB(TPR ~ rep.n+rep.number+bias.level, data = data_combined, family = beta_family())
summary(model_tpr_full)







# Visualization 1

# Fit full model
model_tsr1_full     <- glmmTMB(TSR1 ~ rep.n+rep.number+bias.level+orig.n+bias.level*orig.n, data = data_combined, family = beta_family())
# Get predictions
preds <- ggpredict(model_tsr1_full, terms = c("rep.n", "rep.number", "orig.n", "bias.level"))

# Define custom labels for facet
bias_labels <- c("low" = "Bias = Low", "medium" = "Bias = Medium", "high" = "Bias = High")
orig_labels <- c("20" = "Orig.n = 20", "50" = "Orig.n = 50", "200" = "Orig.n = 200")

# Create the plot
ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_point(size = 2) +
  geom_line(aes(group = group), linewidth = 1) +
  facet_grid(facet ~ panel, labeller = labeller(facet = bias_labels, panel = orig_labels)) +
  labs(
    x = "Replication Sample Size (rep.n)",
    y = "Predicted TSR1",
    color = "Number of Replications",
    title = "Predicted TSR1 by Replication Design, Bias, and Original N"
  ) +
  theme_minimal()


# Visualization 2
# Fit full model
model_tpr_full     <- glmmTMB(TPR ~ rep.n+rep.number+rep.n*rep.number, data = data_combined, family = beta_family())
# Get predictions
preds <- ggpredict(model_tpr_full, terms = c("rep.n", "rep.number"))



# Create the plot
ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_point(size = 2) +
  geom_line(aes(group = group), linewidth = 1) +
  labs(
    x = "Replication Sample Size (rep.n)",
    y = "Predicted TPR",
    color = "Number of Replications",
    title = "Predicted TPR by Replication Design"
  ) +
  theme_minimal()

