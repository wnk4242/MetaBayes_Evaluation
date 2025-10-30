# Load necessary packages
library(ggplot2)
library(tidyverse)
library(purrr)
rm(list = ls()) # housekeeping

# Set the path to the directory containing RDS files
folder_path <- "./MABFanalyses/matrix-wise/rates4RepCurve/"


# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}


#Rep curve function
plot_ROC <- function(scenario_number) {
  condition <- scenario_df %>% filter(scenario == scenario_number)
  # Ensure we extract single values for rep_number and rep_n
  scenario_number <- condition$scenario[1]
  true_es <- condition$true_es[1]
  orig_n <- condition$orig.n[1]
  rep_number <- condition$rep.number[1]
  rep_n <- condition$rep.n[1]
  bias <- condition$bias[1]
  cutoff_p <- condition$threshold_p[1]
  
  # Debug info
  print(paste("Scenario:", scenario_number,
              "Original Sample Size", orig_n,
              "Replications:", rep_number,
              "Sample Size:", rep_n,
              "Bias Level:", bias,
              "Original p Cutoff:", cutoff_p))
  
  title_text <- paste0(
    'Rep Curves when the number of replications is ', rep_number,
    ', replication sample size is ', rep_n,
    ', true effect is ', true_es,
    ', and bias level is ', bias,
    ', and original p-cutoff is ', cutoff_p)
  
  roc_plot <- rates_0.2null %>%
    filter(scenario == as.character(scenario_number)) %>%
    group_by(method) %>%
    ggplot(aes(x = Xrate, y = Yrate, color = method)) +
    geom_step() +
    labs(title = title_text) +
    theme_minimal()
  
  # Save the static plot
  file_name <- paste0("./plots/rep curve/0.2null/",
                      scenario_number, "_", orig_n, "_",
                      rep_number, "_", rep_n, "_", bias, "_", cutoff_p, ".jpg")
  ggsave(file_name, plot = roc_plot, width = 12, height = 8)
  
  # Print static ggplot
  #print(roc_plot)
  
  # ---- NEW interactive zoomable version ----
  p_int <- plotly::ggplotly(roc_plot, tooltip = c("method", "Xrate", "Yrate"))
  print(p_int)
}


# simulation scenario
scenario_df <- data.frame(
  scenario = c(1:243),
  true_es = 0.2,
  orig.n = rep(c(20, 50, 200), each = 27, times = 3),
  rep.number = rep(c(2, 5, 10), each = 3, times = 27),
  rep.n = rep(c(40,100,400), each = 9, times = 9),
  bias = rep(c('none','medium','high'), each=81, times = 1),
  threshold_p = rep(c(0.01, 0.05, 0.1), each = 1, times = 81),
  threshold_ES = rep(c(0), each = 243, times = 1)
)
# 
# scenario_df <- data.frame(
#   scenario = c(1:81),
#   true_es = 0.2,
#   orig.n = rep(c(20, 50, 200), each = 9, times = 3),
#   rep.number = rep(c(2, 5, 10), each = 1, times = 9),
#   rep.n = rep(c(40,100,400), each = 3, times = 9),
#   bias = rep(c('none','medium','high'), each=27, times = 1),
#   threshold_p = rep(c(0.05), each = 81, times = 1),
#   threshold_ES = rep(c(0), each = 81, times = 1)
# )

# Combine datasets of BFMA methods
#rates_0.2null <- rbind(EUBF_TSFS4RepCurve_0.2null)
rates_0.2null <- rbind(BFbMA_TSFS4RepCurve_0.2null, EUBF_TSFS4RepCurve_0.2null, FEMABF_TSFS4RepCurve_0.2null, iBF_TSFS4RepCurve_0.2null)
# Change threshold_p_r to threshold
FEMA_TSFS4RepCurve_0.2null <- FEMA_TSFS4RepCurve_0.2null %>% 
  rename(threshold = threshold_p_r)
# Change threshold_BF to threshold
rates_0.2null <- rates_0.2null  %>% 
  rename(threshold = threshold_BF)
# Combine rows of MABF methods and FEMA
rates_0.2null <- rbind(FEMA_TSFS4RepCurve_0.2null, rates_0.2null)

#When creating TPRvsFPR, original p and original ES do not affect results so no need to include them.
rates_0.2null <- rates_0.2null %>%
  mutate(across(c(`rep number`,rep.n), as.numeric)) %>% 
  mutate(true_es = 0.2) %>% 
  relocate(true_es, .after = `true effect`) %>% 
  relocate(group, .before = `method`) %>% 
  rename(scenario = group) %>% 
  #filter(threshold_p==0.01, threshold_ES==0) %>% 
  mutate(Yrate = TS1/Successes1, Xrate = FS2/Successes2) %>% #TSR1 vs FSR2 #this works #the original ES doesn't affect results
  relocate(Xrate, .after = Yrate)


walk(c(163), plot_ROC)









