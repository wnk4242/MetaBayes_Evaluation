# This script creates pie plots for MABF methods and FE meta-analysis (PieChart1 does not include FEMA)
# However, the pie charts created by this script does not include anecdotal data
# We use pcutoff_o=0.05 for the orignal study, c0 for MABF (BFcutoff = 1) and c1 (pcutoff_r = 0.05) for FEMA
library(ggplot2)
library(dplyr)
library(tidyr)
source("./helper functions_XXR.R")
#######Prepare data for visualization#######
# Import data sets
# Set the path to the directory containing RDS files
folder_path <- "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.05, EScutoff_o=0/c0/" #These dataset contain original study p and ES for calculating TSR, FSR, etc.
# List all RDS files in the directory
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files) {
    # Extract the base name without the extension
    file_name <- tools::file_path_sans_ext(basename(file_path))
    # Create a variable with the name of the file and assign the data from readRDS
    assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}
remove(rates_FEMA_0.2null_c0);remove(rates_FEMA_0.5null_c0);
# Import FEMA rates when pcutoff_r = 0.05 is used (c1)
rates_FEMA_0.2null_c1 <- readRDS("./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.05, EScutoff_o=0/c1/rates_FEMA_0.2null_c1.RDS")
rates_FEMA_0.5null_c1 <- readRDS("./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.05, EScutoff_o=0/c1/rates_FEMA_0.5null_c1.RDS")

# Add a method column containing the method name in each dataset as identifier
rates_BFbMA_0.2null_c0 <- rates_BFbMA_0.2null_c0 %>% 
    mutate(method = rep("BFbMA", 81)) %>% 
    relocate(method, .before = everything()) %>% 
    select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_EUBF_0.2null_c0 <- rates_EUBF_0.2null_c0 %>% 
    mutate(method = rep("EUBF", 81)) %>% 
    relocate(method, .before = everything()) %>% 
    select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_FEMABF_0.2null_c0 <- rates_FEMABF_0.2null_c0 %>% 
    mutate(method = rep("FEMABF", 81)) %>% 
    relocate(method, .before = everything()) %>% 
    select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_iBF_0.2null_c0 <- rates_iBF_0.2null_c0 %>% 
    mutate(method = rep("iBF", 81)) %>% 
    relocate(method, .before = everything()) %>% 
    select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_FEMA_0.2null_c1 <- rates_FEMA_0.2null_c1 %>% 
    mutate(method = rep("FEMA", 81)) %>% 
    relocate(method, .before = everything())
rates_BFbMA_0.5null_c0 <- rates_BFbMA_0.5null_c0 %>% 
    mutate(method = rep("BFbMA", 81)) %>% 
    relocate(method, .before = everything()) %>% 
    select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_EUBF_0.5null_c0 <- rates_EUBF_0.5null_c0 %>% 
    mutate(method = rep("EUBF", 81)) %>% 
    relocate(method, .before = everything()) %>% 
    select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_FEMABF_0.5null_c0 <- rates_FEMABF_0.5null_c0 %>% 
    mutate(method = rep("FEMABF", 81)) %>% 
    relocate(method, .before = everything()) %>% 
    select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_iBF_0.5null_c0 <- rates_iBF_0.5null_c0 %>% 
    mutate(method = rep("iBF", 81)) %>% 
    relocate(method, .before = everything()) %>% 
    select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_FEMA_0.5null_c1 <- rates_FEMA_0.5null_c1 %>% 
    mutate(method = rep("FEMA", 81)) %>% 
    relocate(method, .before = everything())

# Combine MABF and FEMA datasets
rates_MABFEMA_0.2null <- rbind(rates_FEMA_0.2null_c1,rates_BFbMA_0.2null_c0,rates_EUBF_0.2null_c0,rates_FEMABF_0.2null_c0,rates_iBF_0.2null_c0)
rates_MABFEMA_0.5null <- rbind(rates_FEMA_0.5null_c1,rates_BFbMA_0.5null_c0,rates_EUBF_0.5null_c0,rates_FEMABF_0.5null_c0,rates_iBF_0.5null_c0)



#######Underlying effect is null#######
# Prepare the data
plot_data_null <- rates_MABFEMA_0.2null %>%
    mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
    mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
    mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
    mutate(bias.level = factor(censorFunc, levels = c('low', 'medium', 'high'))) %>% 
    relocate(bias.level, .after = orig.n) %>% 
    select(method, rep.n, rep.number, orig.n, bias.level, TS2, FS2, FF2, TF2) %>%
    gather(key = "metric", value = "value", TS2, FS2, FF2, TF2) %>%
    mutate(metric = factor(metric, levels = c("TS2", "FS2", "FF2", "TF2"))) %>%
    mutate(metric_label = factor(metric, 
                                 levels = c("TS2", "FS2", "FF2", "TF2"),
                                 labels = c("true success", "false success", "false failure",  "true failure"))) %>%
    mutate(method = factor(method, levels = c("FEMA","BFbMA","EUBF","FEMABF","iBF"), labels = c("FEMA","BFbMA","EUBF","FEMABF","iBF"))) %>% 
    group_by(method, rep.n, rep.number) %>%
    mutate(proportion = value / sum(value))

# Create custom gradient colors
#null_colors <- c("#0000CD", "red", "pink", "#90EE90")
null_colors <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2")
# Create a named vector for scale_fill_manual
color_mapping <- c("true success" = null_colors[1],
                   "false success" = null_colors[2],
                   "false failure" = null_colors[3],
                   "true failure" = null_colors[4])

# Custom labeller for metrics
metric_labeller <- c(TS2 = "true success",
                     FS2 = "false success",
                     FF2 = "false failure",
                     TF2 = "true failure")

# Create the plot
plot <- ggplot(plot_data_null, aes(x = "", y = proportion, fill = metric_label)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = color_mapping) +  # Modern high contrast colors
    labs(x = "", y = "Proportion", fill = "Evidence Category", title = "Proportion of Meta-Analytic Bayes Factor and Fixed-Effects Effect Sizes Categorized into Success and Failure When the Underlying Effect is Null") +
    facet_grid(method ~ rep.n + rep.number, scales = "free", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "bottom")
print(plot)

# Save the plot
ggsave(filename = "./plots/pie chart/piechart_null_TSTF,pcutoff_o=0.05,BFcutoff=1,pcutoff_r=0.05.jpg", plot = plot, width = 15, height = 10, units = "in", dpi = 300)




#######Underlying effect is 0.2#######
# Prepare the data
plot_data_0.2 <-  rates_MABFEMA_0.2null %>%
    mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
    mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
    mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
    mutate(bias.level = factor(censorFunc, levels = c('low', 'medium', 'high'))) %>% 
    relocate(bias.level, .after = orig.n) %>% 
    select(method, rep.n, rep.number, orig.n, bias.level, TS1, FS1, FF1, TF1) %>%
    gather(key = "metric", value = "value", TS1, FS1, FF1, TF1) %>%
    mutate(metric = factor(metric, levels = c("TS1", "FS1", "FF1", "TF1"))) %>%
    mutate(metric_label = factor(metric, 
                                 levels = c("TS1", "FS1", "FF1", "TF1"),
                                 labels = c("true success", "false success", "false failure", "true failure"))) %>%
    mutate(method = factor(method, levels = c("FEMA","BFbMA","EUBF","FEMABF","iBF"), labels = c("FEMA","BFbMA","EUBF","FEMABF","iBF"))) %>% 
    group_by(method, rep.n, rep.number) %>%
    mutate(proportion = value / sum(value))

# Create custom gradient colors
#true_colors <- c("#0000CD", "red", "pink", "#90EE90")
true_colors <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2")
# Create a named vector for scale_fill_manual
color_mapping <- c("true success" = true_colors[1],
                   "false success" = true_colors[2],
                   "false failure" = true_colors[3],
                   "true failure" = true_colors[4])

# Custom labeller for metrics
metric_labeller <- c(
    TS1 = "true success",
    FS1 = "false success",
    FF1 = "false failure",
    TF1 = "true failure"
)

# Create the plot
plot <- ggplot(plot_data_0.2, aes(x = "", y = proportion, fill = metric_label)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = color_mapping) +  # Modern high contrast colors
    labs(x = "", y = "Proportion", fill = "Evidence Category", title = "Proportion of Meta-Analytic Bayes Factor and Fixed-Effects Effect Sizes Categorized into Success and Failure When the Underlying Effect is 0.2") +
    facet_grid(method ~ rep.n + rep.number, scales = "free", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "bottom")
print(plot)

# Save the plot
ggsave(filename = "./plots/pie chart/piechart_0.2_TSTF,pcutoff_o=0.05,BFcutoff=1,pcutoff_r=0.05.jpg", plot = plot, width = 15, height = 10, units = "in", dpi = 300)




#######Underlying effect is 0.5#######
# Prepare the data
plot_data_0.5 <-  rates_MABFEMA_0.5null %>%
    mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
    mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
    mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
    mutate(bias.level = factor(censorFunc, levels = c('low', 'medium', 'high'))) %>% 
    relocate(bias.level, .after = orig.n) %>% 
    select(method, rep.n, rep.number, orig.n, bias.level, TS1, FS1, FF1, TF1) %>%
    gather(key = "metric", value = "value", TS1, FS1, FF1, TF1) %>%
    mutate(metric = factor(metric, levels = c("TS1", "FS1", "FF1", "TF1"))) %>%
    mutate(metric_label = factor(metric, 
                                 levels = c("TS1", "FS1", "FF1", "TF1"),
                                 labels = c("true success", "false success", "false failure", "true failure"))) %>%
    mutate(method = factor(method, levels = c("FEMA","BFbMA","EUBF","FEMABF","iBF"), labels = c("FEMA","BFbMA","EUBF","FEMABF","iBF"))) %>% 
    group_by(method, rep.n, rep.number) %>%
    mutate(proportion = value / sum(value))

# Create custom gradient colors
#true_colors <- c("#0000CD", "red", "pink", "#90EE90")
true_colors <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2")
# Create a named vector for scale_fill_manual
color_mapping <- c("true success" = true_colors[1],
                   "false success" = true_colors[2],
                   "false failure" = true_colors[3],
                   "true failure" = true_colors[4])

# Custom labeller for metrics
metric_labeller <- c(
    TS1 = "true success",
    FS1 = "false success",
    FF1 = "false failure",
    TF1 = "true failure"
)

# Create the plot
plot <- ggplot(plot_data_0.5, aes(x = "", y = proportion, fill = metric_label)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = color_mapping) +  # Modern high contrast colors
    labs(x = "", y = "Proportion", fill = "Evidence Category", title = "Proportion of Meta-Analytic Bayes Factor and Fixed-Effects Effect Sizes Categorized into Success and Failure When the Underlying Effect is 0.5") +
    facet_grid(method ~ rep.n + rep.number, scales = "free", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "bottom")
print(plot)

# Save the plot
ggsave(filename = "./plots/pie chart/piechart_0.5_TSTF,pcutoff_o=0.05,BFcutoff=1,pcutoff_r=0.05.jpg", plot = plot, width = 15, height = 10, units = "in", dpi = 300)
