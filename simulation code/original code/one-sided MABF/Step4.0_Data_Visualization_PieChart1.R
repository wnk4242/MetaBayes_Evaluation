# This stacked bar plot shows the proportion of BF falling into success/failure, which explains why TSR is high (90%+)
# You can select to use stacked bar chart or pie chart. I recommend pie chart.
# This script only compares MABF methods
# The original p cutoff is 0.05
# You don't want to see many false successes or true failures
library(ggplot2)
library(dplyr)
library(tidyr)
source("./helper functions_XXR.R")
#######When underlying effect is null#########
# Stacked bar chart
# Prepare the data
# MABF_rates_null_rowise dataset is a wide format dataset. This is just for null original study.
# It contains 4 MABF methods * 500 reps * (3 orig samp size * 3 research envir) combinations (original study) * (3 rep samp size * 3 num of reps) combinations (replication study) = 162,000 observations 
# plot_data_null is a long format dataset.
# It contains 4 MABF methods * 500 reps * 9 combinations (original study) * 9 combinations (replication study) * 5 metrics (AD_null, TS2, FS2, FF2, TF2) = 810,000 observations
# AD_null is anecdotal evidence when underlying effect is null
MABF_rates_null_rowise_pcutoff0.05_BFcutoff3 <- readRDS("./MABFanalyses/row-wise/MABF_rates_null_rowise_pcutoff0.05_BFcutoff3.RDS")
plot_data_null <- MABF_rates_null_rowise_pcutoff0.05_BFcutoff3 %>%
    mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
    mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
    mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
    mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>%
    select(method, rep.n, rep.number, AD_null, TS2, FS2, FF2, TF2) %>%
    gather(key = "metric", value = "value", AD_null, TS2, FS2, FF2, TF2) %>%
    mutate(metric = factor(metric, levels = c("AD_null", "TS2", "FS2", "FF2",  "TF2"))) %>%
    mutate(metric_label = factor(metric, 
                                 levels = c("AD_null", "TS2", "FS2", "FF2", "TF2"),
                                 labels = c("anecdotal","true success", "false success", "false failure",  "true failure"))) %>%
    group_by(method, rep.n, rep.number) %>%
    mutate(proportion = value / sum(value))

# Create custom gradient colors for the AD_null metrics
#null_colors <- c("grey",  "#0000CD", "red", "pink", "#90EE90")
null_colors <- c("#E69F00", "#56B4E9", "coral2", "#F0E442", "darkorchid")
# Create a named vector for scale_fill_manual
color_mapping <- c(anecdotal = null_colors[1],
                   "true success" = null_colors[2],
                   "false success" = null_colors[3],
                   "false failure" = null_colors[4],
                   "true failure" = null_colors[5])

# Custom labeller for metrics
metric_labeller <- c(
    AD_null = "anecdotal",
    TS2 = "true success",
    FS2 = "false success",
    FF2 = "false failure",
    TF2 = "true failure"
)

# Create the plot
plot <- ggplot(plot_data_null, aes(x = factor(method), y = proportion, fill = metric_label)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color_mapping) +  # Custom colors
    labs(x = "MABF Method", y = "Proportion", fill = "Evidence Strength", title = "Proportion of Bayes Factors Categoried into Varying Levels of Evidence When the Underlying Effect is Null") +
    facet_grid(rep.n ~ rep.number, scales = "free_x", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")
print(plot)

############
############
# Pie chart
# Prepare the data
MABF_rates_null_rowise_pcutoff0.05_BFcutoff3 <- readRDS("./MABFanalyses/row-wise/MABF_rates_null_rowise_pcutoff0.05_BFcutoff3.RDS")
plot_data_null <- MABF_rates_null_rowise_pcutoff0.05_BFcutoff3 %>%
    mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
    mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
    mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
    mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>%
    select(method, rep.n, rep.number, AD_null, TS2, FS2, TF2, FF2) %>%
    gather(key = "metric", value = "value", AD_null, TS2, FS2, TF2, FF2) %>%
    mutate(metric = factor(metric, levels = c("AD_null", "TS2", "FS2", "TF2", "FF2"))) %>%
    mutate(metric_label = factor(metric, 
                                 levels = c("AD_null", "TS2", "FS2", "TF2", "FF2"),
                                 labels = c("anecdotal", "true success", "false success", "true failure", "false failure"))) %>%
    group_by(method, rep.n, rep.number) %>%
    mutate(proportion = value / sum(value))

# Create modern high contrast colors for the AD_null metrics
null_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")

# Create a named vector for scale_fill_manual
color_mapping <- c(anecdotal = null_colors[1],
                   "true success" = null_colors[2],
                   "false success" = null_colors[3],
                   "true failure" = null_colors[4],
                   "false failure" = null_colors[5])

# Custom labeller for metrics
metric_labeller <- c(
    AD_null = "anecdotal",
    TS2 = "true success",
    FS2 = "false success",
    TF2 = "true failure",
    FF2 = "false failure"
)

# Create the pie chart
plot <- ggplot(plot_data_null, aes(x = "", y = proportion, fill = metric_label)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = color_mapping) +  # Modern high contrast colors
    labs(x = "", y = "Proportion", fill = "Evidence Category", title = "Proportion of Bayes Factors Categorized into Success and Failure When the Underlying Effect is Null") +
    facet_grid(method ~ rep.n + rep.number, scales = "free", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "bottom")

# Save the plot
ggsave(filename = "./plots/pie chart/piechart_null_TSTF,pcutoff=0.05,BFcutoff=3.jpg", plot = plot, width = 15, height = 10, units = "in", dpi = 300)


#######When underlying effect is 0.2#########
# Stacked bar chart
# Prepare the data
MABF_rates_0.2_rowise_pcutoff0.05_BFcutoff3 <- readRDS("./MABFanalyses/row-wise/MABF_rates_0.2_rowise_pcutoff0.05_BFcutoff3.RDS")
plot_data_0.2 <- MABF_rates_0.2_rowise_pcutoff0.05_BFcutoff3 %>%
    mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
    mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
    mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
    mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>%
    select(method, rep.n, rep.number, AD_true, TS1, FS1, FF1, TF1) %>%
    gather(key = "metric", value = "value", AD_true, TS1, FS1, FF1, TF1) %>%
    mutate(metric = factor(metric, levels = c("AD_true", "TS1", "FS1", "FF1",  "TF1"))) %>%
    mutate(metric_label = factor(metric, 
                                 levels = c("AD_true", "TS1", "FS1", "FF1", "TF1"),
                                 labels = c("anecdotal","true success", "false success", "false failure",  "true failure"))) %>%
    group_by(method, rep.n, rep.number) %>%
    mutate(proportion = value / sum(value))

# Create custom gradient colors for the AD_true metrics
true_colors <- c("#E69F00", "#56B4E9", "coral2", "#F0E442", "darkorchid")
# Create a named vector for scale_fill_manual
color_mapping <- c(anecdotal = true_colors[1],
                   "true success" = true_colors[2],
                   "false success" = true_colors[3],
                   "false failure" = true_colors[4],
                   "true failure" = true_colors[5])

# Custom labeller for metrics
metric_labeller <- c(
    AD_true = "anecdotal",
    TS1 = "true success",
    FS1 = "false success",
    FF1 = "false failure",
    TF1 = "true failure"
)

# Create the plot
plot <- ggplot(plot_data_0.2, aes(x = factor(method), y = proportion, fill = metric_label)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color_mapping) +  # Custom colors
    labs(x = "MABF Method", y = "Proportion", fill = "Evidence Strength", title = "Proportion of Bayes Factors Categoried into Varying Levels of Evidence When the Underlying Effect is 0.2") +
    facet_grid(rep.n ~ rep.number, scales = "free_x", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")
print(plot)

############
############
# Pie chart
# Prepare the data
MABF_rates_0.2_rowise_pcutoff0.05_BFcutoff3 <- readRDS("./MABFanalyses/row-wise/MABF_rates_0.2_rowise_pcutoff0.05_BFcutoff3.RDS")
plot_data_0.2 <- MABF_rates_0.2_rowise_pcutoff0.05_BFcutoff3 %>%
    mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
    mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
    mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
    mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>%
    select(method, rep.n, rep.number, AD_true, TS1, FS1, TF1, FF1) %>%
    gather(key = "metric", value = "value", AD_true, TS1, FS1, TF1, FF1) %>%
    mutate(metric = factor(metric, levels = c("AD_true", "TS1", "FS1", "TF1", "FF1"))) %>%
    mutate(metric_label = factor(metric, 
                                 levels = c("AD_true", "TS1", "FS1", "TF1", "FF1"),
                                 labels = c("anecdotal", "true success", "false success", "true failure", "false failure"))) %>%
    group_by(method, rep.n, rep.number) %>%
    mutate(proportion = value / sum(value))

# Create modern high contrast colors for the AD_true metrics
true_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")

# Create a named vector for scale_fill_manual
color_mapping <- c(anecdotal = true_colors[1],
                   "true success" = true_colors[2],
                   "false success" = true_colors[3],
                   "true failure" = true_colors[4],
                   "false failure" = true_colors[5])

# Custom labeller for metrics
metric_labeller <- c(
    AD_true = "anecdotal",
    TS1 = "true success",
    FS1 = "false success",
    TF1 = "true failure",
    FF1 = "false failure"
)

# Create the pie chart
plot <- ggplot(plot_data_0.2, aes(x = "", y = proportion, fill = metric_label)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = color_mapping) +  # Modern high contrast colors
    labs(x = "", y = "Proportion", fill = "Evidence Category", title = "Proportion of Bayes Factors Categorized into Success and Failure When the Underlying Effect is 0.2") +
    facet_grid(method ~ rep.n + rep.number, scales = "free", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "bottom")

# Save the plot
ggsave(filename = "./plots/pie chart/piechart_0.2_TSTF,pcutoff=0.05,BFcutoff=3.jpg", plot = plot, width = 15, height = 10, units = "in", dpi = 300)


#######When underlying effect is 0.5#########
# Stacked bar chart
# Prepare the data
MABF_rates_0.5_rowise_pcutoff0.05_BFcutoff3 <- readRDS("./MABFanalyses/row-wise/MABF_rates_0.5_rowise_pcutoff0.05_BFcutoff3.RDS")
plot_data_0.5 <- MABF_rates_0.5_rowise_pcutoff0.05_BFcutoff3 %>%
    mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
    mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
    mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
    mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>%
    select(method, rep.n, rep.number, AD_true, TS1, FS1, FF1, TF1) %>%
    gather(key = "metric", value = "value", AD_true, TS1, FS1, FF1, TF1) %>%
    mutate(metric = factor(metric, levels = c("AD_true", "TS1", "FS1", "FF1",  "TF1"))) %>%
    mutate(metric_label = factor(metric, 
                                 levels = c("AD_true", "TS1", "FS1", "FF1", "TF1"),
                                 labels = c("anecdotal","true success", "false success", "false failure",  "true failure"))) %>%
    group_by(method, rep.n, rep.number) %>%
    mutate(proportion = value / sum(value))

# Create custom gradient colors for the AD_true metrics
true_colors <- c("#E69F00", "#56B4E9", "coral2", "#F0E442", "darkorchid")
# Create a named vector for scale_fill_manual
color_mapping <- c(anecdotal = true_colors[1],
                   "true success" = true_colors[2],
                   "false success" = true_colors[3],
                   "false failure" = true_colors[4],
                   "true failure" = true_colors[5])

# Custom labeller for metrics
metric_labeller <- c(
    AD_true = "anecdotal",
    TS1 = "true success",
    FS1 = "false success",
    FF1 = "false failure",
    TF1 = "true failure"
)

# Create the plot
plot <- ggplot(plot_data_0.5, aes(x = factor(method), y = proportion, fill = metric_label)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color_mapping) +  # Custom colors
    labs(x = "MABF Method", y = "Proportion", fill = "Evidence Strength", title = "Proportion of Bayes Factors Categoried into Varying Levels of Evidence When the Underlying Effect is 0.5") +
    facet_grid(rep.n ~ rep.number, scales = "free_x", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")
print(plot)

############
############
# Pie chart
# Prepare the data
MABF_rates_0.5_rowise_pcutoff0.05_BFcutoff3 <- readRDS("./MABFanalyses/row-wise/MABF_rates_0.5_rowise_pcutoff0.05_BFcutoff3.RDS")
plot_data_0.5 <- MABF_rates_0.5_rowise_pcutoff0.05_BFcutoff3 %>%
    mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
    mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
    mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
    mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>%
    select(method, rep.n, rep.number, AD_true, TS1, FS1, TF1, FF1) %>%
    gather(key = "metric", value = "value", AD_true, TS1, FS1, TF1, FF1) %>%
    mutate(metric = factor(metric, levels = c("AD_true", "TS1", "FS1", "TF1", "FF1"))) %>%
    mutate(metric_label = factor(metric, 
                                 levels = c("AD_true", "TS1", "FS1", "TF1", "FF1"),
                                 labels = c("anecdotal", "true success", "false success", "true failure", "false failure"))) %>%
    group_by(method, rep.n, rep.number) %>%
    mutate(proportion = value / sum(value))

# Create modern high contrast colors for the AD_true metrics
true_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")

# Create a named vector for scale_fill_manual
color_mapping <- c(anecdotal = true_colors[1],
                   "true success" = true_colors[2],
                   "false success" = true_colors[3],
                   "true failure" = true_colors[4],
                   "false failure" = true_colors[5])

# Custom labeller for metrics
metric_labeller <- c(
    AD_true = "anecdotal",
    TS1 = "true success",
    FS1 = "false success",
    TF1 = "true failure",
    FF1 = "false failure"
)

# Create the pie chart
plot <- ggplot(plot_data_0.5, aes(x = "", y = proportion, fill = metric_label)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = color_mapping) +  # Modern high contrast colors
    labs(x = "", y = "Proportion", fill = "Evidence Category", title = "Proportion of Bayes Factors Categorized into Success and Failure When the Underlying Effect is 0.5") +
    facet_grid(method ~ rep.n + rep.number, scales = "free", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "bottom")

# Save the plot
ggsave(filename = "./plots/pie chart/piechart_0.5_TSTF,pcutoff=0.05,BFcutoff=3.jpg", plot = plot, width = 15, height = 10, units = "in", dpi = 300)
