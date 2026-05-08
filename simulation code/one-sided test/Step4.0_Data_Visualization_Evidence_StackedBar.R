# This stacked bar plot shows the proportion of BF falling into different categories of evidence
library(ggplot2)
library(dplyr)
library(tidyr)
##Underlying effect is null
# Prepare the data
plot_data_null <- MABF_rates_null_rowise_pcutoff0.05_BFcutoff30 %>%
  mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
  mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
  mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
  mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>%
  select(method, rep.n, rep.number, AD_null_anecdotal, AD_null_moderate, AD_null_strong, AD_null_vstrong, AD_null_estrong) %>%
  gather(key = "metric", value = "value", AD_null_anecdotal, AD_null_moderate, AD_null_strong, AD_null_vstrong, AD_null_estrong) %>%
  mutate(metric = factor(metric, levels = c("AD_null_anecdotal", "AD_null_moderate", "AD_null_strong", "AD_null_vstrong", "AD_null_estrong"))) %>%
  mutate(metric_label = factor(metric, 
                               levels = c("AD_null_anecdotal", "AD_null_moderate", "AD_null_strong", "AD_null_vstrong", "AD_null_estrong"),
                               labels = c("anecdotal", "moderate", "strong", "very strong", "extremely strong"))) %>%
  group_by(method, rep.n, rep.number) %>%
  mutate(proportion = value / sum(value))

# Create custom gradient colors for the AD_null metrics
AD_null_colors <- c("#00008B", "#0000CD", "#1E90FF", "#87CEEB", "#B0E0E6")

# Create a named vector for scale_fill_manual
color_mapping <- c(anecdotal = AD_null_colors[1],
                   moderate = AD_null_colors[2],
                   strong = AD_null_colors[3],
                   "very strong" = AD_null_colors[4],
                   "extremely strong" = AD_null_colors[5])

# Custom labeller for metrics
metric_labeller <- c(
  AD_null_anecdotal = "anecdotal",
  AD_null_moderate = "moderate",
  AD_null_strong = "strong",
  AD_null_vstrong = "very strong",
  AD_null_estrong = "extremely strong"
)

# Create the plot
plot <- ggplot(plot_data_null, aes(x = factor(method), y = proportion, fill = metric_label)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_mapping) +  # Custom colors
  labs(x = "MABF Method", y = "Proportion", fill = "Evidence Strength", title = "Proportion of Bayes Factors Categoried into Varying Levels of Evidence When the Underlying Effect is Null") +
  facet_grid(rep.n ~ rep.number, scales = "free_x", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")

# Save the plot with increased width and height
ggsave(filename = "./plots/stacked bar plot/null_stacked.jpg",
       plot = plot, width = 15, height = 12)



##Underlying effect is 0.2
library(ggplot2)
library(dplyr)
library(tidyr)
##Underlying effect is true
# Prepare the data
plot_data_true <- MABF_rates_0.2_rowise_pcutoff0.05_BFcutoff30 %>%
  mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
  mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
  mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
  mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>%
  select(method, rep.n, rep.number, AD_true_anecdotal, AD_true_moderate, AD_true_strong, AD_true_vstrong, AD_true_estrong) %>%
  gather(key = "metric", value = "value", AD_true_anecdotal, AD_true_moderate, AD_true_strong, AD_true_vstrong, AD_true_estrong) %>%
  mutate(metric = factor(metric, levels = c("AD_true_anecdotal", "AD_true_moderate", "AD_true_strong", "AD_true_vstrong", "AD_true_estrong"))) %>%
  mutate(metric_label = factor(metric, 
                               levels = c("AD_true_anecdotal", "AD_true_moderate", "AD_true_strong", "AD_true_vstrong", "AD_true_estrong"),
                               labels = c("anecdotal", "moderate", "strong", "very strong", "extremely strong"))) %>%
  group_by(method, rep.n, rep.number) %>%
  mutate(proportion = value / sum(value))

# Create custom gradient colors for the AD_true metrics
AD_true_colors <- c("#00008B", "#0000CD", "#1E90FF", "#87CEEB", "#B0E0E6")

# Create a named vector for scale_fill_manual
color_mapping <- c(anecdotal = AD_true_colors[1],
                   moderate = AD_true_colors[2],
                   strong = AD_true_colors[3],
                   "very strong" = AD_true_colors[4],
                   "extremely strong" = AD_true_colors[5])

# Custom labeller for metrics
metric_labeller <- c(
  AD_true_anecdotal = "anecdotal",
  AD_true_moderate = "moderate",
  AD_true_strong = "strong",
  AD_true_vstrong = "very strong",
  AD_true_estrong = "extremely strong"
)

# Create the plot
plot <- ggplot(plot_data_true, aes(x = factor(method), y = proportion, fill = metric_label)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_mapping) +  # Custom colors
  labs(x = "MABF Method", y = "Proportion", fill = "Evidence Strength", title = "Proportion of Bayes Factors Categoried into Varying Levels of Evidence When the Underlying Effect is true") +
  facet_grid(rep.n ~ rep.number, scales = "free_x", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")

# Save the plot with increased width and height
ggsave(filename = "./plots/stacked bar plot/true_0.2_stacked.jpg",
       plot = plot, width = 15, height = 12)



##Underlying effect is 0.5
library(ggplot2)
library(dplyr)
library(tidyr)
##Underlying effect is true
# Prepare the data
plot_data_true <- MABF_rates_0.5_rowise_pcutoff0.05_BFcutoff30 %>%
  mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>%
  mutate(rep.number = factor(rep.number, levels = c(2, 5, 10), labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10"))) %>%
  mutate(rep.n = factor(rep.n, levels = c(40, 100, 400), labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400"))) %>%
  mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>%
  select(method, rep.n, rep.number, AD_true_anecdotal, AD_true_moderate, AD_true_strong, AD_true_vstrong, AD_true_estrong) %>%
  gather(key = "metric", value = "value", AD_true_anecdotal, AD_true_moderate, AD_true_strong, AD_true_vstrong, AD_true_estrong) %>%
  mutate(metric = factor(metric, levels = c("AD_true_anecdotal", "AD_true_moderate", "AD_true_strong", "AD_true_vstrong", "AD_true_estrong"))) %>%
  mutate(metric_label = factor(metric, 
                               levels = c("AD_true_anecdotal", "AD_true_moderate", "AD_true_strong", "AD_true_vstrong", "AD_true_estrong"),
                               labels = c("anecdotal", "moderate", "strong", "very strong", "extremely strong"))) %>%
  group_by(method, rep.n, rep.number) %>%
  mutate(proportion = value / sum(value))

# Create custom gradient colors for the AD_true metrics
AD_true_colors <- c("#00008B", "#0000CD", "#1E90FF", "#87CEEB", "#B0E0E6")

# Create a named vector for scale_fill_manual
color_mapping <- c(anecdotal = AD_true_colors[1],
                   moderate = AD_true_colors[2],
                   strong = AD_true_colors[3],
                   "very strong" = AD_true_colors[4],
                   "extremely strong" = AD_true_colors[5])

# Custom labeller for metrics
metric_labeller <- c(
  AD_true_anecdotal = "anecdotal",
  AD_true_moderate = "moderate",
  AD_true_strong = "strong",
  AD_true_vstrong = "very strong",
  AD_true_estrong = "extremely strong"
)

# Create the plot
plot <- ggplot(plot_data_true, aes(x = factor(method), y = proportion, fill = metric_label)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_mapping) +  # Custom colors
  labs(x = "MABF Method", y = "Proportion", fill = "Evidence Strength", title = "Proportion of Bayes Factors Categoried into Varying Levels of Evidence When the Underlying Effect is true") +
  facet_grid(rep.n ~ rep.number, scales = "free_x", labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")

# Save the plot with increased width and height
ggsave(filename = "./plots/stacked bar plot/true_0.5_stacked.jpg",
       plot = plot, width = 15, height = 12)
