# Row-wise data to create bar plot
# Load necessary libraries
library(ggplot2)
library(dplyr)

# One factor as moderator
# Summarize the data to get the average SAR2 for each method and PB.level
plot_data_null <- MABF_rates_null_rowise %>%
  mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>% 
  mutate(rep.number = factor(rep.number, levels = c(2, 5, 10))) %>%
  mutate(rep.n = factor(rep.n, levels = c(40, 100, 400))) %>%
  mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>% 
  filter(!is.na(FSR2)) %>% 
  filter(ADR_null<0.2) %>% 
  group_by(method, rep.n)

avg_SAR2 <- plot_data_null %>%
  summarise(avg_SAR2 = mean(FSR2), .groups = 'drop') # Override grouping

# Create the bar plot with facets
ggplot(avg_SAR2, aes(x = method, y = avg_SAR2, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = sprintf("%.1f%%", avg_SAR2 * 100)), vjust = -0.5, position = position_dodge(width = 0.9)) +
  labs(title = "Average SAR2 for Each Method at Different Levels of PB.level",
       x = "Method",
       y = "Average SAR2") +
  facet_wrap(~ rep.n) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")



# Two factors as moderator
# Summarize the data to get the average SAR2 for each method, PB.level, and orig.n
plot_data_null <- MABF_rates_null_rowise %>%
  mutate(orig.n = factor(orig.n, levels = c(20, 50, 200))) %>% 
  mutate(rep.number = factor(rep.number, levels = c(2, 5, 10))) %>%
  mutate(rep.n = factor(rep.n, levels = c(40, 100, 400))) %>%
  mutate(PB.level = factor(PB.level, levels = c('low', 'medium', 'high'))) %>% 
  filter(!is.na(SAR2)) %>% 
  group_by(method, rep.n, rep.number)

avg_SAR2 <- plot_data_null %>%
  summarise(avg_SAR2 = mean(SAR2), .groups = 'drop') # Override grouping

# Create the bar plot with facets
ggplot(avg_SAR2, aes(x = method, y = avg_SAR2 * 100, fill = method)) + # Scale SAR2 to percentage
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = sprintf("%.2f%%", avg_SAR2 * 100)), vjust = -0.5, position = position_dodge(width = 0.9)) + # Add percentage labels
  scale_y_continuous(labels = scales::percent_format(scale = 1, accuracy = 1)) + # Correct y-axis scale
  labs(title = "Average SAR2 for Each Method at Different Levels of PB.level and orig.n",
       x = "Method",
       y = "Average SAR2 (%)") +
  facet_grid(rep.n ~ rep.number) + # Facet by PB.level and orig.n
  theme_minimal() +
  theme(legend.position = "none")


