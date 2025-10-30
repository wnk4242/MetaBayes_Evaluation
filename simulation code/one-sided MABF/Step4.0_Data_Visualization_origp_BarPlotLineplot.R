# Run Step3.0_OriginalStudy_Data_Analysis first
# This script will create bar plots of TPR and FPR for the original studies' p value.
# Combine rates from datasets of null effect, 0.2, and 0.5
# The bar plot is used in Step5.1_Shiny_origp_Histogram&HeatMap&Bar (only used the idea of bar plot, not the exact code)
rates_orig_0.20.5null <- bind_rows(rates_orig_0.2null, rates_orig_0.5null) %>% 
    mutate(nullES = rep("null", 18)) %>% 
    mutate(trueES = ifelse(delta == "0,0.2", 0.2, 0.5)) %>% 
    mutate(nullES = factor(nullES)) %>% 
    mutate(trueES = factor(trueES)) %>%
    mutate(orig.n = factor(orig.n)) %>%
    mutate(bias.level = censorFunc) %>% 
    relocate(nullES, .after = delta) %>% 
    relocate(trueES, .after = delta) %>% 
    relocate(bias.level, .before = qrpEnv)

# Calculate true positive rate
TPR_0.20.5 <- rates_orig_0.20.5null %>%
    mutate(true_positive_rate = TP / (TP + FN))
# Custom labels for orig.n
trueES_labels <- c(`0.2` = "True ES = 0.2", `0.5` = "True ES = 0.5")
# Generate the bar plot
plot <- ggplot(TPR_0.20.5, aes(x = orig.n, y = true_positive_rate, fill = bias.level)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    #facet_grid(trueES ~ bias.level, scales = "free_y", space = "free", labeller = as_labeller(trueES_labels)) +  # Custom labels for orig.n
    facet_grid(trueES ~ ., scales = "free_y", space = "free", labeller = as_labeller(trueES_labels)) + #another way of displaying bars
    labs(x = "Original Sample Size", y = "True Positive Rate", fill = "Bias level", title = "Average True Positive Rate of the Original Study Influenced by Sample Size and Bias Level When the Underlying Effect is True") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent_format(accuracy = 1)) +  # Set y-axis intervals to 5%
    theme_minimal() +
    theme(strip.text.x = element_blank(),  # Hide facet labels for trueES
          strip.text.y = element_text(angle = 0),  # Keep labels for orig.n
          strip.background = element_blank(),
          legend.position = "right")
ggsave(filename = "./plots/bar plot/TPR_orig_0.20.5_biasfilter.jpg", plot = plot, width = 12, height = 10)
# Calculate false positive rate
FPR_null <- rates_orig_0.20.5null %>%
    mutate(false_positive_rate = FP / (FP + TN))
# Generate the bar plot
plot <- ggplot(FPR_null, aes(x = orig.n, y = false_positive_rate, fill = bias.level)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    #facet_grid(nullES ~ bias.level) +
    facet_grid(nullES ~ .) + #another way of displaying bars
    labs(x = "Original Sample Size", y = "False Positive Rate", fill = "Bias level", title = "Average False Positive Rate of the Original Study Influenced by Sample Size and Bias Level When the Underlying Effect is Null") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.05), labels = scales::percent_format(accuracy = 1)) +  # Set y-axis intervals to 5%
    theme_minimal() +
    theme(strip.text = element_blank(), # Remove facet labels
          strip.background = element_blank(),
          legend.position = "right")
ggsave(filename = "./plots/bar plot/FPR_orig_0.20.5_biasfilter.jpg", plot = plot, width = 12, height = 10)

##################################
# The following script will create trend plot for original study ES based on true ES = 0, 0.2, and 0.5
# Trend plot for true ES = 0.2 and 0.5
# Combine ES from datasets of null effect, 0.2, and 0.5
ES_0.20.5 <- bind_rows(es_orig_0.2null, es_orig_0.5null) %>% 
    select(-avg_d_null, -med_d_null) %>% 
    mutate(trueES = ifelse(delta == "0,0.2", 0.2, 0.5)) %>% 
    mutate(trueES = factor(trueES)) %>%
    mutate(orig.n = factor(orig.n)) %>%
    mutate(bias.level = censorFunc) %>% 
    relocate(trueES, .after = delta) %>% 
    relocate(bias.level, .before = qrpEnv)
# Create labels
trueES_labels <- c(`0.2` = "True ES = 0.2", `0.5` = "True ES = 0.5")
# Generate trend plot
plot <- ggplot(ES_0.20.5, aes(x=orig.n, y=avg_d_true, color=bias.level)) +
    geom_point() +
    geom_line(aes(group=bias.level)) +
    facet_wrap(~ trueES, labeller = labeller(trueES = trueES_labels), strip.position = "top")+
    labs(title="Variation of Average Original Study Effect Size Influenced by Sample Size and Bias Level When the Underlying Effect is True",
         x="Original Sample Size",
         y="Average Effect Size",
         color="Bias Level") +
    scale_y_continuous(breaks = seq(0, 1.2, by = 0.1))+
    theme_minimal() +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y.right = element_text(angle = 90),
          panel.border = element_rect(color = "black", fill = NA))
ggsave(filename = "./plots/line plot/filtered original study ES/ES_orig_0.20.5_biasfilter.jpg", plot = plot, width = 12, height = 10)

# Trend plot for true ES = 0
# Data manipulation
ES_null <- bind_rows(es_orig_0.2null) %>% 
    select(-delta, -avg_d_true, -med_d_true) %>% 
    mutate(trueES = rep("null", 9)) %>% 
    mutate(trueES = factor(trueES)) %>% 
    mutate(orig.n = factor(orig.n)) %>%
    mutate(bias.level = censorFunc) %>% 
    relocate(trueES, .before = orig.n) %>% 
    relocate(bias.level, .before = qrpEnv)
# Create labels
trueES_labels <- c(`null` = "True ES = 0")
# Generate trend plot
plot <- ggplot(ES_null, aes(x=orig.n, y=avg_d_null, color=bias.level)) +
    geom_point() +
    geom_line(aes(group=bias.level)) +
    facet_wrap(~ trueES, labeller = labeller(trueES = trueES_labels), strip.position = "top")+
    labs(title="Variation of Average Original Study Effect Size Influenced by Sample Size and Bias Level When the Underlying Effect is Null",
         x="Original Sample Size",
         y="Average Effect Size",
         color="Bias Level") +
    scale_y_continuous(breaks = seq(0, 1.2, by = 0.1))+
    theme_minimal() +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y.right = element_text(angle = 90),
          panel.border = element_rect(color = "black", fill = NA))
ggsave(filename = "./plots/line plot/filtered original study ES/ES_null_0.20.5_biasfilter.jpg", plot = plot, width = 12, height = 10)

#####
# If you want to compare average original study ES across three true ES level in one plot
# Again, first run Step3.0_OriginalStudy_Data_Analysis
# We need to wrangle some data first to make the format of the datasets (0,0.2,0.5) the same
es_orig_.2 <- cbind(params_orig_0.2null, es_orig_0.2)
row.names(es_orig_.2) <- NULL
es_orig_.5 <- cbind(params_orig_0.5null, es_orig_0.5)
row.names(es_orig_.5) <- NULL
es_orig_Null <- cbind(params_orig_0.2null, es_orig_null)
row.names(es_orig_Null) <- NULL
es_orig_Null <- es_orig_Null %>% 
  mutate(avg_d = avg_d_null) %>% 
  mutate(med_d = med_d_null) %>% 
  select(-avg_d_null, -med_d_null) %>% 
  mutate(trueES = ifelse(delta == "0,0.2", 0)) %>% 
  select(-delta) %>% 
  relocate(trueES, .before = orig.n)
es_orig_.2 <- es_orig_.2 %>% 
  mutate(avg_d = avg_d_true) %>% 
  mutate(med_d = med_d_true) %>% 
  select(-avg_d_true, -med_d_true) %>% 
  mutate(trueES = ifelse(delta == "0,0.2", 0.2)) %>% 
  select(-delta) %>% 
  relocate(trueES, .before = orig.n)
es_orig_.5 <- es_orig_.5 %>% 
  mutate(avg_d = avg_d_true) %>% 
  mutate(med_d = med_d_true) %>% 
  select(-avg_d_true, -med_d_true) %>% 
  mutate(trueES = ifelse(delta == "0,0.5", 0.5)) %>% 
  select(-delta) %>% 
  relocate(trueES, .before = orig.n)
# Combine ES from datasets of null effect, 0.2, and 0.5
ES_.2.5Null <- bind_rows(es_orig_.2, es_orig_.5, es_orig_Null) %>% 
  mutate(trueES = factor(trueES)) %>%
  mutate(orig.n = factor(orig.n)) %>%
  mutate(bias.level = censorFunc) %>% 
  relocate(bias.level, .before = qrpEnv)
ES_.2.5Null <- ES_.2.5Null %>%
  mutate(bias.level = factor(bias.level,
                             levels = c("low", "medium", "high"),
                             labels = c("Low", "Medium", "High")))
# Create labels
trueES_labels <- c(`0` = "Population effect size = 0", `0.2` = "Population effect size = 0.2", `0.5` = "Population effect size = 0.5")
# Generate trend plot
# If want to add a title, it could be "Variation of average original study effect size influenced by sample size and bias level"
plot <- ggplot(ES_.2.5Null, aes(x=orig.n, y=avg_d, color=bias.level)) +
  geom_point() +
  geom_line(aes(group=bias.level)) +
  facet_wrap(~ trueES, labeller = labeller(trueES = trueES_labels), strip.position = "top")+
  labs(
       x="Original Sample Size",
       y="Average Effect Size",
       color="Bias Level") +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.1))+
  theme_minimal() +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.right = element_text(angle = 90),
        panel.border = element_rect(color = "black", fill = NA))
ggsave(filename = "./plots/line plot/filtered original study ES/ES_orig_.2.5Null_biasfilter.jpg", plot = plot, width = 12, height = 10)
