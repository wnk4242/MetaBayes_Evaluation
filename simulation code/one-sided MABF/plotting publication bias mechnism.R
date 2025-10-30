# This script was used in the heterogeneity and ES project. You can find it in that folder.
library(ggplot2)

# The provided publish_decision function
publish_decision <- function(p_value, d, prob) {
  # Ensure that prob is one of the specified values
  if (!(prob %in% c("low", "medium", "high"))) {
    stop("prob must be 'low', 'medium', or 'high'")
  }
  
  # Define probabilities based on 'prob' argument
  prob_values <- list(
    low = 0.9,       # High probability to assign 1
    medium = 0.6,    # Medium probability to assign 1
    high = 0.5       # Lower probability to assign 1
  )
  
  prob_values2 <- list(
    low = 0.9,       # High probability to assign 1
    medium = 0.2,    # Lower probability to assign 1
    high = 0.1       # Lowest probability to assign 1
  )
  
  prob_values3 <- list(
    low = 0.9,       # High probability to assign 1
    medium = 0.1,    # Lower probability to assign 1
    high = 0.05      # Lowest probability to assign 1
  )
  
  # Inverse sigmoid function, as described
  inverse_sigmoid <- function(x, startX, endX, x0, maxY, k, L, prob_assign) {
    return((L) / (1 + exp(-k * (x - x0))) + prob_assign)
  }
  
  # Function to decide publication
  assign_publish <- function(p_value, d, prob_values, prob) {
    if (p_value >= 0.025 & p_value < 0.05 & d > 0) {
      prob_assign <- prob_values[[prob]]
      prob_assigned <- inverse_sigmoid(p_value, 0.05, 0.025, (0.025 + 0.05)/2, 0.9, -1000, 0.9 - prob_assign, prob_assign)
      return(prob_assigned)
    } else if (p_value < 0.025 & d > 0) {
      return(1) # 100% probability to assign 1
    } else if (p_value >= 0.05 & d > 0) {
      prob_assign <- prob_values2[[prob]]
      return(prob_assign)
    } else if (d < 0) {
      prob_assign <- prob_values3[[prob]]
      return(prob_assign)
    }
  }
  # My original function has the following code:
  # publish <- assign_publish(p_value, d, prob_values, prob)
  # return(list(p_value = p_value, publish = publish))
  # For plotting, use the following instead
  prob <- assign_publish(p_value, d, prob_values, prob)
  return(prob)
  
}

# Generate a sequence of p-values between 0 and 0.1
p_values <- seq(0, 0.1, by = 0.0001)
bias_levels <- c("low", "medium", "high")
results <- expand.grid(p_value = p_values, bias = bias_levels)
results$probability <- NA

# Simulate publication decisions for each p-value and bias level
for (i in 1:nrow(results)) {
  results$probability[i] <- publish_decision(results$p_value[i], 1, results$bias[i])
}

# Plot the results with custom colors for each bias level
# Plotting Correct direction
ggplot(results, aes(x = p_value, y = probability, color = bias)) +
  geom_line() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.05)) + # Detailed y-axis
  scale_color_manual(values = c("low" = "green", "medium" = "blue", "high" = "red")) + # Custom colors
  labs(title = "Theoretical Probability of Publication by P-Value and Bias Level",
       x = "P-Value",
       y = "Probability of Publication") +
  theme_minimal()

# Plotting incorrect direction
# Adjust the simulation for an incorrect direction by setting d to a negative value
for (i in 1:nrow(results)) {
  results$probability[i] <- publish_decision(results$p_value[i], -1, results$bias[i])
}

# Continue with plotting as before, using custom colors
ggplot(results, aes(x = p_value, y = probability, color = bias)) +
  geom_line() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.05)) + # Detailed y-axis
  scale_color_manual(values = c("low" = "green", "medium" = "blue", "high" = "red")) + # Custom colors
  labs(title = "Theoretical Probability of Publication by P-Value and Bias Level (Incorrect Direction)",
       x = "P-Value",
       y = "Probability of Publication") +
  theme_minimal()


####
# Modern plot for correct direction
ggplot(results, aes(x = p_value, y = probability, color = bias)) +
  geom_line(size = 0.8) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 0.1, by = 0.01)) +
  scale_color_manual(
    values = c("low" = "#1b9e77", "medium" = "#7570b3", "high" = "#d95f02"),
    labels = c("Low", "Medium", "High")
  ) +
  labs(
    x = "p-Value",
    y = "Probability of Publication",
    color = "Bias Level"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray30"),
    legend.position = "inside",                 # new syntax
    legend.position.inside = c(0.1, 0.15),      # coordinates inside plot
    legend.justification = c(0, 0),             # anchor legend’s bottom-left corner
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    panel.grid.minor = element_blank()
  )


# Modern plot for incorrect direction
ggplot(results, aes(x = p_value, y = probability, color = bias)) +
  geom_line(size = 1.2) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 0.1, by = 0.01)) +
  scale_color_manual(values = c("low" = "#1b9e77", 
                                "medium" = "#7570b3", 
                                "high" = "#d95f02")) + 
  labs(title = "Probability of Publication by p-Value and Bias Level",
       subtitle = "Incorrect Effect Direction",
       x = "p-Value",
       y = "Probability of Publication",
       color = "Bias Level") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray30"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )
