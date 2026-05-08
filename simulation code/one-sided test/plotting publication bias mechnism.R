# This script was used in the heterogeneity and ES project. You can find it in that folder.
library(ggplot2)
library(scales)

# Publication probability function
publish_decision <- function(p_value, d, prob) {
  if (!(prob %in% c("low", "medium", "high"))) {
    stop("prob must be 'low', 'medium', or 'high'")
  }
  
  prob_values <- list(
    low = 0.9,
    medium = 0.6,
    high = 0.5
  )
  
  prob_values2 <- list(
    low = 0.9,
    medium = 0.2,
    high = 0.1
  )
  
  prob_values3 <- list(
    low = 0.9,
    medium = 0.1,
    high = 0.05
  )
  
  inverse_sigmoid <- function(x, startX, endX, x0, maxY, k, L, prob_assign) {
    return((L) / (1 + exp(-k * (x - x0))) + prob_assign)
  }
  
  assign_publish <- function(p_value, d, prob_values, prob) {
    if (p_value >= 0.025 & p_value < 0.05 & d > 0) {
      prob_assign <- prob_values[[prob]]
      prob_assigned <- inverse_sigmoid(
        p_value,
        0.05,
        0.025,
        (0.025 + 0.05) / 2,
        0.9,
        -1000,
        0.9 - prob_assign,
        prob_assign
      )
      return(prob_assigned)
      
    } else if (p_value < 0.025 & d > 0) {
      return(1)
      
    } else if (p_value >= 0.05 & d > 0) {
      prob_assign <- prob_values2[[prob]]
      return(prob_assign)
      
    } else if (d < 0) {
      prob_assign <- prob_values3[[prob]]
      return(prob_assign)
    }
  }
  
  prob <- assign_publish(p_value, d, prob_values, prob)
  return(prob)
}

# Generate data for the plot
p_values <- seq(0, 0.1, by = 0.0001)
bias_levels <- c("low", "medium", "high")

results <- expand.grid(
  p_value = p_values,
  bias = bias_levels
)

results$probability <- NA

# Correct effect direction
for (i in 1:nrow(results)) {
  results$probability[i] <- publish_decision(
    p_value = results$p_value[i],
    d = 1,
    prob = results$bias[i]
  )
}

pub_bias_plot <- ggplot(
  results,
  aes(
    x = p_value,
    y = probability,
    color = bias,
    linetype = bias
  )
) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    breaks = seq(0, 0.1, by = 0.01),
    limits = c(0, 0.1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_color_manual(
    values = c(
      "low" = "#1b9e77",
      "medium" = "#7570b3",
      "high" = "#d95f02"
    ),
    labels = c("Low", "Medium", "High")
  ) +
  scale_linetype_manual(
    values = c(
      "low" = "solid",   # green line
      "medium" = "twodash",
      "high" = "dotted"   # red/orange line
    ),
    labels = c("Low", "Medium", "High")
  ) +
  labs(
    x = expression(italic(p)*"-value"),
    y = "Probability of Publication",
    color = "Bias level",
    linetype = "Bias level"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    
    panel.grid.minor = element_blank(),
    
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_line(color = "black", linewidth = 0.6),
    panel.border = element_blank(),
    
    plot.margin = margin(t = 10, r = 25, b = 10, l = 10),
    
    legend.position = "inside",
    legend.position.inside = c(0.1, 0.15),
    legend.justification = c(0, 0)
  )

pub_bias_plot

ggsave(
  filename = "publication bias.png",
  plot = pub_bias_plot,
  width = 8,
  height = 5.5,
  dpi = 300,
  bg = "white"
)