# Note: I2 is set to 10% in this stacked bar shiny app. I2 in this current shiny app is set to 30%.
library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(scales)

# ============================================================
# Analytic functions
# ============================================================

clip01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

# ============================================================
# Optional stopping function using stacked-bar Shiny logic
# ============================================================

calc_optional_stopping_prob <- function(theta,
                                        n_orig,
                                        alpha_orig,
                                        q,
                                        max_mult = 2,
                                        step = 5) {
  
  n_max <- n_orig + q * (n_orig * max_mult - n_orig)
  n_seq <- seq(n_orig, floor(n_max), by = step)
  
  zcrit <- qnorm(1 - alpha_orig)
  
  if (theta == 0) {
    reject_prob_each_look <- rep(alpha_orig, length(n_seq))
  } else {
    reject_prob_each_look <- 1 - pnorm(
      zcrit - theta * sqrt(n_seq / 2)
    )
  }
  
  pi_os <- 1 - prod(1 - reject_prob_each_look)
  
  clip01(pi_os)
}

calc_outcome_probs <- function(theta,
                               b,
                               alpha_orig,
                               n_orig,
                               n_rep,
                               N_rep,
                               alpha_meta = 0.05,
                               I2 = 0.10) {
  
  if (theta == 0) {
    
    FP_orig_nominal <- alpha_orig
    
    pi0_OS <- calc_optional_stopping_prob(
      theta = 0,
      n_orig = n_orig,
      alpha_orig = alpha_orig,
      q = b
    )
    
    FP_orig <- (1 - b) * FP_orig_nominal + b * pi0_OS
    FP_orig <- clip01(FP_orig)
    TN_orig <- 1 - FP_orig
    
    FP_rep <- alpha_meta
    TN_rep <- 1 - FP_rep
    
    TS <- TN_orig * TN_rep
    FS <- FP_orig * FP_rep
    FF <- FP_orig * TN_rep
    TF <- TN_orig * FP_rep
    
  } else {
    
    zcrit_orig <- qnorm(1 - alpha_orig)
    
    TP_orig_nominal <- 1 - pnorm(
      zcrit_orig - theta * sqrt(n_orig / 2)
    )
    
    pi1_OS <- calc_optional_stopping_prob(
      theta = theta,
      n_orig = n_orig,
      alpha_orig = alpha_orig,
      q = b
    )
    
    TP_orig <- (1 - b) * TP_orig_nominal + b * pi1_OS
    TP_orig <- clip01(TP_orig)
    FN_orig <- 1 - TP_orig
    
    v <- 2 / n_rep + theta^2 / (4 * n_rep)
    tau2 <- I2 * v / (1 - I2)
    vbar <- (v + tau2) / N_rep
    
    lambda <- theta / sqrt(vbar)
    zcrit_meta <- qnorm(1 - alpha_meta / 2)
    
    TP_rep <- 1 - pnorm(zcrit_meta - lambda) +
      pnorm(-zcrit_meta - lambda)
    
    TP_rep <- clip01(TP_rep)
    FN_rep <- 1 - TP_rep
    
    TS <- TP_orig * TP_rep
    FS <- FN_orig * FN_rep
    FF <- FN_orig * TP_rep
    TF <- TP_orig * FN_rep
  }
  
  total <- TS + FS + FF + TF
  
  tibble(
    TS = TS / total,
    FS = FS / total,
    TF = TF / total,
    FF = FF / total
  )
}

# ============================================================
# Helper functions for x-axis values
# ============================================================

get_x_limits <- function(x_var) {
  switch(
    x_var,
    theta = c(0, 0.5),
    b = c(0, 0.8),
    alpha_orig = c(0.01, 0.05),
    n_orig = c(20, 200)
  )
}

get_x_breaks <- function(x_var) {
  switch(
    x_var,
    theta = seq(0, 0.5, by = 0.1),
    b = seq(0, 0.8, by = 0.2),
    alpha_orig = seq(0.01, 0.05, by = 0.01),
    n_orig = seq(20, 200, by = 30)
  )
}

clip_x_value <- function(x, x_var) {
  limits <- get_x_limits(x_var)
  pmin(pmax(x, limits[1]), limits[2])
}

get_default_x <- function(x_var,
                          fixed_theta,
                          fixed_b,
                          fixed_alpha_orig,
                          fixed_n_orig) {
  
  x <- switch(
    x_var,
    theta = {
      if (fixed_theta == 0) 0 else fixed_theta
    },
    b = fixed_b,
    alpha_orig = fixed_alpha_orig,
    n_orig = fixed_n_orig
  )
  
  clip_x_value(x, x_var)
}

format_x_value <- function(x, x_var) {
  if (x_var %in% c("theta", "b", "alpha_orig")) {
    sprintf("%.3f", x)
  } else {
    sprintf("%.0f", x)
  }
}

get_x_name <- function(x_var) {
  switch(
    x_var,
    theta = "true effect size",
    b = "Bias level",
    alpha_orig = "original alpha level",
    n_orig = "original sample size"
  )
}

get_x_label_for_bar <- function(x_var) {
  switch(
    x_var,
    theta = "\u03b8",
    b = "Bias",
    alpha_orig = "\u03b1_orig",
    n_orig = "n_orig"
  )
}

# ============================================================
# Create one-panel curve data with selectable x-axis
# ============================================================

make_curve_data <- function(x_var,
                            fixed_theta,
                            fixed_b,
                            fixed_alpha_orig,
                            fixed_n_orig,
                            n_rep,
                            N_rep,
                            alpha_meta,
                            I2) {
  
  if (x_var == "theta") {
    
    if (fixed_theta == 0) {
      x_vals <- 0
    } else {
      x_vals <- seq(0.001, 0.50, length.out = 500)
    }
    
    x_label_text <- "True effect size"
    x_axis_label_gg <- expression(italic(theta))
    x_axis_label_plotly <- "<i>\u03b8</i>"
  }
  
  if (x_var == "b") {
    x_vals <- seq(0, 0.80, length.out = 500)
    x_label_text <- "Bias level"
    x_axis_label_gg <- "Bias Level"
    x_axis_label_plotly <- "Bias Level"
  }
  
  if (x_var == "alpha_orig") {
    x_vals <- seq(0.01, 0.05, length.out = 500)
    x_label_text <- "Original alpha level"
    x_axis_label_gg <- expression(alpha[orig])
    x_axis_label_plotly <- "<i>\u03b1</i><sub>orig</sub>"
  }
  
  if (x_var == "n_orig") {
    x_vals <- seq(20, 200, length.out = 500)
    x_label_text <- "Original sample size"
    x_axis_label_gg <- expression(n[orig])
    x_axis_label_plotly <- "n<sub>orig</sub>"
  }
  
  plot_data <- lapply(x_vals, function(x) {
    
    theta_val <- fixed_theta
    b_val <- fixed_b
    alpha_orig_val <- fixed_alpha_orig
    n_orig_val <- fixed_n_orig
    
    if (x_var == "theta") {
      theta_val <- x
    }
    
    if (x_var == "b") {
      b_val <- x
    }
    
    if (x_var == "alpha_orig") {
      alpha_orig_val <- x
    }
    
    if (x_var == "n_orig") {
      n_orig_val <- x
    }
    
    probs <- calc_outcome_probs(
      theta = theta_val,
      b = b_val,
      alpha_orig = alpha_orig_val,
      n_orig = n_orig_val,
      n_rep = n_rep,
      N_rep = N_rep,
      alpha_meta = alpha_meta,
      I2 = I2
    )
    
    probs |>
      mutate(
        x = x,
        theta = theta_val,
        b = b_val,
        alpha_orig = alpha_orig_val,
        n_orig = n_orig_val
      ) |>
      pivot_longer(
        cols = c(TS, FS, TF, FF),
        names_to = "category",
        values_to = "proportion"
      ) |>
      mutate(
        category = factor(
          category,
          levels = c("TS", "FS", "TF", "FF")
        ),
        definition = ifelse(
          theta_val == 0,
          "Null-effect classification",
          "Non-null-effect classification"
        ),
        hover_text = paste0(
          "<span style='text-align:left; display:block;'>",
          "<b>",
          recode(
            as.character(category),
            "TS" = "True success",
            "FS" = "False success",
            "TF" = "True failure",
            "FF" = "False failure"
          ),
          " (", category, ")</b>",
          "<br><b>Proportion:</b> ", percent(proportion, accuracy = 0.01),
          "<br><br><b>True effect size (θ):</b> ", round(theta, 4),
          "<br><b>Bias level:</b> ", round(b, 4),
          "<br><b>α<sub>orig</sub>:</b> ", round(alpha_orig, 4),
          "<br><b>n<sub>orig</sub>:</b> ", round(n_orig, 2),
          "<br><b>n<sub>rep</sub>:</b> ", n_rep,
          "<br><b>N<sub>rep</sub>:</b> ", N_rep,
          "</span>"
        )
      )
    
  }) |>
    bind_rows()
  
  list(
    data = plot_data,
    x_axis_label_gg = x_axis_label_gg,
    x_axis_label_plotly = x_axis_label_plotly,
    x_label_text = x_label_text
  )
}

# ============================================================
# Make ggplot object
# ============================================================

make_curve_ggplot <- function(curve_obj, x_var) {
  
  data <- curve_obj$data |>
    mutate(
      category_full = factor(
        dplyr::recode(
          as.character(category),
          "TS" = "True success (TS)",
          "FS" = "False success (FS)",
          "TF" = "True failure (TF)",
          "FF" = "False failure (FF)"
        ),
        levels = c(
          "True success (TS)",
          "False success (FS)",
          "True failure (TF)",
          "False failure (FF)"
        )
      )
    )
  
  title_text <- switch(
    x_var,
    theta = "Replication-Success Classification Across True Effect Size",
    b = "Replication-Success Classification Across Bias Level",
    alpha_orig = "Replication-Success Classification Across Original Alpha Level",
    n_orig = "Replication-Success Classification Across Original Sample Size"
  )
  
  x_breaks <- get_x_breaks(x_var)
  x_limits <- get_x_limits(x_var)
  
  p <- ggplot(
    data,
    aes(
      x = x,
      y = proportion,
      color = category_full,
      linetype = category_full,
      text = hover_text,
      group = category_full
    )
  )
  
  if (x_var == "theta") {
    
    if (all(data$theta == 0)) {
      p <- p +
        geom_point(
          size = 3,
          show.legend = TRUE
        )
    } else {
      p <- p +
        geom_line(linewidth = 0.8)
    }
    
  } else {
    p <- p +
      geom_line(linewidth = 0.8)
  }
  
  p +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_x_continuous(
      limits = x_limits,
      breaks = x_breaks,
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_color_manual(
      values = c(
        "True success (TS)" = "black",
        "False success (FS)" = "#1f77b4",
        "True failure (TF)" = "#2ca02c",
        "False failure (FF)" = "#d62728"
      ),
      breaks = c(
        "True success (TS)",
        "False success (FS)",
        "True failure (TF)",
        "False failure (FF)"
      )
    ) +
    scale_linetype_manual(
      values = c(
        "True success (TS)" = "solid",
        "False success (FS)" = "longdash",
        "True failure (TF)" = "dotted",
        "False failure (FF)" = "solid"
      ),
      breaks = c(
        "True success (TS)",
        "False success (FS)",
        "True failure (TF)",
        "False failure (FF)"
      )
    ) +
    labs(
      title = title_text,
      x = curve_obj$x_axis_label_gg,
      y = "Proportion",
      color = NULL,
      linetype = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title.x = element_text(size = 13, margin = margin(t = 12)),
      axis.title.y = element_text(size = 13, margin = margin(r = 8)),
      axis.text = element_text(size = 11, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 11),
      legend.key.width = unit(1.3, "cm"),
      legend.key.height = unit(0.25, "cm"),
      plot.margin = margin(10, 12, 35, 12)
    ) +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE),
      linetype = guide_legend(nrow = 1, byrow = TRUE)
    )
}

# ============================================================
# Convert ggplot to plotly and add vertical line
# ============================================================

make_curve_plotly <- function(curve_obj, x_var, selected_x_value) {
  
  p <- make_curve_ggplot(curve_obj, x_var)
  
  x_limits <- get_x_limits(x_var)
  selected_x_value <- clip_x_value(selected_x_value, x_var)
  
  ggplotly(
    p,
    tooltip = "text",
    source = "curve_source"
  ) |>
    add_trace(
      x = c(selected_x_value, selected_x_value),
      y = c(0, 1),
      type = "scatter",
      mode = "lines",
      line = list(
        color = "grey",
        width = 2,
        dash = "dot"
      ),
      name = "",
      showlegend = FALSE,
      hoverinfo = "skip",
      inherit = FALSE
    ) |>
    layout(
      hovermode = "closest",
      dragmode = "pan",
      xaxis = list(
        title = list(
          text = curve_obj$x_axis_label_plotly,
          standoff = 18
        ),
        range = x_limits
      ),
      yaxis = list(
        range = c(0, 1),
        tickformat = ".0%"
      ),
      annotations = list(
        list(
          x = selected_x_value,
          y = 1.03,
          xref = "x",
          yref = "y",
          text = paste0(
            get_x_label_for_bar(x_var),
            " = ",
            format_x_value(selected_x_value, x_var)
          ),
          showarrow = FALSE,
          font = list(size = 12, color = "black")
        )
      ),
      legend = list(
        orientation = "h",
        x = 0.5,
        xanchor = "center",
        y = -0.22
      ),
      hoverlabel = list(
        align = "left"
      ),
      margin = list(l = 75, r = 25, t = 70, b = 130)
    ) |>
    config(
      displayModeBar = TRUE,
      scrollZoom = TRUE
    )
}

# ============================================================
# Make cross-sectional stacked bar
# ============================================================

make_cross_section_bar <- function(probs,
                                   selected_params,
                                   x_var,
                                   selected_x_value) {
  
  bar_data <- probs |>
    pivot_longer(
      cols = c(TS, FS, TF, FF),
      names_to = "category",
      values_to = "proportion"
    ) |>
    mutate(
      category = factor(
        category,
        levels = c("TF", "FF", "FS", "TS")
      ),
      category_label = recode(
        as.character(category),
        "TS" = "True success",
        "FS" = "False success",
        "TF" = "True failure",
        "FF" = "False failure"
      ),
      proportion_label = percent(proportion, accuracy = 0.01),
      text_label = ifelse(
        proportion >= 0.035,
        percent(proportion, accuracy = 0.1),
        ""
      ),
      hover_text = paste0(
        "<span style='text-align:left; display:block;'>",
        "<b>", category_label, " (", category, ")</b>",
        "<br><b>Proportion:</b> ", proportion_label,
        "<br><br><b>True effect size (θ):</b> ", round(selected_params$theta, 4),
        "<br><b>Bias level:</b> ", round(selected_params$b, 4),
        "<br><b>α<sub>orig</sub>:</b> ", round(selected_params$alpha_orig, 4),
        "<br><b>n<sub>orig</sub>:</b> ", round(selected_params$n_orig, 2),
        "<br><b>n<sub>rep</sub>:</b> ", selected_params$n_rep,
        "<br><b>N<sub>rep</sub>:</b> ", selected_params$N_rep,
        "</span>"
      )
    ) |>
    arrange(category)
  
  colors <- c(
    "TS" = "black",
    "FS" = "#1f77b4",
    "TF" = "#2ca02c",
    "FF" = "#d62728"
  )
  
  p <- plot_ly()
  
  for (i in seq_len(nrow(bar_data))) {
    
    this_cat <- as.character(bar_data$category[i])
    
    p <- p |>
      add_trace(
        type = "bar",
        x = "Replication-success classification",
        y = bar_data$proportion[i],
        name = this_cat,
        marker = list(
          color = colors[[this_cat]],
          line = list(color = "black", width = 0.5)
        ),
        text = bar_data$text_label[i],
        textposition = "inside",
        insidetextfont = list(color = "white", size = 12),
        hoverinfo = "text",
        hovertext = bar_data$hover_text[i],
        showlegend = TRUE
      )
  }
  
  p |>
    layout(
      title = list(
        text = paste0(
          "Cross-sectional View<br>",
          "<sup>",
          get_x_label_for_bar(x_var),
          " = ",
          format_x_value(selected_x_value, x_var),
          "</sup>"
        ),
        x = 0.02
      ),
      barmode = "stack",
      xaxis = list(
        title = "",
        showgrid = FALSE,
        zeroline = FALSE
      ),
      yaxis = list(
        title = "Proportion",
        range = c(0, 1),
        tickformat = ".0%",
        gridcolor = "rgba(0,0,0,0.12)"
      ),
      legend = list(
        orientation = "h",
        x = 0.5,
        xanchor = "center",
        y = -0.18
      ),
      hoverlabel = list(
        align = "left"
      ),
      margin = list(l = 65, r = 20, t = 70, b = 90),
      plot_bgcolor = "white",
      paper_bgcolor = "white"
    ) |>
    config(displayModeBar = TRUE)
}

# ============================================================
# Shiny UI
# ============================================================

ui <- fluidPage(
  
  titlePanel("Replication-Success Classification Curves"),
  
  sidebarLayout(
    
    sidebarPanel(
      width = 3,
      
      h4("Choose x-axis"),
      
      selectInput(
        inputId = "x_var",
        label = "X-axis variable",
        choices = c(
          "True effect size" = "theta",
          "Bias level" = "b",
          "Original alpha level" = "alpha_orig",
          "Original sample size" = "n_orig"
        ),
        selected = "theta"
      ),
      
      helpText(
        "Click on the curve to move the grey reference line and update the stacked bar."
      ),
      
      tags$hr(),
      
      h4("Original-study settings"),
      
      numericInput(
        inputId = "fixed_theta",
        label = tags$span(
          title = "when the x-axis is not true effect size",
          style = "cursor: help;",
          HTML("True effect size (&theta;)")
        ),
        value = 0.20,
        min = 0,
        max = 0.50,
        step = 0.05
      ),
      
      numericInput(
        inputId = "fixed_b",
        label = tags$span(
          title = "when the x-axis is not bias level",
          style = "cursor: help;",
          HTML("Bias level (b)")
        ),
        value = 0.50,
        min = 0,
        max = 0.80,
        step = 0.05
      ),
      
      numericInput(
        inputId = "fixed_alpha_orig",
        label = tags$span(
          title = "when the x-axis is not original alpha level",
          style = "cursor: help;",
          HTML("Original study nominal sig. level (&alpha;<sub>orig</sub>)")
        ),
        value = 0.05,
        min = 0.01,
        max = 0.05,
        step = 0.001
      ),
      
      numericInput(
        inputId = "fixed_n_orig",
        label = tags$span(
          title = "when the x-axis is not original sample size",
          style = "cursor: help;",
          HTML("Original study sample size (n<sub>orig</sub>)")
        ),
        value = 50,
        min = 20,
        max = 200,
        step = 10
      ),
      
      tags$hr(),
      
      h4("Replication design"),
      
      sliderInput(
        inputId = "n_rep",
        label = HTML("Replication sample size (n<sub>rep</sub>)"),
        min = 40,
        max = 400,
        value = 100,
        step = 10
      ),
      
      sliderInput(
        inputId = "N_rep",
        label = HTML("Number of replications (N<sub>rep</sub>)"),
        min = 2,
        max = 10,
        value = 5,
        step = 1
      ),
      
      tags$hr(),
      
      h4("Model settings"),
      
      sliderInput(
        inputId = "alpha_meta",
        label = HTML("Meta-analytic alpha (&alpha;<sub>meta</sub>)"),
        min = 0.001,
        max = 0.10,
        value = 0.05,
        step = 0.001
      ),
      
      sliderInput(
        inputId = "I2",
        label = HTML("Heterogeneity level (I<sup>2</sup>)"),
        min = 0,
        max = 0.80,
        value = 0.10,
        step = 0.01
      ),
      
      tags$hr(),
      
      downloadButton(
        outputId = "download_plot",
        label = "Download PNG curve"
      )
    ),
    
    mainPanel(
      width = 9,
      
      fluidRow(
        column(
          width = 8,
          plotlyOutput(
            outputId = "curve_plot",
            height = "700px"
          )
        ),
        column(
          width = 4,
          plotlyOutput(
            outputId = "cross_section_bar",
            height = "700px"
          )
        )
      ),
      
      tags$hr(),
      
      htmlOutput("selected_value_text"),
      
      tags$hr(),
      
      htmlOutput("interpretation")
    )
  )
)

# ============================================================
# Shiny server
# ============================================================

server <- function(input, output, session) {
  
  # Stores the selected x-value for the vertical line and stacked bar.
  selected_x <- reactiveVal(NULL)
  
  default_x <- reactive({
    get_default_x(
      x_var = input$x_var,
      fixed_theta = input$fixed_theta,
      fixed_b = input$fixed_b,
      fixed_alpha_orig = input$fixed_alpha_orig,
      fixed_n_orig = input$fixed_n_orig
    )
  })
  
  selected_x_safe <- reactive({
    x <- selected_x()
    
    if (is.null(x) || is.na(x)) {
      x <- default_x()
    }
    
    clip_x_value(x, input$x_var)
  })
  
  # Reset the selected vertical line when key settings change.
  observeEvent(
    list(
      input$x_var,
      input$fixed_theta,
      input$fixed_b,
      input$fixed_alpha_orig,
      input$fixed_n_orig
    ),
    {
      selected_x(default_x())
    },
    ignoreInit = FALSE
  )
  
  curve_obj <- reactive({
    make_curve_data(
      x_var = input$x_var,
      fixed_theta = input$fixed_theta,
      fixed_b = input$fixed_b,
      fixed_alpha_orig = input$fixed_alpha_orig,
      fixed_n_orig = input$fixed_n_orig,
      n_rep = input$n_rep,
      N_rep = input$N_rep,
      alpha_meta = input$alpha_meta,
      I2 = input$I2
    )
  })
  
  current_ggplot <- reactive({
    make_curve_ggplot(
      curve_obj = curve_obj(),
      x_var = input$x_var
    )
  })
  
  current_plotly <- reactive({
    make_curve_plotly(
      curve_obj = curve_obj(),
      x_var = input$x_var,
      selected_x_value = selected_x_safe()
    )
  })
  
  
  
  # Click mode: click on the curve to move the vertical line.
  observeEvent(event_data("plotly_click", source = "curve_source"), {
    
    click_data <- event_data("plotly_click", source = "curve_source")
    
    if (!is.null(click_data$x) && !is.na(click_data$x[1])) {
      selected_x(clip_x_value(click_data$x[1], input$x_var))
    }
    
  }, ignoreInit = TRUE)
  
  # Convert the selected x-value into the actual model parameters.
  selected_params <- reactive({
    
    x <- selected_x_safe()
    
    theta_val <- input$fixed_theta
    b_val <- input$fixed_b
    alpha_orig_val <- input$fixed_alpha_orig
    n_orig_val <- input$fixed_n_orig
    
    if (input$x_var == "theta") {
      theta_val <- x
    }
    
    if (input$x_var == "b") {
      b_val <- x
    }
    
    if (input$x_var == "alpha_orig") {
      alpha_orig_val <- x
    }
    
    if (input$x_var == "n_orig") {
      n_orig_val <- x
    }
    
    list(
      theta = theta_val,
      b = b_val,
      alpha_orig = alpha_orig_val,
      n_orig = n_orig_val,
      n_rep = input$n_rep,
      N_rep = input$N_rep,
      alpha_meta = input$alpha_meta,
      I2 = input$I2
    )
  })
  
  selected_probs <- reactive({
    
    p <- selected_params()
    
    calc_outcome_probs(
      theta = p$theta,
      b = p$b,
      alpha_orig = p$alpha_orig,
      n_orig = p$n_orig,
      n_rep = p$n_rep,
      N_rep = p$N_rep,
      alpha_meta = p$alpha_meta,
      I2 = p$I2
    )
  })
  
  output$curve_plot <- renderPlotly({
    current_plotly()
  })
  
  output$cross_section_bar <- renderPlotly({
    make_cross_section_bar(
      probs = selected_probs(),
      selected_params = selected_params(),
      x_var = input$x_var,
      selected_x_value = selected_x_safe()
    )
  })
  
  output$selected_value_text <- renderUI({
    
    p <- selected_params()
    probs <- selected_probs()
    
    HTML(
      paste0(
        "<p style='font-size:14px;'>",
        "<b>Selected cross-section:</b> ",
        get_x_name(input$x_var),
        " = ",
        format_x_value(selected_x_safe(), input$x_var),
        "</p>",
        
        "<p style='font-size:14px;'>",
        "<b>Current model values:</b> ",
        "&theta; = ", round(p$theta, 4), "; ",
        "bias = ", round(p$b, 4), "; ",
        "&alpha;<sub>orig</sub> = ", round(p$alpha_orig, 4), "; ",
        "n<sub>orig</sub> = ", round(p$n_orig, 2), "; ",
        "n<sub>rep</sub> = ", p$n_rep, "; ",
        "N<sub>rep</sub> = ", p$N_rep, "; ",
        "&alpha;<sub>meta</sub> = ", round(p$alpha_meta, 4), "; ",
        "I<sup>2</sup> = ", round(p$I2, 4),
        "</p>",
        
        "<p style='font-size:14px;'>",
        "<b>Stacked-bar proportions:</b> ",
        "TS = ", percent(probs$TS, accuracy = 0.1), "; ",
        "FS = ", percent(probs$FS, accuracy = 0.1), "; ",
        "TF = ", percent(probs$TF, accuracy = 0.1), "; ",
        "FF = ", percent(probs$FF, accuracy = 0.1),
        "</p>"
      )
    )
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0(
        "replication_success_curve_",
        input$x_var,
        "_theta", input$fixed_theta,
        "_b", input$fixed_b,
        "_alphaorig", input$fixed_alpha_orig,
        "_norig", input$fixed_n_orig,
        "_nrep", input$n_rep,
        "_Nrep", input$N_rep,
        ".png"
      )
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = current_ggplot(),
        width = 10,
        height = 6.5,
        dpi = 300,
        bg = "white"
      )
    },
    contentType = "image/png"
  )
  
  output$interpretation <- renderUI({
    
    x_name <- get_x_name(input$x_var)
    
    if (input$x_var == "theta" && input$fixed_theta == 0) {
      
      theta_text <- paste0(
        "Because the selected fixed value of &theta; is 0, this plot shows only the four classification proportions at &theta; = 0. ",
        "The null-effect classification definitions TS<sub>2</sub>, FS<sub>2</sub>, TF<sub>2</sub>, and FF<sub>2</sub> are used."
      )
      
    } else if (input$x_var == "theta" && input$fixed_theta > 0) {
      
      theta_text <- paste0(
        "Because the selected fixed value of &theta; is greater than 0, this plot shows curves across &theta; &gt; 0. ",
        "The non-null-effect classification definitions TS<sub>1</sub>, FS<sub>1</sub>, TF<sub>1</sub>, and FF<sub>1</sub> are used."
      )
      
    } else {
      
      theta_text <- paste0(
        "When &theta; = 0, the null-effect classification definitions TS<sub>2</sub>, FS<sub>2</sub>, TF<sub>2</sub>, and FF<sub>2</sub> are used. ",
        "When &theta; &gt; 0, the non-null-effect classification definitions TS<sub>1</sub>, FS<sub>1</sub>, TF<sub>1</sub>, and FF<sub>1</sub> are used."
      )
    }
    
    HTML(
      paste0(
        "<p>",
        "This figure shows how the proportions of replication-success classifications change as a function of ",
        x_name,
        ". The other methodological factors are held fixed using the values selected in the sidebar.",
        "</p>",
        "<p>",
        "The vertical dotted line marks the selected x-axis value. The stacked bar on the right is the cross-sectional view of the curve at that selected value. ",
        "The cross-sectional view updates when the user clicks a point on the curve.",
        "</p>",
        "<p>",
        theta_text,
        " The legend uses the general labels TS, FS, TF, and FF to keep the graph simple.",
        "</p>"
      )
    )
  })
}

shinyApp(ui = ui, server = server)
#####
#The following code generate graphs based on the shiny app
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(scales)
# 
# # ============================================================
# # Analytic functions
# # ============================================================
# 
# clip01 <- function(x) {
#   pmin(pmax(x, 0), 1)
# }
# 
# # Same optional-stopping logic as the stacked-bar Shiny app
# calc_optional_stopping_prob <- function(theta,
#                                         n_orig,
#                                         alpha_orig,
#                                         q,
#                                         max_mult = 2,
#                                         step = 5) {
#   
#   n_max <- n_orig + q * (n_orig * max_mult - n_orig)
#   n_seq <- seq(n_orig, floor(n_max), by = step)
#   
#   zcrit <- qnorm(1 - alpha_orig)
#   
#   if (theta == 0) {
#     reject_prob_each_look <- rep(alpha_orig, length(n_seq))
#   } else {
#     reject_prob_each_look <- 1 - pnorm(
#       zcrit - theta * sqrt(n_seq / 2)
#     )
#   }
#   
#   pi_os <- 1 - prod(1 - reject_prob_each_look)
#   
#   clip01(pi_os)
# }
# 
# calc_outcome_probs <- function(theta,
#                                b,
#                                alpha_orig,
#                                n_orig = 50,
#                                n_rep = 100,
#                                N_rep = 5,
#                                alpha_meta = 0.05,
#                                I2 = 0.10) {
#   
#   if (theta == 0) {
#     
#     FP_orig_nominal <- alpha_orig
#     
#     pi0_OS <- calc_optional_stopping_prob(
#       theta = 0,
#       n_orig = n_orig,
#       alpha_orig = alpha_orig,
#       q = b
#     )
#     
#     FP_orig <- (1 - b) * FP_orig_nominal + b * pi0_OS
#     FP_orig <- clip01(FP_orig)
#     TN_orig <- 1 - FP_orig
#     
#     FP_rep <- alpha_meta
#     TN_rep <- 1 - FP_rep
#     
#     TS <- TN_orig * TN_rep
#     FS <- FP_orig * FP_rep
#     FF <- FP_orig * TN_rep
#     TF <- TN_orig * FP_rep
#     
#   } else {
#     
#     zcrit_orig <- qnorm(1 - alpha_orig)
#     
#     TP_orig_nominal <- 1 - pnorm(
#       zcrit_orig - theta * sqrt(n_orig / 2)
#     )
#     
#     pi1_OS <- calc_optional_stopping_prob(
#       theta = theta,
#       n_orig = n_orig,
#       alpha_orig = alpha_orig,
#       q = b
#     )
#     
#     TP_orig <- (1 - b) * TP_orig_nominal + b * pi1_OS
#     TP_orig <- clip01(TP_orig)
#     FN_orig <- 1 - TP_orig
#     
#     v <- 2 / n_rep + theta^2 / (4 * n_rep)
#     tau2 <- I2 * v / (1 - I2)
#     vbar <- (v + tau2) / N_rep
#     
#     lambda <- theta / sqrt(vbar)
#     zcrit_meta <- qnorm(1 - alpha_meta / 2)
#     
#     TP_rep <- 1 - pnorm(zcrit_meta - lambda) +
#       pnorm(-zcrit_meta - lambda)
#     
#     TP_rep <- clip01(TP_rep)
#     FN_rep <- 1 - TP_rep
#     
#     TS <- TP_orig * TP_rep
#     FS <- FN_orig * FN_rep
#     FF <- FN_orig * TP_rep
#     TF <- TP_orig * FN_rep
#   }
#   
#   total <- TS + FS + FF + TF
#   
#   tibble(
#     TS = TS / total,
#     FS = FS / total,
#     TF = TF / total,
#     FF = FF / total
#   )
# }
# 
# # ============================================================
# # Create plotting data
# # ============================================================
# 
# make_theta_data <- function() {
#   
#   theta_vals <- seq(0.001, 0.50, length.out = 500)
#   
#   design_grid <- expand.grid(
#     theta = theta_vals,
#     b = c(0.2, 0.5, 0.8),
#     alpha_orig = c(0.01, 0.05)
#   )
#   
#   plot_data <- lapply(seq_len(nrow(design_grid)), function(i) {
#     
#     theta_i <- design_grid$theta[i]
#     b_i <- design_grid$b[i]
#     alpha_i <- design_grid$alpha_orig[i]
#     
#     calc_outcome_probs(
#       theta = theta_i,
#       b = b_i,
#       alpha_orig = alpha_i,
#       n_orig = 50,
#       n_rep = 100,
#       N_rep = 5,
#       alpha_meta = 0.05,
#       I2 = 0.10
#     ) |>
#       mutate(
#         theta = theta_i,
#         b = b_i,
#         alpha_orig = alpha_i
#       )
#     
#   }) |>
#     bind_rows() |>
#     pivot_longer(
#       cols = c(TS, FS, TF, FF),
#       names_to = "category",
#       values_to = "proportion"
#     ) |>
#     mutate(
#       category = factor(
#         category,
#         levels = c("TS", "FS", "TF", "FF"),
#         labels = c("TS", "FS", "TF", "FF")
#       ),
#       bias_label = factor(
#         b,
#         levels = c(0.2, 0.5, 0.8),
#         labels = c("Bias = Low", "Bias = Medium", "Bias = High")
#       ),
#       alpha_label = factor(
#         alpha_orig,
#         levels = c(0.01, 0.05),
#         labels = c("alpha[orig] == .01", "alpha[orig] == .05")
#       )
#     )
#   
#   plot_data
# }
# 
# plot_data <- make_theta_data()
# 
# # ============================================================
# # Manuscript-style faceted plot
# # ============================================================
# 
# final_plot <- ggplot(
#   plot_data,
#   aes(
#     x = theta,
#     y = proportion,
#     color = category,
#     linetype = category,
#     group = category
#   )
# ) +
#   # Draw left and bottom axis lines inside each panel
#   geom_hline(
#     yintercept = 0,
#     color = "black",
#     linewidth = 0.45
#   ) +
#   geom_vline(
#     xintercept = 0,
#     color = "black",
#     linewidth = 0.45
#   ) +
#   geom_line(linewidth = 0.65) +
#   facet_grid(
#     rows = vars(alpha_label),
#     cols = vars(bias_label),
#     switch = "both",
#     labeller = labeller(
#       alpha_label = label_parsed,
#       bias_label = label_value
#     )
#   ) +
#   scale_y_continuous(
#     limits = c(0, 1),
#     breaks = seq(0, 1, by = 0.2),
#     labels = percent_format(accuracy = 1),
#     expand = expansion(mult = c(0, 0.02))
#   ) +
#   scale_x_continuous(
#     limits = c(0, 0.5),
#     breaks = seq(0, 0.5, by = 0.1),
#     expand = expansion(mult = c(0.01, 0.02))
#   ) +
#   scale_color_manual(
#     name = "Classification",
#     values = c(
#       "TS" = "black",
#       "FS" = "#1f77b4",
#       "TF" = "#2ca02c",
#       "FF" = "#d62728"
#     )
#   ) +
#   scale_linetype_manual(
#     name = "Classification",
#     values = c(
#       "TS" = "solid",
#       "FS" = "longdash",
#       "TF" = "dotted",
#       "FF" = "solid"
#     )
#   ) +
#   labs(
#     x = expression(theta),
#     y = "Proportion"
#   ) +
#   theme_bw(base_size = 12) +
#   theme(
#     axis.title.x = element_text(
#       size = 13,
#       margin = margin(t = 8)
#     ),
#     axis.title.y = element_text(
#       size = 13,
#       margin = margin(r = 8)
#     ),
#     axis.text = element_text(
#       size = 10,
#       color = "black"
#     ),
#     axis.ticks = element_line(
#       color = "black",
#       linewidth = 0.35
#     ),
#     
#     # Remove default panel border so only hline/vline show as axes
#     panel.border = element_blank(),
#     axis.line = element_blank(),
#     
#     panel.grid.major = element_line(
#       color = "grey85",
#       linewidth = 0.35
#     ),
#     panel.grid.minor = element_line(
#       color = "grey93",
#       linewidth = 0.25
#     ),
#     panel.spacing.x = unit(0.65, "lines"),
#     panel.spacing.y = unit(0.65, "lines"),
#     
#     strip.placement = "outside",
#     strip.background = element_blank(),
#     strip.text.x.bottom = element_text(
#       size = 12,
#       face = "bold",
#       margin = margin(t = 6)
#     ),
#     strip.text.y.left = element_text(
#       size = 12,
#       face = "bold",
#       angle = 90,
#       margin = margin(r = 6)
#     ),
#     
#     legend.position = "right",
#     legend.direction = "vertical",
#     legend.title = element_text(
#       size = 11,
#       face = "bold"
#     ),
#     legend.text = element_text(size = 10),
#     legend.key.width = unit(1.2, "cm"),
#     legend.key.height = unit(0.45, "cm"),
#     
#     plot.margin = margin(10, 12, 10, 10)
#   ) +
#   guides(
#     color = guide_legend(
#       ncol = 1,
#       byrow = TRUE,
#       override.aes = list(linewidth = 0.9)
#     ),
#     linetype = guide_legend(
#       ncol = 1,
#       byrow = TRUE,
#       override.aes = list(linewidth = 0.9)
#     )
#   )
# 
# final_plot
# 
# # ============================================================
# # Save figure
# # ============================================================
# 
# ggsave(
#   filename = "replication_success_theta_curves_faceted_by_alpha_bias.png",
#   plot = final_plot,
#   width = 11,
#   height = 6.5,
#   dpi = 300
# )
