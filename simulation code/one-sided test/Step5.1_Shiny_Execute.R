# All my Shiny data-viz apps, together in one place.

ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}
invisible(lapply(c("shiny","bslib","callr","httpuv"), ensure_pkg))

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(callr)
  library(httpuv)
})

# --------- Helper source_lines ----------
source_lines <- function(file, from, to, local = .GlobalEnv, ...) {
  lines <- readLines(file, warn = FALSE)
  txt <- paste(lines[from:to], collapse = "\n")
  source(textConnection(txt), local = local, ...)
}

# --------- Data preload helpers ----------
# These are passed to the child via callr arguments, then redefined there.

preload_registry <- list(
  shiny_only = function(env, root) {
    suppressPackageStartupMessages({
      library(shiny)
    })
    invisible(TRUE)
  },
  tsfs = function(env, root) {
    # For ROC / Donut / Heatmap (model metrics)
    suppressPackageStartupMessages({
      library(shiny); library(ggplot2); library(scales)
    })
    # Try common TS/FS data locations
    dirs <- file.path(root, c("MABF4ROC/4TSFS"))
    load_all_into_env(dirs, env)
  },
  rates4viz = function(env, root) {
    # For Replication curve apps
    suppressPackageStartupMessages({
      library(shiny); library(ggplot2)
    })
    dirs <- file.path(root, c("matrix-wise/rates4Viz"))
    load_all_into_env(dirs, env)
  },
  rowwise = function(env, root) {
    # For Beta regression (row-wise data)
    suppressPackageStartupMessages({
      library(shiny); library(ggplot2); library(dplyr)
    })
    dirs <- file.path(root, c("row-wise"))
    load_all_into_env(dirs, env)
  }
)

# This helper will be defined again INSIDE the child
# It loads all .RData via base::load() and .rds via readRDS() into env
# (silently skips missing paths)
# NOTE: we put the definition in launch_app_bg so it exists in the child.
NULL

# --------- Configure apps (and which preload they need) ----------
apps <- data.frame(
  id    = c("beta","roc","repgrid","rep","heatmap","origplots","pie","donut",
            "stacked1","stacked_fp_tp","stacked2","lineplot"),   # new id
  title = c(
    "Beta regression model",
    "ROC curve",
    "Replication curve (Grid, v1)",
    "Replication curve (Grid, v2)",
    "Heatmap (model metrics)",
    "Original study: hist/heat/bar",
    "Pie chart (model metrics)",
    "Donut: FP/TP vs FS/TS",
    "Stacked bars: metrics: TS, FS",
    "Stacked bars: metrics: FP, TN",
    "Stacked bars: evidence strength",
    "Line plots (rates across design factors)"   # new title
  ),
  file  = c(
    "Step5.1_Shiny_BetaRegModel.R",
    "Step5.1_Shiny_ROC.R",
    "Step5.1_Shiny_RepCurve_GridPlot.R",
    "Step5.1_Shiny_RepCurve_GridPlot.R",
    "Step5.1_Shiny_HeatMap.R",
    "Step5.1_Shiny_origp_Histogram&HeatMap&Bar.R",
    "Step5.1_Shiny_PieChart.R",
    "Step5.1_Shiny_RepRates_Donut.R",
    "Step5.1_Shiny_StackedBar.R",
    "Step5.1_Shiny_StackedBar.R",   # <- FP/TP
    "Step5.1_Shiny_StackedBar.R",
    "Step5.1_Shiny_LinePlot.R"      # <- new file
  ),
  from  = c(  1,   1, 2481,   3145,   1,   1,   1,   1,   1,  176,  529, 1), # <- new from
  to    = c(733, 119, 2802, 3626, 125, 275, 322, 295, 268,  321,  616, 265), # <- new to
  preload = c(
    "rowwise",   # beta
    "tsfs",      # roc
    "rates4viz", # repgrid
    "rates4viz", # rep
    "tsfs",      # heatmap
    "shiny_only",# original study plots
    "shiny_only",# pie
    "tsfs",      # donut
    "shiny_only",# stacked1
    "shiny_only",# stacked2
    "rowwise",   # stacked3
    "rowwise"    # <- new app
  ),
  stringsAsFactors = FALSE
)


# Special case: Pie app needs a second chunk
pie_extra <- list(file = "Step5.1_Shiny_PieChart.R", from = 650, to = 805)

# --------- Background launcher ----------
launch_app_bg <- function(file, from, to, working_dir = getwd(), pie_extra = NULL, preload_key = "shiny_only") {
  port <- httpuv::randomPort()
  url  <- sprintf("http://127.0.0.1:%d", port)
  
  # Select preload function by key (falls back to shiny_only)
  preload_fun <- preload_registry[[preload_key]]
  if (is.null(preload_fun)) preload_fun <- preload_registry[["shiny_only"]]
  
  proc <- callr::r_bg(
    function(file, from, to, working_dir, pie_extra, port, preload_key) {
      # Child process
      setwd(working_dir)
      
      # Define loader helpers INSIDE child
      load_all_into_env <- function(dirs, env) {
        for (d in dirs) {
          if (!dir.exists(d)) next
          # .RData
          rdata_files <- list.files(d, pattern = "\\.RData$", full.names = TRUE, recursive = TRUE)
          for (f in rdata_files) {
            try(load(f, envir = env), silent = TRUE)
          }
          # .rds
          rds_files <- list.files(d, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
          for (f in rds_files) {
            obj <- NULL
            ok <- TRUE
            obj <- tryCatch(readRDS(f), error = function(e) { ok <<- FALSE; NULL })
            if (ok && !is.null(obj)) {
              nm <- tools::file_path_sans_ext(basename(f))
              assign(nm, obj, envir = env)
            }
          }
        }
        invisible(TRUE)
      }
      
      # Recreate the same preload_registry in child
      preload_registry <- list(
        shiny_only = function(env, root) {
          suppressPackageStartupMessages({ library(shiny) })
          invisible(TRUE)
        },
        tsfs = function(env, root) {
          suppressPackageStartupMessages({ library(shiny); library(ggplot2); library(scales) })
          dirs <- file.path(root, c("MABF4ROC/4TSFS"))
          load_all_into_env(dirs, env)
        },
        rates4viz = function(env, root) {
          suppressPackageStartupMessages({ library(shiny); library(ggplot2) })
          dirs <- file.path(root, c("matrix-wise/rates4Viz"))
          load_all_into_env(dirs, env)
        },
        rowwise = function(env, root) {
          suppressPackageStartupMessages({ library(shiny); library(ggplot2); library(dplyr) })
          dirs <- file.path(root, c("row-wise"))
          load_all_into_env(dirs, env)
        }
      )
      
      # Run preload
      env <- new.env(parent = globalenv())
      root <- getwd()
      pre <- preload_registry[[preload_key]]
      if (is.null(pre)) pre <- preload_registry[["shiny_only"]]
      try(pre(env, root), silent = TRUE)  # best-effort preload
      
      # Now source the app chunks into env
      source_lines <- function(file, from, to, local = .GlobalEnv, ...) {
        lines <- readLines(file, warn = FALSE)
        txt <- paste(lines[from:to], collapse = "\n")
        source(textConnection(txt), local = local, ...)
      }
      source_lines(file, from, to, local = env)
      if (!is.null(pie_extra) && basename(file) == basename(pie_extra$file)) {
        source_lines(pie_extra$file, pie_extra$from, pie_extra$to, local = env)
      }
      
      if (!exists("ui", envir = env, inherits = FALSE) || !exists("server", envir = env, inherits = FALSE)) {
        stop("The sourced code did not define objects named 'ui' and 'server'.")
      }
      
      shiny::runApp(
        shiny::shinyApp(ui = env$ui, server = env$server),
        port = port, host = "127.0.0.1", launch.browser = FALSE, quiet = TRUE
      )
    },
    args = list(file = file, from = from, to = to,
                working_dir = working_dir, pie_extra = pie_extra,
                port = port, preload_key = preload_key)
  )
  
  # Open from parent (more reliable on some systems)
  try(utils::browseURL(url), silent = TRUE)
  
  list(proc = proc, url = url)
}

# --------- UI ---------
theme <- bslib::bs_theme(
  bootswatch = "flatly",
  base_font  = bslib::font_google("Inter")
)

ui <- fluidPage(
  theme = theme,
  tags$head(tags$style(HTML("
    .app-card { border: 1px solid #e6e6e6; border-radius: 14px; padding: 16px; margin-bottom: 12px; }
    .app-title { font-weight: 600; font-size: 1.05rem; margin-bottom: 6px; }
    .app-note { color: #666; font-size: 0.9rem; margin-bottom: 10px; }
    .btn-launch { border-radius: 999px; margin-right: 6px; }
    .status-pill { display:inline-block; padding:2px 8px; border-radius:999px; font-size:0.8rem; margin-left:8px; }
    .running { background:#e8f7ee; color:#166534; border:1px solid #bbf7d0; }
    .idle { background:#f3f4f6; color:#374151; border:1px solid #e5e7eb; }
  "))),
  titlePanel("Shiny App Hub"),
  p("Click a button to launch an app in a new tab. The hub stays open so you can launch others."),
  fluidRow(column(width = 12, uiOutput("cards")))
)

# --------- Server ---------
server <- function(input, output, session) {
  procs <- reactiveValues()  # each is list(proc=, url=)
  
  output$cards <- renderUI({
    procs_list <- reactiveValuesToList(procs)
    tagList(lapply(seq_len(nrow(apps)), function(i) {
      a <- apps[i, ]
      
      running <- FALSE
      url_txt <- NULL
      if (!is.null(procs_list[[a$id]])) {
        x <- procs_list[[a$id]]
        if (!is.null(x$proc) && x$proc$is_alive()) {
          running <- TRUE; url_txt <- x$url
        }
      }
      
      pill <- span(class = paste("status-pill", if (running) "running" else "idle"),
                   if (running) "running" else "idle")
      
      launch_id <- paste0("launch_", a$id)
      stop_id   <- paste0("stop_", a$id)
      
      div(class = "app-card",
          div(class = "app-title", a$title, pill),
          div(class = "app-note",
              tags$code(a$file), "  lines ", a$from, "–", a$to,
              " · preload:", a$preload),
          div(
            actionButton(launch_id, "Launch", class = "btn btn-primary btn-launch"),
            actionButton(stop_id,   "Stop",   class = "btn btn-outline-secondary btn-launch",
                         title = "Stop the background R process if running"),
            if (running && !is.null(url_txt)) tags$a(href = url_txt, target = "_blank", "Open app")
          )
      )
    }))
  })
  
  lapply(seq_len(nrow(apps)), function(i) {
    a <- apps[i, ]
    launch_id <- paste0("launch_", a$id)
    stop_id   <- paste0("stop_", a$id)
    
    observeEvent(input[[launch_id]], {
      cur <- isolate(procs[[a$id]])
      if (!is.null(cur) && !is.null(cur$proc) && cur$proc$is_alive()) {
        showNotification(paste(a$title, "is already running."), type = "message")
        return(invisible(NULL))
      }
      showNotification(paste("Starting", a$title, "…"), type = "message")
      res <- tryCatch(
        launch_app_bg(a$file, a$from, a$to,
                      working_dir = getwd(),
                      pie_extra = pie_extra,
                      preload_key = a$preload),
        error = function(e) {
          showNotification(paste("Failed to start", a$title, ":", e$message),
                           type = "error", duration = 8)
          NULL
        }
      )
      if (!is.null(res)) {
        procs[[a$id]] <- res
        showNotification(paste(a$title, "launched at", res$url),
                         type = "message", duration = 6)
      }
    })
    
    observeEvent(input[[stop_id]], {
      cur <- isolate(procs[[a$id]])
      if (!is.null(cur) && !is.null(cur$proc)) {
        try(cur$proc$kill(), silent = TRUE)
        procs[[a$id]] <- NULL
        showNotification(paste("Stopped", a$title), type = "message")
      } else {
        showNotification("Nothing to stop.", type = "warning")
      }
    })
  })
  
  session$onSessionEnded(function() {
    procs_list <- reactiveValuesToList(procs)
    lapply(procs_list, function(x) {
      if (!is.null(x) && !is.null(x$proc)) try(x$proc$kill(), silent = TRUE)
    })
  })
}

shinyApp(ui, server)




# # Function controlling which lines to run
# source_lines <- function(file, from, to, local = .GlobalEnv, ...) {
#   lines <- readLines(file, warn = FALSE)
#   txt <- paste(lines[from:to], collapse = "\n")
#   source(textConnection(txt), local = local, ...)
# }
# 
# 
# # Beta regression model and visualization
# # Need row-wise data
# source_lines("Step5.1_Shiny_BetaRegModel.R", from = 1, to = 733)
# shinyApp(ui, server)
# 
# # ROC curve
# # Need matrix-wise/rates4ROC data (done)
# source_lines("Step5.1_Shiny_ROC.R", from = 1, to = 119)
# shinyApp(ui, server)
# 
# # Replication curve (Grid plot ver.)
# # Need matrix-wise/rates4Viz data (done)
# source_lines("Step5.1_Shiny_RepCurve_GridPlot.R", from = 1856, to = 2133)
# shinyApp(ui, server)
# 
# # Replication curve
# # Need matrix-wise/rates4Viz data (done)
# source_lines("Step5.1_Shiny_RepCurve.R", from = 1, to = 203)
# shinyApp(ui, server)
# 
# # Heatmap for model metrics
# # Need MABF4ROC/4TSFS data (done)
# source_lines("Step5.1_Shiny_HeatMap.R", from = 1, to = 103)
# shinyApp(ui, server)
# 
# # Histogram, heatmap, and bar plot for original study data
# # Need OGDG/10000 data (done)
# source_lines("Step5.1_Shiny_origp_Histogram&HeatMap&Bar.R", from = 1, to = 275)
# shinyApp(ui, server)
# 
# # Pie chart for model metrics
# # Need matrix-wise/rates4Plot/fixed original cutoff data (done)
# source_lines("Step5.1_Shiny_PieChart.R", from = 1, to = 322)
# source_lines("Step5.1_Shiny_PieChart.R", from = 650, to = 805)
# shinyApp(ui, server)
# 
# # Donut chart showing relationship between FP/TP and FS/TS
# # Need MABF4ROC/4TSFS data (done)
# source_lines("Step5.1_Shiny_RepRates_Donut.R", from = 1, to = 295)
# shinyApp(ui, server)
# 
# # Stacked bar plot showing relationship between model metrics (FS/TS)
# # Need matrix-wise/rates4Plot/fixed original cutoff data (done)
# source_lines("Step5.1_Shiny_StackedBar.R", from = 1, to = 168)
# shinyApp(ui, server)

# # Stacked bar plot showing relationship between model metrics (FP/TP)
# # Need matrix-wise/rates4Plot/fixed original cutoff data (done)
# source_lines("Step5.1_Shiny_StackedBar.R", from = 177, to = 318)
# shinyApp(ui, server)

# 
# # Stacked bar plot showing relationship between evidence strength
# # Need matrix-wise/rates4Plot/fixed original cutoff data (done)
# source_lines("Step5.1_Shiny_StackedBar.R", from = 381, to = 605)
# shinyApp(ui, server)

# # Line plot showing a metrics at different conditions of an independent variable (e.g., bias.level)
# # Need row-wise data
# source_lines("Step5.1_Shiny_LinePlot.R", from = 1, to = 265)
# shinyApp(ui, server)
