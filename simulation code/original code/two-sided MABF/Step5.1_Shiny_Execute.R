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

# --------- Configure apps ----------
apps <- data.frame(
  id    = c("roc","heatmap","stacked1","stacked_fp_tp","stacked2"),
  title = c(
    "ROC curve",
    "Heatmap (model metrics)",
    "Stacked bars: metrics: TS, FS",
    "Stacked bars: metrics: FP, TN",
    "Stacked bars: evidence strength"
  ),
  file  = c(
    "Step5.1_Shiny_ROC.R",
    "Step5.1_Shiny_HeatMap.R",
    "Step5.1_Shiny_StackedBar.R",
    "Step5.1_Shiny_StackedBar.R",   # <- new FP/TP tile
    "Step5.1_Shiny_StackedBar.R"
  ),
  from  = c(1, 1, 1, 176, 529),
  to    = c(119, 103, 168, 321, 616),
  preload = c(
    "tsfs",      # roc
    "tsfs",      # heatmap
    "shiny_only",# stacked1 (metrics: TS, FS)
    "shiny_only",# stacked2 (metrics: FP, TN) — same as stacked1
    "rowwise"    # stacked3 (evidence strength) – needs row2matrix.RData
  ),
  stringsAsFactors = FALSE
)



# --------- Background launcher ----------
launch_app_bg <- function(file, from, to, working_dir = getwd(), preload_key = "shiny_only") {
  port <- httpuv::randomPort()
  url  <- sprintf("http://127.0.0.1:%d", port)
  
  proc <- callr::r_bg(
    function(file, from, to, working_dir, port, preload_key) {
      setwd(working_dir)
      env <- new.env(parent = globalenv())
      
      # child-side helpers
      load_all_into_env <- function(dirs, env) {
        for (d in dirs) {
          if (!dir.exists(d)) next
          rdata_files <- list.files(d, pattern = "\\.RData$", full.names = TRUE, recursive = TRUE)
          for (f in rdata_files) try(load(f, envir = env), silent = TRUE)
          rds_files <- list.files(d, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
          for (f in rds_files) {
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
      
      preload_registry <- list(
        shiny_only = function(env, root) { suppressPackageStartupMessages(library(shiny)); invisible(TRUE) },
        tsfs = function(env, root) { suppressPackageStartupMessages({ library(shiny); library(ggplot2); library(scales) }); dirs <- file.path(root, c("MABF4ROC/4TSFS", "MABF4ROC")); load_all_into_env(dirs, env) },
        rates4viz = function(env, root) { suppressPackageStartupMessages({ library(shiny); library(ggplot2) }); dirs <- file.path(root, c("matrix-wise/rates4Viz", "rates4Viz", "matrix-wise")); load_all_into_env(dirs, env) },
        rowwise = function(env, root) { suppressPackageStartupMessages({ library(shiny); library(ggplot2); library(dplyr) }); dirs <- file.path(root, c("row-wise", "MABFanalyses/row-wise")); load_all_into_env(dirs, env) }
      )
      
      # run preload (libs + directory scans)
      pre <- preload_registry[[preload_key]]
      if (is.null(pre)) pre <- preload_registry[["shiny_only"]]
      try(pre(env, getwd()), silent = TRUE)
      
      # source the app lines
      source_lines <- function(file, from, to, local = .GlobalEnv, ...) {
        lines <- readLines(file, warn = FALSE)
        txt <- paste(lines[from:to], collapse = "\n")
        source(textConnection(txt), local = local, ...)
      }
      source_lines(file, from, to, local = env)
      
      if (!exists("ui", envir = env, inherits = FALSE) || !exists("server", envir = env, inherits = FALSE)) {
        stop("The sourced code did not define objects named 'ui' and 'server'.")
      }
      
      shiny::runApp(
        shiny::shinyApp(ui = env$ui, server = env$server),
        port = port, host = "127.0.0.1", launch.browser = FALSE, quiet = TRUE
      )
    },
    args = list(file = file, from = from, to = to,
                working_dir = working_dir, port = port, preload_key = preload_key)
  )
  
  # retry-open browser a few times so you rarely need to refresh
  for (i in 1:25) { # ~5s total
    Sys.sleep(0.2)
    if (!proc$is_alive()) break
    tc <- try(socketConnection("127.0.0.1", port, open = "r", timeout = 0.2), silent = TRUE)
    if (!inherits(tc, "try-error")) {
      try(close(tc), silent = TRUE)
      try(utils::browseURL(url), silent = TRUE)
      break
    }
  }
  
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


