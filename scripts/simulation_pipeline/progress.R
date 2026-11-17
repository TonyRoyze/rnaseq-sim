#!/usr/bin/env Rscript
# Read-only monitor; uses base R and never loads the analysis environment.
read_progress_csv <- function(path) {
  tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
           error = function(e) NULL)
}
latest_progress_rows <- function(paths, keys) {
  rows <- lapply(paths, function(path) {
    x <- read_progress_csv(path)
    if (is.null(x) || !all(c(keys, 'status') %in% names(x))) return(NULL)
    x[!duplicated(x[keys], fromLast = TRUE), , drop = FALSE]
  })
  # Keep files separate: manifests may have different schemas across versions.
  rows
}
progress_counts <- function(rows) {
  status <- unlist(lapply(rows, function(x) x$status), use.names = FALSE)
  c(complete = sum(status == 'complete', na.rm = TRUE),
    failed = sum(status == 'failed', na.rm = TRUE),
    running = sum(status == 'running', na.rm = TRUE))
}
progress_snapshot <- function(root) {
  workers <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  workers <- workers[grepl('^worker-[0-9]+$', basename(workers))]
  lapply(workers, function(worker) {
    batches <- list.dirs(worker, recursive = FALSE, full.names = TRUE)
    batches <- batches[grepl('^jobs-[0-9]+-[0-9]+$', basename(batches))]
    complete <- file.exists(file.path(batches, 'complete.rds'))
    ranges <- strsplit(sub('^jobs-', '', basename(batches)), '-', fixed = TRUE)
    sizes <- vapply(ranges, function(r) as.numeric(r[2]) - as.numeric(r[1]) + 1, numeric(1))
    files <- list.files(worker, recursive = TRUE, full.names = TRUE)
    generation <- latest_progress_rows(files[basename(files) == 'status.csv'], 'output_file')
    analysis <- latest_progress_rows(files[basename(files) == 'manifest.csv'],
                                     c('dataset_id', 'method', 'setting_id'))
    logs <- files[basename(files) %in% c('generation.log', 'analysis.log')]
    latest <- ''
    if (length(logs)) {
      log <- logs[which.max(file.info(logs)$mtime)]
      # Batches bound log size; only show the last nonempty line.
      lines <- tryCatch(readLines(log, warn = FALSE), error = function(e) character())
      lines <- lines[nzchar(trimws(lines))]
      latest <- paste(basename(dirname(log)), basename(log),
                      if (length(lines)) tail(lines, 1) else '(empty)', sep = ' | ')
    }
    list(worker = basename(worker), finished = sum(sizes[complete]),
         batches = sum(complete), generation = progress_counts(generation),
         analysis = progress_counts(analysis), latest = latest,
         unreadable = sum(vapply(c(generation, analysis), is.null, logical(1))))
  })
}
show_progress <- function(root, expected = NA_real_) {
  cat('\n', format(Sys.time(), '%Y-%m-%d %H:%M:%S %z'), '\n', root, '\n', sep = '')
  if (!dir.exists(root)) {
    cat('Output directory does not exist yet.\n'); return(invisible(NULL))
  }
  snapshot <- progress_snapshot(root)
  total <- sum(vapply(snapshot, function(x) x$finished, numeric(1)))
  if (is.finite(expected)) {
    cat(sprintf('Finished datasets (whole batches): %.0f / %.0f (%.1f%%)\n',
                total, expected, 100 * total / expected))
  } else cat(sprintf('Finished datasets (whole batches): %.0f\n', total))
  if (!length(snapshot)) cat('No worker folders found yet.\n')
  cat('Worker       Finished  Batches  Generated  GenFail  AnalysisOK  AnalysisFail\n')
  for (x in snapshot) {
    cat(sprintf('%-12s %8.0f %8d %10d %8d %11d %13d\n',
      x$worker, x$finished, x$batches, x$generation['complete'],
      x$generation['failed'], x$analysis['complete'], x$analysis['failed']))
    if (nzchar(x$latest)) cat('  Latest: ', x$latest, '\n', sep = '')
    if (x$unreadable) cat('  Note: ', x$unreadable, ' unreadable/incomplete status files; retry next refresh.\n', sep = '')
  }
  cat('Generated = latest status per dataset. Analysis counts are method/setting jobs,\n',
      'not datasets, and lag until manifest checkpoints. Failed counts can clear on retry.\n',
      'Finished uses complete.rds batch markers (generation only if Stage: generate).\n',
      'Latest log activity does not prove a worker is still running.\n', sep = '')
  invisible(snapshot)
}
progress_main <- function(args = commandArgs(TRUE)) {
  config <- 'config/machine.dcf'; output <- NULL; live <- FALSE; interval <- 15
  i <- 1L
  while (i <= length(args)) {
    arg <- args[i]
    if (arg %in% c('--help', '-h')) {
      cat('Usage: Rscript scripts/simulation_pipeline/progress.R [--config PATH | --output PATH] [--live] [--interval SECONDS]\n')
      return(invisible(NULL))
    }
    if (arg == '--live') { live <- TRUE; i <- i + 1L; next }
    if (!arg %in% c('--config', '--output', '--interval') || i == length(args))
      stop('Unknown or incomplete option: ', arg)
    value <- args[i + 1L]
    if (arg == '--config') config <- value
    if (arg == '--output') output <- value
    if (arg == '--interval') interval <- suppressWarnings(as.numeric(value))
    i <- i + 2L
  }
  if (length(interval) != 1 || !is.finite(interval) || interval < 1)
    stop('--interval must be at least 1 second.')
  expected <- NA_real_
  if (is.null(output)) {
    config_lines <- readLines(config, warn = FALSE)
    con <- textConnection(config_lines[!grepl('^\\s*#', config_lines)])
    cfg <- tryCatch(read.dcf(con), finally = close(con))
    output <- cfg[1, 'Output']
    expected <- as.numeric(cfg[1, 'Job-End']) - as.numeric(cfg[1, 'Job-Start']) + 1
  }
  output <- normalizePath(output, winslash = '/', mustWork = FALSE)
  repeat {
    show_progress(output, expected)
    if (!live) break
    cat(sprintf('Refreshing in %g seconds; Ctrl+C to stop.\n', interval))
    Sys.sleep(interval)
  }
}
if (sys.nframe() == 0L) progress_main()
