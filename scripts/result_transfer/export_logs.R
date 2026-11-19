#!/usr/bin/env Rscript
# Copy run logs, including logs from incomplete batches.
args <- commandArgs(TRUE)
if (!length(args) || any(args %in% c('--help', '-h'))) {
  cat('Usage: Rscript scripts/result_transfer/export_logs.R --destination PATH [--source PATH]\n',
      'Source defaults to results/assignments_100_v5 from the working directory.\n',
      'Pass the assignment root or a single config folder. Existing log copies are refreshed.\n',
      sep='')
  quit(status=0L)
}

opts <- list(source='results/assignments_100_v5', destination=NULL)
if (length(args) %% 2L) stop('Each option requires a value.')
for (i in seq.int(1L, length(args), by=2L)) {
  if (!args[i] %in% c('--source', '--destination') ||
      !nzchar(args[i+1L]) || startsWith(args[i+1L], '--'))
    stop('Unknown option or missing value: ', args[i])
  opts[[substring(args[i], 3L)]] <- args[i+1L]
}
if (is.null(opts$destination)) stop('--destination is required.')
root <- normalizePath(opts$source, winslash='/', mustWork=TRUE)
if (!dir.exists(root)) stop('Source must be a directory.')

logs <- list.files(root, pattern='[.]log$', recursive=TRUE)
if (!length(logs)) stop('No .log files found in ', root)
prefix <- if (grepl('^config-[0-9]+(-.*)?$', basename(root))) basename(root) else ''
dir.create(opts$destination, recursive=TRUE, showWarnings=FALSE)
destination <- normalizePath(opts$destination, winslash='/', mustWork=TRUE)
if (!dir.exists(destination)) stop('Destination must be a directory.')
canonical <- function(path) if (.Platform$OS.type == 'windows') tolower(path) else path
inside <- function(path, parent) startsWith(canonical(path), paste0(sub('/$', '', canonical(parent)), '/'))
if (canonical(root) == canonical(destination) ||
    inside(destination, root) || inside(root, destination))
  stop('Choose a separate destination outside the source directory tree.')

targets <- file.path(destination, logs)
if (nzchar(prefix)) targets <- file.path(destination, prefix, logs)
for (folder in unique(dirname(targets))) {
  dir.create(folder, recursive=TRUE, showWarnings=FALSE)
  if (!dir.exists(folder)) stop('Could not create destination directory: ', folder)
}
for (i in seq_along(logs)) {
  if (!file.copy(file.path(root, logs[i]), targets[i], overwrite=TRUE, copy.date=TRUE))
    stop('Could not copy log: ', logs[i])
}
cat('Copied', length(logs), 'log files to', destination, '\n')
