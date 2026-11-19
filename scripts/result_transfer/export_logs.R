#!/usr/bin/env Rscript
# Copy run logs, including logs from incomplete batches.
args <- commandArgs(TRUE)
if (!length(args) || any(args %in% c('--help', '-h'))) {
  cat('Usage: Rscript scripts/result_transfer/export_logs.R --source PATH [--destination PATH]\n',
      'Source defaults to results/assignments_100_v5 from the working directory.\n',
      'Destination defaults to a sibling folder named <source>_logs.\n',
      'Recursively copies .log files from any source directory and refreshes existing copies.\n',
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
root <- normalizePath(opts$source, winslash='/', mustWork=TRUE)
if (!dir.exists(root)) stop('Source must be a directory.')
if (is.null(opts$destination)) {
  if (dirname(root) == root) stop('Specify --destination when the source is a filesystem root.')
  opts$destination <- file.path(dirname(root), paste0(basename(root), '_logs'))
}

logs <- list.files(root, pattern='[.]log$', recursive=TRUE, all.files=TRUE, ignore.case=TRUE)
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
