#!/usr/bin/env Rscript
# Copy run logs, including logs from incomplete batches.
args <- commandArgs(TRUE)
if (!length(args) || any(args %in% c('--help', '-h'))) {
  cat('Usage: Rscript scripts/result_transfer/export_logs.R --source PATH [--destination PATH]\n',
      'Source defaults to the gamma-count assignment folder when present.\n',
      'Destination defaults to a sibling folder named <source>_logs.\n',
      'Recursively copies .log files from any source directory and refreshes existing copies.\n',
      'Exported filenames join the original path components with __.\n',
      sep='')
  quit(status=0L)
}

new_root <- 'results/assignments_100_v5_gamma_count'
opts <- list(source=if (dir.exists(new_root)) new_root else 'results/assignments_100_v5', destination=NULL)
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
relative_paths <- if (nzchar(prefix)) paste(prefix, logs, sep='/') else logs
export_names <- gsub('/', '__', relative_paths, fixed=TRUE)
name_keys <- tolower(export_names)
collision <- anyDuplicated(name_keys)
if (collision) {
  stop('Multiple logs would export as ', export_names[collision], ': ',
       paste(relative_paths[name_keys == name_keys[collision]], collapse=', '))
}
dir.create(opts$destination, recursive=TRUE, showWarnings=FALSE)
destination <- normalizePath(opts$destination, winslash='/', mustWork=TRUE)
if (!dir.exists(destination)) stop('Destination must be a directory.')
canonical <- function(path) if (.Platform$OS.type == 'windows') tolower(path) else path
inside <- function(path, parent) startsWith(canonical(path), paste0(sub('/$', '', canonical(parent)), '/'))
if (canonical(root) == canonical(destination) ||
    inside(destination, root) || inside(root, destination))
  stop('Choose a separate destination outside the source directory tree.')

targets <- file.path(destination, export_names)
for (i in seq_along(logs)) {
  if (!file.copy(file.path(root, logs[i]), targets[i], overwrite=TRUE, copy.date=TRUE))
    stop('Could not copy log: ', logs[i])
}
cat('Copied', length(logs), 'log files to', destination, '\n')
