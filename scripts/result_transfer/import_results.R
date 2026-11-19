#!/usr/bin/env Rscript
# Import a verified export created by export_results.R without altering its files.
a <- commandArgs(FALSE)
p <- normalizePath(sub('^--file=', '', grep('^--file=', a, value=TRUE)[1]))
setwd(normalizePath(file.path(dirname(p), '../..')))

args <- commandArgs(TRUE)
opts <- list(source=NULL, destination='results/exported_analysis_v5.1', layout='original')
if (!length(args) || '--help' %in% args) {
  cat('Usage: Rscript scripts/result_transfer/import_results.R --source "E:/RNA-results" [--destination PATH] [--layout original|short]\n',
      'Imports files listed in export_index.csv, verifies MD5 hashes, and preserves existing files.\n',
      'Files already present with matching hashes are left unchanged; conflicting files are preserved and reported.\n', sep='')
  quit(status=0L)
}
if (length(args) %% 2L) stop('Each option requires a value.')
for (i in seq.int(1L, length(args), by=2L)) {
  key <- sub('^--', '', args[i])
  if (!key %in% names(opts)) stop('Unknown option: ', args[i])
  opts[[key]] <- args[i+1L]
}
if (is.null(opts$source)) stop('--source is required and must contain export_index.csv.')
if (!opts$layout %in% c('original', 'short')) stop('--layout must be original or short.')
source_root <- normalizePath(opts$source, winslash='/', mustWork=TRUE)
if (!dir.exists(source_root)) stop('--source must be a directory.')
index_path <- file.path(source_root, 'export_index.csv')
if (!file.exists(index_path)) stop('Missing export_index.csv in source: ', source_root)
index_hash <- unname(tools::md5sum(index_path))
dir.create(opts$destination, recursive=TRUE, showWarnings=FALSE)
destination <- normalizePath(opts$destination, winslash='/', mustWork=TRUE)
canonical <- function(x) if (.Platform$OS.type == 'windows') tolower(x) else x
if (canonical(source_root) == canonical(destination) ||
    startsWith(canonical(destination), paste0(canonical(source_root), '/')) ||
    startsWith(canonical(source_root), paste0(canonical(destination), '/')))
  stop('Choose source and destination outside one another.')

index <- read.csv(index_path, stringsAsFactors=FALSE, check.names=FALSE)
required <- c('source_relative', 'export_relative', 'md5', 'type')
if (!all(required %in% names(index))) stop('export_index.csv is missing required columns: ', paste(setdiff(required, names(index)), collapse=', '))
if (!nrow(index)) stop('export_index.csv contains no files.')
if (anyNA(index[, required]) || any(!nzchar(index$source_relative)) || any(!nzchar(index$export_relative)))
  stop('Index contains empty or missing required values.')
safe_relative <- function(x) {
  x <- gsub('\\\\', '/', x)
  !grepl('^(~|/|[A-Za-z]:|//)', x) && !any(strsplit(x, '/', fixed=TRUE)[[1]] %in% c('', '.', '..'))
}
if (!all(vapply(index$source_relative, safe_relative, logical(1))) ||
    !all(vapply(index$export_relative, safe_relative, logical(1))))
  stop('Index contains an absolute path or path traversal.')
if (anyDuplicated(index$source_relative)) stop('Index contains duplicate destination paths.')
if (any(!grepl('^[[:xdigit:]]{32}$', index$md5))) stop('Index contains an invalid MD5 hash.')

copy_verified <- function(src, dst, expected) {
  actual <- unname(tools::md5sum(src))
  if (is.na(actual) || !identical(tolower(actual), tolower(expected)))
    stop('Source hash mismatch; export may be incomplete or changed: ', src)
  if (file.exists(dst)) {
    existing <- unname(tools::md5sum(dst))
    if (!is.na(existing) && identical(tolower(existing), tolower(expected))) return('unchanged')
    return('conflict')
  }
  dir.create(dirname(dst), recursive=TRUE, showWarnings=FALSE)
  if (!dir.exists(dirname(dst))) stop('Cannot create destination folder: ', dirname(dst))
  tmp <- paste0(dst, '.', Sys.getpid(), '.tmp')
  on.exit(unlink(tmp), add=TRUE)
  if (!file.copy(src, tmp, overwrite=FALSE) ||
      !identical(tolower(unname(tools::md5sum(tmp))), tolower(expected)) ||
      !identical(tolower(unname(tools::md5sum(src))), tolower(expected)))
    stop('Copy verification failed or source changed during import: ', src)
  if (!file.rename(tmp, dst)) stop('Cannot finalize copy: ', dst)
  'copied'
}

actions <- character(nrow(index))
destination_paths <- if (opts$layout == 'short') index$export_relative else index$source_relative
for (i in seq_len(nrow(index))) {
  src <- file.path(source_root, gsub('/', .Platform$file.sep, index$export_relative[i], fixed=TRUE))
  dst <- file.path(destination, gsub('/', .Platform$file.sep, destination_paths[i], fixed=TRUE))
  if (!file.exists(src) || dir.exists(src)) stop('Indexed source file is missing: ', src)
  actions[i] <- copy_verified(src, dst, index$md5[i])
  if (i %% 100L == 0L) { cat('Verified', i, 'of', nrow(index), 'files\n'); flush.console() }
}
cat('Import complete:', sum(actions == 'copied'), 'copied,', sum(actions == 'unchanged'),
    'already present,', sum(actions == 'conflict'), 'conflicts left untouched. Destination:', destination, '\n')
if (any(actions == 'conflict')) {
  cat('Conflicting destination files (not overwritten):\n')
  cat(paste0('  ', destination_paths[actions == 'conflict']), sep='\n')
  cat('\n')
}
if (opts$layout == 'short' && !any(actions == 'conflict')) {
  index_action <- copy_verified(index_path, file.path(destination, 'export_index.csv'), index_hash)
  if (index_action == 'conflict') stop('Destination export index differs. Use a fresh destination to retain the mapping.')
}
