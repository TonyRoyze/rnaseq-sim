#!/usr/bin/env Rscript
# Generate, analyze and release 100 settings at a time.
a <- commandArgs(FALSE)
p <- normalizePath(sub('^--file=', '', grep('^--file=', a, value=TRUE)[1]))
setwd(normalizePath(file.path(dirname(p), '../..')))
source('scripts/shared/version.R')
source('scripts/shared/check_paths.R')
source('scripts/shared/pipeline_logging.R')
source('scripts/simulation_pipeline/resume.R')
args <- commandArgs(TRUE)
if ('--version' %in% args) {
  cat('Custom simulation pipeline version', simulation_version, '\n')
  quit(status=0L)
}
if (any(args %in% c('--output', '--setting-batch', '--settings-per-batch', '--job-index', '--plan-only'))) {
  stop('Pipeline owns output and selection options. Pass simulation options such as --replicates 100 and --no-restore.')
}
option_value <- function(name, default) {
  at <- match(name, args)
  if (is.na(at)) default else {
    if (at == length(args)) stop('Missing value for ', name)
    args[at + 1L]
  }
}
shard_index <- as.integer(option_value('--shard-index', '1'))
shard_count <- as.integer(option_value('--shard-count', '1'))
if (anyNA(c(shard_index, shard_count)) || shard_index < 1L || shard_index > shard_count)
  stop('Invalid shard settings.')
model <- option_value('--underdispersion-model', 'gamma_count')
root <- file.path(paste0('results/custom_simulation_pipeline_v', simulation_version,
                         if (model == 'gamma_count') '_gamma_count' else ''),
                  sprintf('shard-%03d-of-%03d', shard_index, shard_count))
check_pipeline_paths(root)
dir.create(root, recursive=TRUE, showWarnings=FALSE)
log_dir <- file.path(root, 'logs')
dir.create(log_dir, recursive=TRUE, showWarnings=FALSE)
log_path <- normalizePath(file.path(log_dir, paste0('pipeline-', format(Sys.time(), '%Y-%m-%d_%H-%M-%S'), '-', Sys.getpid(), '.log')), winslash='/', mustWork=FALSE)
Sys.setenv(RNA_PIPELINE_LOG=log_path)
start_pipeline_log(log_path)
pipeline_status('START', paste('Version', simulation_version, '| worker', shard_index, 'of', shard_count))
pipeline_status('LOG', log_path)
rscript <- file.path(R.home('bin'), 'Rscript')
run <- function(script, args) {
  pipeline_status('PROCESS', basename(script))
  started_at <- proc.time()[["elapsed"]]
  code <- system2(rscript, c(shQuote(script), shQuote(args)))
  elapsed <- proc.time()[["elapsed"]] - started_at
  pipeline_status('PROCESS COMPLETE', sprintf('%s | elapsed %.3f sec', basename(script), elapsed))
  if (code != 0L) stop('Stopped: ', script, '. Restart with the same options to resume.')
}
# Persist arguments so a resumed pipeline cannot silently change its experiment.
config <- file.path(root, 'options.rds')
if (file.exists(config)) {
  with_model <- function(values, default) {
    pos <- match('--underdispersion-model', values)
    if (is.na(pos)) return(c(values, '--underdispersion-model', default))
    flag <- c(pos, pos + 1L)
    c(values[-flag], values[flag])
  }
  # Before the gamma-count output folder existed, an omitted flag meant GP.
  saved_default <- if (model == 'gamma_count') 'gamma_count' else 'generalized_poisson'
  if (!identical(with_model(readRDS(config), saved_default), with_model(args, 'gamma_count')))
    stop('Pipeline options differ from the existing run.')
}
if (!file.exists(config)) saveRDS(args, config)
index <- 1L
skipped_batches <- 0L
incomplete <- 0L
repeat {
  batch <- file.path(root, sprintf('batch-%06d', index))
  marker <- file.path(batch, 'complete.rds')
  if (file.exists(marker)) { skipped_batches <- skipped_batches + 1L; index <- index + 1L; next }
  if (skipped_batches > 0L) {
    pipeline_status("RESUME", paste(skipped_batches, "completed batches skipped"))
    skipped_batches <- 0L
  }
  cat("\n\n", strrep("=", 72), "\n", sep="")
  pipeline_status("BATCH", sprintf("%06d | worker %d of %d", index, shard_index, shard_count))
  cat(strrep("=", 72), "\n\n", sep="")
  input <- file.path(batch, 'datasets')
  output <- file.path(batch, 'analysis')
  pipeline_status('GENERATE', 'Generate missing datasets; reuse existing datasets')
  run('scripts/simulation_counts/run_batch.R', c(args, '--output', input,
      '--setting-batch', as.character(index), '--settings-per-batch', '100'))
  datasets <- list.files(input, pattern='-rep[0-9]+-n[0-9]+-seed[0-9]+-ref[[:xdigit:]]{8}\\.rds$', recursive=TRUE, full.names=TRUE)
  # Use the planner's full number of settings to distinguish the end from errors.
  source('scripts/simulation_design/factorial_design.R')
  nvar <- 10000L
  pos <- match('--n-var', args); if (!is.na(pos)) nvar <- as.integer(args[pos+1L])
  nde <- 7L
  pos <- match('--de-proportions', args); if (!is.na(pos)) nde <- length(strsplit(args[pos+1L], ',', fixed=TRUE)[[1]])
  pos <- match('--de-mode', args); if (!is.na(pos) && !args[pos+1L] %in% c('3','fixed')) nde <- 1L
  total <- nrow(build_scenarios(nvar, !'--no-imbalance' %in% args, model)) * nde
  if ((index-1L)*100 >= total) break
  if (!length(datasets)) stop('No datasets generated in batch ', index, '; inspect simulation errors.')
  status <- list.files(input, pattern='^status.csv$', recursive=TRUE, full.names=TRUE)
  for (path in status) {
    rows <- read.csv(path)
    rows <- rows[!duplicated(rows$output_file, fromLast=TRUE), ]
    if (any(rows$status != 'complete')) stop('Simulation failures in batch ', index, '; keeping datasets for retry.')
  }
  pipeline_status('ANALYZE', paste(length(datasets), 'datasets; run missing method results'))
  run('scripts/de_analysis/analyze_batch.R', c('--input', input, '--output', output, '--no-restore'))
  manifests <- list.files(output, pattern='^manifest.csv$', recursive=TRUE, full.names=TRUE)
  if (!length(manifests) || !all(vapply(manifests, function(path) {
    analysis_manifest_complete(read.csv(path))
  }, logical(1)))) {
    incomplete <- incomplete + 1L
    pipeline_status('ANALYSIS RETRY NEEDED', batch)
    index <- index + 1L
    next
  }
  files <- list.files(file.path(output, 'results'), pattern='\\.rds$', recursive=TRUE, full.names=TRUE)
  pipeline_status('COMPACT', paste(length(files), 'analysis result files'))
  for (path in files) {
    result <- readRDS(path)
    result$gene_results <- NULL
    tmp <- paste0(path, '.compact.tmp')
    saveRDS(result, tmp, compress='xz')
    if (!file.rename(tmp, path)) stop('Cannot finalize compact result: ', path)
  }
  saveRDS(list(batch=index, settings_per_batch=100L, simulation_options=args), marker)
  pipeline_status('CLEANUP', paste('Remove', length(datasets), 'analyzed dataset files'))
  unlink(datasets)
  pipeline_status('BATCH COMPLETE', sprintf('%06d', index))
  index <- index + 1L
}

if (skipped_batches > 0L) pipeline_status('RESUME', paste(skipped_batches, 'completed batches skipped'))
if (incomplete > 0L) {
  pipeline_status('RETRY NEEDED', sprintf('%d analysis batches need retry', incomplete))
  quit(status = 1L)
}
pipeline_status('FINISHED', paste('Worker', shard_index, 'finished all assigned batches'))
