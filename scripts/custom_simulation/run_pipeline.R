#!/usr/bin/env Rscript
# Generate, analyze and release 100 settings at a time.
a <- commandArgs(FALSE)
p <- normalizePath(sub('^--file=', '', grep('^--file=', a, value=TRUE)[1]))
setwd(normalizePath(file.path(dirname(p), '../..')))
source('scripts/custom_simulation/version.R')
source('scripts/custom_simulation/check_paths.R')
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
root <- file.path(paste0('results/custom_simulation_pipeline_v', simulation_version),
                  sprintf('shard-%03d-of-%03d', shard_index, shard_count))
check_pipeline_paths(root)
dir.create(root, recursive=TRUE, showWarnings=FALSE)
rscript <- file.path(R.home('bin'), 'Rscript')
run <- function(script, args) {
  code <- system2(rscript, c(shQuote(script), shQuote(args)))
  if (code != 0L) stop('Stopped: ', script, '. Restart with the same options to resume.')
}
# Persist arguments so a resumed pipeline cannot silently change its experiment.
config <- file.path(root, 'options.rds')
if (file.exists(config) && !identical(readRDS(config), args)) stop('Pipeline options differ from the existing run.')
if (!file.exists(config)) saveRDS(args, config)
index <- 1L
repeat {
  batch <- file.path(root, sprintf('batch-%06d', index))
  marker <- file.path(batch, 'complete.rds')
  if (file.exists(marker)) { index <- index + 1L; next }
  input <- file.path(batch, 'datasets')
  output <- file.path(batch, 'analysis')
  run('scripts/custom_simulation/run_batch.R', c(args, '--output', input,
      '--setting-batch', as.character(index), '--settings-per-batch', '100'))
  datasets <- list.files(input, pattern='-rep[0-9]+-n[0-9]+-seed[0-9]+-ref[[:xdigit:]]{8}\\.rds$', recursive=TRUE, full.names=TRUE)
  # Use the planner's full number of settings to distinguish the end from errors.
  source('scripts/custom_simulation/factorial_design.R')
  nvar <- 10000L
  pos <- match('--n-var', args); if (!is.na(pos)) nvar <- as.integer(args[pos+1L])
  model <- 'generalized_poisson'
  pos <- match('--underdispersion-model', args); if (!is.na(pos)) model <- args[pos+1L]
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
  run('scripts/custom_simulation/analyze_batch.R', c('--input', input, '--output', output, '--no-restore'))
  manifests <- list.files(output, pattern='^manifest.csv$', recursive=TRUE, full.names=TRUE)
  if (!length(manifests) || any(vapply(manifests, function(path) any(read.csv(path)$status != 'complete'), logical(1))))
    stop('Analysis incomplete in batch ', index, '; keeping datasets for retry.')
  files <- list.files(file.path(output, 'results'), pattern='\\.rds$', recursive=TRUE, full.names=TRUE)
  for (path in files) {
    result <- readRDS(path)
    result$gene_results <- NULL
    tmp <- paste0(path, '.compact.tmp')
    saveRDS(result, tmp, compress='xz')
    if (!file.rename(tmp, path)) stop('Cannot finalize compact result: ', path)
  }
  saveRDS(list(batch=index, settings_per_batch=100L, simulation_options=args), marker)
  unlink(datasets)
  index <- index + 1L
}
