#!/usr/bin/env Rscript
# Reorganize existing simulation datasets without regenerating counts.
a <- commandArgs(FALSE)
p <- normalizePath(sub('^--file=', '', grep('^--file=', a, value=TRUE)[1]))
project_root <- normalizePath(file.path(dirname(p), '../..'))
setwd(project_root)
source('scripts/custom_simulation/version.R')
source('scripts/custom_simulation/factorial_design.R')
source('scripts/custom_simulation/batch_jobs.R')
argv <- commandArgs(TRUE)
separator <- match('--', argv)
if (is.na(separator)) stop('Usage: migrate_pipeline.R --source PATH [--destination PATH] [--apply] -- ORIGINAL_SIMULATION_OPTIONS')
control <- argv[seq_len(separator-1L)]
args <- argv[seq.int(separator+1L, length(argv))]
source_root <- NULL
destination <- paste0('results/custom_simulation_pipeline_v', simulation_version)
apply <- FALSE
i <- 1L
while (i <= length(control)) {
  key <- control[i]
  if (key == '--apply') { apply <- TRUE; i <- i+1L; next }
  if (!key %in% c('--source','--destination') || i == length(control)) stop('Invalid migration option: ', key)
  if (key == '--source') source_root <- control[i+1L] else destination <- control[i+1L]
  i <- i+2L
}
if (is.null(source_root) || !dir.exists(source_root)) stop('--source must name an existing directory.')
if (any(args %in% c('--output','--job-index','--setting-batch','--settings-per-batch','--plan-only','--help','--version')))
  stop('Pass the original experiment flags without output or selection options.')
# Reuse the runner parser without executing its simulation or dependency setup.
expressions <- parse('scripts/custom_simulation/run_batch.R')
for (expr in expressions) {
  if (is.call(expr) && identical(expr[[1]], as.name('<-')) && identical(expr[[2]], as.name('parse_args'))) eval(expr)
}
opts <- parse_args(args)
if (opts$shard_unit != 'setting') stop('Migration requires --shard-unit setting (the default).')
scenarios <- build_scenarios(opts$n_var, opts$include_imbalance, opts$underdispersion_model)
n_base <- nrow(scenarios)
n_de <- if (opts$de_mode == 'fixed') length(opts$de_proportions) else 1L
order <- unlist(lapply(sort(unique(scenarios$n_samples), decreasing=TRUE), function(w) {
  base <- which(scenarios$n_samples == w)
  unlist(lapply(seq_len(n_de)-1L, function(de) de*n_base+base), use.names=FALSE)
}), use.names=FALSE)
owner <- integer(length(order))
owner[order] <- (seq_along(order)-1L) %% opts$shard_count + 1L
refs <- if (opts$reference_mode == 'brca') opts$reference else opts$references
refs <- normalizePath(refs, mustWork=TRUE)
hashes <- unname(tools::md5sum(refs))
nref <- if (opts$reference_mode == 'all') length(refs) else 1L
files <- list.files(source_root, pattern='-rep[0-9]+-n[0-9]+-seed[0-9]+-ref[[:xdigit:]]{8}\\.rds$', recursive=TRUE, full.names=TRUE)
files <- normalizePath(files, mustWork=TRUE)
if (!length(files)) stop('No simulation datasets found.')
plan <- vector('list', length(files))
for (k in seq_along(files)) {
  object <- readRDS(files[k]); job <- object$scenario
  if (is.null(object$simulated$dge$counts) || is.null(job$setting_job_index)) stop('Missing dataset/job metadata: ', files[k])
  setting <- floor((job$setting_job_index-1)/opts$replicates)+1
  base <- (setting-1) %% n_base+1
  de <- floor((setting-1)/n_base)+1
  if (length(setting)!=1L || !is.finite(setting) || setting<1 || setting>length(owner)) stop('Invalid setting index: ', files[k])
  expected_rep <- (job$setting_job_index-1) %% opts$replicates+1
  if (job$replicate != expected_rep || job$replicate != object$replicate || job$n_var != opts$n_var ||
      nrow(object$simulated$dge$counts) != opts$n_var || job$de_mode != opts$de_mode ||
      job$dispersion_estimator != opts$dispersion_estimator) stop('Original options do not match dataset: ', files[k])
  for (column in intersect(names(scenarios), names(job))) {
    if (column %in% c('scenario_id','simulation_version','de_fraction','n_diffexp')) next
    if (!isTRUE(all.equal(job[[column]], scenarios[[column]][base], check.attributes=FALSE)))
      stop('Scenario mismatch (',column,'): ', files[k])
  }
  if (opts$de_mode == 'fixed' && job$de_fraction != opts$de_proportions[de]) stop('DE proportion mismatch: ', files[k])
  # Random mode records the chosen reference (1..N), but has only one
  # logical reference slot per job. Only all mode expands reference slots.
  slot <- if (opts$reference_mode == 'all') job$reference_index else 1L
  index <- (job$setting_job_index-1)*nref + slot
  if (index != job$job_index || job$seed != as.integer((as.double(opts$base_seed)+index*1009+job$replicate) %% .Machine$integer.max)) stop('Seed/job index mismatch: ', files[k])
  ref <- if (opts$reference_mode == 'all') job$reference_index else 1L
  if (opts$reference_mode == 'random') {
    set.seed(as.integer((as.double(opts$base_seed)+job$setting_job_index*104729) %% .Machine$integer.max))
    ref <- sample.int(length(refs),1L)
  }
  if (job$reference_hash != hashes[ref]) stop('Reference checksum mismatch: ', files[k])
  shard <- sprintf('shard-%03d-of-%03d', owner[setting], opts$shard_count)
  batch <- sprintf('batch-%06d', ceiling(setting/100))
  target <- file.path(destination, shard, batch, 'datasets', shard, basename(files[k]))
  if (file.exists(file.path(destination,shard,batch,'complete.rds')) && normalizePath(target,mustWork=FALSE) != files[k]) stop('Target batch is already complete: ', target)
  if (file.exists(target) && normalizePath(target) != files[k]) stop('Destination already exists; no files changed: ', target)
  plan[[k]] <- data.frame(source=files[k], destination=target, worker=owner[setting], batch=batch, job_index=index)
  if (k %% 100L == 0L) {cat('Validated', k, 'of', length(files), 'datasets\n'); flush.console()}
}
plan <- do.call(rbind,plan)
if (anyDuplicated(plan$destination)) stop('Duplicate datasets map to the same target; no files changed.')
cat('Validated', nrow(plan), 'datasets. Planned batches/workers:\n')
print(with(plan, table(worker,batch)))
if (!apply) {cat('Preview only. Repeat with --apply to move files. Stop ALL simulation/pipeline workers first.\n'); quit(status=0L)}
dir.create(destination,recursive=TRUE,showWarnings=FALSE)
log <- file.path(destination,paste0('migration-',format(Sys.time(),'%Y%m%d-%H%M%S'),'.csv'))
write.csv(plan,log,row.names=FALSE)
for (k in seq_len(nrow(plan))) {
  src <- plan$source[k]; dst <- plan$destination[k]
  if (normalizePath(dst,mustWork=FALSE)==src) next
  dir.create(dirname(dst),recursive=TRUE,showWarnings=FALSE)
  if (!file.rename(src,dst)) stop('Move failed: ',src,'. No copy/delete fallback used. Consult ',log)
}
cat('Migration complete. Move log:',log,'\nRestart run_pipeline.R with the SAME experiment flags and each worker index.\n')
