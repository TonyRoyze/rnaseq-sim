# Isolated fixtures only: never read or move existing simulation outputs.
source('scripts/result_transfer/migrate_assignment.R')
source('scripts/simulation_pipeline/machine_config.R')
source('scripts/shared/version.R')
source('scripts/simulation_design/factorial_design.R')
source('scripts/simulation_design/batch_jobs.R')
for (expr in parse('scripts/simulation_counts/run_batch.R')) {
  if (is.call(expr) && identical(expr[[1]],as.name('<-')) &&
      identical(expr[[2]],as.name('parse_args'))) eval(expr)
}
root <- tempfile('assignment-transfer-test-'); dir.create(root)
config_path <- file.path(root,'config.dcf')
lines <- readLines('config/assignments_100/config-085.dcf')
lines <- sub('^N-Var:.*','N-Var: 10',lines)
lines <- sub('^Job-End:.*','Job-End: 2',lines)
writeLines(lines,config_path)
cfg <- read_machine_config(config_path)
src <- file.path(root,'source'); dir.create(src)
snapshot <- list(version=simulation_version,config=cfg)
saveRDS(snapshot,file.path(src,'machine_options.rds'))
opts <- parse_args(c(machine_simulation_args(cfg),'--job-start','1','--job-end','2'))
refs <- data.frame(reference_id=tools::file_path_sans_ext(basename(opts$references)),
 reference=opts$references, reference_hash=unname(tools::md5sum(opts$references)))
jobs <- make_job_iterator(build_scenarios(10),opts,refs)$next_batch()
job <- jobs[1,,drop=FALSE]
relative <- file.path('worker-01','jobs-1-2','datasets','shard-001-of-001',job_output_name(job))
dataset_path <- file.path(src,relative); dir.create(dirname(dataset_path),recursive=TRUE)
dataset <- list(simulation_version=simulation_version,scenario=job,replicate=job$replicate,
 seed=job$seed,simulated=list(dge=list(counts=matrix(1,10,job$n_samples))))
saveRDS(dataset,dataset_path)
spec <- select_method_settings(make_method_settings(cfg$Profile),cfg$Methods)[[1]]
result <- list(dataset_id=sub('[.]rds$','',job_output_name(job)),
 method=spec$method,setting_id=spec$setting_id,setting_args=spec$args,
 simulation_seed=job$seed,analysis_seed=as.integer(job$seed+17L),simulation_scenario=job)
result_path <- file.path(src,'worker-01','jobs-1-2','analysis','results','example.rds')
dir.create(dirname(result_path),recursive=TRUE); saveRDS(result,result_path)
writeLines('RUN example',file.path(src,'worker-01','jobs-1-2','analysis.log'))
writeLines('incomplete',paste0(dataset_path,'.tmp'))
dst <- file.path(root,'destination')
args <- c('--config',config_path,'--source',src,'--destination',dst,'--include-datasets')
run <- function(extra=character()) invisible(capture.output(migrate_assignment(c(args,extra))))
fails <- function(expr, pattern) {
  error <- tryCatch({force(expr); NULL},error=function(e)e)
  stopifnot(inherits(error,'error'),grepl(pattern,conditionMessage(error),fixed=TRUE))
}
run(); stopifnot(!dir.exists(dst))
run('--apply')
stopifnot(file.exists(dataset_path),file.exists(file.path(dst,relative)),
 !file.exists(paste0(file.path(dst,relative),'.tmp')),
 identical(unname(tools::md5sum(dataset_path)),unname(tools::md5sum(file.path(dst,relative)))),
 identical(readRDS(file.path(dst,'machine_options.rds')),snapshot),
 !dir.exists(file.path(src,'.machine-running')), !dir.exists(file.path(dst,'.machine-running')))
run('--apply') # Idempotent retry, including a previous migration receipt.
# The copied folder is itself a valid source for the receiving machine.
receiving <- file.path(root,'receiving')
invisible(capture.output(migrate_assignment(c('--config',config_path,'--source',dst,
 '--destination',receiving,'--include-datasets','--apply'))))
stopifnot(file.exists(file.path(receiving,relative)))
fails(migrate_assignment(args[args!='--include-datasets']),
      'Destination already contains generated datasets')
writeLines('conflict',file.path(dst,'machine_options.rds'))
fails(run('--apply'),'conflicts')
saveRDS(snapshot,file.path(dst,'machine_options.rds'))
dir.create(file.path(src,'.machine-running'))
fails(run(),'locked'); unlink(file.path(src,'.machine-running'),recursive=TRUE)
dir.create(file.path(dst,'.machine-running'))
fails(run(),'locked'); unlink(file.path(dst,'.machine-running'),recursive=TRUE)
snapshot$version <- '5.1'; saveRDS(snapshot,file.path(src,'machine_options.rds'))
fails(run(),'version/config differs')
snapshot$version <- simulation_version; saveRDS(snapshot,file.path(src,'machine_options.rds'))
dataset$seed <- 1L; saveRDS(dataset,dataset_path)
fails(run(),'Dataset version, seed')
dataset$seed <- job$seed; dataset$simulation_version <- '5.1'; saveRDS(dataset,dataset_path)
fails(run(),'Dataset version, seed')
dataset$simulation_version <- simulation_version; saveRDS(dataset,dataset_path)
result$method <- 'ROTS'; saveRDS(result,result_path)
fails(run(),'Analysis method')
saveRDS(within(result,method <- spec$method),result_path)
# Results-only is the default. Do not read/copy count matrices or parameter caches.
status_path <- file.path(dirname(dataset_path),'status.csv')
writeLines('status\ncomplete',status_path)
parameter_path <- file.path(src,'parameters','cache.rds')
dir.create(dirname(parameter_path)); writeLines('not an RDS',parameter_path)
writeLines('not an RDS',dataset_path)
archive <- file.path(root,'results-only')
archive_args <- c('--config',config_path,'--source',src,'--destination',archive)
invisible(capture.output(migrate_assignment(archive_args)))
stopifnot(!dir.exists(archive))
invisible(capture.output(migrate_assignment(c(archive_args,'--apply'))))
stopifnot(!file.exists(file.path(archive,relative)),
          !dir.exists(file.path(archive,'parameters')),
          file.exists(file.path(archive,substring(result_path,nchar(src)+2L))),
          file.exists(file.path(archive,substring(status_path,nchar(src)+2L))),
          file.exists(file.path(archive,'worker-01','jobs-1-2','analysis.log')),
          file.exists(file.path(archive,'.results-only-archive.rds')),
          identical(readLines(dataset_path),'not an RDS'))
invisible(capture.output(migrate_assignment(c(archive_args,'--apply'))))
fails(migrate_assignment(c('--config',config_path,'--source',archive,
 '--destination',file.path(root,'invalid-full'),'--include-datasets')),'results-only archive')
fails(migrate_assignment(c('--config',config_path,'--source',src,
 '--destination',file.path(src,'nested'))),'non-nested')
stopifnot(!dir.exists(file.path(src,'nested')))
unlink(root,recursive=TRUE)
cat('PASS: preview, verified copy, source preservation, idempotence, receiving-machine copy,',
    'results-only exclusion, temporary-file exclusion, conflict/lock/version/seed/method rejection and nested-path protection\n')
