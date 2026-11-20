# End-to-end restart with completed fixture batches: no DE packages required.
source('scripts/simulation_pipeline/machine_config.R')
cfg <- read_machine_config('config/assignments_100/config-030.dcf')
root <- tempfile('restart-summary-')
dir.create(root)
cfg$Output <- root
cfg[['Underdispersion-Model']] <- 'generalized_poisson'
cfg[['Job-Start']] <- 1; cfg[['Job-End']] <- 2
cfg$Workers <- 2; cfg[['Jobs-Per-Batch']] <- 1
config <- tempfile(fileext='.dcf')
write.dcf(as.data.frame(lapply(cfg, function(x) paste(x,collapse=',')),check.names=FALSE),config)
# Preserve a legacy snapshot to exercise the real version/path resume flow.
original <- cfg; original$Output <- 'results/assignments_100_v5.2/config-030-DESeq2'
saveRDS(list(version='5.2',config=original),file.path(root,'machine_options.rds'))
for (i in 1:2) {
  batch <- file.path(root,sprintf('worker-%02d',i),sprintf('jobs-%d-%d',i,i))
  result <- file.path(batch,'analysis','results','result.rds')
  manifest <- file.path(batch,'analysis','shard-001-of-001','manifest.csv')
  dir.create(dirname(result),recursive=TRUE)
  dir.create(dirname(manifest),recursive=TRUE)
  saveRDS(list(),result)
  write.csv(data.frame(dataset_id=paste0('d',i), method='DESeq2',setting_id='default',
    status='complete',output_file=result),manifest,row.names=FALSE)
  saveRDS(list(start=i,end=i,stage='both'),file.path(batch,'complete.rds'))
}
output <- system2(file.path(R.home('bin'),'Rscript'),
  c('--vanilla','scripts/simulation_pipeline/run_machine.R','--config',shQuote(config)),
  stdout=TRUE,stderr=TRUE)
if (!is.null(attr(output,'status'))) cat(paste(output, collapse='\n'), '\n')
stopifnot(is.null(attr(output,'status')),
  any(grepl('STARTUP SUMMARY',output,fixed=TRUE)),
  any(grepl('FINAL SUMMARY',output,fixed=TRUE)),
  any(grepl('2 / 2 (100.0%)',output,fixed=TRUE)),
  any(grepl('Remaining datasets (including partial batches): 0',output,fixed=TRUE)),
  any(grepl('preserving uncapped v5.2',output,fixed=TRUE)),
  !length(list.files(root,pattern='generation.log',recursive=TRUE)),
  identical(readRDS(file.path(root,'machine_options.rds'))$version,'5.2'),
  length(readRDS(file.path(root,'resume_history.rds'))) == 1)
unlink(root,recursive=TRUE); unlink(config)
cat('PASS: startup/final summaries, legacy resume, completed work skipped without generation, history retained\n')
