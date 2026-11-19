#!/usr/bin/env Rscript
# Cross-platform local process pool driven by an editable machine config.
a <- commandArgs(FALSE)
p <- normalizePath(sub('^--file=', '', grep('^--file=',a,value=TRUE)[1]))
setwd(normalizePath(file.path(dirname(p),'../..')))
source('scripts/simulation_pipeline/machine_config.R')
source('scripts/shared/version.R')
source('scripts/simulation_pipeline/resume.R')
source('scripts/simulation_pipeline/progress.R')
args <- commandArgs(TRUE)
if (any(args %in% c('--help','-h'))) {
  cat('Usage: Rscript scripts/simulation_pipeline/run_machine.R [--config config/machine.dcf] [--plan-only]\n')
  quit(status=0)
}
config_path <- 'config/machine.dcf'
plan_only <- FALSE
i <- 1L
while(i<=length(args)) {
  if(args[i]=='--plan-only') { plan_only <- TRUE; i <- i+1L
  } else if(args[i]=='--config' && i<length(args)) { config_path <- args[i+1L]; i <- i+2L
  } else stop('Unknown or incomplete option: ',args[i])
}
cfg <- read_machine_config(config_path)
ranges <- machine_worker_ranges(cfg)
simulation_args <- machine_simulation_args(cfg)
# Reuse the generator's validation without executing its entry point.
for(expr in parse('scripts/simulation_counts/run_batch.R')) {
  if(is.call(expr) && identical(expr[[1]],as.name('<-')) && identical(expr[[2]],as.name('parse_args'))) eval(expr)
}
opts <- parse_args(c(simulation_args,'--job-start',as.character(cfg[['Job-Start']]),
                     '--job-end',as.character(cfg[['Job-End']])))
source('scripts/simulation_design/factorial_design.R')
scenarios <- build_scenarios(opts$n_var,opts$include_imbalance,opts$underdispersion_model)
total <- nrow(scenarios)*as.double(opts$replicates)*if(opts$de_mode=='fixed') length(opts$de_proportions) else 1
if(cfg[['Job-End']]>total) stop('Job-End exceeds experiment size: ',total)
cat('Original job IDs (inclusive); experiment total:',total,'\n')
print(ranges,row.names=FALSE)
cat('Methods:',paste(cfg$Methods,collapse=', '),'\n')
cat('Stage:',cfg$Stage,'| jobs per worker batch:',cfg[['Jobs-Per-Batch']],
    '| output:',cfg$Output,'\n')
if(plan_only) quit(status=0)
# Apply before spawning processes: each worker owns one CPU thread.
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
           VECLIB_MAXIMUM_THREADS='1',BLIS_NUM_THREADS='1',RCPP_PARALLEL_NUM_THREADS='1')
options(mc.cores=1L)
root <- normalizePath(cfg$Output,winslash='/',mustWork=FALSE)
dir.create(root,recursive=TRUE,showWarnings=FALSE)
if(!dir.exists(root)) stop('Cannot create output directory.')
snapshot <- list(version=simulation_version,config=cfg)
state <- file.path(root,'machine_options.rds')
legacy_resume <- FALSE
if(file.exists(state)) legacy_resume <- check_assignment_resume(readRDS(state), cfg, simulation_version)
if(!file.exists(state)) saveRDS(snapshot,state)
cat('\nSTARTUP SUMMARY: recorded progress before workers resume\n')
show_progress(root, cfg[['Job-End']] - cfg[['Job-Start']] + 1)
cat('Resume plan: verify and skip completed batches; reuse existing files in incomplete\n',
    'batches and retry missing/failed work, then continue unstarted batches.\n', sep='')
flush.console()
# Prevent two launchers from writing the same machine output concurrently.
lock <- file.path(root,'.machine-running')
if(!dir.create(lock,showWarnings=FALSE)) stop('Output is locked by another launcher. After a crash, verify all workers stopped before removing ',lock)
if (legacy_resume) {
  cat('RESUME: preserving uncapped v5.2 results; newly generated datasets use the v5.3 baseline cap.\n')
}
tryCatch({
  history_path <- file.path(root, 'resume_history.rds')
  history <- if (file.exists(history_path)) readRDS(history_path) else list()
  history[[length(history)+1L]] <- list(at=Sys.time(), version=simulation_version,
    config=cfg, original=readRDS(state), legacy_resume=legacy_resume)
  saveRDS(history, history_path)
}, error=function(e) {
  unlink(lock, recursive=TRUE)
  stop(e)
})
worker_main <- function(row,cfg,simulation_args,root,project_root) {
  setwd(project_root)
  source('scripts/simulation_pipeline/resume.R', local=TRUE)
  options(mc.cores=1L)
  rscript <- file.path(R.home('bin'),'Rscript')
  run <- function(script,args,log) {
    code <- system2(rscript,c(shQuote(script),shQuote(args)),stdout=log,stderr=log)
    if(code!=0) stop('Process failed; inspect ',log)
  }
  tryCatch({
    worker_root <- file.path(root,sprintf('worker-%02d',row$worker))
    dir.create(worker_root,recursive=TRUE,showWarnings=FALSE)
    start <- row$start
    while(start<=row$end) {
      end <- min(start+cfg[['Jobs-Per-Batch']]-1,row$end)
      batch <- file.path(worker_root,sprintf('jobs-%.0f-%.0f',start,end))
      if (assignment_batch_complete(batch, start, end, cfg$Stage)) {
        cat('SKIP completed batch:', batch, '\n')
        start <- end+1
        next
      }
      dir.create(batch,recursive=TRUE,showWarnings=FALSE)
      input <- file.path(batch,'datasets'); output <- file.path(batch,'analysis')
      run('scripts/simulation_counts/run_batch.R',c(simulation_args,'--parameter-cache-dir',file.path(root,'parameters'),'--output',input,
          '--job-start',as.character(start),'--job-end',as.character(end)),file.path(batch,'generation.log'))
      files <- list.files(input,pattern='-rep[0-9]+-n[0-9]+-seed[0-9]+-ref[[:xdigit:]]{8}\\.rds$',recursive=TRUE,full.names=TRUE)
      if(length(files)!=end-start+1) stop('Missing generated datasets; inspect ',batch)
      if(cfg$Stage=='both') {
        run('scripts/de_analysis/analyze_batch.R',c('--input',input,'--output',output,
             '--profile',cfg$Profile,'--methods',paste(cfg$Methods,collapse=','),'--no-restore'),file.path(batch,'analysis.log'))
        manifests <- list.files(output,pattern='^manifest.csv$',recursive=TRUE,full.names=TRUE)
        if(!length(manifests) || any(vapply(manifests,function(p) {
          rows <- read.csv(p); !nrow(rows) || any(is.na(rows$status) | rows$status!='complete')
        },logical(1)))) stop('Analysis incomplete; inspect ',batch)
      }
      saveRDS(list(start=start,end=end,stage=cfg$Stage),file.path(batch,'complete.rds'))
      start <- end+1
    }
    list(worker=row$worker,success=TRUE,error='')
  },error=function(e) list(worker=row$worker,success=FALSE,error=conditionMessage(e)))
}
run_pool <- function() {
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  cl <- parallel::makePSOCKcluster(nrow(ranges),outfile=file.path(root,'workers.log'))
  on.exit(parallel::stopCluster(cl),add=TRUE)
  parallel::parLapply(cl,lapply(seq_len(nrow(ranges)),function(i) ranges[i,,drop=FALSE]),
    worker_main,cfg=cfg,simulation_args=simulation_args,root=root,project_root=getwd())
}
results <- run_pool()
cat('\nFINAL SUMMARY: recorded progress after workers stopped\n')
show_progress(root, cfg[['Job-End']] - cfg[['Job-Start']] + 1)
for(result in results) cat('Worker',result$worker,if(result$success) 'complete' else result$error,'\n')
if(any(!vapply(results,`[[`,logical(1),'success'))) quit(status=1L)
