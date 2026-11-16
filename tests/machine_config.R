# Pure planner checks: no DE packages, dataset files or large output required.
source('scripts/simulation_pipeline/machine_config.R')
source('scripts/simulation_design/batch_jobs.R')
source('scripts/shared/version.R')
for(e in parse('scripts/simulation_counts/run_batch.R')) {
  if(is.call(e) && identical(e[[1]],as.name('<-')) && identical(e[[2]],as.name('parse_args'))) eval(e)
}
cfg <- read_machine_config('config/machine.dcf')
ranges <- machine_worker_ranges(cfg)
stopifnot(ranges$start[1]==cfg[['Job-Start']],tail(ranges$end,1)==cfg[['Job-End']],
  all(tail(ranges$start,-1)==head(ranges$end,-1)+1))
cfg[['Job-End']] <- cfg[['Job-Start']]+1
stopifnot(nrow(machine_worker_ranges(cfg))==2)
sc <- data.frame(scenario_id=c('a','b','c'),de_setting_id=c('x','x','y'),n_samples=c(6,10,6))
refs <- data.frame(reference_id=c('ref1','ref2'),reference=c('a','b'),reference_hash=c('12345678','87654321'))
collect <- function(opts) {
  it <- make_job_iterator(sc,opts,refs); jobs <- list()
  repeat { x <- it$next_batch(); if(is.null(x)) break; jobs[[length(jobs)+1]] <- x }
  do.call(rbind,jobs)
}
for(mode in c('fixed','setting','combination')) for(ref in c('brca','all','random')) {
  opts <- parse_args(c('--replicates','3','--de-mode',mode,'--reference-mode',ref,'--batch-size','2'))
  expected <- collect(opts)
  opts$job_start <- 2; opts$job_end <- max(expected$job_index)-1
  actual <- collect(opts)
  expected <- expected[expected$job_index>=opts$job_start & expected$job_index<=opts$job_end,]
  expected <- expected[order(expected$job_index),]; rownames(expected) <- NULL
  rownames(actual) <- NULL
  stopifnot(identical(actual,expected))
  opts$batch_size <- 7L
  stopifnot(identical(collect(opts),actual))
}
for(args in list(c('--job-start','0','--job-end','4'),c('--job-start','1.5','--job-end','4'),
  c('--job-start','4'),c('--job-start','5','--job-end','4'),
  c('--job-start','1','--job-end','4','--shard-count','2'))) {
  stopifnot(inherits(try(parse_args(args),silent=TRUE),'try-error'))
}
opts$job_end <- 1e10
stopifnot(inherits(try(make_job_iterator(sc,opts,refs),silent=TRUE),'try-error'))
cat('PASS: machine splits, exact original job IDs/seeds/truth for all DE/reference modes, batch independence and invalid ranges\n')
