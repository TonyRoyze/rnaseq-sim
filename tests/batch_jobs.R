# Run from the project root: Rscript --vanilla tests/batch_jobs.R
simulation_version <- '3.0'
for (e in parse('scripts/custom_simulation/run_batch.R')) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) &&
      as.character(e[[2]])[1] %in% c('parse_args', 'build_scenarios', 'run_job')) eval(e)
}
source('scripts/custom_simulation/batch_jobs.R')
source('scripts/custom_simulation/de_proportion_modes.R')
opts <- parse_args(c('--replicates','3','--batch-size','7'))
stopifnot(inherits(try(parse_args(c('--batch-size','0')),silent=TRUE),'try-error'))
sc <- build_scenarios(100)
# Several weights and a DE setting shared across size levels and count models.
sc <- sc[c(1:12, nrow(sc)-11:0), ]
refs <- data.frame(reference_id=c('a','b'),reference=c('a.rds','b.rds'),reference_hash=c('12345678','87654321'))
collect <- function(iterator) {
  rows <- list()
  repeat {
    batch <- iterator$next_batch()
    if(is.null(batch)) break
    stopifnot(nrow(batch)<=opts$batch_size)
    rows[[length(rows)+1]] <- batch
  }
  if(!length(rows)) return(NULL)
  do.call(rbind,rows)
}
set.seed(123); previous_rng <- .Random.seed
for (mode in c('fixed','combination','setting')) for (reference_mode in c('brca','all','random')) {
  opts$de_mode <- mode
  opts$reference_mode <- reference_mode
  legacy <- sc[rep(seq_len(nrow(sc)),each=opts$replicates),]
  legacy$replicate <- rep(seq_len(opts$replicates),nrow(sc))
  legacy <- assign_de_proportions(legacy,mode,opts$de_proportions,opts$base_seed,opts$n_var)
  legacy$setting_job_index <- seq_len(nrow(legacy))
  if(reference_mode=='all') {
    legacy <- legacy[rep(seq_len(nrow(legacy)),each=nrow(refs)),]
    legacy$reference_index <- rep(seq_len(nrow(refs)),nrow(legacy)/nrow(refs))
  } else if(reference_mode=='random') {
    legacy$reference_index <- preserve_rng(function() vapply(legacy$setting_job_index,function(i) {
      set.seed(as.integer((as.double(opts$base_seed)+i*104729)%%.Machine$integer.max))
      sample.int(nrow(refs),1L)
    },integer(1)))
  } else legacy$reference_index <- 1L
  legacy$job_index <- seq_len(nrow(legacy))
  legacy$seed <- as.integer((as.double(opts$base_seed)+legacy$job_index*1009+legacy$replicate)%%.Machine$integer.max)
  shards <- list()
  opts$shard_count <- 4L
  for(shard in 1:4) {
    opts$shard_index <- shard
    iterator <- make_job_iterator(sc,opts,refs)
    part <- collect(iterator)
    stopifnot(nrow(part)==iterator$assigned)
    shards[[shard]] <- part
  }
  actual <- do.call(rbind,shards)
  actual <- actual[order(actual$job_index),]
  stopifnot(nrow(actual)==nrow(legacy),!anyDuplicated(actual$job_index))
  for(column in c('job_index','scenario_id','replicate','de_fraction','seed','reference_index','setting_job_index'))
    stopifnot(isTRUE(all.equal(actual[[column]],legacy[[column]],check.attributes=FALSE)))
  # Changing batch size does not change ownership, RNG draws, or filenames.
  opts$batch_size <- 1L
  one <- collect(make_job_iterator(sc,opts,refs))
  expected <- shards[[4]]
  stopifnot(identical(job_output_name(one),job_output_name(expected)))
  opts$batch_size <- 7L
}
stopifnot(identical(previous_rng,.Random.seed))
opts$replicates <- 100L; opts$de_mode <- 'fixed'; opts$reference_mode <- 'brca'
full <- make_job_iterator(build_scenarios(100),opts,refs)
stopifnot(full$total==16329600,full$assigned==4082400,nrow(full$next_batch())==7)
opts$shard_count <- 1000000L; opts$shard_index <- 1000000L; opts$replicates <- 1L
empty <- make_job_iterator(sc,opts,refs)
stopifnot(empty$assigned==0,is.null(empty$next_batch()))

# Exercise real output writes, failed-job retry, and resume from old shards,
# using a tiny simulated payload so no edgeR installation is needed.
output_root <- tempfile('batch-jobs-'); dir.create(output_root)
shard_dir <- file.path(output_root,'shard-001-of-004'); dir.create(shard_dir)
old_dir <- file.path(output_root,'shard-002-of-004'); dir.create(old_dir)
output_dirs <- c(shard_dir,old_dir)
manifest_path <- file.path(shard_dir,'status.csv')
iterator <- list(assigned=2)
job <- actual[1,,drop=FALSE]
job$reference <- 'unused'; job$n_var <- 100
get_reference_parameters <- function(...) list()
simulate_custom_counts <- function(...) list(counts=matrix(1L,2,2))
run_job(job,0)
stopifnot(file.exists(file.path(shard_dir,job_output_name(job))))
log <- read.csv(manifest_path)
stopifnot(identical(log$status,c('running','complete')))
payload <- readRDS(log$output_file[2]); stopifnot(payload$seed==job$seed)
# Previous sharding layout: find the identical completed file and skip it.
stopifnot(file.rename(log$output_file[2],file.path(old_dir,job_output_name(job))))
simulate_custom_counts <- function(...) stop('should be skipped')
run_job(job,0)
stopifnot(tail(read.csv(manifest_path)$output_file,1)==file.path(old_dir,job_output_name(job)))
job$scenario_id <- paste0(job$scenario_id,'_retry')
simulate_custom_counts <- function(...) stop('failure with comma, and\nnewline')
run_job(job,0)
stopifnot(tail(read.csv(manifest_path)$status,1)=='failed')
simulate_custom_counts <- function(...) list(counts=matrix(2L,2,2))
run_job(job,0)
stopifnot(tail(read.csv(manifest_path)$status,1)=='complete')
unlink(output_root,recursive=TRUE)
cat('PASS: all modes match legacy IDs/seeds/draws; disjoint shards; batch independence; 100-replicate counts; append log; resume and retry\n')
