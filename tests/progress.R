source('scripts/simulation_pipeline/progress.R')
root <- tempfile('progress-')
batch <- file.path(root, 'worker-01', 'jobs-1-20')
gen <- file.path(batch, 'datasets', 'shard-001-of-001')
analysis <- file.path(batch, 'analysis', 'shard-001-of-001')
dir.create(gen, recursive = TRUE)
dir.create(analysis, recursive = TRUE)
write.csv(data.frame(output_file = c('a', 'a', 'b', 'b', 'c'),
  status = c('failed', 'complete', 'running', 'failed', 'complete')),
  file.path(gen, 'status.csv'), row.names = FALSE)
write.csv(data.frame(dataset_id = c('a', 'a'), method = 'test',
  setting_id = 'default', status = c('failed', 'complete')),
  file.path(analysis, 'manifest.csv'), row.names = FALSE)
writeLines('FAILED example', file.path(batch, 'generation.log'))
x <- progress_snapshot(root)[[1]]
stopifnot(x$finished == 0, x$generation['complete'] == 2,
  x$generation['failed'] == 1, x$analysis['complete'] == 1,
  x$analysis['failed'] == 0, grepl('FAILED example', x$latest))
saveRDS(list(start = 1, end = 20), file.path(batch, 'complete.rds'))
x <- progress_snapshot(root)[[1]]
stopifnot(x$finished == 20, x$batches == 1)
stopifnot(length(progress_snapshot(tempfile())) == 0)
text <- capture.output(show_progress(root, 40))
stopifnot(any(grepl('20 / 40 (50.0%)', text, fixed = TRUE)),
  any(grepl('Remaining datasets (including partial batches): 20', text, fixed = TRUE)))
unlink(root, recursive = TRUE)
cat('PASS: monitor retry deduplication, batch totals, analysis settings, latest logs, and empty outputs\n')
config_text <- capture.output(progress_main(c('--config', 'config/assignments_100/config-007.dcf')))
stopifnot(any(grepl('config-007-NBPSeq', config_text, fixed = TRUE)))
