source('scripts/shared/pipeline_logging.R')
rscript <- file.path(R.home('bin'), 'Rscript')
script <- tempfile(fileext='.R')
log <- tempfile(fileext='.log')
writeLines(c(
  'source("scripts/shared/pipeline_logging.R")',
  'start_pipeline_log("")',
  'pipeline_context("analysis job 3/20")',
  'broken_fit <- function() stop("fit failed")',
  'tryCatch(log_pipeline_errors(broken_fit()), error=function(e) NULL)',
  'cat("continued after caught error\\n")',
  'pipeline_context("analysis final manifest")',
  'stop("manifest write failed")'
), script)

failure <- tryCatch(
  run_logged_process(rscript, c('--vanilla', shQuote(script)), log, 'worker-01 test'),
  error=identity
)
lines <- readLines(log)
stopifnot(
  inherits(failure, 'error'),
  grepl('exit status 1', conditionMessage(failure), fixed=TRUE),
  any(grepl('analysis job 3/20 | fit failed', lines, fixed=TRUE)),
  any(grepl('at broken_fit()', lines, fixed=TRUE)),
  any(lines == 'continued after caught error'),
  any(grepl('FATAL pid=', lines, fixed=TRUE) &
        grepl('analysis final manifest', lines, fixed=TRUE)),
  any(grepl('PROCESS EXIT worker-01 test', lines, fixed=TRUE))
)
source('scripts/simulation_counts/generalized_poisson.R')
gp_error <- tryCatch(gp_probabilities(120000, -.1), error=identity)
stopifnot(grepl('upper=1199999, limit=1000000', conditionMessage(gp_error), fixed=TRUE))
unlink(c(script, log))
cat('PASS: caught and fatal R diagnostics, child exit status, and GP support details\n')
