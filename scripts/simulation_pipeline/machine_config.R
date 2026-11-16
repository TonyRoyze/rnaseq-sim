# Plain-text DCF configuration with full-line # comments.
source("scripts/de_analysis/method_settings.R")
read_machine_config <- function(path) {
  lines <- readLines(path, warn=FALSE)
  lines <- lines[!grepl("^[[:space:]]*#", lines)]
  con <- textConnection(lines)
  on.exit(close(con), add=TRUE)
  rows <- read.dcf(con)
  if (nrow(rows) != 1L) stop('Config must contain exactly one DCF record.')
  cfg <- as.list(rows[1, ])
  required <- c('Job-Start','Job-End','Workers','Jobs-Per-Batch','Output','Stage',
    'Methods','Profile','Replicates','N-Var','Base-Seed','Reference-Mode','Reference',
    'Dispersion-Estimator','Underdispersion-Model','Include-Imbalance','DE-Mode','DE-Proportions')
  if (!setequal(names(cfg), required)) stop('Config fields missing or unknown: ',
    paste(c(setdiff(required,names(cfg)),setdiff(names(cfg),required)),collapse=', '))
  for (key in c('Job-Start','Job-End','Workers','Jobs-Per-Batch','Replicates','N-Var','Base-Seed')) {
    value <- suppressWarnings(as.numeric(cfg[[key]]))
    if (length(value)!=1 || !is.finite(value) || value<1 || value!=floor(value))
      stop('Expected positive integer: ',key)
    cfg[[key]] <- value
  }
  if (cfg[['Job-End']] < cfg[['Job-Start']]) stop('Job-End must be >= Job-Start.')
  if (!cfg$Stage %in% c('generate','both')) stop('Stage must be generate or both.')
  if (!cfg$Profile %in% c('core','quick')) stop('Profile must be core or quick.')
  if (!cfg[['Include-Imbalance']] %in% c('true','false')) stop('Include-Imbalance must be true or false.')
  if (!nzchar(cfg$Output)) stop('Output cannot be empty.')
  methods <- cfg$Methods
  if (is.na(methods) || !nzchar(trimws(methods))) stop('Uncomment at least one method under Methods.')
  cfg$Methods <- strsplit(trimws(methods), '[,[:space:]]+')[[1]]
  invisible(select_method_settings(make_method_settings(cfg$Profile), cfg$Methods))
  cfg
}
machine_simulation_args <- function(cfg) {
  fields <- c('Replicates','N-Var','Base-Seed','Reference-Mode','Reference',
              'Dispersion-Estimator','Underdispersion-Model','DE-Mode','DE-Proportions')
  c(unlist(lapply(fields,function(key) c(paste0('--',tolower(key)),as.character(cfg[[key]]))),use.names=FALSE),
    if(cfg[['Include-Imbalance']]=='true') '--include-imbalance' else '--no-imbalance', '--no-restore')
}
machine_worker_ranges <- function(cfg) {
  n <- cfg[['Job-End']]-cfg[['Job-Start']]+1
  workers <- min(cfg$Workers,n)
  sizes <- rep(floor(n/workers),workers)
  if(n%%workers>0) sizes[seq_len(n%%workers)] <- sizes[seq_len(n%%workers)]+1
  ends <- cfg[['Job-Start']]-1+cumsum(sizes)
  data.frame(worker=seq_len(workers),start=c(cfg[['Job-Start']],head(ends,-1)+1),end=ends)
}
