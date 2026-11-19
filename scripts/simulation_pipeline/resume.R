# Resume only the same assignment; moving its output folder is allowed.
check_assignment_resume <- function(saved, cfg, version) {
  comparable <- cfg
  comparable$Output <- saved$config$Output
  if (!identical(saved$config, comparable))
    stop('Config differs from existing output. Use a new Output directory for a different range or experiment.')
  if (!identical(saved$version, version) &&
      !(identical(version, '5.3') && identical(saved$version, '5.2')))
    stop('Unsupported saved simulation version: ', saved$version)
  !identical(saved$version, version)
}

# A marker alone cannot certify that result files still exist. Rebase saved
# absolute manifest paths after the assignment folder has been renamed/moved.
assignment_batch_complete <- function(batch, start, end, stage) {
  marker <- file.path(batch, 'complete.rds')
  if (!file.exists(marker)) return(FALSE)
  saved <- readRDS(marker)
  if (saved$start != start || saved$end != end || !identical(saved$stage, stage))
    stop('Completion marker does not match this batch: ', batch)
  if (stage == 'both') {
    manifests <- list.files(file.path(batch, 'analysis'), pattern='^manifest.csv$',
                            recursive=TRUE, full.names=TRUE)
    if (!length(manifests)) stop('Completed batch is missing analysis manifests: ', batch)
    for (path in manifests) {
      rows <- read.csv(path, stringsAsFactors=FALSE)
      if (!nrow(rows) || !all(c('status','output_file') %in% names(rows)) ||
          anyNA(rows$status) || any(rows$status != 'complete'))
        stop('Completed batch has an incomplete analysis manifest: ', path)
      old <- chartr("\\", "/", rows$output_file)
      relative <- sub('^.*/analysis/', '', old)
      rebased <- ifelse(relative != old, file.path(batch, 'analysis', relative), rows$output_file)
      present <- file.exists(rebased)
      if (any(!present)) stop('Completed batch is missing analysis result files; restore them before resuming: ', batch)
    }
  } else {
    files <- list.files(file.path(batch, 'datasets'), pattern='[.]rds$', recursive=TRUE)
    if (length(files) != end-start+1) stop('Completed batch is missing generated datasets: ', batch)
  }
  TRUE
}
