#!/usr/bin/env Rscript
# Copy pipeline analysis outputs to portable storage; never alter source files.
a <- commandArgs(FALSE)
p <- normalizePath(sub('^--file=', '', grep('^--file=', a, value=TRUE)[1]))
setwd(normalizePath(file.path(dirname(p), '../..')))
source('scripts/shared/version.R')
args <- commandArgs(TRUE)
assignment_root <- 'results/assignments_100_v5'
opts <- list(source=if (dir.exists(assignment_root)) assignment_root else paste0('results/custom_simulation_pipeline_v',simulation_version), destination=NULL)
if (!length(args) || '--help' %in% args) {
  cat('Usage: Rscript scripts/result_transfer/export_results.R --destination "E:/RNA-results" [--source PATH]\nExports finalized analysis files from all config folders or a single config folder. Also supports the older shard/batch pipeline.\n')
  quit(status=0L)
}
if (length(args) %% 2L) stop('Each option requires a value.')
for(i in seq.int(1L,length(args),by=2L)) {
  key <- sub('^--','',args[i])
  if (!key %in% names(opts)) stop('Unknown option: ',args[i])
  opts[[key]] <- args[i+1L]
}
if (is.null(opts$destination)) stop('--destination is required, for example E:/RNA-results')
root <- normalizePath(opts$source,winslash='/',mustWork=TRUE)
if (!dir.exists(root)) stop('Source must be a directory.')
dir.create(opts$destination,recursive=TRUE,showWarnings=FALSE)
destination <- normalizePath(opts$destination,winslash='/',mustWork=TRUE)
canonical <- function(x) if (.Platform$OS.type=='windows') tolower(x) else x
if (canonical(root)==canonical(destination) || startsWith(canonical(destination),paste0(canonical(root),'/')) || startsWith(canonical(root),paste0(canonical(destination),'/')))
  stop('Choose a destination outside the source directory tree.')
config_pattern <- '^config-([0-9]+)(-.*)?$'
single_config <- grepl(config_pattern, basename(root))
configs <- if (single_config) root else list.dirs(root, recursive=FALSE)
configs <- configs[grepl(config_pattern, basename(configs))]
assignment <- length(configs) > 0L
if (assignment) {
  ids <- as.integer(sub(config_pattern, '\\1', basename(configs)))
  if (anyNA(ids) || anyDuplicated(ids)) stop('Expected unique numeric config indexes in the source.')
  result_files <- function(config) {
    workers <- list.dirs(config, recursive=FALSE)
    workers <- workers[grepl('^worker-[0-9]+$', basename(workers))]
    batches <- unlist(lapply(workers, list.dirs, recursive=FALSE), use.names=FALSE)
    batches <- batches[grepl('^jobs-[0-9]+-[0-9]+$', basename(batches))]
    if (!length(batches)) return(character())
    paths <- list.files(file.path(batches, 'analysis', 'results'), pattern='[.]rds$',
                        recursive=TRUE, full.names=TRUE)
    paths[!grepl('[.]tmp([.]|$)', basename(paths))]
  }
  by_config <- lapply(configs, result_files)
  files <- as.character(unlist(by_config, use.names=FALSE))
  folders <- rep(sprintf('%03d', ids), lengths(by_config))
  relative_root <- if (single_config) dirname(root) else root
} else {
  files <- list.files(root,recursive=TRUE,full.names=TRUE)
  relative_root <- root
}
relative <- substring(chartr('\\', '/', files), nchar(relative_root)+2L)
result <- grepl('^shard-[0-9]+-of-[0-9]+/batch-[0-9]+/analysis/results/.+\\.rds$',relative)
metadata <- grepl('^shard-[0-9]+-of-[0-9]+/(options\\.rds|batch-[0-9]+/(complete\\.rds|analysis/shard-[0-9]+-of-[0-9]+/manifest\\.csv))$',relative)
logs <- grepl('^shard-[0-9]+-of-[0-9]+/logs/pipeline-[^/]+\\.log$', relative)
if (assignment) result <- rep(TRUE, length(files))
keep <- which(result | metadata | logs)
if (!any(result)) stop('No finalized analysis files found. Check --source and the analysis/results folders.')
rows <- vector('list',length(keep))
verified_copy <- function(src,dst,hash) {
  if (file.exists(dst) && identical(unname(tools::md5sum(dst)),hash)) return('unchanged')
  dir.create(dirname(dst),recursive=TRUE,showWarnings=FALSE)
  if (!dir.exists(dirname(dst))) stop('Cannot create destination folder: ',dirname(dst))
  tmp <- paste0(dst,'.',Sys.getpid(),'.tmp')
  on.exit(unlink(tmp),add=TRUE)
  if (!file.copy(src,tmp,overwrite=TRUE) || !identical(unname(tools::md5sum(tmp)),hash) || !identical(unname(tools::md5sum(src)),hash))
    stop('Copy verification failed or source changed during export: ',src,'. Rerun after the writer finishes.')
  backup <- paste0(dst,'.previous')
  if (file.exists(backup)) stop('Previous interrupted replacement found: ',backup)
  had <- file.exists(dst)
  if (had && !file.rename(dst,backup)) stop('Cannot replace destination: ',dst)
  if (!file.rename(tmp,dst)) {
    if(had) file.rename(backup,dst)
    stop('Cannot finalize copy: ',dst)
  }
  if(had) unlink(backup)
  'copied'
}
index_path <- file.path(destination,'export_index.csv')
index_already_there <- file.exists(index_path)
previous <- if (index_already_there)
  read.csv(index_path, colClasses='character', check.names=FALSE) else NULL
if (!is.null(previous)) {
  required <- c('source_relative','export_relative','md5','type')
  if (!all(required %in% names(previous)) || anyNA(previous[required]))
    stop('Destination export index is incomplete.')
  indexed <- grepl('^[0-9]{3,}/[0-9]{6,}[.]rds$', previous$export_relative)
  if (any(indexed != assignment))
    stop('Destination has an incompatible export index. Choose a separate export folder.')
  if (assignment && (anyDuplicated(previous$source_relative) || anyDuplicated(previous$export_relative)))
    stop('Destination export index has duplicate paths.')
}
if (assignment) {
  next_number <- vapply(unique(folders), function(folder) {
    existing <- list.files(file.path(destination, folder), pattern='^[0-9]+[.]rds$')
    indexed <- previous$export_relative[startsWith(as.character(previous$export_relative), paste0(folder, '/'))]
    max(c(0, as.numeric(sub('[.]rds$', '', c(existing, basename(as.character(indexed))))))) + 1
  }, numeric(1))
  previous_row <- if (is.null(previous)) rep(NA_integer_, length(relative)) else
    match(relative, previous$source_relative)
}
for(k in seq_along(keep)) {
  i <- keep[k]; src <- files[i]; rel <- relative[i]
  hash <- unname(tools::md5sum(src))
  if (is.na(hash)) stop('Cannot read: ',src)
  # Short names prevent the repeated dataset ID from exceeding Windows limits.
  target <- if(result[i]) paste0(sub('/analysis/results/.*$','',rel),'/analysis/results/',hash,'.rds') else rel
  if (assignment) {
    folder <- folders[i]
    old <- previous_row[i]
    if (!is.na(old) && identical(previous$md5[old], hash)) {
      target <- previous$export_relative[old]
    } else {
      target <- sprintf('%s/%06.0f.rds', folder, next_number[[folder]])
      next_number[[folder]] <- next_number[[folder]] + 1
    }
  }
  destination_file <- file.path(destination,target)
  already_there <- file.exists(destination_file)
  action <- verified_copy(src,destination_file,hash)
  rows[[k]] <- data.frame(source_relative=rel,export_relative=target,md5=hash,type=if(result[i])'result' else if(logs[i])'log' else 'metadata',already_there=already_there,action=action)
  if(k %% 100L==0L) {cat('Verified',k,'of',length(keep),'files\n');flush.console()}
}
index <- do.call(rbind,rows)
if (assignment && !is.null(previous)) {
  # Keep other configs indexed when exporting one config at a time.
  retained <- previous[!previous$source_relative %in% index$source_relative, , drop=FALSE]
  if (nrow(retained)) {
    retained$already_there <- TRUE
    retained$action <- 'unchanged'
    index <- rbind(index, retained[names(index)])
  }
}
tmp_index <- tempfile('export-index-',tmpdir=destination,fileext='.csv')
write.csv(index,tmp_index,row.names=FALSE)
invisible(verified_copy(tmp_index,index_path,unname(tools::md5sum(tmp_index))))
unlink(tmp_index)
cat('Export complete:',sum(result),'analysis results verified;',sum(vapply(rows, function(row) row$type == 'result' && row$already_there, logical(1))),'already present in destination. Destination:',destination,'\nIndex:',index_path,'(already present:',index_already_there,')\n')
