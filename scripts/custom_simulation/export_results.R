#!/usr/bin/env Rscript
# Copy pipeline analysis outputs to portable storage; never alter source files.
a <- commandArgs(FALSE)
p <- normalizePath(sub('^--file=', '', grep('^--file=', a, value=TRUE)[1]))
setwd(normalizePath(file.path(dirname(p), '../..')))
source('scripts/custom_simulation/version.R')
args <- commandArgs(TRUE)
opts <- list(source=paste0('results/custom_simulation_pipeline_v',simulation_version), destination=NULL)
if (!length(args) || '--help' %in% args) {
  cat('Usage: Rscript scripts/custom_simulation/export_results.R --destination "E:/RNA-results" [--source PATH]\nCopies analysis RDS files, manifests, pipeline logs and batch/configuration metadata; excludes datasets.\n')
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
files <- list.files(root,recursive=TRUE,full.names=TRUE)
relative <- substring(gsub('\\\\','/',files), nchar(root)+2L)
result <- grepl('^shard-[0-9]+-of-[0-9]+/batch-[0-9]+/analysis/results/.+\\.rds$',relative)
metadata <- grepl('^shard-[0-9]+-of-[0-9]+/(options\\.rds|batch-[0-9]+/(complete\\.rds|analysis/shard-[0-9]+-of-[0-9]+/manifest\\.csv))$',relative)
logs <- grepl('^shard-[0-9]+-of-[0-9]+/logs/pipeline-[^/]+\\.log$', relative)
keep <- which(result | metadata | logs)
if (!any(result)) stop('No pipeline analysis result files found. Point --source at the pipeline root containing shard folders.')
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
for(k in seq_along(keep)) {
  i <- keep[k]; src <- files[i]; rel <- relative[i]
  hash <- unname(tools::md5sum(src))
  if (is.na(hash)) stop('Cannot read: ',src)
  # Short names prevent the repeated dataset ID from exceeding Windows limits.
  target <- if(result[i]) paste0(sub('/analysis/results/.*$','',rel),'/analysis/results/',hash,'.rds') else rel
  action <- verified_copy(src,file.path(destination,target),hash)
  rows[[k]] <- data.frame(source_relative=rel,export_relative=target,md5=hash,type=if(result[i])'result' else if(logs[i])'log' else 'metadata',action=action)
  if(k %% 100L==0L) {cat('Verified',k,'of',length(keep),'files\n');flush.console()}
}
index <- do.call(rbind,rows)
index_path <- file.path(destination,'export_index.csv')
tmp_index <- tempfile('export-index-',tmpdir=destination,fileext='.csv')
write.csv(index,tmp_index,row.names=FALSE)
invisible(verified_copy(tmp_index,index_path,unname(tools::md5sum(tmp_index))))
unlink(tmp_index)
cat('Export complete:',sum(result),'analysis results verified. Destination:',destination,'\nIndex:',index_path,'\n')
