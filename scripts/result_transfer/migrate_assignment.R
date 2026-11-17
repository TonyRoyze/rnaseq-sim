#!/usr/bin/env Rscript
# Archive assignment results by default; optionally include datasets for resume.
assignment_path <- function(path) {
  path <- path.expand(path)
  if (file.exists(path) || dir.exists(path))
    return(normalizePath(path, winslash='/', mustWork=TRUE))
  file.path(assignment_path(dirname(path)), basename(path))
}

migrate_assignment <- function(args) {
  if (!length(args) || any(args %in% c('--help', '-h'))) {
    cat('Usage: Rscript scripts/result_transfer/migrate_assignment.R',
        '--config PATH [--source PATH] --destination PATH [--include-datasets] [--apply]\n',
        'Preview by default. Copies results, status, logs and configuration; excludes generated datasets.\n',
        'Use --include-datasets only for a full resume transfer.\n',
        'Source defaults to the config Output. Destination names the exact assignment folder.\n')
    return(invisible(NULL))
  }
  opts <- list(config=NULL, source=NULL, destination=NULL, apply=FALSE, include_datasets=FALSE)
  i <- 1L
  while (i <= length(args)) {
    if (args[i] == '--apply') { opts$apply <- TRUE; i <- i+1L; next }
    if (args[i] %in% c('--include-datasets','--results-only')) {
      opts$include_datasets <- args[i] == '--include-datasets'; i <- i+1L; next
    }
    key <- sub('^--', '', args[i])
    if (!key %in% c('config','source','destination') || i == length(args))
      stop('Unknown option or missing value: ', args[i])
    opts[[key]] <- args[i+1L]; i <- i+2L
  }
  if (is.null(opts$config) || is.null(opts$destination))
    stop('--config and --destination are required.')
  source('scripts/simulation_pipeline/machine_config.R', local=TRUE)
  source('scripts/shared/version.R', local=TRUE)
  source('scripts/simulation_design/factorial_design.R', local=TRUE)
  source('scripts/simulation_design/batch_jobs.R', local=TRUE)
  for (expr in parse('scripts/simulation_counts/run_batch.R')) {
    if (is.call(expr) && identical(expr[[1]], as.name('<-')) &&
        identical(expr[[2]], as.name('parse_args'))) eval(expr)
  }
  cfg <- read_machine_config(opts$config)
  src <- assignment_path(if (is.null(opts$source)) cfg$Output else opts$source)
  dst <- assignment_path(opts$destination)
  if (!dir.exists(src)) stop('Source directory not found: ', src)
  canonical <- function(path) if (.Platform$OS.type == 'windows') tolower(path) else path
  if (canonical(src) == canonical(dst) ||
      startsWith(canonical(src), paste0(canonical(dst), '/')) ||
      startsWith(canonical(dst), paste0(canonical(src), '/')))
    stop('Source and destination must be separate, non-nested directories.')
  locks <- file.path(c(src, dst), '.machine-running')
  archive_marker <- '.results-only-archive.rds'
  if (opts$include_datasets && any(file.exists(file.path(c(src,dst),archive_marker))))
    stop('A results-only archive cannot be used as a full resume transfer. Omit --include-datasets.')
  if (!opts$include_datasets && dir.exists(dst)) {
    old_files <- list.files(dst,recursive=TRUE,all.files=TRUE)
    if (any(grepl('(^|/)datasets/.+[.]rds$|(^|/)parameters(/|$)|^parameters-',old_files)))
      stop('Destination already contains generated datasets or caches. Choose a new results-only folder; existing files are preserved.')
  }
  if (any(file.exists(locks)))
    stop('Assignment is locked. Stop all workers first. If a crash left a lock,',
         ' verify all workers stopped before manually removing it.')
  inventory <- function(root) {
    entries <- list.files(root, recursive=TRUE, full.names=TRUE,
                         all.files=TRUE, no..=TRUE, include.dirs=TRUE)
    if (any(nzchar(Sys.readlink(entries)), na.rm=TRUE))
      stop('Symbolic links are not supported in an assignment transfer.')
    files <- entries[!dir.exists(entries)]
    relative <- substring(files, nchar(root)+2L)
    keep <- !grepl('(^|/)([.]machine-running|[.]migration-receipts)(/|$)', relative) &
      !grepl('[.]tmp($|[.])', basename(files)) & relative != archive_marker
    if (!opts$include_datasets) {
      keep <- keep & (relative == 'machine_options.rds' |
        basename(relative) == 'complete.rds' |
        grepl('(^|/)analysis/results/.+[.]rds$', relative) |
        grepl('[.](log|csv)$', relative)) &
        !grepl('(^|/)parameters(/|$)|^parameters-',relative)
    }
    sort(relative[keep])
  }
  snapshot_file <- file.path(src, 'machine_options.rds')
  if (!file.exists(snapshot_file)) stop('Missing machine_options.rds; this is not an assignment output.')
  snapshot <- readRDS(snapshot_file)
  if (!identical(snapshot, list(version=simulation_version, config=cfg)))
    stop('Saved assignment version/config differs. Use its original config and matching code version.')
  sim_opts <- parse_args(c(machine_simulation_args(cfg),
    '--job-start', as.character(cfg[['Job-Start']]),
    '--job-end', as.character(cfg[['Job-End']])))
  ref_paths <- if (sim_opts$reference_mode == 'brca') sim_opts$reference else sim_opts$references
  if (any(!file.exists(ref_paths))) stop('All configured reference files must be available locally.')
  refs <- data.frame(reference_id=tools::file_path_sans_ext(basename(ref_paths)),
                     reference=ref_paths, reference_hash=unname(tools::md5sum(ref_paths)))
  iterator <- make_job_iterator(build_scenarios(sim_opts$n_var, sim_opts$include_imbalance,
    sim_opts$underdispersion_model), sim_opts, refs)
  expected <- new.env(parent=emptyenv())
  repeat {
    jobs <- iterator$next_batch()
    if (is.null(jobs)) break
    for (j in seq_len(nrow(jobs))) assign(job_output_name(jobs[j,,drop=FALSE]),
                                         jobs[j,,drop=FALSE], envir=expected)
  }
  settings <- select_method_settings(make_method_settings(cfg$Profile), cfg$Methods)
  files <- inventory(src)
  if (!length(files)) stop('No assignment files found.')
  hashes <- unname(tools::md5sum(file.path(src, files)))
  if (anyNA(hashes)) stop('Cannot checksum all source files.')
  dataset_pattern <- '-rep[0-9]+-n[0-9]+-seed[0-9]+-ref[[:xdigit:]]{8}[.]rds$'
  lookup <- function(name) {
    if (!exists(name, envir=expected, inherits=FALSE))
      stop('Dataset is outside this assignment or has an unexpected identity: ', name)
    get(name, envir=expected, inherits=FALSE)
  }
  same <- function(a,b) isTRUE(all.equal(a,b,check.attributes=FALSE))
  for (k in seq_along(files)) {
    rel <- files[k]; path <- file.path(src,rel)
    if (grepl('(^|/)datasets/', rel) && grepl(dataset_pattern, basename(rel))) {
      job <- lookup(basename(rel)); object <- readRDS(path)
      columns <- c('scenario_id','job_index','setting_job_index','replicate','seed',
        'n_var','n_group1','n_group2','count_model','de_mode','de_fraction','n_diffexp',
        'reference_id','reference_index','reference_hash','fixed_fold_change')
      if (!identical(object$simulation_version, simulation_version) ||
          !same(object$seed,job$seed) || !same(object$replicate,job$replicate) ||
          any(!vapply(columns, function(column) same(object$scenario[[column]],job[[column]]), logical(1))) ||
          !same(dim(object$simulated$dge$counts), c(cfg[['N-Var']],job$n_samples)))
        stop('Dataset version, seed, dimensions or scenario mismatch: ', rel)
    } else if (grepl('(^|/)analysis/results/.+[.]rds$', rel)) {
      result <- readRDS(path)
      if (!is.character(result$dataset_id) || length(result$dataset_id) != 1L)
        stop('Missing analysis dataset identity: ', rel)
      job <- lookup(paste0(result$dataset_id,'.rds'))
      selected <- which(vapply(settings,function(spec)
        identical(spec$method,result$method) && identical(spec$setting_id,result$setting_id),logical(1)))
      if (length(selected)!=1L || !same(result$simulation_seed,job$seed) ||
          !identical(result$setting_args,settings[[selected]]$args) ||
          !same(result$analysis_seed, as.integer(job$seed + selected*17L)) ||
          !same(result$simulation_scenario$scenario_id, job$scenario_id))
        stop('Analysis method, parameters, seed or scenario mismatch: ', rel)
    }
    if (k %% 100L == 0L) { cat('Validated',k,'of',length(files),'files\n'); flush.console() }
  }
  existing <- if (dir.exists(dst)) inventory(dst) else character()
  extra <- setdiff(existing,files)
  if (length(extra)) stop('Destination contains extra files; refusing to merge assignments: ',extra[1])
  destinations <- file.path(dst, files)
  present <- file.exists(destinations)
  if (any(dir.exists(destinations))) stop('A destination directory conflicts with a source file.')
  if (any(present) && any(unname(tools::md5sum(destinations[present])) != hashes[present]))
    stop('Destination file conflicts with source; nothing overwritten.')
  cat('Validated assignment:', paste(cfg$Methods,collapse=', '), '| jobs',cfg[['Job-Start']],
      'through',cfg[['Job-End']], '\n')
  cat(length(files),'files;',sum(!present),'to copy;',sum(present),'identical files to reuse.\n')
  cat('Mode:',if(opts$include_datasets) 'full resume transfer' else
      'results/status/log archive (generated datasets and parameter caches excluded)','\n')
  cat('Source:',src,'\nDestination:',dst,'\n')
  if (!opts$apply) {
    cat('Preview only; no files changed. Stop all workers, then repeat with --apply.\n')
    return(invisible(NULL))
  }
  if (!dir.exists(dst) && !dir.create(dst, recursive=TRUE)) stop('Cannot create destination.')
  # Use the launcher's own lock name, preventing concurrent launch or migration.
  owned_locks <- character()
  on.exit(unlink(owned_locks,recursive=TRUE),add=TRUE)
  for (lock in locks) {
    if (!dir.create(lock, showWarnings=FALSE)) stop('Cannot acquire assignment lock: ',lock)
    owned_locks <- c(owned_locks,lock)
  }
  if (!identical(inventory(src),files) ||
      !identical(unname(tools::md5sum(file.path(src,files))),hashes))
    stop('Source changed during validation; stop all writers and retry.')
  for (k in seq_along(files)) {
    source_file <- file.path(src, files[k]); target <- destinations[k]
    if (!file.exists(target)) {
      dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
      tmp <- tempfile('assignment-copy-',tmpdir=dirname(target),fileext='.tmp')
      ok <- file.copy(source_file,tmp,overwrite=FALSE,copy.date=TRUE)
      if (!ok || !identical(unname(tools::md5sum(tmp)),hashes[k])) {
        unlink(tmp); stop('Copy verification failed: ',files[k])
      }
      if (file.exists(target) || !file.rename(tmp,target)) {
        unlink(tmp); stop('Cannot finalize destination file: ',files[k])
      }
    }
    if (!identical(unname(tools::md5sum(target)),hashes[k]))
      stop('Destination changed during copy: ',files[k])
    if (k %% 100L == 0L) { cat('Copied/verified',k,'of',length(files),'files\n'); flush.console() }
  }
  if (!identical(inventory(src),files) ||
      !identical(unname(tools::md5sum(file.path(src,files))),hashes))
    stop('Source changed during copy. Do not resume destination until transfer is verified again.')
  receipt_dir <- file.path(dst,'.migration-receipts')
  dir.create(receipt_dir,showWarnings=FALSE)
  receipt <- tempfile(paste0(format(Sys.time(),'%Y%m%d-%H%M%S'),'-'),
                      tmpdir=receipt_dir,fileext='.csv')
  write.csv(data.frame(relative_path=files,md5=hashes),receipt,row.names=FALSE)
  if (!opts$include_datasets) {
    saveRDS(list(version=simulation_version,config=cfg,includes_datasets=FALSE),
            file.path(dst,archive_marker))
    cat('Results archive verified. Source preserved. Receipt:',receipt,'\n',
        'Generated datasets excluded. This archive is for collection, not a full resume transfer.\n')
  } else {
    cat('Migration copy verified. Source preserved. Receipt:',receipt,'\n',
        'At the receiving repository, place this folder at the config Output path:\n',
        cfg$Output,'\nResume with the SAME config; do not run both copies.\n')
  }
  invisible(receipt)
}

if (sys.nframe() == 0L) {
  script <- grep('^--file=',commandArgs(FALSE),value=TRUE)[1]
  root <- normalizePath(file.path(dirname(sub('^--file=','',script)),'../..'))
  setwd(root)
  migrate_assignment(commandArgs(TRUE))
}
