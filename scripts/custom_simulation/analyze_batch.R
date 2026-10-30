#!/usr/bin/env Rscript
# Evaluate supported DE methods and method-specific settings on generated RDS
# datasets. Each dataset/method/setting result is independently resumable.

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
setwd(project_root)

parse_args <- function(args) {
  opts <- list(input = "results/custom_simulation",
               output = "results/custom_simulation_analysis",
               shard_index = 1L, shard_count = 1L,
               profile = "core", restore = TRUE, plan_only = FALSE)
  if ("--help" %in% args || "-h" %in% args) {
    cat(paste(
      "Usage: Rscript scripts/custom_simulation/analyze_batch.R [options]",
      "Options:",
      "  --input PATH        Generated simulation directory",
      "  --output PATH       Analysis result directory",
      "  --profile NAME      core or quick (default: core)",
      "  --shard-index N     Worker shard, starting at 1",
      "  --shard-count N     Number of workers",
      "  --no-restore        Skip renv::restore()",
      "  --plan-only         Print assigned analysis jobs and exit",
      sep = "\n"
    ), "\n")
    return(NULL)
  }
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--no-restore") {
      opts$restore <- FALSE
      i <- i + 1L
    } else if (key == "--plan-only") {
      opts$plan_only <- TRUE
      i <- i + 1L
    } else if (key %in% c("--input", "--output", "--profile", "--shard-index", "--shard-count")) {
      if (i == length(args)) stop("Missing value for ", key)
      value <- args[i + 1L]
      name <- switch(key, "--input" = "input", "--output" = "output",
                     "--profile" = "profile", "--shard-index" = "shard_index",
                     "--shard-count" = "shard_count")
      opts[[name]] <- if (name %in% c("shard_index", "shard_count")) {
        suppressWarnings(as.integer(value))
      } else value
      i <- i + 2L
    } else stop("Unknown option: ", key)
  }
  if (is.na(opts$shard_index) || is.na(opts$shard_count) ||
      opts$shard_index < 1L || opts$shard_count < 1L ||
      opts$shard_index > opts$shard_count) stop("Invalid shard settings.")
  if (!opts$profile %in% c("core", "quick")) stop("profile must be 'core' or 'quick'.")
  opts
}

options <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(options)) quit(status = 0L)
source("renv/activate.R")
if (options$restore) renv::restore(prompt = FALSE)
source("scripts/custom_simulation/method_settings.R")

resolve_path <- function(path) {
  is_absolute <- grepl("^(~|/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
  if (is_absolute) path.expand(path) else file.path(project_root, path)
}
input_root <- resolve_path(options$input)
output_root <- resolve_path(options$output)
if (!dir.exists(input_root)) stop("Input directory not found: ", input_root)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

input_files <- sort(list.files(
  input_root,
  pattern = "-rep[0-9]+-n[0-9]+-seed[0-9]+-ref[[:xdigit:]]{8}\\.rds$",
  recursive = TRUE, full.names = TRUE
))
if (!length(input_files)) stop("No generated simulation datasets found under: ", input_root)
settings <- make_method_settings(options$profile)

dataset_info <- lapply(input_files, function(path) {
  object <- readRDS(path)
  if (is.null(object$simulated$dge) || is.null(object$simulated$truth$is_de)) {
    stop("Not a generated custom-simulation dataset: ", path)
  }
  list(
    path = path,
    dataset_id = sub("\\.rds$", "", basename(path)),
    genes = nrow(object$simulated$dge$counts),
    n_group1 = sum(object$simulated$dge$samples$group == levels(object$simulated$dge$samples$group)[1]),
    n_group2 = sum(object$simulated$dge$samples$group == levels(object$simulated$dge$samples$group)[2]),
    n_samples = ncol(object$simulated$dge$counts),
    samples_per_group = object$simulated$parameters$s
  )
})

jobs <- list()
for (dataset_index in seq_along(dataset_info)) {
  for (setting_index in seq_along(settings)) {
    setting <- settings[[setting_index]]
    weight <- if (setting$method == "ROTS") 10 else if (setting$method == "NOISeq") 4 else
      if (setting$method %in% c("DESeq2", "DESeq", "ABSSeq", "DSS")) 2 else 1
    jobs[[length(jobs) + 1L]] <- data.frame(
      dataset_index = dataset_index,
      dataset_id = dataset_info[[dataset_index]]$dataset_id,
      input_file = dataset_info[[dataset_index]]$path,
      setting_index = setting_index,
      method = setting$method,
      setting_id = setting$setting_id,
      samples_per_group = dataset_info[[dataset_index]]$samples_per_group,
      n_group1 = dataset_info[[dataset_index]]$n_group1,
      n_group2 = dataset_info[[dataset_index]]$n_group2,
      n_samples = dataset_info[[dataset_index]]$n_samples,
      genes = dataset_info[[dataset_index]]$genes,
      estimated_work = dataset_info[[dataset_index]]$n_samples *
        dataset_info[[dataset_index]]$genes * weight,
      stringsAsFactors = FALSE
    )
  }
}
jobs <- do.call(rbind, jobs)
jobs$job_index <- seq_len(nrow(jobs))
jobs$shard <- NA_integer_
shard_load <- rep(0, options$shard_count)
for (job_index in order(-jobs$estimated_work, jobs$job_index)) {
  target <- which(shard_load == min(shard_load))[1]
  jobs$shard[job_index] <- target
  shard_load[target] <- shard_load[target] + jobs$estimated_work[job_index]
}
jobs <- jobs[jobs$shard == options$shard_index, , drop = FALSE]
jobs <- jobs[order(jobs$dataset_index, jobs$setting_index), , drop = FALSE]

if (options$plan_only) {
  print(jobs[, c("job_index", "dataset_id", "method", "setting_id",
                 "samples_per_group", "estimated_work", "shard")])
  cat("\n", nrow(jobs), " analysis jobs assigned to shard ",
      options$shard_index, " of ", options$shard_count, "\n", sep = "")
  quit(status = 0L)
}
if (nrow(jobs) == 0L) {
  cat("No analysis jobs assigned to shard ", options$shard_index, ".\n", sep = "")
  quit(status = 0L)
}

shard_dir <- file.path(output_root,
  sprintf("shard-%03d-of-%03d", options$shard_index, options$shard_count))
dir.create(shard_dir, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(shard_dir, "manifest.csv")
manifest <- if (file.exists(manifest_path)) {
  utils::read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
} else NULL
rows <- lapply(seq_len(nrow(jobs)), function(i) {
  job <- jobs[i, , drop = FALSE]
  spec <- settings[[job$setting_index]]
  file_name <- paste(job$dataset_id, job$method, spec$setting_id, sep = "--")
  output_dir <- file.path(output_root, "results", job$dataset_id)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  record <- job
  record$setting_args <- paste(names(spec$args), vapply(spec$args, function(x) {
    if (is.null(x)) "NULL" else paste(x, collapse = ",")
  }, character(1)), sep = "=", collapse = ";")
  record$output_file <- file.path(output_dir, paste0(file_name, ".rds"))
  record$status <- "pending"
  record$error <- ""
  if (!is.null(manifest)) {
    previous <- manifest[manifest$dataset_id == job$dataset_id &
      manifest$method == job$method & manifest$setting_id == job$setting_id, , drop = FALSE]
    if (nrow(previous)) {
      record$status <- previous$status[1]
      if ("error" %in% names(previous)) record$error <- previous$error[1]
    }
  }
  record
})
manifest_extras <- if (!is.null(manifest)) {
  manifest[!paste(manifest$dataset_id, manifest$method, manifest$setting_id) %in%
    paste(jobs$dataset_id, jobs$method, jobs$setting_id), , drop = FALSE]
} else NULL
write_manifest <- function() {
  current <- rows
  if (!is.null(manifest_extras) && nrow(manifest_extras)) current <- c(current, list(manifest_extras))
  path_tmp <- paste0(manifest_path, ".tmp")
  utils::write.csv(do.call(rbind, current), path_tmp, row.names = FALSE, na = "")
  if (!file.rename(path_tmp, manifest_path)) stop("Could not update manifest: ", manifest_path)
}

last_dataset <- NULL
dataset <- NULL
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, , drop = FALSE]
  record <- rows[[i]]
  output_path <- record$output_file
  spec <- settings[[job$setting_index]]
  if (file.exists(output_path)) {
    record$status <- "complete"
    record$error <- ""
    rows[[i]] <- record
    write_manifest()
    next
  }
  record$status <- "running"
  rows[[i]] <- record
  write_manifest()
  cat(sprintf("RUN [%d/%d] %s / %s / %s\n", i, nrow(jobs),
              job$dataset_id, job$method, job$setting_id))

  tryCatch({
    if (!identical(last_dataset, job$dataset_id)) {
      dataset <- readRDS(job$input_file)
      last_dataset <- job$dataset_id
    }
    simulated <- dataset$simulated
    counts <- simulated$dge$counts
    truth <- simulated$truth$is_de[match(rownames(counts), simulated$truth$gene)]
    keep <- rowSums(counts) > 0
    truth <- truth[keep]
    analysis_dge <- simulated$dge[keep, , keep.lib.sizes = FALSE]
    analysis_seed <- as.integer(dataset$seed + job$setting_index * 17L)
    p <- run_method_setting(analysis_dge, spec, seed = analysis_seed)
    gene_ids <- rownames(analysis_dge$counts)
    truth <- truth[match(names(p), gene_ids)]
    metrics <- summarize_method_result(p, truth)
    result <- list(
      dataset_id = job$dataset_id,
      input_file = job$input_file,
      simulation_scenario = dataset$scenario,
      simulation_seed = dataset$seed,
      analysis_seed = analysis_seed,
      method = job$method,
      setting_id = spec$setting_id,
      setting_args = spec$args,
      score_type = if (job$method == "NOISeq") "1 - NOISeq probability (ranking score)" else "p-value",
      gene_results = data.frame(gene = names(p), p_value = p,
                                is_de = truth, adjusted_p_value = stats::p.adjust(p, "BH")),
      metrics = metrics
    )
    tmp <- paste0(output_path, ".", Sys.getpid(), ".tmp")
    saveRDS(result, tmp, compress = "xz")
    if (!file.rename(tmp, output_path)) stop("Could not finalize analysis result.")
    record$status <- "complete"
    cat("DONE ", basename(output_path), "\n", sep = "")
  }, error = function(e) {
    record$status <<- "failed"
    record$error <<- conditionMessage(e)
    cat("ERROR ", job$method, "/", spec$setting_id, ": ",
        conditionMessage(e), "\n", sep = "")
  })
  rows[[i]] <- record
  write_manifest()
}
write_manifest()
cat("Finished analysis shard ", options$shard_index, " of ",
    options$shard_count, ". Manifest: ", manifest_path, "\n", sep = "")
