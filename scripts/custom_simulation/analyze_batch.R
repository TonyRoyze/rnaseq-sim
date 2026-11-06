#!/usr/bin/env Rscript
# Evaluate supported DE methods and method-specific settings on generated RDS
# datasets. Each dataset/method/setting result is independently resumable.

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
setwd(project_root)
source("scripts/custom_simulation/pipeline_logging.R")
start_pipeline_log()

source("scripts/custom_simulation/version.R")

parse_args <- function(args) {
  opts <- list(input = paste0("results/custom_simulation_v", simulation_version),
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

dataset_ids <- sub("\\.rds$", "", basename(input_files))
extract_capture <- function(value, pattern, label, capture = 2L) {
  match <- regmatches(value, regexec(pattern, value, perl = TRUE))[[1]]
  if (length(match) < capture) stop("Cannot parse ", label, " from dataset filename: ", value)
  match[capture]
}
group1 <- as.integer(vapply(dataset_ids, extract_capture, character(1),
                            pattern = "_g([0-9]+)-([0-9]+)_", label = "group sizes", capture = 2L))
group2 <- as.integer(vapply(dataset_ids, extract_capture, character(1),
                            pattern = "_g([0-9]+)-([0-9]+)_", label = "group sizes", capture = 3L))
n_genes <- as.integer(vapply(dataset_ids, extract_capture, character(1),
                             pattern = "-n([0-9]+)-seed", label = "gene count"))
dataset_info_table <- data.frame(
  dataset_index = seq_along(input_files), dataset_id = dataset_ids,
  input_file = input_files, genes = n_genes,
  n_group1 = group1, n_group2 = group2, n_samples = group1 + group2,
  samples_per_group = ifelse(group1 == group2, group1, NA_integer_),
  stringsAsFactors = FALSE
)
setting_methods <- vapply(settings, `[[`, character(1), "method")
setting_ids <- vapply(settings, `[[`, character(1), "setting_id")
setting_weights <- ifelse(setting_methods == "ROTS", 10,
  ifelse(setting_methods == "NOISeq", 4,
    ifelse(setting_methods %in% c("DESeq2", "DESeq", "ABSSeq", "DSS"), 2, 1)))
dataset_indices <- rep(seq_len(nrow(dataset_info_table)), each = length(settings))
setting_indices <- rep(seq_along(settings), times = nrow(dataset_info_table))
jobs <- dataset_info_table[dataset_indices, , drop = FALSE]
jobs$setting_index <- setting_indices
jobs$method <- setting_methods[setting_indices]
jobs$setting_id <- setting_ids[setting_indices]
jobs$estimated_work <- jobs$n_samples * jobs$genes * setting_weights[setting_indices]
jobs <- jobs[c("dataset_index", "dataset_id", "input_file", "setting_index", "method",
               "setting_id", "samples_per_group", "n_group1", "n_group2", "n_samples",
               "genes", "estimated_work")]
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
setting_args <- vapply(settings, function(spec) paste(
  names(spec$args), vapply(spec$args, function(x) {
    if (is.null(x)) "NULL" else paste(x, collapse = ",")
  }, character(1)), sep = "=", collapse = ";"
), character(1))
job_keys <- paste(jobs$dataset_id, jobs$method, setting_ids[jobs$setting_index], sep = "\r")
previous_index <- rep(NA_integer_, nrow(jobs))
if (!is.null(manifest)) {
  manifest_keys <- paste(manifest$dataset_id, manifest$method, manifest$setting_id, sep = "\r")
  previous_index <- match(job_keys, manifest_keys)
}
previous_status <- rep("pending", nrow(jobs))
previous_error <- rep("", nrow(jobs))
found_previous <- !is.na(previous_index)
if (any(found_previous)) {
  previous_status[found_previous] <- manifest$status[previous_index[found_previous]]
  if ("error" %in% names(manifest)) {
    previous_error[found_previous] <- manifest$error[previous_index[found_previous]]
  }
}

file_names <- paste(jobs$dataset_id, jobs$method, setting_ids[jobs$setting_index], sep = "--")
jobs$setting_args <- setting_args[jobs$setting_index]
jobs$output_file <- file.path(output_root, "results", jobs$dataset_id, paste0(file_names, ".rds"))
jobs$status <- previous_status
jobs$error <- previous_error
manifest_extras <- if (!is.null(manifest)) {
  manifest[!manifest_keys %in% job_keys, , drop = FALSE]
} else NULL
write_manifest <- function() {
  current <- jobs
  if (!is.null(manifest_extras) && nrow(manifest_extras)) current <- rbind(current, manifest_extras)
  path_tmp <- paste0(manifest_path, ".tmp")
  utils::write.csv(current, path_tmp, row.names = FALSE, na = "")
  if (!file.rename(path_tmp, manifest_path)) stop("Could not update manifest: ", manifest_path)
}

manifest_checkpoint_every <- 10000L
completed_since_checkpoint <- 0L
analysis_counts <- c(complete=0L, failed=0L, skipped=0L)
last_dataset <- NULL
dataset <- NULL
cat(sprintf("Prepared %d analysis jobs; manifest checkpoints every %d jobs.\n",
            nrow(jobs), manifest_checkpoint_every))
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, , drop = FALSE]
  output_path <- jobs$output_file[i]
  spec <- settings[[job$setting_index]]
  if (file.exists(output_path)) {
    jobs$status[i] <- "complete"
    jobs$error[i] <- ""
    analysis_counts["skipped"] <- analysis_counts["skipped"] + 1L
  } else {
    jobs$status[i] <- "running"
    jobs$error[i] <- ""
    cat(sprintf("RUN [%d/%d] %s / %s / %s\n", i, nrow(jobs),
                job$dataset_id, job$method, job$setting_id))

    run_error <- tryCatch({
      if (!identical(last_dataset, job$dataset_id)) {
        dataset <- readRDS(job$input_file)
        if (is.null(dataset$simulated$dge) || is.null(dataset$simulated$truth$is_de)) {
          stop("Not a generated custom-simulation dataset: ", job$input_file)
        }
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
      dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
      tmp <- paste0(output_path, ".", Sys.getpid(), ".tmp")
      saveRDS(result, tmp, compress = "xz")
      if (!file.rename(tmp, output_path)) stop("Could not finalize analysis result.")
      NULL
    }, error = function(e) e)

    if (inherits(run_error, "error")) {
      jobs$status[i] <- "failed"
      jobs$error[i] <- conditionMessage(run_error)
      analysis_counts["failed"] <- analysis_counts["failed"] + 1L
      cat("ERROR ", job$method, "/", spec$setting_id, ": ",
          conditionMessage(run_error), "\n", sep = "")
    } else {
      jobs$status[i] <- "complete"
      analysis_counts["complete"] <- analysis_counts["complete"] + 1L
      cat("DONE ", basename(output_path), "\n", sep = "")
    }
  }

  completed_since_checkpoint <- completed_since_checkpoint + 1L
  if (completed_since_checkpoint >= manifest_checkpoint_every) {
    write_manifest()
    cat(sprintf("CHECKPOINT [%d/%d] manifest saved; %d new results calculated, %d failed, %d existing skipped.\n",
                i, nrow(jobs), analysis_counts["complete"], analysis_counts["failed"],
                analysis_counts["skipped"]))
    completed_since_checkpoint <- 0L
  }
}
write_manifest()
cat(sprintf("ANALYSIS SUMMARY: %d calculated, %d failed, %d existing results skipped.\n",
            analysis_counts["complete"], analysis_counts["failed"], analysis_counts["skipped"]))
cat("Finished analysis shard ", options$shard_index, " of ",
    options$shard_count, ". Manifest: ", manifest_path, "\n", sep = "")
