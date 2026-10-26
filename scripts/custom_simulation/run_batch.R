#!/usr/bin/env Rscript
# Reproducibly generate and checkpoint a manageable custom-simulation suite.
# Run from any working directory with Rscript scripts/custom_simulation/run_batch.R.

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
setwd(project_root)

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(args) {
  opts <- list(
    reference = "data/DGElist_BRCA.rds",
    output = "results/custom_simulation",
    shard_index = 1L,
    shard_count = 1L,
    replicates = 1L,
    n_var = 10000L,
    base_seed = 20261024L,
    restore = TRUE,
    plan_only = FALSE
  )
  if ("--help" %in% args || "-h" %in% args) {
    cat(paste(
      "Usage: Rscript scripts/custom_simulation/run_batch.R [options]",
      "",
      "Options:",
      "  --reference PATH      Real edgeR DGEList RDS (default: data/DGElist_BRCA.rds)",
      "  --output PATH         Output directory (default: results/custom_simulation)",
      "  --shard-index N       This shard, starting at 1 (default: 1)",
      "  --shard-count N       Total shards/workers (default: 1)",
      "  --replicates N        Replicates per setting (default: 1)",
      "  --n-var N             Genes per simulated dataset (default: 10000)",
      "  --base-seed N         Base random seed (default: 20261024)",
      "  --no-restore          Do not run renv::restore() before generating data",
      "  --plan-only           Print the assigned settings and exit",
      "",
      "Settings use one-factor-at-a-time comparisons around a baseline, with a",
      "separate Poisson baseline and Poisson sensitivity runs. Shards are",
      "assigned by estimated workload so large sample-size scenarios spread out.",
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
    } else if (key %in% c("--reference", "--output", "--shard-index", "--shard-count",
                          "--replicates", "--n-var", "--base-seed")) {
      if (i == length(args)) stop("Missing value for ", key)
      value <- args[i + 1L]
      name <- switch(key,
        "--reference" = "reference", "--output" = "output",
        "--shard-index" = "shard_index", "--shard-count" = "shard_count",
        "--replicates" = "replicates", "--n-var" = "n_var",
        "--base-seed" = "base_seed"
      )
      opts[[name]] <- if (name %in% c("shard_index", "shard_count", "replicates", "n_var", "base_seed")) {
        suppressWarnings(as.integer(value))
      } else value
      i <- i + 2L
    } else {
      stop("Unknown option: ", key)
    }
  }
  if (anyNA(unlist(opts[c("shard_index", "shard_count", "replicates", "n_var", "base_seed")]))) {
    stop("Numeric options must be whole numbers.")
  }
  if (opts$shard_count < 1L || opts$shard_index < 1L ||
      opts$shard_index > opts$shard_count || opts$replicates < 1L || opts$n_var < 1L) {
    stop("Shard values, replicates, and n-var must be positive; shard-index cannot exceed shard-count.")
  }
  opts
}

options <- parse_args(args)
if (is.null(options)) quit(status = 0L)

if (!file.exists("renv/activate.R")) stop("renv/activate.R is missing from the project.")
source("renv/activate.R")
if (options$restore) renv::restore(prompt = FALSE)
if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop("edgeR is unavailable. Restore the renv environment or install the locked packages.")
}
source("scripts/shared/simulate_custom.R")

build_scenarios <- function(n_var) {
  # One-factor-at-a-time design: each scenario changes one setting from the
  # model-specific baseline. This keeps the first run tractable and interpretable.
  baseline <- list(
    s = 10L, dispersion_scale = 1, library_size_sdlog = 0,
    de_fraction = 0.30, fold_change_min = 2,
    upregulated_fraction = 0.50, outlier_fraction = 0,
    sparsity = 0, batch_fraction = 0.50, batch_effect = 1
  )
  rows <- list()
  add_row <- function(model, id, factor_name, values) {
    row <- c(list(scenario_id = id, count_model = model, varied_factor = factor_name),
             utils::modifyList(baseline, values))
    row$n_diffexp <- as.integer(round(n_var * row$de_fraction))
    rows[[length(rows) + 1L]] <<- as.data.frame(row, stringsAsFactors = FALSE)
  }
  add_family <- function(model) {
    add_row(model, paste0(model, "_baseline"), "baseline", list())
    variations <- list(
      sample_size = list(s = c(3L, 5L, 30L, 50L)),
      library_size_sdlog = list(library_size_sdlog = c(0.25, 0.50)),
      de_fraction = list(de_fraction = c(0.10, 0.50)),
      fold_change_min = list(fold_change_min = c(1.25, 4)),
      upregulated_fraction = list(upregulated_fraction = c(0.25, 0.75)),
      outlier_fraction = list(outlier_fraction = c(0.001, 0.01)),
      sparsity = list(sparsity = c(0.10, 0.30)),
      batch_effect = list(batch_effect = c(1.25, 1.50))
    )
    if (model == "negative_binomial") {
      variations$dispersion_scale <- list(dispersion_scale = c(0.5, 2))
    }
    for (factor_name in names(variations)) {
      setting <- names(variations[[factor_name]])[1]
      for (value in variations[[factor_name]][[1]]) {
        add_row(model, paste(model, setting, format(value, trim = TRUE), sep = "_"),
                factor_name, setNames(list(value), setting))
      }
    }
  }
  add_family("negative_binomial")
  add_family("poisson")

  scenarios <- do.call(rbind, rows)
  rownames(scenarios) <- NULL
  scenarios$scenario_index <- seq_len(nrow(scenarios))
  scenarios$n_var <- n_var
  scenarios
}

assign_shards <- function(scenarios, replicates, shard_count, n_var) {
  jobs <- scenarios[rep(seq_len(nrow(scenarios)), each = replicates), , drop = FALSE]
  jobs$replicate <- rep(seq_len(replicates), times = nrow(scenarios))
  jobs$job_index <- seq_len(nrow(jobs))
  jobs$estimated_work <- as.numeric(jobs$s) * n_var
  jobs$shard <- NA_integer_
  load <- rep(0, shard_count)
  # Longest jobs first; tie-breaking by job index makes assignments repeatable.
  order_jobs <- order(-jobs$estimated_work, jobs$job_index)
  for (job in order_jobs) {
    target <- which(load == min(load))[1]
    jobs$shard[job] <- target
    load[target] <- load[target] + jobs$estimated_work[job]
  }
  jobs[order(jobs$job_index), , drop = FALSE]
}

atomic_write_csv <- function(x, path) {
  temporary <- paste0(path, ".tmp")
  utils::write.csv(x, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) stop("Could not replace manifest: ", path)
}

resolve_path <- function(path) {
  is_absolute <- grepl("^(~|/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
  if (is_absolute) path.expand(path) else file.path(project_root, path)
}

reference_path <- resolve_path(options$reference)
if (!file.exists(reference_path)) stop("Reference RDS not found: ", reference_path)
output_root <- resolve_path(options$output)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
shard_label <- sprintf("shard-%03d-of-%03d", options$shard_index, options$shard_count)
shard_dir <- file.path(output_root, shard_label)
dir.create(shard_dir, recursive = TRUE, showWarnings = FALSE)

jobs <- assign_shards(build_scenarios(options$n_var), options$replicates,
                      options$shard_count, options$n_var)
jobs <- jobs[jobs$shard == options$shard_index, , drop = FALSE]
manifest_path <- file.path(shard_dir, "manifest.csv")
if (options$plan_only) {
  print(jobs[, c("job_index", "scenario_id", "replicate", "count_model", "s",
                 "dispersion_scale", "library_size_sdlog", "de_fraction",
                 "fold_change_min", "upregulated_fraction", "outlier_fraction",
                 "sparsity", "batch_effect", "estimated_work", "shard")])
  cat("\n", nrow(jobs), " jobs assigned to ", shard_label, "\n", sep = "")
  quit(status = 0L)
}
if (nrow(jobs) == 0L) {
  cat("No jobs assigned to ", shard_label, ".\n", sep = "")
  quit(status = 0L)
}

reference_hash <- unname(tools::md5sum(reference_path))
cache_path <- file.path(output_root, paste0("parameters-", reference_hash, ".rds"))
load_or_estimate_parameters <- function(cache_path, reference_path, reference_hash) {
  if (file.exists(cache_path)) return(readRDS(cache_path)$parameters)
  lock_dir <- paste0(cache_path, ".lock")
  owns_lock <- FALSE
  for (attempt in seq_len(1800L)) {
    if (file.exists(cache_path)) return(readRDS(cache_path)$parameters)
    owns_lock <- dir.create(lock_dir, showWarnings = FALSE)
    if (owns_lock) break
    lock_info <- file.info(lock_dir)
    lock_age <- as.numeric(difftime(Sys.time(), lock_info$mtime, units = "secs"))
    if (is.finite(lock_age) && lock_age > 3600) {
      unlink(lock_dir, recursive = TRUE, force = TRUE)
    } else {
      Sys.sleep(2)
    }
  }
  if (!owns_lock) stop("Timed out waiting for the shared parameter cache lock.")
  on.exit(unlink(lock_dir, recursive = TRUE, force = TRUE), add = TRUE)
  if (file.exists(cache_path)) return(readRDS(cache_path)$parameters)

  reference_dge <- readRDS(reference_path)
  parameters <- estimate_simulation_parameters(reference_dge)
  temporary_cache <- paste0(cache_path, ".", Sys.getpid(), ".tmp")
  saveRDS(list(reference = reference_path, reference_hash = reference_hash,
               parameters = parameters), temporary_cache)
  if (!file.rename(temporary_cache, cache_path)) stop("Could not write parameter cache.")
  parameters
}
simulation_params <- load_or_estimate_parameters(cache_path, reference_path, reference_hash)

if (file.exists(manifest_path)) {
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE,
                              check.names = FALSE)
} else {
  manifest <- NULL
}
manifest_rows <- lapply(seq_len(nrow(jobs)), function(i) {
  job <- jobs[i, , drop = FALSE]
  seed <- as.integer((as.double(options$base_seed) + job$job_index * 1009 + job$replicate) %% .Machine$integer.max)
  output_name <- sprintf("%s-rep%03d-n%d-seed%d-ref%s.rds", job$scenario_id,
                         job$replicate, options$n_var, seed, substr(reference_hash, 1, 8))
  record <- job
  record$seed <- seed
  record$reference <- reference_path
  record$output_file <- file.path(shard_dir, output_name)
  record$status <- "pending"
  record$error <- ""
  if (!is.null(manifest)) {
    previous <- manifest[manifest$scenario_id == job$scenario_id &
      manifest$replicate == job$replicate, , drop = FALSE]
    if (nrow(previous)) {
      record$status <- previous$status[1]
      if ("error" %in% names(previous)) record$error <- previous$error[1]
    }
  }
  record
})
manifest_extras <- if (!is.null(manifest)) {
  manifest[!manifest$scenario_id %in% jobs$scenario_id, , drop = FALSE]
} else NULL
write_manifest <- function() {
  rows <- manifest_rows
  if (!is.null(manifest_extras) && nrow(manifest_extras)) {
    rows <- c(rows, list(manifest_extras))
  }
  atomic_write_csv(do.call(rbind, rows), manifest_path)
}

for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, , drop = FALSE]
  record <- manifest_rows[[i]]
  output_path <- record$output_file
  output_name <- basename(output_path)
  seed <- record$seed
  if (file.exists(output_path)) {
    record$status <- "complete"
    record$error <- ""
    manifest_rows[[i]] <- record
    write_manifest()
    cat("SKIP ", output_name, " (already complete)\n", sep = "")
    next
  }

  record$status <- "running"
  record$error <- ""
  manifest_rows[[i]] <- record
  write_manifest()

  cat(sprintf("RUN  [%d/%d] %s (%s, s=%d, rep=%d)\n",
              i, nrow(jobs), job$scenario_id, job$count_model, job$s, job$replicate))
  tryCatch({
    simulation <- simulate_custom_counts(
      params = simulation_params,
      s = job$s,
      n_var = options$n_var,
      n_diffexp = job$n_diffexp,
      count_model = job$count_model,
      dispersion_scale = if (job$count_model == "negative_binomial") job$dispersion_scale else 1,
      library_size_sdlog = job$library_size_sdlog,
      fold_change_min = job$fold_change_min,
      upregulated_fraction = job$upregulated_fraction,
      outlier_fraction = job$outlier_fraction,
      sparsity = job$sparsity,
      batch_fraction = job$batch_fraction,
      batch_effect = job$batch_effect,
      seed = seed
    )
    temporary_output <- paste0(output_path, ".tmp")
    saveRDS(list(
      scenario = job,
      replicate = job$replicate,
      seed = seed,
      reference = reference_path,
      reference_hash = reference_hash,
      simulated = simulation
    ), temporary_output, compress = "xz")
    if (!file.rename(temporary_output, output_path)) stop("Could not finalize dataset file.")
    record$status <- "complete"
    cat("DONE ", output_name, "\n", sep = "")
  }, error = function(e) {
    record$status <<- "failed"
    record$error <<- conditionMessage(e)
    cat("ERROR ", job$scenario_id, ": ", conditionMessage(e), "\n", sep = "")
  })
  manifest_rows[[i]] <- record
  write_manifest()
}

write_manifest()
cat("Finished ", shard_label, ". Manifest: ", manifest_path, "\n", sep = "")
