#!/usr/bin/env Rscript
# Reproducibly generate and checkpoint the full-factorial simulation suite.
# Run from any working directory with Rscript scripts/custom_simulation/run_batch.R.

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
setwd(project_root)

simulation_version <- "3.0"

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(args) {
  opts <- list(
    reference = "data/DGElist_BRCA.rds",
    reference_mode = "brca",
    references = c(
      "data/DGElist_BRCA.rds",
      "data/DGElist_E-ENAD-34.rds",
      "data/DGElist_GSE150910.rds",
      "data/DGElist_LUSC.rds"
    ),
    output = paste0("results/custom_simulation_v", simulation_version),
    shard_index = 1L,
    shard_count = 1L,
    replicates = 1L,
    de_mode = "fixed",
    de_proportions = c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75),
    n_var = 10000L,
    base_seed = 20261024L,
    dispersion_estimator = "method_of_moments",
    include_imbalance = TRUE,
    restore = TRUE,
    plan_only = FALSE
  )
  if ("--version" %in% args) {
    cat("Custom simulation version", simulation_version, "\n")
    return(NULL)
  }
  if ("--help" %in% args || "-h" %in% args) {
    cat(paste(
      "Usage: Rscript scripts/custom_simulation/run_batch.R [options]",
      "",
      "Options:",
      "  --reference PATH      Real edgeR DGEList RDS (default: data/DGElist_BRCA.rds)",
      "  --reference-mode MODE brca (default), all, or random",
      "  --references PATHS    Comma-separated RDS paths for all/random modes",
      paste0("  --output PATH         Output directory (default: results/custom_simulation_v", simulation_version, ")"),
      "  --version             Print simulation version and exit",
      "  --shard-index N       This shard, starting at 1 (default: 1)",
      "  --shard-count N       Total shards/workers (default: 1)",
      "  --replicates N        Replicates per setting (default: 1)",
      "  --de-mode MODE        1/combination, 2/setting, 3/fixed (default)",
      "  --de-proportions LIST Comma-separated DE proportions (default: 0.1,0.2,0.3,0.4,0.5,0.6,0.75)",
      "  --n-var N             Genes per simulated dataset (default: 10000)",
      "  --base-seed N         Base random seed (default: 20261024)",
      "  --dispersion-estimator METHOD  method_of_moments (default) or edgeR",
      "  --include-imbalance   Include 1:3 and 1:7 splits at each size (default)",
      "  --no-imbalance        Generate balanced groups only",
      "  --no-restore          Do not run renv::restore() before generating data",
      "  --plan-only           Print the assigned settings and exit",
      "",
      "Settings use a full factorial design across both count models. Shards are",
      "assigned by estimated workload so large sample-size scenarios spread out.",
      sep = "\n"
    ), "\n")
    return(NULL)
  }
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--no-imbalance") {
      opts$include_imbalance <- FALSE
      i <- i + 1L
    } else if (key == "--include-imbalance") {
      opts$include_imbalance <- TRUE
      i <- i + 1L
    } else if (key == "--no-restore") {
      opts$restore <- FALSE
      i <- i + 1L
    } else if (key == "--plan-only") {
      opts$plan_only <- TRUE
      i <- i + 1L
    } else if (key %in% c("--reference", "--reference-mode", "--references", "--output", "--shard-index", "--shard-count",
                          "--replicates", "--n-var", "--base-seed", "--dispersion-estimator",
                          "--de-mode", "--de-proportions")) {
      if (i == length(args)) stop("Missing value for ", key)
      value <- args[i + 1L]
      name <- switch(key,
        "--reference" = "reference", "--output" = "output",
        "--reference-mode" = "reference_mode", "--references" = "references",
        "--shard-index" = "shard_index", "--shard-count" = "shard_count",
        "--de-mode" = "de_mode", "--de-proportions" = "de_proportions",
        "--replicates" = "replicates", "--n-var" = "n_var",
        "--base-seed" = "base_seed", "--dispersion-estimator" = "dispersion_estimator"
      )
      opts[[name]] <- if (name == "references") {
        trimws(strsplit(value, ",", fixed = TRUE)[[1]])
      } else if (name == "de_proportions") {
        suppressWarnings(as.numeric(trimws(strsplit(value, ",", fixed = TRUE)[[1]])))
      } else if (name %in% c("shard_index", "shard_count", "replicates", "n_var", "base_seed")) {
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
  if (!opts$reference_mode %in% c("brca", "all", "random")) {
    stop("--reference-mode must be one of: brca, all, random.")
  }
  if (!opts$dispersion_estimator %in% c("method_of_moments", "edgeR")) {
    stop("--dispersion-estimator must be method_of_moments or edgeR.")
  }
  if (!length(opts$references) || any(!nzchar(opts$references))) {
    stop("--references must contain one or more comma-separated paths.")
  }
  aliases <- c("1" = "combination", "2" = "setting", "3" = "fixed")
  if (opts$de_mode %in% names(aliases)) opts$de_mode <- unname(aliases[opts$de_mode])
  if (!opts$de_mode %in% c("combination", "setting", "fixed")) {
    stop("--de-mode must be 1/combination, 2/setting, or 3/fixed.")
  }
  if (!length(opts$de_proportions) || any(!is.finite(opts$de_proportions)) ||
      any(opts$de_proportions <= 0 | opts$de_proportions >= 1) ||
      anyDuplicated(opts$de_proportions)) {
    stop("--de-proportions must contain unique finite proportions strictly between 0 and 1.")
  }
  opts
}

options <- parse_args(args)
if (is.null(options)) quit(status = 0L)

build_scenarios <- function(n_var, include_imbalance = TRUE) {
  # s denotes the balanced size; each split preserves the total of 2*s.
  sizes <- do.call(rbind, lapply(c(3L, 5L, 10L, 30L, 50L), function(s) {
    fractions <- if (include_imbalance) c(0.5, 0.25, 0.125) else 0.5
    n1 <- as.integer(floor(2 * s * fractions + 0.5))
    result <- data.frame(s = s, n_group1 = n1, n_group2 = 2L * s - n1,
                         minority_fraction_target = fractions)
    result <- result[result$n_group1 >= 3L & result$n_group2 >= 3L, , drop = FALSE]
    result[!duplicated(result[c("n_group1", "n_group2")]), , drop = FALSE]
  }))
  make_family <- function(model) {
    grid <- expand.grid(
      size_index = seq_len(nrow(sizes)),
      fold_change_min = c(1, 2, 4),
      upregulated_fraction = c(0.25, 0.50, 0.75),
      library_size_sdlog = c(0, 0.25, 0.50),
      outlier_fraction = c(0, 0.001, 0.01),
      sparsity = c(0, 0.10, 0.30),
      batch_effect = c(1, 1.5),
      dispersion_scale = if (model == "negative_binomial") c(0.5, 1, 2) else 1,
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    grid <- cbind(grid, sizes[grid$size_index, , drop = FALSE])
    grid$size_index <- NULL
    grid$n_samples <- grid$n_group1 + grid$n_group2
    grid$count_model <- model
    grid$design <- "full_factorial"
    grid$sample_size_baseline <- grid$s
    grid$s[grid$n_group1 != grid$n_group2] <- NA_integer_
    grid$varied_factor <- "full_factorial"
    grid$batch_fraction <- 0.5
    grid$de_fraction <- 0.3 # Replaced by DE selection before simulation.
    grid$n_diffexp <- as.integer(round(n_var * grid$de_fraction))
    # A setting fixes all non-sample factors; mode 2 shares DE across sizes.
    grid$de_setting_id <- do.call(paste, c(grid[c(
      "fold_change_min", "upregulated_fraction", "library_size_sdlog",
      "outlier_fraction", "sparsity", "batch_effect", "dispersion_scale"
    )], sep = "_"))
    grid$scenario_id <- paste0(model, "_factorial_g", grid$n_group1, "-", grid$n_group2,
      "_fc", grid$fold_change_min, "_up", grid$upregulated_fraction,
      "_lib", grid$library_size_sdlog, "_out", grid$outlier_fraction,
      "_zero", grid$sparsity, "_batch", grid$batch_effect, "_disp", grid$dispersion_scale)
    grid
  }
  scenarios <- rbind(make_family("negative_binomial"), make_family("poisson"))
  rownames(scenarios) <- NULL
  scenarios$scenario_index <- seq_len(nrow(scenarios))
  scenarios$simulation_version <- simulation_version
  scenarios$n_var <- n_var
  scenarios
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

reference_specs <- if (options$reference_mode == "brca") {
  options$reference
} else {
  options$references
}
reference_paths <- vapply(reference_specs, resolve_path, character(1))
if (anyDuplicated(reference_paths)) stop("Reference paths must be unique.")
missing_references <- reference_paths[!file.exists(reference_paths)]
if (length(missing_references)) {
  stop("Reference RDS not found: ", paste(missing_references, collapse = ", "))
}
reference_ids <- tools::file_path_sans_ext(basename(reference_paths))
if (any(tolower(reference_ids) == "sclc_cell2024")) {
  stop("DGElist_SCLC_Cell2024.rds is reserved as the independent external test set and cannot be used as a simulation reference.")
}
if (anyDuplicated(reference_ids)) stop("Reference filenames must be unique across --references.")
reference_table <- data.frame(
  reference_id = reference_ids,
  reference = reference_paths,
  reference_hash = unname(tools::md5sum(reference_paths)),
  stringsAsFactors = FALSE
)
output_root <- resolve_path(options$output)
shard_label <- sprintf("shard-%03d-of-%03d", options$shard_index, options$shard_count)
shard_dir <- file.path(output_root, shard_label)

scenarios <- build_scenarios(options$n_var, options$include_imbalance)
jobs <- scenarios[rep(seq_len(nrow(scenarios)), each = options$replicates), , drop = FALSE]
jobs$replicate <- rep(seq_len(options$replicates), times = nrow(scenarios))
source("scripts/custom_simulation/de_proportion_modes.R")
jobs <- assign_de_proportions(jobs, options$de_mode, options$de_proportions,
                              options$base_seed, options$n_var)
jobs$setting_job_index <- seq_len(nrow(jobs))
if (options$reference_mode == "all") {
  jobs <- jobs[rep(seq_len(nrow(jobs)), each = nrow(reference_table)), , drop = FALSE]
  jobs$reference_index <- rep(seq_len(nrow(reference_table)), times = nrow(jobs) / nrow(reference_table))
} else if (options$reference_mode == "random") {
  jobs$reference_index <- vapply(jobs$setting_job_index, function(index) {
    selection_seed <- as.integer((as.double(options$base_seed) + index * 104729) %% .Machine$integer.max)
    set.seed(selection_seed)
    sample.int(nrow(reference_table), 1L)
  }, integer(1))
} else {
  jobs$reference_index <- 1L
}
jobs$reference_id <- reference_table$reference_id[jobs$reference_index]
jobs$reference <- reference_table$reference[jobs$reference_index]
jobs$reference_hash <- reference_table$reference_hash[jobs$reference_index]
jobs$job_index <- seq_len(nrow(jobs))
jobs$estimated_work <- as.numeric(jobs$n_samples) * options$n_var
jobs$shard <- NA_integer_
shard_load <- rep(0, options$shard_count)
for (job_index in order(-jobs$estimated_work, jobs$job_index)) {
  target <- which(shard_load == min(shard_load))[1]
  jobs$shard[job_index] <- target
  shard_load[target] <- shard_load[target] + jobs$estimated_work[job_index]
}
jobs <- jobs[order(jobs$job_index), , drop = FALSE]
jobs <- jobs[jobs$shard == options$shard_index, , drop = FALSE]
manifest_path <- file.path(shard_dir, "manifest.csv")
if (options$plan_only) {
  cat("Simulation version:", simulation_version, "\n")
  cat("DE proportion mode:", options$de_mode, "\n")
  cat("Dispersion estimator:", options$dispersion_estimator, "\n")
  print(jobs[, c("job_index", "scenario_id", "replicate", "count_model", "sample_size_baseline", "s", "n_group1", "n_group2", "n_samples",
                 "reference_id",
                 "dispersion_scale", "library_size_sdlog", "de_mode", "de_fraction",
                 "fold_change_min", "upregulated_fraction", "outlier_fraction",
                 "sparsity", "batch_effect", "estimated_work", "shard")])
  cat("\n", nrow(jobs), " jobs assigned to ", shard_label, "\n", sep = "")
  quit(status = 0L)
}
if (nrow(jobs) == 0L) {
  cat("No jobs assigned to ", shard_label, ".\n", sep = "")
  quit(status = 0L)
}

if (!file.exists("renv/activate.R")) stop("renv/activate.R is missing from the project.")
source("renv/activate.R")
if (options$restore) renv::restore(prompt = FALSE)
if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop("edgeR is unavailable. Restore the renv environment or install the locked packages.")
}
source("scripts/shared/simulate_custom.R")

dir.create(shard_dir, recursive = TRUE, showWarnings = FALSE)

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
parameter_cache <- new.env(parent = emptyenv())
get_reference_parameters <- function(reference_path, reference_hash) {
  if (exists(reference_hash, envir = parameter_cache, inherits = FALSE)) {
    return(get(reference_hash, envir = parameter_cache, inherits = FALSE))
  }
  cache_path <- file.path(output_root, paste0("parameters-dispersion-v2-", reference_hash, ".rds"))
  parameters <- load_or_estimate_parameters(cache_path, reference_path, reference_hash)
  assign(reference_hash, parameters, envir = parameter_cache)
  parameters
}

if (file.exists(manifest_path)) {
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE,
                              check.names = FALSE)
} else {
  manifest <- NULL
}
manifest_rows <- lapply(seq_len(nrow(jobs)), function(i) {
  job <- jobs[i, , drop = FALSE]
  seed <- as.integer((as.double(options$base_seed) + job$job_index * 1009 + job$replicate) %% .Machine$integer.max)
  reference_id <- job$reference_id
  reference_path <- job$reference
  reference_hash <- job$reference_hash
  output_name <- sprintf("%s-%s-disp-%s-rep%03d-n%d-seed%d-ref%s.rds", job$scenario_id,
                         reference_id, ifelse(options$dispersion_estimator == "method_of_moments", "mom", "edgeR"),
                         job$replicate, options$n_var, seed,
                         substr(reference_hash, 1, 8))
  record <- job
  record$seed <- seed
  record$dispersion_estimator <- options$dispersion_estimator
  record$reference <- reference_path
  record$reference_id <- reference_id
  record$reference_hash <- reference_hash
  record$output_file <- file.path(shard_dir, output_name)
  record$status <- "pending"
  record$error <- ""
  if (!is.null(manifest)) {
    if ("reference_id" %in% names(manifest)) {
      previous <- manifest[manifest$scenario_id == job$scenario_id &
        manifest$replicate == job$replicate &
        manifest$reference_id == reference_id, , drop = FALSE]
    } else {
      previous <- manifest[FALSE, , drop = FALSE]
    }
    if (nrow(previous)) {
      record$status <- previous$status[1]
      if ("error" %in% names(previous)) record$error <- previous$error[1]
    }
  }
  record
})
write_manifest <- function() {
  atomic_write_csv(do.call(rbind, manifest_rows), manifest_path)
}

for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, , drop = FALSE]
  record <- manifest_rows[[i]]
  output_path <- record$output_file
  output_name <- basename(output_path)
  seed <- record$seed
  reference_id <- job$reference_id
  reference_path <- job$reference
  reference_hash <- job$reference_hash
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

  cat(sprintf("RUN  [%d/%d] %s (%s, groups=%d/%d, rep=%d, reference=%s)\n",
              i, nrow(jobs), job$scenario_id, job$count_model, job$n_group1, job$n_group2,
              job$replicate, reference_id))
  tryCatch({
    simulation_params <- get_reference_parameters(reference_path, reference_hash)
    simulation <- simulate_custom_counts(
      params = simulation_params,
      group_sizes = c(job$n_group1, job$n_group2),
      n_var = options$n_var,
      n_diffexp = job$n_diffexp,
      count_model = job$count_model,
      dispersion_estimator = options$dispersion_estimator,
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
      simulation_version = simulation_version,
      scenario = job,
      replicate = job$replicate,
      seed = seed,
      reference = reference_path,
      reference_id = reference_id,
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
