#!/usr/bin/env Rscript
# Generate ten small simulation datasets and run every registered DE method.
# Baseline methods run on all ten datasets; quick one-factor method-setting
# alternatives run on the first dataset. Results are checkpointed as they finish.

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
setwd(project_root)

parse_args <- function(args) {
  opts <- list(
    reference = "data/DGElist_E-ENAD-34.rds",
    output = "results/custom_simulation_smoke",
    n_var = 500L,
    base_seed = 20261024L,
    restore = TRUE
  )
  if ("--help" %in% args || "-h" %in% args) {
    cat(paste(
      "Usage: Rscript scripts/custom_simulation/smoke_test_10.R [options]",
      "  --reference PATH   Two-group edgeR DGEList RDS (default: E-ENAD-34)",
      "  --output PATH      Test output directory",
      "  --n-var N          Genes per dataset (default: 500)",
      "  --base-seed N      Reproducible base seed",
      "  --no-restore       Skip renv::restore()",
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
    } else if (key %in% c("--reference", "--output", "--n-var", "--base-seed")) {
      if (i == length(args)) stop("Missing value for ", key)
      value <- args[i + 1L]
      name <- switch(key, "--reference" = "reference", "--output" = "output",
                     "--n-var" = "n_var", "--base-seed" = "base_seed")
      opts[[name]] <- if (name %in% c("n_var", "base_seed")) {
        suppressWarnings(as.integer(value))
      } else value
      i <- i + 2L
    } else stop("Unknown option: ", key)
  }
  if (is.na(opts$n_var) || is.na(opts$base_seed) || opts$n_var < 1L) {
    stop("n-var and base-seed must be whole numbers; n-var must be positive.")
  }
  opts
}

config <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(config)) quit(status = 0L)
source("renv/activate.R")
if (config$restore) renv::restore(prompt = FALSE)
source("scripts/shared/simulate_custom.R")
source("scripts/custom_simulation/method_settings.R")

resolve_path <- function(path) {
  is_absolute <- grepl("^(~|/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
  if (is_absolute) path.expand(path) else file.path(project_root, path)
}
reference_path <- resolve_path(config$reference)
output_root <- resolve_path(config$output)
if (!file.exists(reference_path)) stop("Reference data not found: ", reference_path)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dataset_dir <- file.path(output_root, "datasets")
result_dir <- file.path(output_root, "method_results")
dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

scenarios <- data.frame(
  scenario_id = c(
    "nb_baseline", "nb_sample5", "nb_high_dispersion", "nb_outlier_sparse",
    "nb_weaker_effect", "poisson_baseline", "poisson_sample5",
    "poisson_library_spread", "poisson_direction_skew", "poisson_batch_sparse"
  ),
  count_model = c(rep("negative_binomial", 5), rep("poisson", 5)),
  s = c(3, 5, 3, 3, 5, 3, 5, 3, 3, 5),
  de_fraction = c(0.30, 0.30, 0.30, 0.20, 0.10, 0.30, 0.30, 0.30, 0.10, 0.20),
  dispersion_scale = c(1, 1, 2, 1, 0.5, 1, 1, 1, 1, 1),
  library_size_sdlog = c(0, 0, 0, 0, 0, 0, 0, 0.3, 0, 0.2),
  fold_change_min = c(2, 2, 2, 2, 1.25, 2, 2, 2, 1.5, 2),
  upregulated_fraction = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.8, 0.5),
  outlier_fraction = c(0, 0, 0, 0.01, 0, 0, 0, 0, 0, 0.005),
  sparsity = c(0, 0, 0, 0.10, 0, 0, 0, 0, 0, 0.10),
  batch_effect = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1.5),
  batch_fraction = rep(0.5, 10),
  stringsAsFactors = FALSE
)
scenarios$n_var <- config$n_var
scenarios$n_diffexp <- as.integer(round(config$n_var * scenarios$de_fraction))

reference_hash <- unname(tools::md5sum(reference_path))
parameter_path <- file.path(output_root, paste0("parameters-", reference_hash, ".rds"))
if (file.exists(parameter_path)) {
  simulation_params <- readRDS(parameter_path)$parameters
} else {
  reference_dge <- readRDS(reference_path)
  simulation_params <- estimate_simulation_parameters(reference_dge)
  saveRDS(list(reference = reference_path, reference_hash = reference_hash,
               parameters = simulation_params), parameter_path)
}

simulation_files <- character(nrow(scenarios))
for (i in seq_len(nrow(scenarios))) {
  scenario <- scenarios[i, , drop = FALSE]
  seed <- as.integer(config$base_seed + i * 101L)
  simulation_path <- file.path(dataset_dir,
    sprintf("%02d-%s-n%d-seed%d.rds", i, scenario$scenario_id, config$n_var, seed))
  simulation_files[i] <- simulation_path
  if (file.exists(simulation_path)) next
  set.seed(seed)
  simulated <- simulate_custom_counts(
    params = simulation_params,
    s = scenario$s,
    n_var = scenario$n_var,
    n_diffexp = scenario$n_diffexp,
    count_model = scenario$count_model,
    dispersion_scale = scenario$dispersion_scale,
    library_size_sdlog = scenario$library_size_sdlog,
    fold_change_min = scenario$fold_change_min,
    upregulated_fraction = scenario$upregulated_fraction,
    outlier_fraction = scenario$outlier_fraction,
    sparsity = scenario$sparsity,
    batch_fraction = scenario$batch_fraction,
    batch_effect = scenario$batch_effect,
    seed = seed
  )
  tmp <- paste0(simulation_path, ".tmp")
  saveRDS(list(scenario = scenario, seed = seed, reference = reference_path,
               reference_hash = reference_hash, simulated = simulated), tmp,
          compress = "xz")
  if (!file.rename(tmp, simulation_path)) stop("Could not save ", simulation_path)
  cat("GENERATED ", basename(simulation_path), "\n", sep = "")
}

all_settings <- make_method_settings("core")
baseline_settings <- Filter(function(x) x$baseline, all_settings)
quick_settings <- make_method_settings("quick")
extra_settings <- Filter(function(x) !x$baseline, quick_settings)
jobs <- list()
for (i in seq_along(simulation_files)) {
  for (setting in baseline_settings) {
    jobs[[length(jobs) + 1L]] <- list(dataset_index = i, setting = setting)
  }
}
for (setting in extra_settings) {
  jobs[[length(jobs) + 1L]] <- list(dataset_index = 1L, setting = setting)
}

status_path <- file.path(output_root, "method_status.csv")
status <- if (file.exists(status_path)) {
  utils::read.csv(status_path, stringsAsFactors = FALSE, check.names = FALSE)
} else data.frame()
write_status <- function(row) {
  key <- paste(row$dataset_id, row$method, row$setting_id, sep = "::")
  if (nrow(status)) {
    existing_key <- paste(status$dataset_id, status$method, status$setting_id, sep = "::")
    status <<- status[existing_key != key, , drop = FALSE]
  }
  status <<- rbind(status, row)
  tmp <- paste0(status_path, ".tmp")
  utils::write.csv(status, tmp, row.names = FALSE, na = "")
  if (!file.rename(tmp, status_path)) stop("Could not update status file.")
}

for (job_index in seq_along(jobs)) {
  job <- jobs[[job_index]]
  dataset_index <- job$dataset_index
  setting <- job$setting
  scenario <- scenarios[dataset_index, , drop = FALSE]
  dataset_id <- sub("\\.rds$", "", basename(simulation_files[dataset_index]))
  setting_id <- setting$setting_id
  target_dir <- file.path(result_dir, dataset_id)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(target_dir,
    paste0(setting$method, "--", setting_id, ".rds"))
  if (file.exists(result_path)) next

  seed <- as.integer(config$base_seed + dataset_index * 101L + job_index * 17L)
  cat(sprintf("TEST [%d/%d] %s / %s / %s\n", job_index, length(jobs),
              scenario$scenario_id, setting$method, setting_id))
  record <- data.frame(
    dataset_id = dataset_id,
    scenario_id = scenario$scenario_id,
    count_model = scenario$count_model,
    method = setting$method,
    setting_id = setting_id,
    setting_args = paste(names(setting$args), vapply(setting$args, function(x) {
      if (is.null(x)) "NULL" else paste(x, collapse = ",")
    }, character(1)), sep = "=", collapse = ";"),
    seed = seed,
    status = "running",
    auc = NA_real_,
    true_positive_rate = NA_real_,
    observed_fdr = NA_real_,
    error = "",
    output_file = result_path,
    stringsAsFactors = FALSE
  )
  write_status(record)

  outcome <- tryCatch({
    dataset <- readRDS(simulation_files[dataset_index])
    dge <- dataset$simulated$dge
    truth <- dataset$simulated$truth$is_de[
      match(rownames(dge$counts), dataset$simulated$truth$gene)
    ]
    set.seed(seed)
    p <- run_method_setting(dge, setting, seed = seed)
    truth <- truth[match(names(p), rownames(dge$counts))]
    metrics <- summarize_method_result(p, truth)
    result <- list(
      dataset_id = dataset_id,
      scenario = dataset$scenario,
      method = setting$method,
      setting_id = setting_id,
      setting_args = setting$args,
      seed = seed,
      score_type = if (setting$method == "NOISeq") "1 - NOISeq probability" else "p-value",
      gene_results = data.frame(
        gene = names(p), score = p, is_de = truth,
        adjusted_score = stats::p.adjust(p, method = "BH")
      ),
      metrics = metrics
    )
    tmp <- paste0(result_path, ".tmp")
    saveRDS(result, tmp, compress = "xz")
    if (!file.rename(tmp, result_path)) stop("Could not finalize method result.")
    record$status <- "complete"
    record$auc <- metrics$auc
    record$true_positive_rate <- metrics$true_positive_rate
    record$observed_fdr <- metrics$observed_fdr
    record
  }, error = function(e) {
    record$status <<- "failed"
    record$error <<- conditionMessage(e)
    record
  })
  write_status(outcome)
}

failures <- if (nrow(status)) sum(status$status == "failed") else 0L
cat("Smoke test finished: ", nrow(status), " recorded runs, ",
    failures, " failed. Summary: ", status_path, "\n", sep = "")
if (failures > 0L) quit(status = 1L)
