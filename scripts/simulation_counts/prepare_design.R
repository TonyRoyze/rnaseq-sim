# Turn user settings and estimated real-data parameters into a simulation
# design: selected genes, expected-count multipliers, DE truth, and dispersions.
prepare_simulation_design <- function(
    params, s, n_var, n_diffexp, count_model, dispersion_scale,
    dispersion_source, dispersion_estimator, library_size_factors, library_size_sdlog,
    fold_change_min, fold_change_excess_rate, upregulated_fraction,
    batch, batch_fraction, batch_effect, group_sizes = NULL) {

  validate_simulation_scalar(s, "s", 1, Inf)
  validate_simulation_scalar(n_var, "n_var", 1, Inf)
  validate_simulation_scalar(n_diffexp, "n_diffexp", 0, n_var)
  if (any(c(s, n_var, n_diffexp) != floor(c(s, n_var, n_diffexp)))) {
    stop("s, n_var, and n_diffexp must be whole numbers.")
  }
  validate_simulation_scalar(dispersion_scale, "dispersion_scale", 0, Inf, FALSE)
  validate_simulation_scalar(library_size_sdlog, "library_size_sdlog", 0, Inf)
  validate_simulation_scalar(fold_change_min, "fold_change_min", 1, Inf)
  validate_simulation_scalar(fold_change_excess_rate, "fold_change_excess_rate", 0, Inf, FALSE)
  validate_simulation_scalar(upregulated_fraction, "upregulated_fraction", 0, 1)
  validate_simulation_scalar(batch_fraction, "batch_fraction", 0, 1)
  validate_simulation_scalar(batch_effect, "batch_effect", 0, Inf, FALSE)

  s <- as.integer(s)
  n_var <- as.integer(n_var)
  n_diffexp <- as.integer(n_diffexp)
  estimator_prefix <- if (dispersion_estimator == "edgeR") "edger" else "mom"
  needed <- "mean_total"
  if (count_model == "negative_binomial") {
    needed <- c(needed, paste0("dispersion_", estimator_prefix, "_total"))
  }
  if (!all(needed %in% names(params))) {
    stop("params must contain mean_total and estimates for the selected dispersion estimator.")
  }
  n_source <- length(params$mean_total)
  if (n_source < n_var) stop("n_var cannot exceed the number of available source genes.")
  if (count_model == "negative_binomial" &&
      length(params[[paste0("dispersion_", estimator_prefix, "_total")]]) != n_source) {
    stop("Simulation mean and dispersion vectors must have the same length.")
  }

  source_index <- sample.int(n_source, n_var)
  baseline_mean <- pmax(as.numeric(params$mean_total[source_index]), 0)
  dispersion_field <- function(group) {
    modern <- paste0("dispersion_", estimator_prefix, "_", group)
    legacy <- paste0("dispersion_", group)
    if (modern %in% names(params)) modern else legacy
  }
  if (dispersion_source == "total") {
    dispersion_unscaled <- params[[dispersion_field("total")]][source_index]
  } else {
    required <- c(dispersion_field("group1"), dispersion_field("group2"))
    if (!all(required %in% names(params))) {
      stop("dispersion_source='group' requires group-specific dispersions in params.")
    }
    dispersion_unscaled <- (params[[required[1]]][source_index] +
      params[[required[2]]][source_index]) / 2
  }
  dispersion_mom <- if (dispersion_source == "total") {
    params$dispersion_mom_total[source_index]
  } else (params$dispersion_mom_group1[source_index] +
    params$dispersion_mom_group2[source_index]) / 2
  dispersion_edger <- if (dispersion_source == "total") {
    params$dispersion_edger_total[source_index]
  } else (params$dispersion_edger_group1[source_index] +
    params$dispersion_edger_group2[source_index]) / 2
  dispersion_unscaled <- as.numeric(dispersion_unscaled)
  dispersion_used <- if (count_model == "negative_binomial") {
    pmax(dispersion_unscaled * dispersion_scale, 1e-8)
  } else rep(1, n_var)
  if (any(!is.finite(baseline_mean)) || any(!is.finite(dispersion_used))) {
    stop("Estimated means and dispersions must be finite.")
  }

  if (is.null(group_sizes)) group_sizes <- c(s, s)
  if (!is.numeric(group_sizes) || length(group_sizes) != 2L ||
      any(!is.finite(group_sizes)) || any(group_sizes < 1) ||
      any(group_sizes != floor(group_sizes))) {
    stop("group_sizes must contain two positive whole numbers.")
  }
  group_sizes <- as.integer(group_sizes)
  n_samples <- sum(group_sizes)
  condition <- factor(rep(c("group1", "group2"), times = group_sizes),
                      levels = c("group1", "group2"))
  if (is.null(library_size_factors)) {
    library_size_factors <- if (library_size_sdlog == 0) rep(1, n_samples) else
      exp(stats::rnorm(n_samples, sd = library_size_sdlog))
  }
  if (length(library_size_factors) != n_samples ||
      any(!is.finite(library_size_factors)) || any(library_size_factors <= 0)) {
    stop("library_size_factors must contain one positive finite value per sample.")
  }
  library_size_factors <- library_size_factors / exp(mean(log(library_size_factors)))

  if (is.null(batch)) {
    batch <- if (batch_effect == 1) rep("batch1", n_samples) else
      ifelse(stats::runif(n_samples) < batch_fraction, "batch2", "batch1")
  }
  if (length(batch) != n_samples || anyNA(batch)) {
    stop("batch must contain one non-missing label per sample.")
  }
  batch <- factor(batch)
  if (nlevels(batch) > 2L) stop("At most two batch levels are supported.")
  batch_multiplier <- ifelse(batch == levels(batch)[1], 1, batch_effect)

  de_index <- if (n_diffexp == 0L) integer(0) else sample.int(n_var, n_diffexp)
  direction <- rep("none", n_var)
  fold_change <- rep(1, n_var)
  n_up <- round(n_diffexp * upregulated_fraction)
  if (n_diffexp > 0L) {
    up_index <- de_index[seq_len(n_up)]
    down_index <- if (n_up < n_diffexp) de_index[seq.int(n_up + 1L, n_diffexp)] else integer(0)
    if (length(up_index)) {
      fold_change[up_index] <- fold_change_min + stats::rexp(length(up_index), fold_change_excess_rate)
      direction[up_index] <- "up"
    }
    if (length(down_index)) {
      fold_change[down_index] <- 1 / (fold_change_min +
        stats::rexp(length(down_index), fold_change_excess_rate))
      direction[down_index] <- "down"
    }
  }

  list(
    s = if (group_sizes[1] == group_sizes[2]) group_sizes[1] else NA_integer_,
    group_sizes = group_sizes, n_var = n_var, n_diffexp = n_diffexp,
    n_samples = n_samples, condition = condition,
    source_index = source_index, baseline_mean = baseline_mean,
    dispersion = dispersion_used, dispersion_unscaled = dispersion_unscaled,
    dispersion_mom = dispersion_mom, dispersion_edger = dispersion_edger,
    library_size_factors = library_size_factors,
    batch = batch, batch_multiplier = batch_multiplier,
    fold_change = fold_change, direction = direction
  )
}
