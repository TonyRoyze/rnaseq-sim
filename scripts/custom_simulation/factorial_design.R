build_scenarios <- function(n_var, include_imbalance = TRUE, underdispersion_model = "generalized_poisson") {
  underdispersion_model <- match.arg(underdispersion_model, c("generalized_poisson", "com_poisson"))
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
    grid$variance_mean_ratio <- if (model == "poisson") 1 else NA_real_
    grid$dispersion_regime <- if (model == "poisson") "equi" else "over"
    grid$design <- "full_factorial"
    grid$sample_size_baseline <- grid$s
    grid$s[grid$n_group1 != grid$n_group2] <- NA_integer_
    grid$varied_factor <- "full_factorial"
    grid$batch_fraction <- 0.5
    grid$de_fraction <- 0.3 # Replaced by DE selection before simulation.
    grid$n_diffexp <- as.integer(round(n_var * grid$de_fraction))
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
  under <- do.call(rbind, lapply(c(0.2, 0.5, 0.81), function(k) {
    g <- make_family(underdispersion_model)
    g$variance_mean_ratio <- k
    g$dispersion_regime <- "under"
    g$scenario_id <- paste0(g$scenario_id, "_k", k)
    g$de_setting_id <- paste0(g$de_setting_id, "_k", k)
    g
  }))
  scenarios <- rbind(scenarios, under)
  rownames(scenarios) <- NULL
  scenarios$scenario_index <- seq_len(nrow(scenarios))
  scenarios$simulation_version <- simulation_version
  scenarios$n_var <- n_var
  scenarios
}
