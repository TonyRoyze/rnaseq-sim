# User-facing simulation function. The other modules handle parameter
# estimation, simulation design, and count generation separately.
simulate_custom_counts <- function(
    params,
    s = 10,
    n_var = 10000,
    n_diffexp = 3000,
    count_model = c("negative_binomial", "poisson"),
    dispersion_scale = 1,
    dispersion_source = c("total", "group"),
    dispersion_estimator = c("method_of_moments", "edgeR"),
    library_size_factors = NULL,
    library_size_sdlog = 0,
    fold_change_min = 2,
    fold_change_excess_rate = 1,
    upregulated_fraction = 0.5,
    outlier_fraction = 0,
    outlier_multiplier = c(5, 10),
    sparsity = 0,
    batch = NULL,
    batch_fraction = 0.5,
    batch_effect = 1,
    seed = NULL,
    group_sizes = NULL) {

  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required to construct the simulated DGEList.")
  }
  count_model <- match.arg(count_model)
  dispersion_source <- match.arg(dispersion_source)
  dispersion_estimator <- match.arg(dispersion_estimator)
  validate_simulation_scalar(outlier_fraction, "outlier_fraction", 0, 1)
  validate_simulation_scalar(sparsity, "sparsity", 0, 1)
  if (length(outlier_multiplier) != 2L || any(!is.finite(outlier_multiplier)) ||
      outlier_multiplier[1] <= 0 || outlier_multiplier[2] < outlier_multiplier[1]) {
    stop("outlier_multiplier must be two positive values in increasing order.")
  }
  if (!is.null(seed)) set.seed(seed)

  design <- prepare_simulation_design(
    params = params, s = s, n_var = n_var, n_diffexp = n_diffexp,
    count_model = count_model, dispersion_scale = dispersion_scale,
    dispersion_source = dispersion_source, dispersion_estimator = dispersion_estimator,
    library_size_factors = library_size_factors,
    library_size_sdlog = library_size_sdlog,
    fold_change_min = fold_change_min,
    fold_change_excess_rate = fold_change_excess_rate,
    upregulated_fraction = upregulated_fraction,
    batch = batch, batch_fraction = batch_fraction, batch_effect = batch_effect,
    group_sizes = group_sizes
  )
  counts <- generate_simulated_count_matrix(
    design, count_model, outlier_fraction, outlier_multiplier, sparsity
  )

  rownames(counts) <- paste0("g", seq_len(design$n_var))
  colnames(counts) <- paste0("sample", seq_len(design$n_samples))
  samples <- data.frame(
    group = design$condition,
    batch = design$batch,
    library_factor = design$library_size_factors,
    row.names = colnames(counts)
  )
  dge <- edgeR::DGEList(counts = counts, group = design$condition)
  dge$samples$batch <- design$batch
  dge$samples$library.factor <- design$library_size_factors

  truth <- data.frame(
    gene = rownames(counts),
    source_gene = (params$gene_id %||% seq_along(params$mean_total))[design$source_index],
    is_de = design$direction != "none",
    direction = design$direction,
    fold_change_group2_vs_group1 = design$fold_change,
    row.names = rownames(counts)
  )
  dispersion_estimates <- data.frame(
    gene = rownames(counts),
    source_gene = truth$source_gene,
    method_of_moments = design$dispersion_mom,
    edgeR = design$dispersion_edger,
    selected_estimator = dispersion_estimator,
    selected_unscaled = design$dispersion_unscaled,
    dispersion_used = if (count_model == "negative_binomial") design$dispersion else NA_real_,
    row.names = rownames(counts)
  )
  list(
    dge = dge,
    truth = truth,
    samples = samples,
    parameters = list(
      s = design$s, group_sizes = design$group_sizes, n_samples = design$n_samples,
      n_var = design$n_var, n_diffexp = design$n_diffexp,
      count_model = count_model, dispersion = design$dispersion,
      dispersion_scale = dispersion_scale, dispersion_source = dispersion_source,
      dispersion_estimator = dispersion_estimator,
      dispersion_estimates = dispersion_estimates,
      library_size_factors = design$library_size_factors,
      fold_change = design$fold_change,
      outlier_fraction = outlier_fraction, sparsity = sparsity,
      batch_effect = batch_effect, batch_multiplier = design$batch_multiplier,
      source_index = design$source_index
    )
  )
}
