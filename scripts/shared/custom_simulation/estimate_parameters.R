# Estimate baseline abundance and gene-wise dispersion from real count data.
# Input should be an edgeR DGEList with exactly two groups.
estimate_simulation_parameters <- function(dge, min_mean = 7) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required to estimate simulation parameters.")
  }
  if (is.null(dge$counts) || is.null(dge$samples$group)) {
    stop("dge must contain counts and samples$group, as in an edgeR DGEList.")
  }

  counts <- as.matrix(dge$counts)
  keep <- rowMeans(counts) > min_mean
  counts <- counts[keep, , drop = FALSE]
  groups <- droplevels(factor(dge$samples$group))
  if (length(groups) != ncol(counts)) stop("Sample groups do not match count columns.")
  if (nlevels(groups) != 2L) stop("The input data must have exactly two groups.")
  if (nrow(counts) < 1L) stop("No genes remain after the mean-count filter.")

  fit_dispersion <- function(group_counts, group) {
    group <- droplevels(factor(group))
    y <- edgeR::DGEList(counts = group_counts, group = group)
    y <- edgeR::calcNormFactors(y)
    design <- if (nlevels(group) == 1L) {
      matrix(1, nrow = ncol(group_counts), ncol = 1L)
    } else {
      stats::model.matrix(~ group)
    }
    y <- edgeR::estimateDisp(y, design = design)
    list(mean = rowMeans(y$counts), dispersion = y$tagwise.dispersion)
  }

  group1 <- groups == levels(groups)[1]
  group2 <- groups == levels(groups)[2]
  fit1 <- fit_dispersion(counts[, group1, drop = FALSE], groups[group1])
  fit2 <- fit_dispersion(counts[, group2, drop = FALSE], groups[group2])
  fit_all <- fit_dispersion(counts, groups)

  # Base-R median-ratio normalization, followed by a within-group
  # method-of-moments estimate for Var(Y) = mu + alpha * mu^2.
  positive <- apply(counts > 0, 1, all)
  if (any(positive)) {
    geometric_mean <- exp(rowMeans(log(counts[positive, , drop = FALSE])))
    ratios <- sweep(counts[positive, , drop = FALSE], 1, geometric_mean, "/")
    size_factors <- apply(ratios, 2, stats::median)
  } else {
    size_factors <- colSums(counts)
  }
  if (any(!is.finite(size_factors)) || any(size_factors <= 0)) {
    size_factors <- colSums(counts)
  }
  if (any(!is.finite(size_factors)) || any(size_factors <= 0)) {
    stop("Cannot estimate method-of-moments dispersion: invalid sample size factors.")
  }
  size_factors <- size_factors / exp(mean(log(size_factors)))

  estimate_moments <- function(sample_index, group_label) {
    n_genes <- nrow(counts)
    numerator <- numeric(n_genes)
    denominator <- numeric(n_genes)
    for (level in unique(group_label)) {
      index <- sample_index[group_label == level]
      n <- length(index)
      if (n < 2L) next
      normalized <- sweep(counts[, index, drop = FALSE], 2, size_factors[index], "/")
      gene_mean <- rowMeans(normalized)
      residual_ss <- rowSums((normalized - gene_mean)^2)
      # Normalized Poisson sampling variance is mu / size_factor.
      numerator <- numerator + residual_ss -
        (n - 1) * gene_mean * mean(1 / size_factors[index])
      denominator <- denominator + (n - 1) * gene_mean^2
    }
    if (!any(denominator > 0)) {
      return(rep(NA_real_, n_genes))
    }
    estimate <- pmax(numerator / denominator, 0)
    estimate[denominator <= 0 | !is.finite(estimate)] <- NA_real_
    estimate
  }
  dispersion_mom_group1 <- estimate_moments(which(group1), groups[group1])
  dispersion_mom_group2 <- estimate_moments(which(group2), groups[group2])
  dispersion_mom_total <- estimate_moments(seq_len(ncol(counts)), groups)
  # If a reference group has too few replicates for its own estimate, retain
  # the pooled within-group estimate for that gene.
  dispersion_mom_group1[!is.finite(dispersion_mom_group1)] <-
    dispersion_mom_total[!is.finite(dispersion_mom_group1)]
  dispersion_mom_group2[!is.finite(dispersion_mom_group2)] <-
    dispersion_mom_total[!is.finite(dispersion_mom_group2)]

  list(
    gene_id = rownames(counts) %||% paste0("gene", seq_len(nrow(counts))),
    mean_group1 = fit1$mean,
    mean_group2 = fit2$mean,
    mean_total = rowMeans(counts),
    dispersion_mom_group1 = dispersion_mom_group1,
    dispersion_mom_group2 = dispersion_mom_group2,
    dispersion_mom_total = dispersion_mom_total,
    dispersion_edger_group1 = fit1$dispersion,
    dispersion_edger_group2 = fit2$dispersion,
    dispersion_edger_total = fit_all$dispersion,
    # Retain the old field as an edgeR alias for existing callers.
    dispersion_group1 = fit1$dispersion,
    dispersion_group2 = fit2$dispersion,
    dispersion_total = fit_all$dispersion,
    size_factors_mom = size_factors
  )
}
