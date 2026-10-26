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

  list(
    gene_id = rownames(counts) %||% paste0("gene", seq_len(nrow(counts))),
    mean_group1 = fit1$mean,
    mean_group2 = fit2$mean,
    mean_total = rowMeans(counts),
    dispersion_group1 = fit1$dispersion,
    dispersion_group2 = fit2$dispersion,
    dispersion_total = fit_all$dispersion
  )
}
