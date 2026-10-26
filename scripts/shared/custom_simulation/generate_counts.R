# Draw the count matrix, then optionally add outlier counts and extra zeros.
generate_simulated_count_matrix <- function(
    design, count_model, outlier_fraction, outlier_multiplier, sparsity) {
  counts <- matrix(0L, nrow = design$n_var, ncol = design$n_samples)
  expected <- outer(design$baseline_mean, rep(1, design$n_samples))
  group2 <- design$condition == "group2"
  expected[, group2] <- expected[, group2, drop = FALSE] * design$fold_change
  expected <- sweep(expected, 2, design$library_size_factors, `*`)
  expected <- sweep(expected, 2, design$batch_multiplier, `*`)

  for (sample_index in seq_len(design$n_samples)) {
    counts[, sample_index] <- if (count_model == "poisson") {
      stats::rpois(design$n_var, lambda = expected[, sample_index])
    } else {
      stats::rnbinom(design$n_var, size = 1 / design$dispersion,
                     mu = expected[, sample_index])
    }
  }

  n_outlier <- floor(length(counts) * outlier_fraction)
  if (n_outlier > 0L) {
    index <- sample.int(length(counts), n_outlier)
    counts[index] <- round(counts[index] *
      stats::runif(n_outlier, outlier_multiplier[1], outlier_multiplier[2]))
  }
  n_zero <- floor(length(counts) * sparsity)
  if (n_zero > 0L) counts[sample.int(length(counts), n_zero)] <- 0L
  counts
}
