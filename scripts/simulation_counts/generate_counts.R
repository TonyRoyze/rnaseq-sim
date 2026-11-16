# Calibrate both CMP moments; never substitute lambda for the mean.
r_underdispersed_cmp <- function(mu, k) {
  if (!requireNamespace("COMPoissonReg", quietly = TRUE))
    stop("Install COMPoissonReg to generate underdispersed counts.")
  vapply(mu, function(m) {
    if (m == 0) return(0)
    minimum_variance <- (m - floor(m)) * (ceiling(m) - m)
    if (k * m <= minimum_variance)
      stop("Requested underdispersion is infeasible for this noninteger mean.")
    moments <- function(theta) {
      lambda <- exp(theta[1]); nu <- 1 + exp(theta[2])
      c(COMPoissonReg::ecmp(lambda, nu), COMPoissonReg::vcmp(lambda, nu))
    }
    objective <- function(theta) {
      z <- tryCatch(moments(theta), error = function(e) c(NA_real_, NA_real_))
      if (any(!is.finite(z)) || any(z <= 0)) return(1e20)
      sum(log(z / c(m, k * m))^2)
    }
    nu0 <- 1 / k
    fit <- stats::optim(c(nu0 * log(m), log(nu0 - 1)), objective,
                        method = "Nelder-Mead", control = list(maxit = 500, reltol = 1e-12))
    if (objective(fit$par) > 1e-6)
      stop("COM-Poisson moment calibration failed; target was not approximated silently.")
    COMPoissonReg::rcmp(1, lambda = exp(fit$par[1]), nu = 1 + exp(fit$par[2]))
  }, numeric(1))
}

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
    counts[, sample_index] <- if (count_model == "generalized_poisson") {
      r_underdispersed_gp(expected[, sample_index], design$variance_mean_ratio)
    } else if (count_model == "com_poisson") {
      r_underdispersed_cmp(expected[, sample_index], design$variance_mean_ratio)
    } else if (count_model == "poisson") {
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
