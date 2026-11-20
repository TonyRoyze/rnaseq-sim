# Counts from a renewal process with Gamma waiting times. The survival
# probability at j is the chance that the jth arrival occurs within one unit.
gamma_count_moments <- function(rate, shape, arrivals) {
  tail <- stats::pgamma(rate, shape = arrivals * shape)
  offset <- arrivals[1] - 1
  if (any(!is.finite(tail)) || tail[length(tail)] > 1e-12 ||
      (offset > 0 && tail[1] < 1 - 1e-12))
    return(c(NA_real_, NA_real_))
  centered_mean <- sum(tail)
  variance <- sum((2 * (arrivals - offset) - 1) * tail) - centered_mean^2
  c(offset + centered_mean, variance)
}

gamma_count_parameters <- function(mu, k, tolerance = 1e-3) {
  fraction <- mu - floor(mu)
  if (mu * k <= fraction * (1 - fraction))
    stop("Requested gamma-count variance is infeasible for this mean.")

  # For rate / shape = L and inverse shape = q, the large-mean renewal
  # moments are L - (1-q)/2 and L*q + (1-q^2)/12. Solve for q and L.
  target_variance <- k * mu
  q <- k
  if (target_variance > 1 / 12) {
    linear <- mu + 0.5
    constant <- target_variance - 1 / 12
    q <- 2 * constant / (linear + sqrt(linear^2 - 7 * constant / 3))
  }
  shape <- 1 / q
  rate <- (mu + (1 - q) / 2) * shape
  if (mu >= 1000 && target_variance >= 200) return(c(rate = rate, shape = shape))

  target <- c(mu, target_variance)
  radius <- 14 * sqrt(k * mu) + 30
  arrivals <- seq.int(max(1, floor(mu - radius)), ceiling(mu + radius))
  error <- function(par) {
    parameters <- exp(par)
    moments <- gamma_count_moments(parameters[1], 1 + parameters[2], arrivals)
    if (any(!is.finite(moments))) return(c(Inf, Inf))
    moments / target - 1
  }
  initial <- log(c(rate, shape - 1))
  if (max(abs(error(initial))) <= tolerance) return(c(rate = rate, shape = shape))

  fit <- stats::optim(initial, function(par) sum(error(par)^2),
                      method = "Nelder-Mead",
                      control = list(maxit = 500L, reltol = 1e-12))
  if (max(abs(error(fit$par))) > tolerance)
    stop(sprintf("Gamma-count moment calibration failed (mean=%.10g, variance/mean=%.10g).", mu, k))
  c(rate = exp(fit$par[1]), shape = 1 + exp(fit$par[2]))
}

r_underdispersed_gamma_count <- function(mu, k) {
  validate_simulation_scalar(k, "k", 0, 1, FALSE)
  if (k >= 1) stop("Underdispersed gamma-count requires k < 1.")
  if (any(!is.finite(mu)) || any(mu < 0)) stop("Means must be finite and nonnegative.")
  result <- numeric(length(mu))
  positive <- mu > 0
  if (!any(positive)) return(result)

  mu <- mu[positive]
  unique_mu <- unique(mu)
  parameters <- vapply(unique_mu, gamma_count_parameters, numeric(2), k = k)
  index <- match(mu, unique_mu)
  rate <- parameters[1, index]
  shape <- parameters[2, index]
  uniform <- stats::runif(length(mu))

  # Start near the target mean, then expand if a draw is in the upper tail.
  # Binary search inverts the CDF without allocating its entire support.
  upper <- ceiling(mu + 12 * sqrt(k * mu) + 16)
  outside <- stats::pgamma(rate, shape = (upper + 1) * shape, lower.tail = FALSE) < uniform
  while (any(outside)) {
    upper[outside] <- 2 * upper[outside]
    outside[outside] <- stats::pgamma(rate[outside], shape = (upper[outside] + 1) * shape[outside],
                                     lower.tail = FALSE) < uniform[outside]
  }
  lower <- numeric(length(mu))
  active <- which(lower < upper)
  while (length(active)) {
    mid <- lower[active] + floor((upper[active] - lower[active]) / 2)
    reached <- stats::pgamma(rate[active], shape = (mid + 1) * shape[active],
                             lower.tail = FALSE) >= uniform[active]
    upper[active[reached]] <- mid[reached]
    lower[active[!reached]] <- mid[!reached] + 1
    active <- active[lower[active] < upper[active]]
  }
  result[positive] <- lower
  result
}
