# Normalized finite-support Consul-Jain GP, allowing any negative shape.
# Normalization changes the usual nominal moment formulas at low means;
# calibrate the actual moments instead of claiming those formulas are exact.
gp_probabilities <- function(theta, shape, max_support = 1000000L) {
  if (!is.finite(theta) || !is.finite(shape) || theta <= 0 || shape >= 0)
    stop("Generalized-Poisson parameters must be finite with theta > 0 and shape < 0.")
  upper <- ceiling(theta / -shape) - 1
  if (!is.finite(upper) || upper > max_support)
    stop("Generalized-Poisson support exceeds the configured numerical limit.")
  y <- seq.int(0, upper)
  z <- theta + shape * y
  logp <- log(theta) + (y - 1) * log(z) - z - lgamma(y + 1)
  p <- exp(logp - max(logp)); p <- p / sum(p)
  mean <- sum(y * p)
  variance <- sum((y - mean)^2 * p)
  list(y = y, p = p, mean = mean, variance = variance)
}

calibrate_underdispersed_gp <- function(mu, k, tolerance = 1e-3) {
  validate_simulation_scalar(mu, "mu", 0, Inf, FALSE)
  validate_simulation_scalar(k, "k", 0, 1, FALSE)
  if (k >= 1) stop("Underdispersed Generalized Poisson requires k < 1.")
  fraction <- mu - floor(mu)
  if (k * mu <= fraction * (1 - fraction))
    stop("Requested variance is infeasible for integer counts at this mean.")
  # Nominal GP moments: mean=theta/(1-shape), variance=theta/(1-shape)^3.
  initial <- c(log(mu / sqrt(k)), log(1 / sqrt(k) - 1))
  evaluate <- function(par) gp_probabilities(exp(par[1]), -exp(par[2]))
  objective <- function(par) {
    d <- tryCatch(evaluate(par), error = function(e) NULL)
    if (is.null(d) || d$mean <= 0 || d$variance <= 0) return(1e20)
    sum(log(c(d$mean / mu, d$variance / (mu * k)))^2)
  }
  par <- initial
  d <- evaluate(par)
  error <- function(d) max(abs(c(d$mean / mu - 1, d$variance / (mu * k) - 1)))
  if (error(d) > tolerance) {
    fit <- stats::optim(par, objective, method = "Nelder-Mead",
                        control = list(maxit = 1000, reltol = 1e-12))
    par <- fit$par; d <- evaluate(par)
  }
  if (error(d) > tolerance)
    stop("Generalized-Poisson moment calibration failed for requested mu and k.")
  d$theta <- exp(par[1]); d$shape <- -exp(par[2]); d
}

r_underdispersed_gp <- function(mu, k) {
  validate_simulation_scalar(k, "k", 0, 1, FALSE)
  if (k >= 1) stop("Underdispersed Generalized Poisson requires k < 1.")
  if (any(!is.finite(mu)) || any(mu < 0)) stop("Means must be finite and nonnegative.")
  # Reuse probability vectors for repeated means within this column only.
  result <- numeric(length(mu))
  for (at in split(seq_along(mu), match(mu, unique(mu)))) {
    m <- mu[at[1]]
    if (m == 0) next
    d <- calibrate_underdispersed_gp(m, k)
    result[at] <- sample.int(length(d$p), length(at), replace = TRUE, prob = d$p) - 1L
  }
  result
}
