# Standalone sampler checks without edgeR or other external packages.
source("scripts/simulation_counts/helpers.R")
source("scripts/simulation_counts/gamma_count.R")
stopifnot(identical(r_underdispersed_gamma_count(c(0, 0), .5), c(0, 0)))
bad <- try(r_underdispersed_gamma_count(.5, .2), silent = TRUE)
stopifnot(inherits(bad, "try-error"))

for (k in c(.2, .5, .81)) {
  for (mu in c(1, 3.5, 7, 100, 999, 1000, 540000, 1000000)) {
    parameters <- gamma_count_parameters(mu, k)
    radius <- 16 * sqrt(k * mu) + 50
    arrivals <- seq.int(max(1, floor(mu - radius)), ceiling(mu + radius))
    moments <- gamma_count_moments(parameters["rate"], parameters["shape"], arrivals)
    stopifnot(max(abs(moments / c(mu, k * mu) - 1)) <= 1e-3)
  }

  mu <- rep(c(7, 180000, 540000, 1000000), each = 100L)
  parameters <- vapply(mu, gamma_count_parameters, numeric(2), k = k)
  set.seed(985)
  u <- runif(length(mu))
  set.seed(985)
  y <- r_underdispersed_gamma_count(mu, k)
  p <- pgamma(parameters[1, ], (y + 1) * parameters[2, ], lower.tail = FALSE)
  previous <- numeric(length(y))
  positive <- y > 0
  previous[positive] <- pgamma(parameters[1, positive], y[positive] * parameters[2, positive],
                              lower.tail = FALSE)
  stopifnot(all(y == floor(y)), all(y >= 0), all(previous < u), all(u <= p))
}
cat("PASS: gamma-count calibration and discrete CDF inversion at small and million-count means\n")
