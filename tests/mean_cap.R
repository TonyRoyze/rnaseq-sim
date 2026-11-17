source('renv/activate.R')
source('scripts/simulation_counts/load.R')

p <- list(gene_id = c('below', 'boundary', 'above', 'extreme'),
          mean_total = c(50, 180000, 180001, 957191.2),
          dispersion_mom_total = rep(.1, 4),
          dispersion_edger_total = rep(.1, 4))
for (model in c('poisson', 'negative_binomial', 'generalized_poisson')) {
  x <- simulate_custom_counts(p, n_var = 4, n_diffexp = 2,
    count_model = model, group_sizes = c(2, 2), variance_mean_ratio = .2,
    fixed_fold_change = TRUE, fold_change_min = 2,
    upregulated_fraction = 1, batch = c('a', 'b', 'a', 'b'),
    batch_effect = 1.5, seed = 42)
  t <- x$truth
  original <- p$mean_total[match(t$source_gene, p$gene_id)]
  stopifnot(identical(t$baseline_mean_original, original),
    identical(t$baseline_mean, pmin(original, 180000)),
    identical(t$baseline_mean_capped, original > 180000),
    x$parameters$baseline_mean_cap == 180000,
    all(t$fold_change_group2_vs_group1[t$is_de] == 2),
    all(is.finite(x$dge$counts)), all(x$dge$counts >= 0))
  path <- tempfile(fileext = '.rds')
  saveRDS(x, path)
  stopifnot(identical(readRDS(path)$truth, t))
  unlink(path)
}
stopifnot(identical(p$mean_total, c(50, 180000, 180001, 957191.2)))
# The cap must not conceal invalid reference means.
p$mean_total[1] <- Inf
stopifnot(inherits(try(simulate_custom_counts(p, n_var = 4, n_diffexp = 0,
  count_model = 'poisson', seed = 42), silent = TRUE), 'try-error'))
cat('PASS: baseline cap, flags, provenance, persistence, and GP with FC2/batch1.5\n')
