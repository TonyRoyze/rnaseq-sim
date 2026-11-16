# Run from the project root: Rscript --vanilla tests/full_factorial.R
simulation_version <- '4.0'
for (e in parse('scripts/simulation_counts/run_batch.R')) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) &&
      as.character(e[[2]])[1] %in% c('parse_args')) eval(e)
}
source('scripts/simulation_design/factorial_design.R')
source('scripts/simulation_design/de_proportion_modes.R')
opts <- parse_args(character())
stopifnot(opts$de_mode == 'fixed', opts$include_imbalance,
          identical(opts$de_proportions, c(.05,.25,.75)))
stopifnot(!parse_args('--no-imbalance')$include_imbalance)
sc <- build_scenarios(10000)
balanced <- build_scenarios(10000, FALSE)
stopifnot(nrow(sc) == 1512, nrow(balanced) == 840,
          !anyDuplicated(sc$scenario_id), all(sc$n_samples == 2 * sc$sample_size_baseline),
          all(sc$n_group1 >= 3), all(sc$n_group2 >= 3))
splits <- unique(sc[c('sample_size_baseline','n_group1','n_group2')])
stopifnot(nrow(splits) == 9,
          identical(as.integer(table(splits$sample_size_baseline)), c(1L,2L,2L,2L,2L)))
factors <- c('fold_change_min','upregulated_fraction','library_size_sdlog',
             'outlier_fraction','sparsity','batch_effect','dispersion_scale',
             'n_group1','n_group2')
for (model in c('poisson', 'negative_binomial')) {
  family <- sc[sc$count_model == model, ]
  stopifnot(!anyDuplicated(family[factors]),
            nrow(family) == if (model == 'poisson') 216 else 648)
}
stopifnot(all(sc$fold_change_min == 2), all(sc$library_size_sdlog == 0),
          setequal(sc$upregulated_fraction, c(.25,.5)),
          setequal(sc$sparsity, c(0,.1)), all(sc$fixed_fold_change))
sc$replicate <- 1L
x <- assign_de_proportions(sc, 'fixed', opts$de_proportions, 42, 10000)
stopifnot(nrow(x) == 4536, !anyDuplicated(x$scenario_id),
          all(table(x$de_fraction) == nrow(sc)),
          all(x$n_diffexp == round(10000 * x$de_fraction)))
stopifnot(sum(x$count_model == 'poisson') == 648,
          sum(x$count_model == 'negative_binomial') == 1944)
y <- assign_de_proportions(sc, 'setting', opts$de_proportions, 42, 10000)
stopifnot(identical(y, assign_de_proportions(sc, 'setting', opts$de_proportions, 42, 10000)),
          all(vapply(split(y$de_fraction,y$de_setting_id),function(z) length(unique(z))==1,logical(1))))
z <- assign_de_proportions(sc, 'combination', opts$de_proportions, 42, 10000)
stopifnot(nrow(z) == nrow(sc), length(unique(z$de_fraction)) == 3)
for (module in c('helpers.R', 'prepare_design.R', 'generate_counts.R')) {
  source(file.path('scripts/simulation_counts', module))
}
p <- list(mean_total=rep(100,100),dispersion_mom_total=rep(.1,100),dispersion_edger_total=rep(.1,100))
for (i in seq_len(nrow(splits))) {
  sizes <- as.integer(splits[i,c('n_group1','n_group2')])
  for (model in c('poisson','negative_binomial')) {
    d <- prepare_simulation_design(p,10,100,30,model,1,'total','method_of_moments',NULL,0,1,1,.5,NULL,.5,1,group_sizes=sizes)
    counts <- generate_simulated_count_matrix(d,model,0,c(5,10),0)
    stopifnot(ncol(counts)==sum(sizes),all(as.integer(table(d$condition))==sizes))
  }
}
d <- prepare_simulation_design(p,10,100,30,'negative_binomial',1,'total',
  'method_of_moments',NULL,0,2,1,.5,NULL,.5,1,fixed_fold_change=TRUE)
stopifnot(all(d$fold_change[d$direction=='up']==2),
          all(d$fold_change[d$direction=='down']==.5), all(d$library_size_factors==1))
cat('PASS: factorial coverage, counts, imbalance sizes, DE modes, and 18 count-generation cases\n')

filenames <- sprintf('%s-DGElist_GSE150910-disp-edgeR-rep100-n10000-seed2147483647-ref12345678.rds', x$scenario_id)
stopifnot(max(nchar(filenames, type = 'bytes')) <= 255L)
cat('PASS: filenames fit filesystem limits; 100-replicate total =', nrow(x) * 100, '\n')
