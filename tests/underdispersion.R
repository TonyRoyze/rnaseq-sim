source('renv/activate.R')
source('scripts/shared/simulate_custom.R')
source('scripts/custom_simulation/factorial_design.R')
simulation_version <- '4.0'
for (e in parse('scripts/custom_simulation/run_batch.R')) {
 if (is.call(e) && identical(e[[1]], as.name('<-')) && identical(e[[2]], as.name('parse_args'))) eval(e)
}
stopifnot(parse_args(character())$underdispersion_model == 'generalized_poisson')
stopifnot(inherits(try(parse_args(c('--underdispersion-model','invalid')),silent=TRUE),'try-error'))
for(model in c('generalized_poisson','com_poisson')) {
 sc <- build_scenarios(100,underdispersion_model=model)
 stopifnot(nrow(sc)==40824, sum(sc$count_model==model)==17496,
           !anyDuplicated(sc$scenario_id))
}
set.seed(42)
for (k in c(.2,.5,.81)) for(m in c(1,2,5,10,50,100,1000)) {
 d <- calibrate_underdispersed_gp(m,k)
 stopifnot(abs(d$mean/m-1)<.001,abs(d$variance/(k*m)-1)<.001,
           abs(sum(d$p)-1)<1e-12,all(d$p>=0))
 if(m %in% c(5,50)) {
  draws <- r_underdispersed_gp(rep(m,50000),k)
  stopifnot(all(draws>=0),all(draws==floor(draws)),abs(mean(draws)/m-1)<.01,
            abs(var(draws)/mean(draws)/k-1)<.04)
 }
}
stopifnot(all(r_underdispersed_gp(rep(0,5),.2)==0))
stopifnot(inherits(try(calibrate_underdispersed_gp(.5,.2),silent=TRUE),'try-error'))
p <- list(mean_total=rep(50,100),dispersion_mom_total=rep(.1,100),dispersion_edger_total=rep(.1,100))
for(model in c('generalized_poisson','com_poisson','poisson','negative_binomial')) {
 x <- simulate_custom_counts(p,n_var=20,n_diffexp=6,count_model=model,
                            group_sizes=c(3,7),variance_mean_ratio=.5,seed=42)
 stopifnot(identical(dim(x$dge$counts),c(20L,10L)),all(x$dge$counts>=0))
}
cat('PASS: model selection, exact GP moments, empirical moments, support, and end-to-end count generation for all four families\n')
