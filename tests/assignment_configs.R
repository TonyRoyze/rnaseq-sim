# Pure validation of the 100 portable, single-method work assignments.
source('scripts/simulation_pipeline/machine_config.R')
source('scripts/simulation_design/factorial_design.R')
source('scripts/shared/version.R')
for (expr in parse('scripts/simulation_counts/run_batch.R')) {
  if (is.call(expr) && identical(expr[[1]], as.name('<-')) &&
      identical(expr[[2]], as.name('parse_args'))) eval(expr)
}
paths <- sprintf('config/assignments_100/config-%03d.dcf', 1:100)
stopifnot(all(file.exists(paths)),
          length(list.files('config/assignments_100', pattern='[.]dcf$')) == 100L)
configs <- lapply(paths, read_machine_config)
expected <- c(NBPSeq=27L, DESeq2=27L, DESeq=10L, edgeR.lrt=4L,
              edgeR.qlf=4L, DSS=4L, voom=4L, Wilcoxon=4L,
              T.test=4L, NOISeq=4L, ABSSeq=4L, ROTS=4L)
methods <- vapply(configs, function(cfg) {
  stopifnot(length(cfg$Methods) == 1L)
  cfg$Methods
}, character(1))
stopifnot(identical(as.integer(table(methods)[names(expected)]), unname(expected)))
scientific_fields <- setdiff(names(configs[[1]]),
                            c('Job-Start', 'Job-End', 'Output', 'Methods'))
outputs <- vapply(configs, `[[`, character(1), 'Output')
stopifnot(!anyDuplicated(outputs))
stopifnot(all(startsWith(outputs, paste0('results/assignments_100_v', simulation_version, '/'))))
opts <- parse_args(machine_simulation_args(configs[[1]]))
total <- nrow(build_scenarios(opts$n_var, opts$include_imbalance,
                             opts$underdispersion_model)) * opts$replicates *
  length(opts$de_proportions)
stopifnot(total == 136080, opts$reference_mode == 'random', opts$replicates == 30L)
for (method in names(expected)) {
  pool <- configs[methods == method]
  starts <- vapply(pool, `[[`, numeric(1), 'Job-Start')
  ends <- vapply(pool, `[[`, numeric(1), 'Job-End')
  stopifnot(starts[1] == 1, tail(ends, 1) == total,
            all(tail(starts, -1) == head(ends, -1) + 1),
            all(ends - starts + 1 == total / expected[[method]]))
}
for (cfg in configs) {
  stopifnot(identical(cfg[scientific_fields], configs[[1]][scientific_fields]))
  parsed <- parse_args(c(machine_simulation_args(cfg),
    '--job-start', as.character(cfg[['Job-Start']]),
    '--job-end', as.character(cfg[['Job-End']])))
  stopifnot(parsed$job_end <= total)
  workers <- machine_worker_ranges(cfg)
  stopifnot(nrow(workers) == 6L,
            workers$start[1] == cfg[['Job-Start']],
            tail(workers$end, 1) == cfg[['Job-End']],
            all(tail(workers$start, -1) == head(workers$end, -1) + 1))
}
# Check the suggested first wave has 30 distinct assignments in the right pools.
initial <- c(1:8, 28:35, 55:57, 65, 69, 73, 77, 81, 85, 89:90, 93, 97:98)
stopifnot(length(initial) == 30L, !anyDuplicated(initial),
          sum(methods[initial] == 'NBPSeq') == 8L,
          sum(methods[initial] == 'DESeq2') == 8L,
          sum(methods[initial] == 'DESeq') == 3L)
cat('PASS: 100 configs; one method each; all 12 methods cover every dataset once;',
    'unique outputs; identical science; six-worker coverage; 30-machine first wave\n')
