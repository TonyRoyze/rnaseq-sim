# Bounded-memory factorial job planning. Logical job IDs retain the original
# DE -> scenario -> replicate -> reference order, independently of sharding.
preserve_rng <- function(fun) {
  had <- exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE)
  if (had) previous <- get('.Random.seed', envir = .GlobalEnv)
  on.exit({
    if (had) assign('.Random.seed', previous, envir = .GlobalEnv) else
      if (exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE))
        rm('.Random.seed', envir = .GlobalEnv)
  }, add = TRUE)
  fun()
}

make_job_iterator <- function(scenarios, opts, references) {
  n_base <- nrow(scenarios)
  n_rep <- opts$replicates
  n_de <- if (opts$de_mode == 'fixed') length(opts$de_proportions) else 1L
  n_ref <- if (opts$reference_mode == 'all') nrow(references) else 1L
  n_total <- as.double(n_base) * n_rep * n_de * n_ref
  # Sorting only the small base table permits balanced cyclic assignment
  # within each sample-count block without sorting millions of job rows.
  weights <- sort(unique(scenarios$n_samples), decreasing = TRUE)
  base_groups <- lapply(weights, function(w) which(scenarios$n_samples == w))
  group_lengths <- lengths(base_groups)
  ends <- cumsum(as.double(group_lengths) * n_rep * n_de * n_ref)
  starts <- c(0, head(ends, -1))
  n_assigned <- if (opts$shard_index > n_total) 0 else
    floor((n_total - opts$shard_index) / opts$shard_count) + 1
  cursor <- 0

  # Store one RNG checkpoint per base scenario/setting, never per replicate.
  # This preserves the previous random DE draws without a full expanded table.
  checkpoints <- NULL
  de_groups <- NULL
  if (opts$de_mode != 'fixed') {
    keys <- if (opts$de_mode == 'setting') scenarios$de_setting_id else scenarios$scenario_id
    de_groups <- match(keys, unique(keys))
    checkpoints <- preserve_rng(function() {
      set.seed(as.integer((as.double(opts$base_seed) + 7919) %% .Machine$integer.max))
      lapply(seq_len(max(de_groups)), function(g) {
        checkpoint <- get('.Random.seed', envir = .GlobalEnv)
        remaining <- as.double(n_rep)
        while (remaining > 0) {
          take <- min(remaining, opts$batch_size)
          invisible(sample.int(length(opts$de_proportions), take, replace = TRUE))
          remaining <- remaining - take
        }
        checkpoint
      })
    })
  }
  random_de <- function(base, replicate) preserve_rng(function() {
    result <- numeric(length(base))
    groups <- de_groups[base]
    for (g in unique(groups)) {
      positions <- which(groups == g)
      wanted <- replicate[positions]
      assign('.Random.seed', checkpoints[[g]], envir = .GlobalEnv)
      offset <- 0
      while (offset < max(wanted)) {
        take <- min(opts$batch_size, max(wanted) - offset)
        draws <- sample.int(length(opts$de_proportions), take, replace = TRUE)
        keep <- which(wanted > offset & wanted <= offset + take)
        result[positions[keep]] <- opts$de_proportions[draws[wanted[keep] - offset]]
        offset <- offset + take
      }
    }
    result
  })

  next_batch <- function(limit = opts$batch_size) {
    count <- min(limit, opts$batch_size, n_assigned - cursor)
    if (count <= 0) return(NULL)
    ranks <- (opts$shard_index - 1) + (cursor + seq_len(count) - 1) * opts$shard_count
    blocks <- findInterval(ranks, starts)
    local <- ranks - starts[blocks]
    block_width <- as.double(group_lengths[blocks]) * n_rep * n_ref
    de_index <- floor(local / block_width)
    local <- local %% block_width
    local_base <- floor(local / (as.double(n_rep) * n_ref)) + 1
    base <- integer(count)
    for (b in unique(blocks)) {
      at <- which(blocks == b)
      base[at] <- base_groups[[b]][local_base[at]]
    }
    replicate <- as.integer(floor((local %% (as.double(n_rep) * n_ref)) / n_ref) + 1)
    reference_index <- as.integer(local %% n_ref + 1)
    setting_job_index <- de_index * n_base * n_rep + (as.double(base) - 1) * n_rep + replicate
    job_index <- (setting_job_index - 1) * n_ref + reference_index
    job <- scenarios[base, , drop = FALSE]
    job$replicate <- replicate
    job$de_fraction <- if (opts$de_mode == 'fixed') opts$de_proportions[de_index + 1] else
      random_de(base, replicate)
    job$de_mode <- opts$de_mode
    job$n_diffexp <- as.integer(round(opts$n_var * job$de_fraction))
    job$scenario_id <- paste0(job$scenario_id, '_de_', opts$de_mode, '_',
                             vapply(job$de_fraction, function(v) sprintf('%.15g', v), character(1)))
    job$setting_job_index <- setting_job_index
    if (opts$reference_mode == 'random') {
      reference_index <- preserve_rng(function() vapply(setting_job_index, function(index) {
        set.seed(as.integer((as.double(opts$base_seed) + index * 104729) %% .Machine$integer.max))
        sample.int(nrow(references), 1L)
      }, integer(1)))
    }
    job$reference_index <- reference_index
    job$reference_id <- references$reference_id[reference_index]
    job$reference <- references$reference[reference_index]
    job$reference_hash <- references$reference_hash[reference_index]
    job$job_index <- job_index
    job$estimated_work <- as.double(job$n_samples) * opts$n_var
    job$shard <- opts$shard_index
    job$seed <- as.integer((as.double(opts$base_seed) + job_index * 1009 + replicate) %% .Machine$integer.max)
    job$dispersion_estimator <- opts$dispersion_estimator
    job$planner <- 'bounded_v1'
    rownames(job) <- NULL
    cursor <<- cursor + count
    job
  }
  list(next_batch = next_batch, total = n_total, assigned = n_assigned,
       base_count = n_base)
}

job_output_name <- function(job) {
  sprintf('%s-%s-disp-%s-rep%03d-n%d-seed%d-ref%s.rds', job$scenario_id,
          job$reference_id, ifelse(job$dispersion_estimator == 'method_of_moments', 'mom', 'edgeR'),
          job$replicate, job$n_var, job$seed, substr(job$reference_hash, 1, 8))
}

append_status <- function(record, path) {
  # This file is an event log: latest row for an output_file is its status.
  record$updated_at <- format(Sys.time(), '%Y-%m-%dT%H:%M:%OS3%z')
  has_header <- file.exists(path) && file.info(path)$size > 0
  utils::write.table(record, file = path, append = has_header, sep = ',',
                     row.names = FALSE, col.names = !has_header, na = '', qmethod = 'double')
}
