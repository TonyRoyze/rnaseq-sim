# Assign before reference expansion and sharding so workers agree on DE truth.
assign_de_proportions <- function(jobs, mode, proportions, base_seed, n_var) {
  stopifnot(mode %in% c("combination", "setting", "fixed"))
  if (mode == "fixed") {
    jobs <- jobs[rep(seq_len(nrow(jobs)), times = length(proportions)), , drop = FALSE]
    jobs$de_fraction <- rep(proportions, each = nrow(jobs) / length(proportions))
  } else {
    # Repeat draws per replicate. Mode 2 shares draws across count models and
    # sample sizes with the other factors fixed; mode 1 draws per scenario.
    setting_id <- if ("de_setting_id" %in% names(jobs)) jobs$de_setting_id else jobs$varied_factor
    keys <- if (mode == "setting") paste(jobs$replicate, setting_id) else
      paste(jobs$replicate, jobs$scenario_id)
    groups <- unique(keys)
    # Keep this stream separate from count generation and reference selection.
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv) else
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer((as.double(base_seed) + 7919) %% .Machine$integer.max))
    draws <- proportions[sample.int(length(proportions), length(groups), replace = TRUE)]
    jobs$de_fraction <- draws[match(keys, groups)]
  }
  jobs$de_mode <- mode
  jobs$n_diffexp <- as.integer(round(n_var * jobs$de_fraction))
  # Keep files distinct across modes and proportion selections on reruns.
  jobs$scenario_id <- paste0(jobs$scenario_id, "_de_", mode, "_",
                             vapply(jobs$de_fraction, function(x) sprintf("%.15g", x), character(1)))
  rownames(jobs) <- NULL
  jobs
}
