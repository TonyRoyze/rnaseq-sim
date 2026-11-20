#!/usr/bin/env Rscript
# Migrate one stopped generalized-Poisson assignment in place.
migrate_gamma_count <- function(config, root, apply = FALSE) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  backup <- file.path(paste0(dirname(root), "_generalized_poisson_backup"), basename(root))
  state <- file.path(root, "machine_options.rds")
  pending <- file.path(root, ".gamma-count-migration.rds")
  lock <- file.path(root, ".machine-running")
  if (apply) {
    if (!dir.create(lock, showWarnings = FALSE)) stop("Assignment is locked: ", root)
    on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  } else if (dir.exists(lock)) stop("Assignment is locked: ", root)

  cfg <- read_machine_config(config)
  if (cfg[["Underdispersion-Model"]] != "gamma_count") stop("Supply a gamma-count config.")
  cfg$Output <- root
  resume_config <- file.path(root, "gamma_count_config.dcf")
  old_state <- file.path(backup, basename(state))
  in_progress <- file.exists(pending)
  if (in_progress) {
    if (!identical(readRDS(pending), cfg)) stop("Resume the migration with the same config and source.")
  } else {
    if (!file.exists(state)) stop("No saved assignment state at ", state)
    if (identical(readRDS(state), list(version = simulation_version, config = cfg)) &&
        file.exists(resume_config)) {
      cat("Already migrated. Resume with ", resume_config, "\n", sep = "")
      return(invisible(NULL))
    }
    if (file.exists(backup) || file.exists(resume_config))
      stop("Existing backup or resume config. Inspect before applying: ", backup)
  }
  saved <- readRDS(if (file.exists(old_state)) old_state else state)
  expected <- cfg
  expected[["Underdispersion-Model"]] <- "generalized_poisson"
  legacy <- check_assignment_resume(saved, expected, simulation_version)

  scenarios <- build_scenarios(cfg[["N-Var"]], cfg[["Include-Imbalance"]] == "true",
                               "generalized_poisson")
  repetitions <- cfg$Replicates * if (cfg[["Reference-Mode"]] == "all") 4 else 1
  is_gp <- function(id) grepl("^generalized_poisson_factorial_", basename(gsub("\\\\", "/", id)))
  ranges <- machine_worker_ranges(cfg)
  batches <- list.dirs(root, recursive = TRUE, full.names = TRUE)
  batches <- batches[grepl("^jobs-[0-9]+-[0-9]+$", basename(batches)) &
                       grepl("^worker-[0-9]+$", basename(dirname(batches))) &
                       dirname(dirname(batches)) == root]
  move <- character()
  reusable <- 0L
  updates <- list()
  reopened <- 0L
  archive_path <- function(path) file.path(backup, substring(path, nchar(root) + 2L))
  for (batch in batches) {
    worker <- as.integer(sub("^worker-", "", basename(dirname(batch))))
    span <- as.numeric(strsplit(sub("^jobs-", "", basename(batch)), "-", fixed = TRUE)[[1]])
    if (worker < 1 || worker > nrow(ranges) ||
        span[1] < ranges$start[worker] || span[2] > ranges$end[worker] ||
        (span[1] - ranges$start[worker]) %% cfg[["Jobs-Per-Batch"]] != 0 ||
        span[2] != min(span[1] + cfg[["Jobs-Per-Batch"]] - 1, ranges$end[worker]))
      stop("Unexpected batch range: ", batch)
    index <- (floor((seq.int(span[1], span[2]) - 1) / repetitions) %% nrow(scenarios)) + 1
    files <- list.files(batch, recursive = TRUE, full.names = TRUE)
    datasets <- files[startsWith(files, paste0(batch, "/datasets/")) & grepl("\\.rds($|\\.)", files)]
    results <- files[startsWith(files, paste0(batch, "/analysis/results/")) & grepl("\\.rds($|\\.)", files)]
    obsolete <- c(datasets[is_gp(datasets)], results[is_gp(results) | is_gp(dirname(results))])
    reusable <- reusable + sum(!results %in% obsolete & endsWith(results, ".rds"))
    affected <- length(obsolete) > 0L || any(scenarios$count_model[index] == "generalized_poisson")
    csv <- files[basename(files) %in% c("status.csv", "manifest.csv")]
    archived <- list.files(archive_path(batch), recursive = TRUE, full.names = TRUE,
                           pattern = "^(status|manifest)\\.csv$")
    csv <- unique(c(csv, file.path(root, substring(archived, nchar(backup) + 2L))))
    for (path in csv) {
      original <- archive_path(path)
      rows <- read.csv(if (file.exists(original)) original else path, check.names = FALSE)
      key <- if (basename(path) == "manifest.csv") "dataset_id" else "output_file"
      if (!key %in% names(rows) || anyNA(rows[[key]])) stop("Invalid status file: ", path)
      old <- is_gp(rows[[key]])
      if (any(old)) {
        updates[[path]] <- rows[!old, , drop = FALSE]
        affected <- TRUE
      }
    }
    marker <- file.path(batch, "complete.rds")
    if (file.exists(marker)) {
      if (!affected) {
        # An invalid marker causes the launcher to stop. Reopen it so absent files can be rebuilt.
        affected <- !isTRUE(tryCatch(assignment_batch_complete(batch, span[1], span[2], cfg$Stage),
                                     error = function(e) FALSE))
      }
      if (affected) {
        obsolete <- c(obsolete, marker)
        reopened <- reopened + 1L
      }
    }
    move <- c(move, obsolete)
  }

  cat(if (apply) "Apply" else "Preview", "gamma-count migration for", root, "\n",
      "Reusable analysis files present:", reusable, "\n",
      "Files to archive:", length(move), "\n", "Batches to reopen:", reopened, "\n",
      "Status files to filter:", length(updates), "\n", "Backup:", backup, "\n",
      "Resume config:", resume_config, "\n")
  if (legacy) cat("Retained files from v5.2 may contain uncapped baseline means.\n")
  if (!apply) return(invisible(NULL))

  publish <- function(path, write) {
    tmp <- tempfile(".gamma-migration-", tmpdir = dirname(path))
    on.exit(unlink(tmp))
    write(tmp)
    if (!file.rename(tmp, path)) stop("Could not publish ", path)
  }
  archive <- function(path) {
    dest <- archive_path(path)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(dest) || !file.rename(path, dest)) stop("Could not archive ", path)
  }
  if (!in_progress) publish(pending, function(path) saveRDS(cfg, path))
  for (path in move) archive(path)
  for (path in names(updates)) {
    if (!file.exists(archive_path(path))) archive(path)
    # An interrupted migration may have already written the filtered file.
    if (!file.exists(path)) publish(path, function(tmp) write.csv(updates[[path]], tmp, row.names = FALSE, na = ""))
  }
  if (!file.exists(resume_config)) {
    fields <- lapply(cfg, function(value) paste(value, collapse = ","))
    publish(resume_config, function(path) write.dcf(as.data.frame(fields, check.names = FALSE), path))
  }
  if (!identical(read_machine_config(resume_config), cfg)) stop("Resume config verification failed.")
  if (!file.exists(old_state)) archive(state)
  snapshot <- list(version = simulation_version, config = cfg)
  if (!file.exists(state)) publish(state, function(path) saveRDS(snapshot, path))
  if (!identical(readRDS(state), snapshot)) stop("Saved assignment verification failed.")
  if (unlink(pending) != 0L) stop("Could not clear migration marker.")
  cat("Migration complete. Resume with ", resume_config, "\n", sep = "")
}

if (sys.nframe() == 0L) {
  script <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  setwd(normalizePath(file.path(dirname(sub("^--file=", "", script)), "../..")))
  args <- commandArgs(TRUE)
  if (!length(args) || "--help" %in% args) {
    cat("Usage: Rscript scripts/result_transfer/migrate_gamma_count.R --config PATH --source PATH [--apply]\n")
  } else {
    apply <- "--apply" %in% args
    args <- args[args != "--apply"]
    if (length(args) != 4L || !setequal(args[c(1L, 3L)], c("--config", "--source")))
      stop("Supply --config and --source. Use --help for usage.")
    source("scripts/simulation_pipeline/machine_config.R")
    source("scripts/simulation_pipeline/resume.R")
    source("scripts/simulation_design/factorial_design.R")
    source("scripts/shared/version.R")
    migrate_gamma_count(args[match("--config", args) + 1L], args[match("--source", args) + 1L], apply)
  }
}
