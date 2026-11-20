source("scripts/simulation_pipeline/machine_config.R")
source("scripts/simulation_pipeline/resume.R")
source("scripts/simulation_design/factorial_design.R")
source("scripts/shared/version.R")
source("scripts/result_transfer/migrate_gamma_count.R")

work <- tempfile("gamma-migration-")
root <- file.path(work, "old", "config-001-NBPSeq")
dir.create(root, recursive = TRUE)
root <- normalizePath(root, winslash = "/")
cfg <- read_machine_config("config/assignments_100/config-001.dcf")
cfg[["Job-Start"]] <- 862
cfg[["Job-End"]] <- 865
cfg[["Jobs-Per-Batch"]] <- 2
cfg$Workers <- 1
cfg$Replicates <- 1
cfg$Output <- root
config <- file.path(work, "config.dcf")
write.dcf(as.data.frame(lapply(cfg, paste, collapse = ","), check.names = FALSE), config)
old <- cfg
old[["Underdispersion-Model"]] <- "generalized_poisson"
state <- file.path(root, "machine_options.rds")
saveRDS(list(version = simulation_version, config = old), state)
backup <- file.path(work, "old_generalized_poisson_backup", basename(root))

make_batch <- function(start, ids) {
  batch <- file.path(root, "worker-01", sprintf("jobs-%d-%d", start, start + 1))
  datasets <- file.path(batch, "datasets", "shard-001-of-001")
  analysis <- file.path(batch, "analysis")
  shard <- file.path(analysis, "shard-001-of-001")
  dir.create(datasets, recursive = TRUE)
  dir.create(shard, recursive = TRUE)
  outputs <- file.path(datasets, paste0(ids, ".rds"))
  results <- file.path(analysis, "results", ids, paste0(ids, "--NBPSeq--default.rds"))
  for (i in seq_along(ids)) {
    dir.create(dirname(results[i]), recursive = TRUE)
    saveRDS(ids[i], outputs[i])
    saveRDS(ids[i], results[i])
  }
  windows <- paste0("C:/original/", basename(outputs))
  status <- file.path(datasets, "status.csv")
  write.csv(data.frame(output_file = windows, status = "complete"), status, row.names = FALSE)
  manifest <- file.path(shard, "manifest.csv")
  write.csv(data.frame(dataset_id = ids, output_file = results, status = "complete"),
            manifest, row.names = FALSE)
  marker <- file.path(batch, "complete.rds")
  saveRDS(list(start = start, end = start + 1, stage = "both"), marker)
  list(marker = marker, outputs = outputs, results = results, status = status, manifest = manifest)
}
pure <- make_batch(862, c("poisson_factorial_a", "poisson_factorial_b"))
mixed <- make_batch(864, c("poisson_factorial_c", "generalized_poisson_factorial_d"))
keep <- c(pure$marker, pure$outputs, pure$results, mixed$outputs[1], mixed$results[1])
hashes <- tools::md5sum(keep)
migrate_gamma_count(config, root)
stopifnot(!dir.exists(backup), file.exists(mixed$marker), identical(readRDS(state)$config, old))
dir.create(file.path(root, ".machine-running"))
stopifnot(inherits(try(migrate_gamma_count(config, root, TRUE), silent = TRUE), "try-error"))
unlink(file.path(root, ".machine-running"), recursive = TRUE)

# Reproduce an interruption after a CSV was archived and before its replacement.
pending <- file.path(root, ".gamma-count-migration.rds")
saveRDS(cfg, pending)
archived_status <- file.path(backup, substring(mixed$status, nchar(root) + 2))
dir.create(dirname(archived_status), recursive = TRUE)
stopifnot(file.rename(mixed$status, archived_status))
migrate_gamma_count(config, root, TRUE)
stopifnot(identical(tools::md5sum(keep), hashes),
          !any(file.exists(c(mixed$marker, mixed$outputs[2], mixed$results[2], pending))),
          identical(read.csv(mixed$manifest)$dataset_id, "poisson_factorial_c"),
          identical(read.csv(mixed$status)$output_file, "C:/original/poisson_factorial_c.rds"),
          identical(readRDS(file.path(backup, "machine_options.rds"))$config, old),
          identical(readRDS(state)$config, cfg),
          identical(read_machine_config(file.path(root, "gamma_count_config.dcf")), cfg))
migrate_gamma_count(config, root, TRUE)
unlink(work, recursive = TRUE)
cat("PASS: preview, locking, restart, mixed-batch reuse, and configuration cutover\n")
