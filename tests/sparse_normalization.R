source("scripts/de_analysis/method_settings.R")
source("scripts/simulation_pipeline/resume.R")

settings <- make_method_settings("core")
sparse <- rbind(c(0, 4, 1), c(2, 0, 3), c(1, 3, 0))
unsupported <- function(counts) vapply(settings, function(setting) {
  !is.null(sparse_normalization_reason(counts, setting))
}, logical(1))

expected <- c("DESeq2/default", "DESeq/default", "edgeR.lrt/norm_method-RLE",
              "edgeR.qlf/norm_method-RLE", "voom/norm_method-RLE",
              "Wilcoxon/norm_method-RLE", "T.test/norm_method-RLE",
              "NBPSeq/default", "ABSSeq/norm_method-geometric")
labels <- vapply(settings, function(setting) {
  paste(setting$method, setting$setting_id, sep = "/")
}, character(1))
stopifnot(setequal(labels[unsupported(sparse)], expected),
          !any(unsupported(rbind(sparse, c(1, 1, 1)))),
          !any(unsupported(sparse[FALSE, , drop = FALSE])))

batch <- tempfile()
dir.create(file.path(batch, "analysis"), recursive = TRUE)
saveRDS(list(start = 1, end = 1, stage = "both"), file.path(batch, "complete.rds"))
manifest <- data.frame(status = "unsupported", output_file = "never-written.rds")
write.csv(manifest, file.path(batch, "analysis", "manifest.csv"), row.names = FALSE)
stopifnot(assignment_batch_complete(batch, 1, 1, "both"),
          !analysis_manifest_complete(data.frame(status = c("unsupported", "failed"))))
unlink(batch, recursive = TRUE)
cat("PASS: sparse setting classification and completion without result artifacts\n")
