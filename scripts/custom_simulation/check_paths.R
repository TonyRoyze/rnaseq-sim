# Verify actual Windows long-path support before a migration or large run.
check_pipeline_paths <- function(root) {
  if (.Platform$OS.type != 'windows') return(invisible(TRUE))
  if (!dir.exists(root) && !dir.create(root, recursive=TRUE)) stop('Cannot create output root: ', root)
  probe_root <- tempfile('path-check-', tmpdir=normalizePath(root, winslash='/'))
  on.exit(unlink(probe_root, recursive=TRUE, force=TRUE), add=TRUE)
  probe <- probe_root
  while (nchar(probe) < 300L) probe <- file.path(probe, 'long-path-probe')
  ok <- suppressWarnings(dir.create(probe, recursive=TRUE))
  if (ok) ok <- suppressWarnings(file.create(file.path(probe, 'check.tmp')))
  if (!isTRUE(ok)) stop(paste(
    'Windows long-path support is unavailable for this output location.',
    'Enable LongPathsEnabled=1 under HKLM\\SYSTEM\\CurrentControlSet\\Control\\FileSystem',
    'using an Administrator PowerShell, restart Windows, and use R >= 4.3.',
    'No datasets were moved by this invocation.', sep='\n'))
  invisible(TRUE)
}
