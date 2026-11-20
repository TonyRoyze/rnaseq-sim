# Archive assignment results or transfer a full assignment

Use `migrate_assignment.R` for output produced by `run_machine.R`, including
the 100 numbered configs. **Results-only is the default**: analysis result RDS
files, status/manifests, logs, completion markers and saved machine config are
copied. Generated count datasets and reference-parameter caches are excluded.
Source files are preserved; the tool never deletes generated datasets.
No simulation or analysis runs during migration.

## Collect results, status and logs (default)

Stop the assignment's workers, then preview:

```powershell
Rscript scripts/result_transfer/migrate_assignment.R --config config/assignments_100/config-001.dcf --destination "E:/assignment-001"
```

Add `--apply` to copy:

```powershell
Rscript scripts/result_transfer/migrate_assignment.R --config config/assignments_100/config-001.dcf --destination "E:/assignment-001" --apply
```

No `--results-only` flag is required (it is accepted as an explicit alias).
The source defaults to the config Output. Change `001` for each assignment.
Generation `status.csv` and `progress.log` may retain their paths inside a
`datasets/` folder in the archive; count-matrix RDS files in that folder are
excluded. Analysis results are copied intact, including any gene-level scores.
This preserves partial progress too; it does not mark incomplete work Done.

Keep one archive folder per assignment. A `.results-only-archive.rds` marker
identifies the archive as excluding generated datasets. Use these folders for
central results collection; they are not full resume transfers. Code, package
libraries and the four reference RDS files are not copied. Reference files
must be available on the exporting computer for validation.
If a previous destination already contains generated datasets or parameter
caches, choose a fresh archive folder; the tool will not delete those files.

## Full transfer for resuming on another computer (optional)

Use `--include-datasets` on **both** copy steps below to transfer generated
datasets along with the results and metadata. This is a separate use case
from collecting finished results. A results-only archive cannot be promoted
to a full transfer by adding the flag; return to the original output instead.

## 1. Stop the assignment on its original computer

Wait for completion, or stop its launcher and all child R workers. The tool
refuses an existing `.machine-running` lock at either end. A lock left by a
crash must be removed manually only after confirming all its workers stopped.
Do not run both source and destination copies of an assignment.

## 2. Copy to portable storage

From the repository root, preview first (PowerShell commands are single lines):

```powershell
Rscript scripts/result_transfer/migrate_assignment.R --config config/assignments_100/config-001.dcf --destination "E:/assignment-001" --include-datasets
```

The source defaults to the config's `Output` directory. Add `--apply` to copy:

```powershell
Rscript scripts/result_transfer/migrate_assignment.R --config config/assignments_100/config-001.dcf --destination "E:/assignment-001" --include-datasets --apply
```

The destination is the exact assignment folder, not its parent. Existing
identical files are reused. Conflicting or extra destination files stop the
transfer; files are not overwritten. An interrupted copy can be rerun.

## 3. Copy into the receiving computer's repository

Use the same code version, original config, pinned packages and identical
four reference RDS files. From that repository root, preview:

```powershell
Rscript scripts/result_transfer/migrate_assignment.R --config config/assignments_100/config-001.dcf --source "E:/assignment-001" --destination "results/assignments_100_v5_gamma_count/config-001-NBPSeq" --include-datasets
```

Repeat that command with `--apply`. Then resume:

```powershell
.\automations\run_machine.ps1 --config config/assignments_100/config-001.dcf
```

Replace config number, source and destination for the actual assignment.
Restore the destination at its config's `Output` path before resuming. A
new, unstarted assignment needs no migration: run its config directly.
The same R migration commands work on macOS/Linux with appropriate drive paths.

## Validation and limits

- Preview validates without creating destination directories or changing files.
- The saved version/config must exactly match the supplied config and current
  code version. Changing workers, job bounds, methods, seed or scientific
  settings is not a migration of the same assignment.
- In full-transfer mode, dataset filenames, scenario metadata, versions, dimensions, seeds, selected
  references and reference checksums are checked against the current planner.
- Analysis identities, methods, settings and RNG seeds are checked against the
  assigned method profile. The tool does not rerun numerical DE analysis.
- All included files are checksummed before copying; copied files and the
  source are verified again. A receipt of relative paths and MD5 checksums is
  saved under `.migration-receipts/` in the destination.
- Existing launcher locks are respected, and migration temporarily acquires
  those same locks while copying. This blocks new local launchers. It cannot
  stop an already orphaned R process or coordinate copies on separate machines;
  stop all writers before using it.
- Temporary `.tmp` files and prior migration receipts are excluded. They are
  not completed results. Symbolic links and nested source/destination trees
  are refused. Full dataset/result names are preserved, so Windows long-path
  support and sufficient disk space are still necessary.
- File contents are copied without rewriting historical absolute paths in
  logs, manifests or result metadata. The current launcher and analysis runner
  rediscover local files when resuming. Do not use historical paths as transfer
  destinations. Existing batch logs may be overwritten by the normal launcher
  on resume; the preserved source copy retains the pre-transfer logs.
- Results-only archives exclude generated datasets and parameter caches without
  reading or hashing them. Full-transfer mode includes and validates them.
- A successful migration means all included files were verified, not that the
  scientific assignment is complete. Partial work is intentionally supported.
  Use the launcher's final success status to mark an assignment Done.
- This tool does not change the experiment, redistribute jobs between configs,
  merge different assignments, or convert old-version datasets to version 5.2.

## Export completed analysis files with short names

From the project folder, run this PowerShell command to export every numbered
assignment to a drive or a separate folder.

```powershell
Rscript --vanilla scripts/result_transfer/export_results.R --source results/assignments_100_v5_gamma_count --destination "E:/RNA-results"
```

For a single assignment, point `--source` at its `config-NNN-*` folder. For
older output, supply its existing source folder and use a
separate destination. Assignment exports contain finalized analysis RDS files,
including successful results from partial batches. Files use paths such as
`030/000001.rds`. Keep `export_index.csv` with them. It records each original
path, short path, and MD5 hash. Rerunning an export reuses indexed names for
unchanged files. Exporting a single config preserves other indexed configs.

To verify and copy an export while retaining its short paths, run this from
the workspace root:

```powershell
Rscript --vanilla import_results.R --source "E:/RNA-results" --destination "C:/RNA-results" --layout short
```

The import default restores the original indexed paths. The exporter also
supports the older shard/batch pipeline and retains its existing hash-based
naming and metadata. Use `migrate_assignment.R --include-datasets` for a full
assignment resume transfer.

## Export logs separately

From the project folder on Windows, run:

```powershell
Rscript --vanilla scripts/result_transfer/export_logs.R --source results/assignments_100_v5_gamma_count
```

This copies all `.log` files (case-insensitive) from the source tree directly
into `results/assignments_100_v5_gamma_count_logs`. Each filename joins the original path
components with `__`, for example
`config-001-NBPSeq__worker-06__jobs-4341-4360__datasets__shard-001-of-001__progress.log`.
Use `--destination "E:/RNA-logs"` to override the sibling folder. Rerunning
refreshes existing flat copies.

Validation: `Rscript --vanilla tests/assignment_migration.R` exercises isolated
fixtures for preview, idempotent copying, receiving-machine transfer, temporary
files, locks, conflicts and mismatched versions/seeds/methods, plus default
results-only exclusion of counts/caches while retaining logs and status. A real one-worker
smoke test generated two datasets and ten core T.test results, copied them
through a portable folder into a restored output, and resumed without rewriting
any completed result. It used a small synthetic reference on macOS. Real
Windows drive and long-path behavior has not been exercised in this session.
