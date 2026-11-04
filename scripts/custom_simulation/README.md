# Run custom simulation batches

`run_batch.R` restores the R packages pinned in this project's `renv.lock`,
estimates expression and dispersion parameters from a real `edgeR::DGEList`,
generates simulated datasets, and checkpoints each result so interrupted work
can resume.

Run from the project root:

```sh
Rscript scripts/custom_simulation/run_batch.R --plan-only
Rscript scripts/custom_simulation/run_batch.R
```

The default reference mode uses `data/DGElist_BRCA.rds`, and outputs go to
`results/custom_simulation_v5.1/`. The runner supports three reference modes:

- `--reference-mode brca` (default): use only `--reference` (defaults to
  `data/DGElist_BRCA.rds`).
- `--reference-mode all`: create a separate job for every reference and
  scenario combination.
- `--reference-mode random`: choose one reference reproducibly for each
  scenario replicate. The chosen reference and checksum are saved in both the
  result RDS and the append-only `status.csv` log.

The default reference list for `all` and `random` includes BRCA, E-ENAD-34,
GSE150910, and LUSC. `DGElist_SCLC_Cell2024.rds` is reserved as an independent
external test set and the runner rejects it as a simulation reference. Override
the eligible reference list with comma-separated paths:

```sh
Rscript scripts/custom_simulation/run_batch.R --reference-mode all
Rscript scripts/custom_simulation/run_batch.R --reference-mode random --replicates 5
Rscript scripts/custom_simulation/run_batch.R --reference-mode all \
  --references data/DGElist_BRCA.rds,data/DGElist_LUSC.rds
```

Every reference must be a two-group edgeR `DGEList` with enough filtered genes
for `--n-var`. The output directory can point to a mounted shared drive.
Reference-specific parameter caches are keyed by the input file checksum.
Both method-of-moments and edgeR gene-wise dispersion estimates are stored in
the parameter cache and alongside each simulated dataset. The batch runner
uses method-of-moments by default; select `--dispersion-estimator edgeR` to
generate a comparison set using the edgeR estimates. The selected estimator and
both per-gene estimates are recorded in each dataset RDS.

## Design used

The runner generates a **full factorial design**: every combination of these
levels, separately for Poisson and negative-binomial count models:

- Samples per group: 3, 5, 10, 30, 50.
- Minimum fold-change magnitude: 1, 2, 4.
- Upregulated fraction among DE genes: 0.25, 0.50, 0.75.
- Library-size variation (sdlog): 0, 0.25, 0.50.
- Outlier fraction: 0, 0.001, 0.01.
- Extra-zero sparsity: 0, 0.10, 0.30.
- Batch-effect multiplier: 1, 1.5 (batch fraction stays 0.50).
- Dispersion scale: 0.5, 1, 2 for negative binomial; fixed at 1 for Poisson
  (dispersion does not affect Poisson generation).

DE mode 3 is now the default and crosses every combination with all seven
DE proportions. There are 2,430 Poisson and 7,290 negative-binomial base
combinations, or **68,040 balanced datasets per replicate** after DE expansion.
With 100 replicates and one reference, the balanced-only subset produces
**6,804,000 datasets**.
Reference mode `all` multiplies this by the number of references. Shards divide
this total; they do not multiply it.

Sample imbalance is **included by default** and crossed with all other
factors. Each nominal size `s` keeps total samples at `2*s`, with balanced,
1:3, and 1:7 splits. Minority sizes round to nearest integer (halves upward).
Splits with fewer than 3 samples in either group are omitted; duplicate
rounded splits are removed. The actual splits are:

- s=3: 3/3.
- s=5: 5/5, 3/7.
- s=10: 10/10, 5/15, 3/17.
- s=30: 30/30, 15/45, 8/52.
- s=50: 50/50, 25/75, 13/87.

These 12 sample splits yield **163,296 datasets per replicate** in default
DE mode 3: 40,824 Poisson and 122,472 negative-binomial. With 100 replicates
and one reference, the total is **16,329,600 datasets**. Use `--no-imbalance`
for the balanced-only subset; `--include-imbalance` explicitly selects the
default. `sample_size_baseline` records nominal `s`; `n_group1`, `n_group2`,
and `n_samples` record actual sizes. For unequal groups the legacy `s` field
is `NA`. The simulator continues to accept arbitrary explicit `group_sizes`.

A minimum fold change of 1 allows weak effects: DE magnitudes are the minimum
plus an exponential draw, with reciprocals for downregulated genes.
Scenario IDs encode all factor levels and the design is recorded as
`full_factorial` in manifests and dataset metadata. Existing output files are
not relabeled. Use a separate `--output` folder if it already contains results
from the previous one-factor-at-a-time design.

The planner retains only the base combinations and at most `--batch-size`
expanded job records (default 1000). It computes each shard's job IDs directly,
cycling through blocks of equal sample counts to distribute workload, without
expanding or sorting the complete experiment. Changing batch size changes
neither job ownership nor dataset IDs, DE draws, references, or seeds.
Random DE modes retain one RNG checkpoint per base scenario/setting rather
than per replicate; this metadata stays bounded by the base design.

`--plan-only` prints exact total/shard counts and up to 10 example jobs. It
requires no simulation packages, creates no output directories, and does not
expand the full experiment. Startup and each processing batch print progress.

## Distribute and resume work

Ask for a plan before computing. For example, split the jobs over four workers:

```sh
Rscript scripts/custom_simulation/run_batch.R --plan-only --shard-index 1 --shard-count 4
Rscript scripts/custom_simulation/run_batch.R --shard-index 1 --shard-count 4
```

On other workers, use the same project version, reference RDS, `--n-var`,
`--replicates`, `--base-seed`, `--shard-count`, and `--output`, changing only
`--shard-index` to 2, 3, or 4. Give each worker a distinct shard index. A
completed dataset is skipped on restart; failed runs are recorded in
`status.csv` and attempted again. Resume uses completed dataset files rather
than loading a manifest. The planner also checks sibling shard directories
under the same output root, so existing v3.0 datasets from the original
planner can be reused with the same experiment options. Seeds, DE selections,
and filenames remain compatible, although shard ownership and execution order
have changed. Stop old workers before restarting with the new planner. Use `--replicates 3` to generate three
independent datasets per setting. `--no-restore` skips the automatic
`renv::restore()` when the environment is already prepared.

### Run 100 replicates with bounded memory

```sh
# Preview worker 1 of 4, without generating datasets
Rscript scripts/custom_simulation/run_batch.R \
  --replicates 100 --shard-count 4 --shard-index 1 --plan-only

# Run worker 1; launch workers 2, 3, and 4 separately
Rscript scripts/custom_simulation/run_batch.R \
  --replicates 100 --shard-count 4 --shard-index 1 --batch-size 1000
```

The default design still contains 16,329,600 datasets per reference with 100
replicates, or 4,082,400 per worker with four shards. Reducing `--batch-size`
reduces job-record memory without changing the experiment. Count matrices,
reference data, and R packages also require memory; this setting bounds job
bookkeeping, not the entire R process.

Each shard appends `running`, `complete`, or `failed` events to `status.csv`.
The latest event for an `output_file` is its current status; repeated complete
events on resume are expected. Existing `manifest.csv` files are left untouched
and are not rewritten. No full pending-job manifest is built. Analysis still
uses its separate analysis manifest; simulation datasets remain compatible.

## Ten-dataset method smoke test

To check that all 12 methods run on a small spread of simulated data, use:

```sh
Rscript scripts/custom_simulation/smoke_test_10.R
```

The script creates five negative-binomial and five Poisson datasets, runs every
method's baseline on all ten, then runs the `quick` method-setting alternatives
on the first dataset. It writes generated data, per-method results, and a
`method_status.csv` with errors under `results/custom_simulation_smoke/`. It
exits with status 1 if any method or setting fails; failed runs stay in the
CSV so the issue is visible.

## Compare differential-expression method settings

After simulation, `analyze_batch.R` runs the project's 12 differential
expression methods on each generated dataset. It compares each method's
baseline with one-factor-at-a-time settings, avoiding the enormous Cartesian
product of all method options.

```sh
Rscript scripts/custom_simulation/analyze_batch.R --plan-only --shard-count 8
Rscript scripts/custom_simulation/analyze_batch.R --shard-index 1 --shard-count 8
```

Run the second command on each worker with a distinct index from 1 to 8. Use
the same input, output, and `--profile` values on every worker. `--profile quick`
runs each method baseline and one representative alternative; `core`
(the default) includes all listed alternatives. Results are saved individually
under `results/custom_simulation_analysis/results/<dataset>/`, and each
analysis shard has a resumable `manifest.csv`.

The core profile covers these options where supported:

- **DESeq2:** dispersion fit (`parametric`, `local`, `mean`), size-factor
  estimation (`ratio`, `poscounts`, `iterate`), independent filtering, and
  Cook's cutoff/outlier handling.
- **DESeq:** dispersion fit type (`local`, `parametric`).
- **edgeR LRT/QLF:** TMM, TMMwsp, RLE, upper-quartile, or no composition
  normalization; robust dispersion estimation.
- **DSS:** log-ratio, quantile, total-count, or median normalization; mean
  versus abundance-trended dispersion; equal versus separate group variance.
- **voom:** edgeR normalization methods, quantile versus no between-array
  normalization, robust empirical Bayes, and an intensity trend.
- **Wilcoxon / t-test:** normalization, log-CPM prior count, exact versus
  asymptotic Wilcoxon p-values, and equal-variance versus Welch t-tests.
- **NBPSeq:** Anders-Huber normalization versus no normalization; NBQ, NBP,
  and NB2 dispersion models.
- **NOISeq:** TMM, upper-quartile, or no normalization; zero replacement,
  permutation count, smoothing, and low-count filtering.
- **ABSSeq:** package default and TMM, qtotal, total, quartile, or geometric
  normalization; outlier replacement and aFold options.
- **ROTS:** TMM versus no composition adjustment, bootstrap count, and top-list
  size.

Each saved result includes the method settings, gene-level scores, BH adjusted
scores, known DE labels, AUC, true-positive rate, observed FDR, and simulation
and analysis seeds. NOISeq produces `1 - probability` scores rather than
classical p-values, so its BH and observed-FDR summaries are descriptive and
should not be interpreted as calibrated FDR control. Normalization remains
package-specific:
DESeq2 varies size-factor algorithms, while edgeR and voom vary library
normalization factors. Failed settings remain marked in the manifest so they
can be inspected and retried.

Each shard writes into its own directory, for example:

```text
results/custom_simulation_v5.1/
├── parameters-<reference-md5>.rds
├── shard-001-of-004/
│   ├── status.csv
│   └── <scenario>-rep001.rds
  └── shard-002-of-004/
    ├── status.csv
    └── <scenario>-rep001.rds
```

Dataset filenames include the scenario, replicate, gene count, seed, and a
short reference-data checksum. Do not run the same shard from multiple workers
at once. Different shard directories can be written concurrently on a shared
mounted filesystem.

## Shared storage recommendation

Keep code and small CSV manifests in Git. Keep generated `.rds` datasets in a
shared drive or object store mounted by the workers; set `--output` to that
location. This works well when datasets are frequently regenerated or large.

If by “GHLF” you mean GitHub LFS, it is suitable for a modest, stable result
set when every collaborator needs versioned files. Use a separate private data
repository, enable Git LFS, and track only the generated result format, for
example:

```sh
git lfs install
git lfs track "*.rds"
git add .gitattributes
```

GitHub currently includes 10 GiB of LFS storage and download bandwidth on Free
and Pro accounts; Team and Enterprise Cloud include 250 GiB each. Downloads
count against bandwidth, and every changed version of a dataset consumes more
storage, so repeated large outputs can outgrow the included quota. Check the
current [Git LFS billing details](https://docs.github.com/en/billing/concepts/product-billing/git-lfs)
before committing a large batch. For bigger or frequently changing simulation
outputs, use shared/object storage and keep a manifest plus checksums in Git.

## DE proportion selection

Use `--de-mode 1`, `--de-mode 2`, or `--de-mode 3` (names below also work).
All modes use `0.1,0.2,0.3,0.4,0.5,0.6,0.75` unless overridden with
`--de-proportions 0.1,0.3,0.5`.

1. **combination**: independently select one proportion for each scenario
   (a combination of parameter values), per replicate. Draws can coincide.
2. **setting**: select one proportion for each combination of non-sample
   factors, per replicate. Hold it fixed across sample-size levels (including
   imbalance splits). Both count models share the draw when their
   non-sample factors match, including dispersion scale 1.
3. **fixed** (default): run the entire base setting set at each listed proportion.
   Each full pass fixes the proportion for every setting and level. Use a
   single value (e.g. `--de-proportions 0.3`) to fix the whole run to one value.

`--replicates N` repeats the experiment N times. Random modes redraw on each
replicate, reproducibly from `--base-seed`. Proportions are assigned before
reference expansion and sharding, so references and workers agree. These modes
select the proportion; gene identities are still sampled by each simulation.
The old separate DE-proportion factor sweep is replaced by this mode selection.
Mode 1/2 creates 23,328 scenarios per replicate; mode 3 creates 163,296 with all seven
proportions (with default imbalance), before multiplying by references in reference mode `all`.
Mode and selected proportion are saved in scenario IDs, manifests, and RDS
scenario metadata. Existing output datasets are left untouched; use a fresh
`--output` directory for the new experiment to keep analyses separate.

```sh
# Optional: one DE proportion shared across sample-size levels
Rscript scripts/custom_simulation/run_batch.R --de-mode 2 --replicates 3

# Entire suite at every requested DE proportion, with three repetitions
Rscript scripts/custom_simulation/run_batch.R --de-mode 3 --replicates 3
```

## Simulation version

The current simulation version is **5.1**. The default generation output and
analysis input are `results/custom_simulation_v5.1`. Use `--output PATH` on the
generator or `--input PATH` on the analysis runner to override these paths.
Each generated dataset RDS and each manifest row record `simulation_version`.
Existing datasets are not moved or relabeled.

```sh
Rscript scripts/custom_simulation/run_batch.R --version
```

## Verify simulation output

Run the checker with the **same experiment options** used by the generator. It
recreates expected jobs one bounded batch at a time, checks every expected RDS
path across the shard folders, and exits nonzero if datasets are missing. The
check does not require `edgeR` or restore the R environment.

```sh
Rscript scripts/custom_simulation/check_batch.R \\
  --replicates 100 --shard-count 4
```

Add `--deep` to read every RDS and validate simulation version, scenario ID,
seed, count dimensions and values, group sizes, and the number of DE genes.
Deep checking reads and decompresses all files, so it can take a while. Missing paths and deep validation failures are written to timestamped
`verification-missing-*.txt` and `verification-invalid-*.txt` reports. If checking an experiment that overrode any design
flags (such as `--reference-mode`, `--n-var`, `--de-mode`, or
`--no-imbalance`), pass the same values to the checker.

## Underdispersion extension

The default factorial includes Generalized Poisson at target variance/mean
ratios 0.2, 0.5, and 0.81, alongside Poisson at 1 and the three NB scales.
Select `--underdispersion-model generalized_poisson` (default) or
`--underdispersion-model com_poisson`. This selects one generator, not both;
the total number of combinations is the same for either choice. The model is
encoded in IDs and saved as `count_model` in datasets and status records.
The Generalized-Poisson implementation uses base R only. The optional
COM-Poisson dependency COMPoissonReg 0.8.2 is locked in renv.lock and restored
and checked by the Windows setup script. Manual setup uses `renv::restore()`.
Both choices are available in package-free `--plan-only` mode.
There are 285,768 combinations per replicate with imbalance, or 28,576,800
with 100 replicates and one reference. Earlier totals above describe the
Poisson/NB subset only. Underdispersion adds 122,472 combinations per replicate.

COM-Poisson parameters are calibrated numerically for each conditional expected
count to match E(Y)=mu and Var(Y)=k*mu; parameters are not interchangeable with
NB dispersion scales. This is substantially slower than Poisson/NB generation.
The target applies before optional outliers and zero insertion, and conditional
on group, library size, and batch. Those modifications can increase the
observed variance/mean ratio. For a clean dispersion comparison, use zero
outlier fraction, zero sparsity, and fixed library/batch multipliers.
For integer-valued counts with noninteger mean mu, variance cannot be less
than frac(mu)*(1-frac(mu)). Infeasible targets and failed numerical calibrations
fail explicitly and are recorded in the status log.

### Generalized-Poisson equations and numerical checks

The sampler normalizes finite-support weights
`w(y) = theta * (theta + b*y)^(y-1) * exp(-theta-b*y) / y!`,
for integers `y >= 0` satisfying `theta+b*y > 0`, with `theta > 0`, `b < 0`.
It samples with `p(y) = w(y) / sum(w)`. Nominal starting parameters are
`b = 1 - 1/sqrt(k)` and `theta = mu/sqrt(k)`. Finite-support normalization can
change the nominal mean and variance at low means, so the code computes actual
moments and calibrates theta and b when necessary. Actual mean and variance
must each be within 0.1% of mu and k*mu, respectively. Infeasible targets or
failed calibration raise errors. Support above one million points is rejected
explicitly to bound numerical work. A zero mean produces zero counts.

```sh
Rscript scripts/custom_simulation/run_batch.R --replicates 100 \
  --shard-count 4 --shard-index 1 --underdispersion-model generalized_poisson

Rscript scripts/custom_simulation/run_batch.R --plan-only \
  --underdispersion-model com_poisson
```

## Process 100 settings at a time

Run `Rscript scripts/custom_simulation/run_pipeline.R --replicates 100`.
The pipeline generates all replicates for 100 settings, runs the core method
analysis, retains per-replicate metrics and reproducibility metadata, and
removes simulated datasets after successful analysis. Gene-level analysis
results are removed to bound storage. Results live under
`results/custom_simulation_pipeline_v5.1/batch-*/analysis/`. Completed batches are
skipped on restart with the same arguments. Failed batches keep their datasets.
Stop existing generation workers before switching to this pipeline.

For generation alone, use `run_batch.R --setting-batch 1 --settings-per-batch 100`.
A setting includes its DE proportion; all selected replicates and references
are included. Original job IDs and seeds remain unchanged. Regenerate one
original job with `run_batch.R --job-index 12345` and the original experiment
options, including replicate count, seed, reference options and gene count.
Use a separate output directory when recovering a deleted dataset. Preserve
reference files, parameter caches, renv.lock and the simulator code version.

For the random-reference, edgeR-dispersion experiment on five workers:

```sh
Rscript scripts/custom_simulation/run_pipeline.R \
  --replicates 100 --reference-mode random \
  --shard-count 5 --shard-index 1 \
  --dispersion-estimator edgeR --no-restore
```

Use shard indices 1 through 5 on separate workers. Each pipeline worker owns
its own configuration, completion markers and analysis directory under
`results/custom_simulation_pipeline_v5.1/shard-001-of-005/` (and equivalent shards).

Settings now belong entirely to one worker (`--shard-unit setting`, the
default): replicas 1 through 100 run consecutively on that worker. The legacy
replicate-level assignment remains available as `--shard-unit job`. Stop all
old workers before changing assignment. Do not reuse completed pipeline batch
markers from the old assignment; start with a fresh pipeline directory.
The progress denominator refers to the selected batch, not the full experiment.

## Import datasets from the old flat shard folders

Stop all generation and pipeline workers first. Preview the migration using
exactly the experiment flags used to generate the datasets:

```sh
Rscript scripts/custom_simulation/migrate_pipeline.R \
  --source results/custom_simulation_pipeline_v5.1 -- \
  --replicates 100 --reference-mode random --shard-count 5 --shard-index 3 \
  --dispersion-estimator edgeR --no-restore
```

Repeat with `--apply` before `--` to perform the moves. The source may be the
whole old results root (recommended, covering every old worker) or a single
shard directory. The default destination is the versioned pipeline root.
Existing completed batches and conflicting files are protected. The script
reads and validates every dataset before moving any, records a CSV move plan,
and preserves original datasets and metadata without recompressing them.
Files are reassigned by setting to the correct worker and batch, regardless
of their old worker. No completion markers are created: the pipeline checks
existing dataset filenames, generates missing jobs, then analyzes the batch.
Restart the pipeline with the same experiment flags, separately for each
worker. Migration does not run analysis itself. Old status logs remain in
place; fresh batch status logs are written when the pipeline resumes.

## Copy results to a portable disk

In PowerShell, replace `E:` with your portable disk drive letter:

```powershell
Rscript scripts/custom_simulation/export_results.R --destination "E:/RNA-results"
```

The source defaults to `results/custom_simulation_pipeline_v5.1`. Override with
`--source "C:/path/to/results/custom_simulation_pipeline_v5.1"` if needed.
The export includes saved analysis result RDS files, analysis manifests,
pipeline run logs, worker options and completed-batch markers. It excludes simulated datasets,
parameter caches and simulation logs. It also includes saved results from
unfinished batches; use complete.rds and manifest statuses to assess completion.
For a stable snapshot, stop workers first. A file that changes while copying
causes an error; rerun to retry. Copies are checksum-verified and unchanged
files are skipped when rerunning. Source files are never removed or modified.
Analysis RDS filenames are shortened to content checksums for Windows path
compatibility; export_index.csv maps original paths to exported paths.
Embedded original input paths are historical metadata and need not exist on
the analysis computer. Each RDS still contains its dataset/method/settings IDs
and metrics. The export does not aggregate replicates. Rerunning adds or
updates results; it does not delete older exports from the portable disk.

## Pipeline console status and run logs

Each pipeline invocation creates a log under its worker directory:
`logs/pipeline-YYYY-MM-DD_HH-MM-SS-PID.log`. The start date/time identifies the
invocation; restarts create a new file. Stage messages show GENERATE, ANALYZE,
COMPACT, CLEANUP and BATCH COMPLETE, with blank lines and a banner between
batches. Generation and analysis output is also saved to the run log.
Completed dataset/result skips are summarized as counts rather than printed
per filename; completed pipeline batches are summarized on resume. Detailed
per-job states remain in status.csv and the analysis manifest.csv. Standalone
generation/analysis commands keep their usual console output; automatic run
logs are created by run_pipeline.R.
