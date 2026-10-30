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
`results/custom_simulation_v3.0/`. The runner supports three reference modes:

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
results/custom_simulation_v3.0/
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

The current simulation version is **3.0**. The default generation output and
analysis input are `results/custom_simulation_v3.0`. Use `--output PATH` on the
generator or `--input PATH` on the analysis runner to override these paths.
Each generated dataset RDS and each manifest row record `simulation_version`.
Existing datasets are not moved or relabeled.

```sh
Rscript scripts/custom_simulation/run_batch.R --version
```
