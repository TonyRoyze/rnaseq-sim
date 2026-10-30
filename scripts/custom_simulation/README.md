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
`results/custom_simulation/`. The runner supports three reference modes:

- `--reference-mode brca` (default): use only `--reference` (defaults to
  `data/DGElist_BRCA.rds`).
- `--reference-mode all`: create a separate job for every reference and
  scenario combination.
- `--reference-mode random`: choose one reference reproducibly for each
  scenario replicate. The chosen reference and checksum are saved in both the
  result RDS and `manifest.csv`.

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

The runner avoids the enormous Cartesian product of every possible parameter.
It creates a one-factor-at-a-time suite: a baseline plus scenarios that change
one parameter at a time. Reference mode `all` multiplies each setting by the
number of selected references; `random` assigns exactly one selected reference
to each setting replicate. There is a negative-binomial and a Poisson family;
dispersion scaling is only varied for negative-binomial simulations. The
settings include balanced sample sizes 3, 5, 10, 20, 30, and 50 per group, dispersion scales
0.25/0.5/1/2/4, library-size spread, DE proportion, fold-change strength,
direction balance, outlier fraction, sparsity, and batch effect. This produces
44 base settings: 24 negative-binomial and 20 Poisson. DE selection is
applied to these settings using the mode below. The
negative-binomial baseline (dispersion scale 1) is included once; its dispersion
sweep adds 0.25, 0.5, 2, and 4. The Poisson family has no dispersion parameter.

The baseline uses 20 samples per group (40 total), a mode-selected DE proportion,
minimum fold change 2, and upregulated fraction 0.50. DE proportions are
0.10/0.20/0.30/0.40/0.50/0.60/0.75 and minimum fold changes are 1/2/4.
A minimum of 1 allows weak effects; DE magnitudes are the minimum plus an
exponential draw, with reciprocals for downregulated genes.
Sample imbalance scenarios use 10/30 and 5/35, compared with the 20/20
baseline, keeping 40 total samples. Other settings stay at baseline in these
comparisons. Group sizes and total sample counts are recorded in manifests.
The simulator accepts `group_sizes = c(10, 30)`; omitting it preserves the
existing `s` samples per group API. For unequal sizes, returned `s` is `NA`.

Jobs are assigned reproducibly to shards using a greedy workload estimate
based on `sample size × gene count`. This spreads the larger sample-size runs
among workers. The complete setting and seed are saved inside each `.rds` file
and in the shard's `manifest.csv`.

## Distribute and resume work

Ask for a plan before computing. For example, split the jobs over four workers:

```sh
Rscript scripts/custom_simulation/run_batch.R --plan-only --shard-index 1 --shard-count 4
Rscript scripts/custom_simulation/run_batch.R --shard-index 1 --shard-count 4
```

On other workers, use the same project version, reference RDS, `--n-var`,
`--replicates`, `--base-seed`, `--shard-count`, and `--output`, changing only
`--shard-index` to 2, 3, or 4. Give each worker a distinct shard index. A
completed dataset is skipped on restart; failed runs remain marked in the
manifest and are attempted again. Use `--replicates 3` to generate three
independent datasets per setting. `--no-restore` skips the automatic
`renv::restore()` when the environment is already prepared.

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
results/custom_simulation/
├── parameters-<reference-md5>.rds
├── shard-001-of-004/
│   ├── manifest.csv
│   └── <scenario>-rep001.rds
  └── shard-002-of-004/
    ├── manifest.csv
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
2. **setting** (default): select one proportion per varied factor, per
   replicate. All levels of a factor, such as sample size, share that value,
   including across the Poisson and negative-binomial families. The global
   baseline has its own draw. For a matched comparison, compare levels within
   the same factor; the global baseline may have a different DE proportion.
3. **fixed**: run the entire base setting set at each listed proportion.
   Each full pass fixes the proportion for every setting and level. Use a
   single value (e.g. `--de-proportions 0.3`) to fix the whole run to one value.

`--replicates N` repeats the experiment N times. Random modes redraw on each
replicate, reproducibly from `--base-seed`. Proportions are assigned before
reference expansion and sharding, so references and workers agree. These modes
select the proportion; gene identities are still sampled by each simulation.
The old separate DE-proportion factor sweep is replaced by this mode selection.
Mode 1/2 creates 44 scenarios per replicate; mode 3 creates 308 with all seven
proportions, before multiplying by references in reference mode `all`.
Mode and selected proportion are saved in scenario IDs, manifests, and RDS
scenario metadata. Existing output datasets are left untouched; use a fresh
`--output` directory for the new experiment to keep analyses separate.

```sh
# Default: one DE proportion per factor, shared across its levels
Rscript scripts/custom_simulation/run_batch.R --de-mode 2 --replicates 3

# Entire suite at every requested DE proportion, with three repetitions
Rscript scripts/custom_simulation/run_batch.R --de-mode 3 --replicates 3
```
