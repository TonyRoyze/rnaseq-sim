# Machine configuration and six-core execution

Edit `config/machine.dcf` on each computer. It is plain text in DCF format:
one `Name: value` per line, with no blank lines between fields. Full-line `#`
comments are supported, including indented comments in the method list.

```text
Job-Start: 1
Job-End: 600
Workers: 6
Jobs-Per-Batch: 20
```

Bounds are inclusive original **dataset job IDs**, not method-result IDs,
setting IDs or positions in the load-balanced shard order. A job is one
simulated dataset. Every selected dataset runs only the uncommented methods, including all their
settings in the chosen profile, when `Stage: both`. `Stage: generate` only generates datasets.
Original IDs, seeds, DE truth and reference draws are preserved.

Within the same method pool, machine two can use `Job-Start: 601`,
`Job-End: 1200`. Use a different Output for each machine on shared storage.
Allocate nonoverlapping ranges within a method pool; different method pools
should cover the same overall dataset range for comparable results.
Keep all experiment fields identical across machines: Replicates, N-Var,
Base-Seed, reference choices, dispersion, imbalance and DE settings. Changing
Replicates changes what a job ID means. The shipped example uses 100 replicates;
IDs 1–100 are the replicates of the first setting. The launcher prints the
experiment total and rejects ranges beyond it.

## Launch

Prepare the pinned packages once, before starting workers: `renv::restore()`
in an R session opened in the project. Workers skip package restoration so
six processes do not try to install into the same library simultaneously.

PowerShell, from the project root:

```powershell
.\automations\run_machine.ps1 --plan-only
.\automations\run_machine.ps1
```

macOS/Linux:

```sh
./automations/run_machine.sh --plan-only
./automations/run_machine.sh
```

Use `--config config/another-machine.dcf` for another config file. Plan-only
validates the config and prints ranges without launching workers or creating
output. The R entry point is `scripts/simulation_pipeline/run_machine.R`.

## CPU, memory and disk

`Workers: 6` starts six independent processes, each handling its own inclusive
subrange, generating a batch and then analyzing it. Numerical-library thread
limits are set to one per worker to avoid six workers each asking for six
threads. Windows uses the same socket process pool as macOS/Linux.
This is concurrency, not a measured guarantee of optimal throughput. Start
with fewer workers if RAM is limited; each active method has its own R memory
and package overhead. Disk throughput can also limit speed. Workers have
static subranges, so the slowest subrange determines finishing time.

`Jobs-Per-Batch` limits the number of datasets passed to one analysis process;
it does not limit total RAM or total stored output. Reference parameter caches
are shared across local workers and batches, with the generator's existing
cache lock preventing duplicate estimation. Datasets and gene-level analysis
results are retained, and disk usage grows over the full configured range.
This launcher does not perform the older pipeline's compaction or deletion.

## Choose a method pool by uncommenting

At the bottom of `config/machine.dcf`, remove `#` from the method lines assigned
to this machine. Keep the indentation. Comment out any active method you do
not want. At least one method must be active; names are case-sensitive.
The shipped example enables **only NBPSeq**.

```text
Methods:
    # NBPSeq
    DESeq2
    # DESeq
    # edgeR.lrt
    # edgeR.qlf
    # DSS
    # voom
    # Wilcoxon
    # T.test
    # NOISeq
    # ABSSeq
    # ROTS
```

This example runs only DESeq2. Uncomment several lines for a combined pool.
`Profile: core` or `quick` selects settings within each selected method; it
never enables another method. `--plan-only` prints the active methods.
Standalone analysis also supports `--methods NBPSeq,DESeq2`.

To implement the proposed 20-machine layout:

- On the five NBPSeq machines, enable only NBPSeq and divide the full dataset
  range into five nonoverlapping ranges.
- On the five DESeq2 machines, enable only DESeq2 and cover that **same full
  dataset range**, again divided into five ranges.
- On the two DESeq machines, enable only DESeq and divide the same full range
  into two ranges.
- On the other eight machines, enable the remaining nine methods (edgeR.lrt,
  edgeR.qlf, DSS, voom, Wilcoxon, T.test, NOISeq, ABSSeq, ROTS), and divide the
  same full range into eight ranges.

Every machine can run six local workers within its assigned range. Identical
experiment settings and reference files preserve dataset IDs/seeds across
method pools. Each pool generates its own datasets locally; there is no shared
cross-machine dataset distributor. The 5/5/2/8 allocation is supported as a
manual deployment layout, but has not been benchmarked for balanced finishing
times. Faster pools may finish earlier.

## Resume and errors

Each worker writes separate batch directories and generation/analysis logs
under Output. Dataset and method-result files are reused on restart. The
launcher checks generation coverage and analysis manifest completion before
marking a batch complete; failures cause a nonzero final exit code.
`machine_options.rds` records the experiment and machine configuration. Resume
with the same file values. Use a fresh Output directory when changing them,
including changing the uncommented methods. Reordering or editing comments
alone does not change the selected list.
The output lock prevents a second launcher using that output simultaneously;
after a crash, stop all remaining workers before manually removing
`.machine-running`. Stop workers before updating source code. Full numerical
DE validation is separate from the launcher/planner checks.

## Validation — 2026-11-16

Planner range tests and the existing factorial, batch/resume, DSS-setting and
underdispersion checks passed. A real six-process generation smoke test using
a small synthetic edgeR reference completed six distinct jobs successfully.
The full BRCA-reference performance run and six-worker DE-method throughput
have not been benchmarked. PowerShell launcher execution was not tested on
Windows in this session.

Method-selection validation: all 12 methods retain their complete core/quick
profiles when selected; commented methods are excluded before job expansion.
Invalid, duplicate and empty method selections are rejected. Standalone
analysis plan-only validation does not restore packages or create output.

A real two-worker `Stage: both` smoke test completed generation plus all quick
T.test settings for two datasets, with no other methods in the manifests.
Standalone T.test/Wilcoxon selection also completed with no unselected methods.
These checks used a small synthetic reference, not a full throughput benchmark.

## Reduced experiment configuration (version 5.2)

The current factor levels and exact group sizes are listed in
[custom_simulation.md](custom_simulation.md#design-used). DE proportions in
`config/machine.dcf` are `0.05,0.25,0.75`; the example uses 30 replicates,
10,000 genes, six workers, jobs 1–600, and NBPSeq only. Uncomment the other
methods to enable them. The new default Output is `results/machine_run_002`;
use a fresh folder for this experiment and stop existing workers before updating.

`Profile: core` now includes only each method's baseline and the requested
normalization variants. Other fitting, robustness, filtering and resampling
variants are excluded. `quick` selects a smaller subset for smoke checks.

| Method | Baseline | Additional normalization labels |
| --- | --- | --- |
| ABSSeq | Package default | TMM, geometric, qtotal, quartile, total |
| DESeq | Package default | none |
| DESeq2 | Package default | iterate, poscounts |
| DSS | lr | median, quantile, total |
| NBPSeq | AH2010 | none |
| NOISeq | TMM (`tmm`) | n, upper quartile (`uqua`), none |
| ROTS | TMM | none |
| T.test | none | RLE, TMM, TMMwsp, upperquartile |
| Wilcoxon | none | RLE, TMM, TMMwsp, upperquartile |
| edgeR.lrt | TMM | RLE, TMMwsp, none, upperquartile |
| edgeR.qlf | TMM | RLE, TMMwsp, none, upperquartile |
| voom | TMM | RLE, TMMwsp, none, upperquartile |

There are 47 labels across all 12 methods. NOISeq's `n` and `none` labels both
pass `norm="n"`; they are retained as requested but represent the same
normalization. ABSSeq's explicit qtotal variant is retained alongside its
package-default baseline. NBPSeq's none variant passes `NULL`.

The machine config now uses `Reference-Mode: random`: each dataset selects one
of BRCA, E-ENAD-34, GSE150910 or LUSC with equal probability using its seeded
reference-selection stream. SCLC is excluded. This does not generate four
copies of every job. At the configured 30 replicates the full experiment has
136,080 dataset jobs; the example machine still runs only IDs 1–600.
The `Reference` field is used only in `brca` mode; random mode uses the four
paths in the generator's `references` defaults. Keep Base-Seed and experiment
settings identical across method pools so reference draws agree.

## 100 numbered assignments

[config/assignments_100/README.md](../config/assignments_100/README.md) contains
100 ready-to-run configs, the method-pool allocation, initial assignments for
30 machines and a manual assignment register. Pass `--config` explicitly;
the default machine config does not automatically select these files.
