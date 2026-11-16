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
