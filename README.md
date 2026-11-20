# MultiDimDE simulation

A configurable RNA-seq simulation and differential-expression benchmark.
Start with the custom workflow below; the original paper workflows remain in
`scripts/` and are described in [the legacy benchmark guide](docs/legacy_benchmark.md).

## Layout

```text
scripts/
  simulation_design/  Factorial scenarios, DE proportions, bounded job planning
  simulation_counts/  Parameter estimation, count generation, simulator loader
  de_analysis/        DE methods, profiles, analysis runner and smoke test
  simulation_pipeline/ Batch orchestration and result-layout migration
  result_transfer/    Portable result import and export
  shared/             Shared DE utilities, version, path checks and logging
  nb_*/               Original negative-binomial workflows
  poisson_*/          Original Poisson workflows
  brca/, sclc/        Real-data benchmark workflows
  plot/               Paper figures and plot styling
  local_data_download/ Reference data preparation
automations/      Shell and PowerShell launchers only (.sh, .bat, .ps1)
docs/             Usage, architecture, research plans and agent handoff
tests/            Regression checks for design, jobs, dispersion and DSS
data/             Reference DGEList inputs
results/          Generated datasets and analysis outputs (existing paths retained)
renv/             Project environment bootstrap
```

## Run

For editable job ranges, method pools and six-core execution, edit
`config/machine.dcf` (uncomment the methods assigned to this machine) and
follow [the machine guide](docs/machine_config.md). Preview with
`./automations/run_machine.sh --plan-only` or
`.\automations\run_machine.ps1 --plan-only`.

From this project directory:

```sh
Rscript scripts/simulation_counts/run_batch.R --plan-only
Rscript scripts/simulation_pipeline/run_pipeline.R --replicates 100 --no-restore
```

Review the plan and [workflow guide](docs/custom_simulation.md) before running a
large experiment. The default design is large. Restore the pinned environment
with `renv::restore()` when needed; `--no-restore` assumes it is already ready.
For Windows, see [setup](docs/setup.md) and `automations/setup_windows.ps1`.

For interactive exploration, source the reusable simulator:

```r
source("scripts/simulation_counts/load.R")
source("scripts/simulation_design/factorial_design.R")
source("scripts/de_analysis/method_settings.R")
```

See [docs](docs/README.md) for the full runbook, module map and future-agent notes.

### Sparse normalization outcomes

If no gene has positive counts in every sample, analysis records the following
settings as `unsupported` with the reason in the analysis manifest:

- RLE for edgeR LRT, edgeR QLF, voom, Wilcoxon, and T.test
- Default AH2010 normalization for NBPSeq
- Geometric normalization for ABSSeq
- Legacy DESeq size-factor estimation and DESeq2 ratio size-factor estimation

These settings count toward batch completion and produce no result RDS.
Other settings still run on the same dataset. Unexpected analysis errors stay
`failed`. The runner continues to later batches and leaves the affected
inputs available for retry. When summarizing the benchmark, report unsupported
counts against all assigned replicates and compute metrics from completed fits.
