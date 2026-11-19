# Architecture and agent handoff

## Start here

Read `README.md`, this page and `custom_simulation.md`. For broader paper and
source-dataset caveats, read `../HANDOFF.md` in the enclosing agent-rick workspace
(the file is outside this nested repository). Use `renv.lock` for dependencies.

## Responsibilities

| Location | Responsibility |
| --- | --- |
| `scripts/simulation_design/factorial_design.R` | Factorial scenarios and group-size design |
| `scripts/simulation_design/de_proportion_modes.R` | Reproducible DE truth assignment |
| `scripts/simulation_design/batch_jobs.R` | Bounded-memory job iteration and sharding |
| `scripts/simulation_counts/load.R` | Load simulator modules in dependency order |
| `scripts/simulation_counts/estimate_parameters.R` | Estimate reference expression and dispersion |
| `scripts/simulation_counts/prepare_design.R` | Prepare simulation parameters and group design |
| `scripts/simulation_counts/generalized_poisson.R` | Underdispersed count sampler |
| `scripts/simulation_counts/generate_counts.R` | Generate and perturb counts |
| `scripts/simulation_counts/simulate_counts.R` | Public simulation function |
| `scripts/simulation_counts/helpers.R` | Simulator validation and helper functions |
| `scripts/de_analysis/method_settings.R` | DE method profiles and execution helpers |
| `scripts/shared/` | Version, path checks and logging |
| `scripts/simulation_pipeline/run_pipeline.R` | Generate, analyze, compact and clean up batches |
| `scripts/simulation_counts/run_batch.R` | Standalone generation and planning |
| `scripts/de_analysis/analyze_batch.R` | Resumable method analysis |
| `scripts/simulation_counts/check_batch.R` | Verify generated batch coverage |
| `scripts/de_analysis/smoke_test_10.R` | Exercise all methods on ten datasets |
| `scripts/simulation_pipeline/migrate_pipeline.R` | Preview/apply old result-layout migration |
| `scripts/result_transfer/export_results.R`, workspace-root `import_results.R` and `import_exports.py` | Portable result transfer |
| `automations/setup_windows.ps1` | Windows environment checks |

## Refactor record — 2026-11-13

The former `scripts/custom_simulation/` files were separated by responsibility
under `scripts/`, using task-based snake_case folder names consistent with the
existing paper workflows. The interim top-level `R/` layout was removed at the
user’s request.
The simulator modules formerly in `scripts/shared/custom_simulation/` now live
in `scripts/simulation_counts/`, with `load.R` replacing `scripts/shared/simulate_custom.R`.
R runners live in the relevant task folders under `scripts/`.
`automations/` is reserved for shell and PowerShell launchers (`.sh`, `.bat`,
`.ps1`); it contains no R files. Update external worker commands
to these paths. No compatibility wrappers are kept at the old paths.
The paper reproduction scripts remain in `scripts/`; their historical setup
requirements are documented in `legacy_benchmark.md`.

Existing uncommitted work was preserved during relocation. This is a source
layout change: simulation version, seeds, job IDs, output directory names,
reference inputs and result formats are retained. Do not rename old results
merely to match source folders. Stop existing workers before updating code.
The batch checker historically defaults to version 4.0 while generation uses
`scripts/shared/version.R`; pass an explicit `--output` when checking v5.1.

## Validation and maintenance

Run these from the project root:

```sh
Rscript --vanilla tests/full_factorial.R
Rscript --vanilla tests/batch_jobs.R
Rscript --vanilla tests/dss_settings.R
Rscript tests/underdispersion.R
Rscript --vanilla scripts/simulation_counts/run_batch.R --plan-only
```

The ten-dataset smoke test requires the restored R packages and runs actual DE
methods. Planning alone does not validate their numerical results. Do not
launch the full factorial experiment merely to test path changes.

Keep reusable functions in task folders under `scripts/`, R execution and filesystem orchestration in their relevant `scripts/` task folders,
shell launchers in `automations/`, and research/handoff material in `docs/`. The root README is the
entry point. Preserve resumability and RNG behavior when editing job planning.

## Editable machine runner — 2026-11-16

`config/machine.dcf` configures inclusive dataset job IDs, worker count and
experiment settings. `scripts/simulation_pipeline/run_machine.R` uses a local
socket process pool; shell launchers remain in `automations/`. See
`machine_config.md` for resume, resource limits and retained outputs.
`run_batch.R --job-start N --job-end M` decodes original IDs directly without
scanning preceding jobs. Range selection cannot be combined with sharding,
setting batches or single-job selection. `tests/machine_config.R` verifies
range identity against the existing iterator across DE/reference modes.

Method pools: `Methods:` in the DCF config is an indented list with full-line
`#` comments. `read_machine_config` removes commented lines and validates the
active list. The machine runner passes it as `analyze_batch.R --methods`;
selection occurs before analysis job expansion. The default example is
NBPSeq only. Keep identical experiment fields across pools, overlap their
overall ranges, and split ranges without overlap within each method pool.
Changing selected methods requires a fresh Output directory. Validation:
`Rscript --vanilla tests/method_selection.R`.
