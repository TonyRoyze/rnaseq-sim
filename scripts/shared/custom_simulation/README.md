# Custom simulation modules

The simulator is split into small files so each stage can be read on its own.
Load them through the public entry point from the project root:

```r
source("scripts/shared/simulate_custom.R")
```

That entry point links the modules in this directory and makes two functions
available:

1. `estimate_simulation_parameters(real_dge)` filters and estimates means and
   gene-wise dispersions from a real two-group `edgeR::DGEList`.
2. `simulate_custom_counts(params, ...)` chooses genes and simulation settings,
   generates counts, then returns the count data and known DE truth.

## Files

- `helpers.R`: shared scalar validation and null fallback helper.
- `estimate_parameters.R`: parameter estimation from real count data.
- `prepare_design.R`: validates settings and builds sample, batch, library,
  dispersion, and DE-effect values.
- `generate_counts.R`: draws Poisson or negative-binomial counts and applies
  outliers and extra zero counts.
- `simulate_counts.R`: public interface that connects the stages and packages
  the result as a `DGEList`.

## Example

```r
source("scripts/shared/simulate_custom.R")
params <- estimate_simulation_parameters(real_dge)

sim <- simulate_custom_counts(
  params,
  s = 10,                       # samples per group
  n_var = 10000,                # genes to simulate
  n_diffexp = 3000,             # DE genes
  count_model = "negative_binomial",
  dispersion_scale = 1.5,       # multiply estimated dispersions
  library_size_sdlog = 0.25,    # sample-to-sample library factor spread
  fold_change_min = 1.5,        # minimum DE fold-change magnitude
  upregulated_fraction = 0.6,   # fraction of DE genes going up in group 2
  outlier_fraction = 0.001,     # fraction of individual counts affected
  sparsity = 0.05,              # fraction of individual counts set to zero
  batch_fraction = 0.5,         # probability a sample is in batch 2
  batch_effect = 1.3,           # expected-count multiplier for batch 2
  seed = 1
)

sim$dge       # edgeR DGEList; simulated counts and sample metadata
sim$truth     # gene-level DE labels, directions, and fold changes
sim$samples   # condition, batch, and library factors by sample
sim$parameters
```

Set `count_model = "poisson"` to remove the negative-binomial dispersion
component. You can pass explicit per-sample `library_size_factors` or `batch`
labels instead of having the simulator generate them. The batch effect is
applied to the second batch level; generated labels are `batch1` and `batch2`.
