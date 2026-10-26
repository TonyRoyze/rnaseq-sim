# benchDE

Code for benchmarking 12 differential-expression (DE) methods in bulk RNA-seq
data under standard, outlier, and noise settings. The project includes
negative-binomial and Poisson simulations together with false-positive-rate
control and ranking-stability analyses in TCGA-BRCA and SCLC data.

This README documents the scripts currently contained in
`data_scripts_cleaned`. No additional analysis scripts are assumed.

## Project Structure

### Public Data Preparation

- **`local_data_download/get_datasets.R`**: Downloads or imports the public
  E-ENAD-34, GSE150910, TCGA-LUSC, and TCGA-BRCA datasets; prepares filtered
  `edgeR::DGEList` objects; and saves them under `<BENCHDE_DIR>/datasets/`.

The SCLC workflow expects the following source files under
`<BENCHDE_DIR>/datasets/`:

- `count_ensg_STAR.txt.gz`
- `Cell2023_SCLC_data_clinical.txt`

### Negative-Binomial Simulations

#### Standard simulation

- **`nb_standard/sim1_function.R`**: Estimates mean and dispersion parameters
  from real count data and generates standard negative-binomial simulations.
- **`nb_standard/sim1-assessment.R`**: Runs the standard negative-binomial
  benchmark across sample sizes and public reference datasets; writes
  simulated counts, DE results, and AUC values to `sim1/`.

#### Outlier simulation

- **`nb_outlier/sim2_function.R`**: Generates negative-binomial counts and
  introduces random outlier observations at the specified proportions.
- **`nb_outlier/sim2-assessment.R`**: Runs the negative-binomial outlier
  benchmark and writes `counts/`, `datp/`, and `auc/` outputs to `sim2/`.

#### Noise simulation

- **`nb_noise/sim3_function.R`**: Generates negative-binomial counts with
  multiplicative noise at the specified noise levels.
- **`nb_noise/sim3_assessment.R`**: Runs the negative-binomial noise benchmark
  and writes outputs to `sim3_restored/`.

### Poisson Simulations

#### Standard simulation

- **`poisson_standard/poisson_function.R`**: Generates standard Poisson counts
  using expression-derived mean parameters.
- **`poisson_standard/poisson_assessment.R`**: Runs the standard Poisson
  benchmark and writes outputs to `poisson/`.

#### Outlier simulation

- **`poisson_outlier/poisson_sim2_function.R`**: Generates Poisson counts and
  applies the random-outlier mechanism.
- **`poisson_outlier/poisson_sim2_assessment.R`**: Runs the resumable Poisson
  outlier benchmark and writes outputs to `poisson_sim2/`.

#### Noise simulation

- **`poisson_noise/poisson_sim3_function.R`**: Generates Poisson counts after
  perturbing the expected count with multiplicative uniform noise.
- **`poisson_noise/poisson_sim3_assessment.R`**: Runs the Poisson noise
  benchmark and writes outputs to `poisson_sim3/`.

### Real-Data Analysis

#### TCGA-BRCA

- **`brca/fpr_assessment.R`**: Creates null comparisons from TCGA-BRCA samples
  and evaluates empirical false-positive-rate control over sample sizes 3, 5,
  10, 30, and 50 per group. Outputs are written to `FPR/`.
- **`brca/CAT_assessment.R`**: Repeatedly samples independent TCGA-BRCA subsets
  and evaluates concordance of method-specific gene rankings. Counts, DE
  results, and concordance outputs are written to `CAT/`.

#### SCLC

- **`sclc/sclc_preprocess.R`**: Matches the SCLC count and clinical tables,
  filters genes, constructs the SCLC `DGEList`, and prepares FPR and stability
  benchmark inputs.
- **`sclc/sclc_process.R`**: Runs the 12 DE methods for SCLC FPR and CAT
  analyses with batched parallel processing, checkpoints, and resumable result
  files. Outputs are written to `FPR_SCLC/` and `CAT_SCLC/`.

### Shared Analysis Utilities

- **`shared/runDE.R`**: Defines the common interface for ABSSeq, DESeq, DESeq2,
  DSS, edgeR likelihood-ratio and quasi-likelihood tests, NBPSeq, NOISeq, ROTS,
  Student's t-test, voom, and Wilcoxon rank-sum test.
- **`shared/sim_DEanalysis_auc.R`**: Runs the configured DE methods on simulated
  count objects, stores p-values, and calculates AUC values against the known
  differential-expression labels.
- **`shared/get_conc_function.R`**: Calculates cross-replicate concordance of
  ranked gene lists for the real-data stability analyses.

### Plotting

- **`plot/fig2_plot.R`**: Figure 2, negative-binomial standard-simulation AUC.
- **`plot/fig3_plot.R`**: Figure 3, negative-binomial FDR-TPR curves.
- **`plot/fig4_plot.R`**: Figure 4, Poisson standard-simulation AUC.
- **`plot/fig5_plot.R`**: Figure 5, Poisson FDR-TPR curves.
- **`plot/fig6_plot.R`**: Figure 6, negative-binomial outlier robustness.
- **`plot/fig7_plot.R`**: Figure 7, Poisson outlier robustness.
- **`plot/fig8_plot.R`**: Figure 8, negative-binomial noise robustness.
- **`plot/fig9_plot.R`**: Figure 9, Poisson noise robustness.
- **`plot/fig10_plot.R`**: Figure 10, TCGA-BRCA ranking stability.
- **`plot/fig11_plot.R`**: Figure 11, SCLC ranking stability.
- **`plot/fig12_plot.R`**: Figure 12, TCGA-BRCA FPR control.
- **`plot/fig13_plot.R`**: Figure 13, SCLC FPR control.
- **`plot/fig14_plot.R`**: Figure 14, combined bubble panels summarizing
  simulation accuracy/robustness and real-data FPR/stability results.
- **`plot/plot_aes.R`**: Shared method ordering, colors, and ggplot styling.

## Software Requirements

- R 4.3.3
- A Unix-like system for scripts that use `parallel::mclapply`
- R packages and tested versions listed in **`R_packages.txt`**

### Isolated R environment

This repository uses `renv` to keep its R packages in a project-local library.
The pinned environment is recorded in `renv.lock` and uses R 4.3.3 with
Bioconductor 3.18. Open R from this directory so `.Rprofile` activates the
project library:

For Windows installation and shard-running instructions, see [setup.md](setup.md).
Windows users can run [setup_windows.ps1](setup_windows.ps1) to check Git LFS,
the reference data, R 4.3.3, the restored packages, and the simulation runner.

```bash
rig run -r 4.3-arm64
```

To restore the lockfile on another setup, run `renv::restore()` from the
project's R session.

The current lockfile uses `dplyr` 1.2.1 because `TCGAbiolinks` 2.30.4 fails to
load with the listed `dplyr` 1.1.4 (`filter_out` is not exported there). It
uses DSS 2.48.0 from the Bioconductor 3.17 source archive and the legacy DESeq
1.36.0 source from Bioconductor 3.9. The `zoo` version written as `1.8.12` in
`R_packages.txt` is installed as CRAN's `1.8-12`.

The benchmark compares the following methods:

`ABSSeq`, `DESeq`, `DESeq2`, `DSS`, `edgeR.lrt`, `edgeR.qlf`, `NBPSeq`,
`NOISeq`, `ROTS`, `T.test`, `voom`, and `Wilcoxon`.

## Project Layout for Running the Scripts

Set the project root before running an analysis:

```bash
export BENCHDE_DIR=/path/to/MultiDimDE
```

The scripts expect this working layout:

```text
<BENCHDE_DIR>/
├── datasets/          # Processed DGEList objects and SCLC source tables
├── utils/             # Shared and simulation helper scripts
├── sim1/
├── sim2/
├── sim3_restored/
├── poisson/
├── poisson_sim2/
├── poisson_sim3/
├── FPR/
├── CAT/
├── FPR_SCLC/
└── CAT_SCLC/
```

Before execution, copy or link the relevant function files and the three
scripts in `shared/` into `<BENCHDE_DIR>/utils/`, because the assessment scripts
source helper files from that location.

## Suggested Run Order

1. Prepare public datasets with `local_data_download/get_datasets.R` and place
   the two SCLC source tables in `<BENCHDE_DIR>/datasets/`.
2. Run the six simulation assessment scripts to generate `counts/`, `datp/`,
   and `auc/` results.
3. Run `brca/fpr_assessment.R` and `brca/CAT_assessment.R`.
4. Run `sclc/sclc_preprocess.R`, followed by `sclc/sclc_process.R`.
5. Run the corresponding scripts under `plot/` after their RDS inputs have
   been generated.

Example:

```bash
Rscript local_data_download/get_datasets.R
Rscript nb_standard/sim1-assessment.R
Rscript brca/fpr_assessment.R
Rscript plot/fig2_plot.R
```

Some plotting scripts expose input locations through environment variables;
others retain an explicit historical server path near the beginning of the
file. Set the documented environment variable or update that input assignment
to the location of the generated results before plotting.

## Default Benchmark Settings

- random seed: 8848 in the revision workflows;
- simulation and FPR replicates: 100;
- CAT replicates: 10;
- sample sizes per group: 3, 5, 10, 30, and 50;
- outlier proportions: 0.1%, 0.5%, 1%, 2%, and 5%;
- noise levels: 0, 0.2, 0.4, 0.6, and 0.8.

## Data and Results

Raw public datasets and generated RDS files are not included in the code
repository. The supplied `gitignore.txt` excludes local datasets, simulation
outputs, real-data benchmark outputs, figures, logs, and R session files.
