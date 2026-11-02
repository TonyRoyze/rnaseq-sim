# Windows setup for the simulation

This guide prepares a Windows computer to generate and analyze datasets in
this repository. The project lockfile pins **R 4.3.3** and **Bioconductor
3.18**. Use the same versions on all workers so the results are reproducible.

## 1. Install prerequisites

Install these before cloning the repository. In an Administrator PowerShell,
install Git, Git LFS, and Posit's R Installation Manager (`rig`) with WinGet:

```powershell
winget install --id Git.Git -e --source winget
winget install --id GitHub.GitLFS -e --source winget
winget install --id Posit.rig -e --source winget
```

Close and reopen PowerShell after installation so the updated `PATH` is loaded.
Confirm `rig --version` identifies Posit's R Installation Manager; another
program named `rig` may be on `PATH` and is not compatible with this guide.

Install the remaining version-specific prerequisites:

1. **R 4.3.3 for Windows (64-bit)**. The setup script can install it with
   Posit rig if it is not already present.
2. **Rtools43**, matching R 4.3. Download and install it from the
   [official Rtools43 page](https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html).
   Some locked packages may need to be compiled from source during restore.
   The setup script checks for its `sh.exe` prerequisite before starting
   `renv::restore()`.

If the lab has restricted internet access, allow R and Git to reach CRAN,
Bioconductor, Posit Package Manager, and GitHub during setup. Package restore
downloads the versions recorded in `renv.lock`.

## 2. Clone the repository and retrieve its data

In PowerShell, choose a local working folder and run:

```powershell
git lfs install
git clone https://github.com/TonyRoyze/rnaseq-sim.git
cd rnaseq-sim
git lfs pull
```

Check that `data/DGElist_BRCA.rds` exists and is a real data file, rather than
a small Git LFS pointer text file.

## 3. Run the automated prerequisite and environment check

From PowerShell in the repository folder, run:

```powershell
powershell -ExecutionPolicy Bypass -File .\setup_windows.ps1
```

The script checks that all tracked Git LFS objects are present and all tracked
data files are downloaded rather than pointer files. It then checks the
reference RDS, searches the standard Windows install locations for R 4.3.3,
restores `renv.lock`, verifies the simulation packages, and previews the
simulation plan without generating datasets. If R 4.3.3 is not available, it
uses WinGet to install `rig`, then installs R 4.3.3. If `rig list` shows the R
version with an architecture suffix, pass that exact label, for example:

```powershell
powershell -ExecutionPolicy Bypass -File .\setup_windows.ps1 -RVersion 4.3-x86_64
```

If setup reports that `sh.exe` is missing, install Rtools43, restart
PowerShell, and rerun the setup script. If it is already installed in a
nonstandard location, add its `usr\bin` directory to `PATH` first. For
example, the standard location is `C:\rtools43\usr\bin`.

## 4. Create the project-local R environment

From the repository root, run:

```powershell
Rscript.exe -e "source('renv/activate.R'); renv::restore(prompt = FALSE)"
```

This creates the project-local package library under `renv/library/` using
`renv.lock`. The library is local to this checkout and is intentionally not
stored in Git. Keep the computer connected to the internet while packages are
being downloaded and installed. To repair or refresh the environment later,
run the same command again.

Confirm that key packages are available:

```powershell
Rscript.exe -e "source('renv/activate.R'); stopifnot(requireNamespace('edgeR', quietly=TRUE), requireNamespace('DESeq2', quietly=TRUE)); cat('Simulation packages are available.\n')"
```

If restore fails with `'sh' is not recognized`, install Rtools43, restart
PowerShell, and rerun `setup_windows.ps1` (or the restore command). This means
R tried to install a package from source but could not find Rtools' shell. If
restore later reports a compiler error, confirm that the matching Rtools43
toolchain is installed. If restore reports that R is not 4.3.3, install/select
R 4.3.3 and run the commands again.

## 5. Preview and run a shard

Preview the jobs assigned to worker 1 of 4:

```powershell
Rscript.exe scripts/custom_simulation/run_batch.R --plan-only --shard-index 1 --shard-count 4
```

Then generate that shard's datasets:

```powershell
Rscript.exe scripts/custom_simulation/run_batch.R --shard-index 1 --shard-count 4 --output "Z:/rnaseq-sim/results/custom_simulation"
```

Replace `Z:/rnaseq-sim/results/custom_simulation` with the lab's shared
storage path. The shared folder must be mounted at the same drive/path on all
workers. Run the same command on each computer, changing the shard index to
`2`, `3`, or `4`. Keep the shard count, reference data, `--n-var`, replicate
count, and base seed identical across workers; do not run the same shard on
two computers at once.

The default run has 40 settings and generates 10,000 genes per dataset. For a
quick single-dataset setup check, use one shard out of 40 and reduce the gene
count:

```powershell
Rscript.exe scripts/custom_simulation/run_batch.R --shard-index 1 --shard-count 40 --n-var 500 --output "Z:/rnaseq-sim/results/setup-check"
```

Use that reduced run only to verify the installation. It is not a substitute
for the full benchmark settings.

## Notes

- The custom simulation scripts use ordinary sequential R execution and are
  designed to run from Windows. Some older scripts elsewhere in the repository
  use `parallel::mclapply` and require a Unix-like system.
- Windows drive paths such as `Z:/...` and UNC paths such as
  `//server/share/...` are accepted for `--output` and `--reference`.
- Generated datasets and manifests are resumable. A completed dataset is
  skipped when its shard is run again. Keep outputs on reliable shared storage
  with enough free space; avoid synchronizing active result directories via
  consumer sync tools while multiple workers are writing.
- For all runner options, setting definitions, and method analysis, see
  [`scripts/custom_simulation/README.md`](scripts/custom_simulation/README.md).

## Underdispersion generators

Generalized Poisson is the default (`--underdispersion-model generalized_poisson`).
Its normalized finite-support sampler uses the base R `stats` package, requiring
no extra package installation. The optional COM-Poisson generator is selected
with `--underdispersion-model com_poisson`. `COMPoissonReg` 0.8.2 and its
Rcpp/numDeriv dependencies are in `renv.lock`; `setup_windows.ps1` restores and
checks the package. Manual setup on any platform:

```r
source("renv/activate.R")
renv::restore(prompt = FALSE)
stopifnot(requireNamespace("COMPoissonReg", quietly = TRUE))
```

Pass the same `--underdispersion-model` flag to all workers and to
`check_batch.R` when checking completion. The model name is saved in dataset
IDs, so results from the two underdispersion generators are distinguishable.
