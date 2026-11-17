# Read-only progress monitor. Arguments are forwarded to the base-R script.
$ErrorActionPreference = 'Stop'
Push-Location -LiteralPath (Split-Path -Parent $PSScriptRoot)
try {
    & Rscript --vanilla 'scripts/simulation_pipeline/progress.R' @args
    $MonitorExitCode = $LASTEXITCODE
} finally {
    Pop-Location
}
exit $MonitorExitCode
