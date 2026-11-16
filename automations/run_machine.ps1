# Arguments are forwarded unchanged, e.g. --config config/machine.dcf --plan-only.
$ErrorActionPreference = 'Stop'
$ProjectRoot = Split-Path -Parent $PSScriptRoot
Push-Location -LiteralPath $ProjectRoot
try {
  & Rscript 'scripts/simulation_pipeline/run_machine.R' @args
  $RunExitCode = $LASTEXITCODE
} finally {
  Pop-Location
}
exit $RunExitCode
