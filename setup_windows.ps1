[CmdletBinding()]
param(
  # Pass the exact version label shown by `rig list` if it includes an
  # architecture suffix, such as 4.3-x86_64.
  [string]$RVersion = "4.3.3"
)

$ErrorActionPreference = "Stop"

if ($env:OS -ne "Windows_NT") {
  throw "This setup script is for Windows PowerShell."
}

$ProjectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location -LiteralPath $ProjectRoot

function Require-Command {
  param([Parameter(Mandatory = $true)][string]$Name)
  $command = Get-Command $Name -ErrorAction SilentlyContinue
  if (-not $command) {
    throw "Required command '$Name' was not found. Install it and open a new PowerShell window."
  }
  return $command
}

$git = Require-Command "git"
& $git.Source lfs version
if ($LASTEXITCODE -ne 0) {
  throw "Git LFS is unavailable. Install Git LFS, then run 'git lfs install' and restart PowerShell."
}

Write-Host "Checking all Git LFS objects and tracked data files..."
& $git.Source lfs fsck --objects --pointers --dry-run
if ($LASTEXITCODE -ne 0) {
  throw "Git LFS found missing or invalid objects/pointers. Run 'git lfs pull' from the repository root, then rerun this script."
}
$lfsPaths = @(& $git.Source lfs ls-files --name-only)
if ($LASTEXITCODE -ne 0 -or $lfsPaths.Count -eq 0) {
  throw "Could not list Git LFS files in this repository. Verify Git LFS is installed and the repository is a Git clone."
}
foreach ($relativePath in $lfsPaths) {
  $trackedFile = Join-Path $ProjectRoot "$relativePath"
  if (-not (Test-Path -LiteralPath $trackedFile -PathType Leaf)) {
    throw "Git LFS data file is missing: $relativePath. Run 'git lfs pull' from the repository root."
  }
  $stream = [System.IO.File]::OpenRead($trackedFile)
  try {
    $headerBytes = New-Object byte[] 80
    $headerLength = $stream.Read($headerBytes, 0, $headerBytes.Length)
  } finally {
    $stream.Dispose()
  }
  $header = [System.Text.Encoding]::UTF8.GetString($headerBytes, 0, $headerLength)
  if ($header.StartsWith("version https://git-lfs.github.com/spec/v1")) {
    throw "Git LFS pointer found instead of downloaded data: $relativePath. Run 'git lfs pull' from the repository root."
  }
}
Write-Host "PASS: checked $($lfsPaths.Count) Git LFS data files"

$reference = Join-Path $ProjectRoot "data/DGElist_BRCA.rds"
if (-not (Test-Path -LiteralPath $reference -PathType Leaf)) {
  throw "Reference data is missing: $reference. Run 'git lfs pull' from the repository root."
}
if ((Get-Item -LiteralPath $reference).Length -lt 1MB) {
  throw "Reference data is unexpectedly small: $reference. It may still be a Git LFS pointer; run 'git lfs pull'."
}
$stream = [System.IO.File]::OpenRead($reference)
try {
  $headerBytes = New-Object byte[] 80
  $headerLength = $stream.Read($headerBytes, 0, $headerBytes.Length)
} finally {
  $stream.Dispose()
}
$header = [System.Text.Encoding]::UTF8.GetString($headerBytes, 0, $headerLength)
if ($header.StartsWith("version https://git-lfs.github.com/spec/v1")) {
  throw "Reference data is still an LFS pointer. Run 'git lfs pull' from the repository root."
}

function Find-InstalledRscript {
  $candidatePaths = @()
  $registryLocations = @(
    "Registry::HKEY_LOCAL_MACHINE\SOFTWARE\R-core\R64\4.3.3",
    "Registry::HKEY_CURRENT_USER\SOFTWARE\R-core\R64\4.3.3",
    "Registry::HKEY_LOCAL_MACHINE\SOFTWARE\R-core\R\4.3.3",
    "Registry::HKEY_CURRENT_USER\SOFTWARE\R-core\R\4.3.3"
  )
  foreach ($location in $registryLocations) {
    if (Test-Path -LiteralPath $location) {
      try {
        $key = Get-Item -LiteralPath $location
        $installPath = $key.GetValue("InstallPath")
        if (-not $installPath) { $installPath = $key.GetValue("") }
        if ($installPath) {
          $candidatePaths += (Join-Path $installPath "bin\Rscript.exe")
          $candidatePaths += (Join-Path $installPath "bin\x64\Rscript.exe")
        }
      } catch {
        # Keep searching other registry locations and standard install paths.
      }
    }
  }
  $candidatePaths += (Join-Path $env:ProgramFiles "R\R-4.3.3\bin\Rscript.exe")
  $candidatePaths += (Join-Path $env:ProgramFiles "R\R-4.3.3\bin\x64\Rscript.exe")
  foreach ($candidate in $candidatePaths) {
    if (Test-Path -LiteralPath $candidate -PathType Leaf) {
      return (Get-Item -LiteralPath $candidate).FullName
    }
  }
  return $null
}

function Install-RigViaWinget {
  $winget = Get-Command "winget.exe" -ErrorAction SilentlyContinue
  if (-not $winget) {
    throw "The selected R is not R 4.3.3 and rig is unavailable. Install rig with 'winget install --id posit.rig', then rerun this script."
  }
  Write-Host "Installing rig with WinGet so the project can use R 4.3.3..."
  & $winget.Source install --id posit.rig --exact --accept-package-agreements --accept-source-agreements
  if ($LASTEXITCODE -ne 0) {
    throw "WinGet could not install rig (exit code $LASTEXITCODE). Check network access and installation permissions, then rerun this script."
  }
  $machinePath = [Environment]::GetEnvironmentVariable("Path", "Machine")
  $userPath = [Environment]::GetEnvironmentVariable("Path", "User")
  $env:Path = "$machinePath;$userPath"
  $installedRig = Get-Command "rig.exe" -ErrorAction SilentlyContinue
  if ($installedRig) { return $installedRig.Source }
  $knownRigPath = Join-Path $env:USERPROFILE ".local\bin\rig.exe"
  if (Test-Path -LiteralPath $knownRigPath -PathType Leaf) { return $knownRigPath }
  throw "WinGet installed rig, but it is not available in this PowerShell session. Close PowerShell, open a new one, and rerun this script."
}

$rig = Get-Command "rig" -ErrorAction SilentlyContinue
if (-not $rig) {
  $knownRigPath = Join-Path $env:USERPROFILE ".local\bin\rig.exe"
  if (Test-Path -LiteralPath $knownRigPath -PathType Leaf) {
    $rig = Get-Item -LiteralPath $knownRigPath
  }
}
$rscript = Get-Command "Rscript.exe" -ErrorAction SilentlyContinue
$rscriptPath = if ($rscript) { $rscript.Source } else { Find-InstalledRscript }
$rigPath = if ($rig -is [System.IO.FileInfo]) { $rig.FullName } elseif ($rig) { $rig.Source } else { $null }

if (-not $rig -and -not $rscriptPath) {
  $rigPath = Install-RigViaWinget
  Write-Host "Installing project R $RVersion with rig..."
  & $rigPath add $RVersion
  if ($LASTEXITCODE -ne 0) {
    throw "rig could not install R $RVersion. If this computer has no administrator rights, run 'rig system user-mode', then rerun this script."
  }
}

function Invoke-ProjectR {
  param([Parameter(Mandatory = $true)][string[]]$RCommandArgs)

  if ($script:RigPath) {
    $invokeArgs = @("run", "-r", $script:RVersion, "--rscript") + $RCommandArgs
    & $script:RigPath @invokeArgs
  } else {
    & $script:RscriptPath @RCommandArgs
  }
  if ($LASTEXITCODE -ne 0) {
    throw "R command failed with exit code $LASTEXITCODE. See the error above."
  }
}

$script:RigPath = $rigPath
$script:RscriptPath = $rscriptPath
$script:RVersion = $RVersion

$versionExpression = "cat('R_VERSION=', as.character(getRversion()), '\n', sep='')"
try {
  $versionOutput = Invoke-ProjectR -RCommandArgs @("-e", $versionExpression)
} catch {
  if (-not $script:RigPath) { throw }
  Write-Host "R $RVersion is not available through rig. Installing it now..."
  & $script:RigPath add $RVersion
  if ($LASTEXITCODE -ne 0) {
    throw "rig could not install R $RVersion. If this computer has no administrator rights, run 'rig system user-mode', then rerun this script."
  }
  try {
    $versionOutput = Invoke-ProjectR -RCommandArgs @("-e", $versionExpression)
  } catch {
    throw "Could not start R version '$RVersion'. Run 'rig list' and retry with the installed 4.3.3 label, for example: .\setup_windows.ps1 -RVersion 4.3-x86_64"
  }
}
$versionLine = @($versionOutput) |
  Where-Object { "$_" -match "^R_VERSION=" } |
  Select-Object -Last 1
$actualRVersion = if ($versionLine) { ("$versionLine" -replace "^R_VERSION=", "").Trim() } else { "" }
if ($actualRVersion -ne "4.3.3") {
  if (-not $script:RigPath) {
    $script:RigPath = Install-RigViaWinget
    Write-Host "Installing project R $RVersion with rig..."
    & $script:RigPath add $RVersion
    if ($LASTEXITCODE -ne 0) {
      throw "rig could not install R $RVersion. If this computer has no administrator rights, run 'rig system user-mode', then rerun this script."
    }
    $versionOutput = Invoke-ProjectR -RCommandArgs @("-e", $versionExpression)
    $versionLine = @($versionOutput) |
      Where-Object { "$_" -match "^R_VERSION=" } |
      Select-Object -Last 1
    $actualRVersion = if ($versionLine) { ("$versionLine" -replace "^R_VERSION=", "").Trim() } else { "" }
  }
  if ($actualRVersion -ne "4.3.3") {
    throw "This project requires R 4.3.3, but the selected R is '$actualRVersion'. Use rig to install/select R 4.3.3. Do not change the R version recorded in renv.lock to match a newer R."
  }
}
Write-Host "PASS: selected R $actualRVersion"

if (-not (Test-Path -LiteralPath (Join-Path $ProjectRoot "renv/activate.R") -PathType Leaf) -or
    -not (Test-Path -LiteralPath (Join-Path $ProjectRoot "renv.lock") -PathType Leaf)) {
  throw "renv/activate.R or renv.lock is missing from the project."
}

Write-Host "Restoring the project library from renv.lock..."
$restoreExpression = "source('renv/activate.R'); renv::restore(prompt = FALSE)"
Invoke-ProjectR -RCommandArgs @("-e", $restoreExpression)

$packageExpression = @"
required <- c('edgeR', 'limma', 'DESeq2', 'DESeq', 'DSS', 'NBPSeq', 'NOISeq', 'ABSSeq', 'ROTS')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop('Missing required packages: ', paste(missing, collapse = ', '))
cat('PACKAGE_CHECK=OK\n')
"@
$packageOutput = Invoke-ProjectR -RCommandArgs @("-e", $packageExpression)
if (-not (@($packageOutput) | Where-Object { "$_" -match "^PACKAGE_CHECK=OK" })) {
  throw "Required package verification did not report success."
}
Write-Host "PASS: required simulation packages are available"

$checkOutput = Join-Path $env:TEMP "rnaseq-sim-setup-check"
$runnerArgs = @(
  "-f", "scripts/custom_simulation/run_batch.R",
  "--plan-only", "--n-var", "500",
  "--shard-index", "1", "--shard-count", "40",
  "--output", $checkOutput
)
Write-Host "Checking the simulation runner and reference data (plan only)..."
Invoke-ProjectR -RCommandArgs $runnerArgs

Write-Host "Setup checks passed. No datasets were generated by this setup check."
Write-Host "For the next step, see setup.md."
