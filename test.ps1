[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot
$app = Join-Path $root 'build\app\QTL-Cartographer.exe'
if (-not (Test-Path -LiteralPath $app)) { throw 'Build the application before testing it.' }

Write-Host 'Testing package metadata and checksum...'
$productInfo = Get-Content -Raw (Join-Path $root 'app\ProductInfo.cs')
if ($productInfo -notmatch 'public const string Version = "([^"]+)"') { throw 'Application version metadata is invalid.' }
$version = $Matches[1]
$zip = Join-Path $root "build\QTL-Cartographer-v$version-Windows-x64.zip"
$checksumFile = $zip + '.sha256'
if (-not (Test-Path -LiteralPath $zip) -or -not (Test-Path -LiteralPath $checksumFile)) { throw 'Versioned package or SHA-256 file is missing.' }
$expectedChecksum = ((Get-Content -Raw $checksumFile).Trim() -split '\s+')[0]
$actualChecksum = (Get-FileHash -Algorithm SHA256 -LiteralPath $zip).Hash.ToLowerInvariant()
if ($expectedChecksum -ne $actualChecksum) { throw 'Package SHA-256 does not match the distributable ZIP.' }
$runtimeConfig = Get-Content -Raw (Join-Path $root 'build\app\QTL-Cartographer.exe.config')
if ($runtimeConfig -notmatch 'Version=v4\.8') { throw 'The packaged application does not require the supported .NET Framework 4.8 runtime.' }

Write-Host 'Testing GUI discovery of every analysis program...'
$guiTest = Start-Process -FilePath $app -ArgumentList '--smoke-test' -Wait -PassThru
if ($guiTest.ExitCode -ne 0) { throw "GUI smoke test failed with exit code $($guiTest.ExitCode)." }

$run = Join-Path $root ("build\tests\" + (Get-Date -Format 'yyyyMMdd-HHmmss'))
New-Item -ItemType Directory -Force -Path $run | Out-Null
Copy-Item -Force (Join-Path $root 'example\sample.*') $run
$tools = Join-Path $root 'build\app\tools'
$steps = @(
    @('Rmap', '-i', 'sample.mps', '-A'), @('Rcross', '-i', 'sample.raw', '-A'),
    @('Qstats', '-A'), @('LRmapqtl', '-A'), @('SRmapqtl', '-A'),
    @('Zmapqtl', '-A'), @('MImapqtl', '-A'), @('Eqtl', '-A'), @('Preplot', '-A')
)
Push-Location $run
try {
    foreach ($step in $steps) {
        $name = $step[0]
        $arguments = $step[1..($step.Count - 1)]
        & (Join-Path $tools "$name.exe") @arguments *> "$name.log"
        if ($LASTEXITCODE -ne 0) { throw "$name failed with exit code $LASTEXITCODE." }
    }
}
finally { Pop-Location }

$expected = 'qtlcart.map', 'qtlcart.cro', 'qtlcart.qst', 'qtlcart.lr', 'qtlcart.sr', 'qtlcart.z', 'qtlcart.mim', 'qtlcart.eqt', 'qtlcart.plt'
foreach ($name in $expected) {
    $path = Join-Path $run $name
    if (-not (Test-Path -LiteralPath $path) -or (Get-Item $path).Length -eq 0) { throw "Missing output: $name" }
}

Write-Host 'Testing results parsing, exports, import validation, and project persistence...'
$featureTest = Start-Process -FilePath $app -ArgumentList @('--feature-smoke-test', "`"$run`"") -Wait -PassThru
if ($featureTest.ExitCode -ne 0) { throw "Feature smoke test failed with exit code $($featureTest.ExitCode)." }

Write-Host 'Testing five-chromosome map continuity and project stems...'
$fiveChrom = Join-Path $root ("build\tests\nuzhdin-" + (Get-Date -Format 'yyyyMMdd-HHmmss'))
New-Item -ItemType Directory -Force -Path $fiveChrom | Out-Null
Copy-Item -Force (Join-Path $root 'example\nuzhdinm.inp'), (Join-Path $root 'example\nuzhdinc.inp') $fiveChrom
Push-Location $fiveChrom
try {
    & (Join-Path $tools 'Rmap.exe') -i nuzhdinm.inp -X nuzhdin -A -V
    if ($LASTEXITCODE -ne 0) { throw "Five-chromosome Rmap failed with exit code $LASTEXITCODE." }
    & (Join-Path $tools 'Rcross.exe') -i nuzhdinc.inp -X nuzhdin -A -V
    if ($LASTEXITCODE -ne 0) { throw "Five-chromosome Rcross failed with exit code $LASTEXITCODE." }
}
finally { Pop-Location }

$mapText = Get-Content -Raw (Join-Path $fiveChrom 'nuzhdin.map')
if ($mapText -notmatch '(?m)^-c\s+5\s') { throw 'Rcross did not retain all five chromosomes in the project map.' }
$resourceText = Get-Content -Raw (Join-Path $fiveChrom 'qtlcart.rc')
if ($resourceText -notmatch '(?m)^-stem\s+nuzhdin\s') { throw 'The native resource file did not retain the nuzhdin stem.' }
if ($resourceText -notmatch '(?m)^-chrom\s+5\s') { throw 'The native resource file did not retain the chromosome count.' }
if (-not (Test-Path -LiteralPath (Join-Path $fiveChrom 'nuzhdin.cro'))) { throw 'Five-chromosome cross output was not created.' }
Write-Host "All tests passed. Results: $run"
