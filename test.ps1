[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot
$app = Join-Path $root 'build\app\QTL-Cartographer.exe'
if (-not (Test-Path -LiteralPath $app)) { throw 'Build the application before testing it.' }

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
Write-Host "All tests passed. Results: $run"
