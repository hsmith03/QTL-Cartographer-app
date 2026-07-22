[CmdletBinding()]
param(
    [ValidateSet('Release', 'Debug')]
    [string] $Configuration = 'Release'
)

$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot
$vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
if (-not (Test-Path -LiteralPath $vswhere)) {
    throw 'Visual Studio Build Tools were not found. Install the Desktop development with C++ workload.'
}
$vsInstall = & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
if (-not $vsInstall) { throw 'The Visual Studio C++ x64 toolchain was not found.' }

$vcvars = Join-Path $vsInstall 'VC\Auxiliary\Build\vcvars64.bat'
$environment = & cmd.exe /d /s /c "`"$vcvars`" >nul && set"
if ($LASTEXITCODE -ne 0) { throw 'Visual Studio could not initialize its x64 build environment.' }
foreach ($line in $environment) {
    if ($line -match '^([^=]+)=(.*)$') { Set-Item -Path "Env:$($matches[1])" -Value $matches[2] }
}

$common = @('params.c', 'Utilities.c', 'Mdatain.c', 'Genome.c', 'Qdatain.c', 'Idatain.c', 'NumRec.c')
$programs = [ordered]@{
    Rmap = @('RMmain.c', 'RMfunc.c'); Rqtl = @('RQmain.c', 'RQfunc.c')
    Eqtl = @('EQmain.c', 'EQfunc.c', 'MIfunc.c', 'MissMark.c'); Rcross = @('RCmain.c', 'RCfunc.c')
    Emap = @('Emain.c', 'Efunc.c', 'MissMark.c', 'Blas.c'); Preplot = @('Preplot.c', 'EQfunc.c')
    Prune = @('Prune.c'); Qstats = @('QSmain.c', 'QSfunc.c', 'MissMark.c')
    LRmapqtl = @('LRmain.c', 'LRfunc.c', 'Blas.c', 'Linpak.c', 'MLnpkws.c', 'Otraits.c', 'MissMark.c', 'QSfunc.c')
    BTmapqtl = @('BTL_LRmain.c', 'BTL_LRfunc.c', 'LRfunc.c', 'Blas.c', 'Linpak.c', 'MLnpkws.c', 'Otraits.c', 'MissMark.c', 'QSfunc.c')
    Zmapqtl = @('Zmain.c', 'Zfunc.c', 'QSfunc.c', 'Otraits.c', 'MissMark.c', 'Linpak.c', 'MLnpkws.c', 'Blas.c')
    SRmapqtl = @('SRmain.c', 'SRfunc.c', 'Otraits.c', 'MissMark.c', 'Linpak.c', 'MLnpkws.c', 'Blas.c')
    JZmapqtl = @('MZmain.c', 'MZfunc.c', 'Otraits.c', 'MissMark.c', 'Linpak.c', 'MLnpkws.c', 'Blas.c')
    MImapqtl = @('MImain.c', 'MIfunc.c', 'Otraits.c', 'MissMark.c', 'MLnpkws.c')
    MultiRegress = @('MRmain.c', 'MRfunc.c', 'MLnpkws.c', 'Linpak.c', 'Blas.c')
    Bmapqtl = @('Bmain.c', 'QTL_MOVE.c', 'birth.c', 'mcmc.c', 'update_parms.c', 'INITVALS.c', 'accept_birth.c', 'death.c', 'read_data.c', 'update_qtl.c', 'LAPACK.c', 'f.c', 'ranlib.c', 'Otraits.c', 'MissMark.c', 'RCfunc.c')
    QTLcart = @('QTLcart.c')
}

$buildRoot = Join-Path $root 'build'
$toolsDir = Join-Path $buildRoot "app\tools"
$objRoot = Join-Path $buildRoot "obj\$Configuration"
$appDir = Join-Path $buildRoot 'app'
New-Item -ItemType Directory -Force -Path $toolsDir, $objRoot, $appDir | Out-Null
$optimization = if ($Configuration -eq 'Release') { '/O2' } else { '/Od', '/Zi' }
$definitions = @('/DUNIX', '/DWIN32', '/DITOA', '/D_CRT_SECURE_NO_WARNINGS', '/DSIXTYFOURBIT')

Push-Location $root
try {
    foreach ($entry in $programs.GetEnumerator()) {
        $name = $entry.Key
        $relativeObjDir = "build\obj\$Configuration\$name"
        New-Item -ItemType Directory -Force -Path (Join-Path $root $relativeObjDir) | Out-Null
        $sources = @($entry.Value + $common | ForEach-Object { Join-Path 'native\src' $_ })
        $arguments = @('/nologo', '/TC', '/MT', '/W3', '/wd4244', '/wd4267', '/wd4996', $optimization, $definitions,
            '/I', 'native\src', "/Fo:$relativeObjDir\", "/Fe:build\app\tools\$name.exe") + $sources
        Write-Host "Building native tool $name.exe"
        $ErrorActionPreference = 'Continue'
        & cl.exe @arguments
        $ErrorActionPreference = 'Stop'
        if ($LASTEXITCODE -ne 0) { throw "Compilation failed for $name.exe." }
    }

    $csc = Join-Path $vsInstall 'MSBuild\Current\Bin\Roslyn\csc.exe'
    if (-not (Test-Path -LiteralPath $csc)) { throw 'The Visual Studio C# compiler was not found.' }
    $guiSources = Get-ChildItem (Join-Path $root 'app\*.cs') | ForEach-Object FullName
    $guiArgs = @('/nologo', '/target:winexe', '/platform:anycpu', '/warn:4',
        '/reference:System.dll', '/reference:System.Core.dll', '/reference:System.Drawing.dll', '/reference:System.Windows.Forms.dll',
        '/out:build\app\QTL-Cartographer.exe') + $guiSources
    if ($Configuration -eq 'Release') { $guiArgs += '/optimize+' } else { $guiArgs += @('/debug+', '/optimize-') }
    Write-Host 'Building QTL-Cartographer.exe'
    & $csc @guiArgs
    if ($LASTEXITCODE -ne 0) { throw 'GUI compilation failed.' }
}
finally { Pop-Location }

Copy-Item -Force (Join-Path $root 'app\QTL-Cartographer.exe.config') $appDir
Copy-Item -Force (Join-Path $root 'LICENSE'), (Join-Path $root 'README.md') $appDir
Copy-Item -Recurse -Force (Join-Path $root 'example'), (Join-Path $root 'doc') $appDir

$zip = Join-Path $buildRoot 'QTL-Cartographer-Windows-x64.zip'
Compress-Archive -Force -Path (Join-Path $appDir '*') -DestinationPath $zip
Write-Host "Application: $(Join-Path $appDir 'QTL-Cartographer.exe')"
Write-Host "Package: $zip"
