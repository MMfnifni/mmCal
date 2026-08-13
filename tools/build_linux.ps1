param(
    [ValidateSet("Debug", "Release")]
    [string]$Configuration = "Release"
)

$ErrorActionPreference = "Stop"

# ------------------------------------------------------------
# Resolve wsl.exe without relying on PATH.
# ------------------------------------------------------------
$wslCandidates = @(
    (Join-Path $env:WINDIR "System32\wsl.exe"),
    (Join-Path $env:WINDIR "Sysnative\wsl.exe")
)

$wslExe = $null
foreach ($candidate in $wslCandidates) {
    if (Test-Path $candidate) {
        $wslExe = $candidate
        break
    }
}

if (-not $wslExe) {
    throw "wsl.exe was not found in System32/Sysnative."
}

Write-Host "[mmCal.Linux] WSL executable : $wslExe"

# ------------------------------------------------------------
# Do NOT use 'wsl --status' here.
# Some WSL/console combinations return a non-zero code or emit text
# that becomes garbled when Visual Studio/MSBuild captures it.
#
# Instead, execute a tiny command in the default Linux distribution.
# ------------------------------------------------------------
$probe = & $wslExe -e sh -lc "printf MMCAL_WSL_READY" 2>$null
$probeExit = $LASTEXITCODE

if ($probeExit -ne 0 -or (($probe -join "") -notmatch "MMCAL_WSL_READY")) {
    throw @"
wsl.exe exists, but no usable default Linux distribution could be started.

Open PowerShell (preferably as Administrator) and run:

    wsl -l -v

If Ubuntu is not listed, install it with:

    wsl --install -d Ubuntu

After installation, start Ubuntu once and complete the initial user setup.
Then install the build tools inside Ubuntu:

    sudo apt update
    sudo apt install -y build-essential cmake ninja-build

Finally verify from PowerShell:

    wsl -e sh -lc "g++ --version"
    wsl -e sh -lc "cmake --version"
"@
}

Write-Host "[mmCal.Linux] WSL distribution: OK"

# ------------------------------------------------------------
# Project root
# ------------------------------------------------------------
$root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path

# Convert Windows path -> Linux path.
$wslRootRaw = & $wslExe -e wslpath -a -- $root 2>$null
if ($LASTEXITCODE -ne 0) {
    throw "wslpath failed while converting project root: $root"
}

$wslRoot = (($wslRootRaw | Select-Object -First 1) -as [string])
if ([string]::IsNullOrWhiteSpace($wslRoot)) {
    throw "wslpath returned an empty project path."
}
$wslRoot = $wslRoot.Trim()

# ------------------------------------------------------------
# Check build tools without printing their native output.
# ------------------------------------------------------------
& $wslExe -e sh -lc "command -v cmake >/dev/null 2>&1 && command -v g++ >/dev/null 2>&1"
if ($LASTEXITCODE -ne 0) {
    throw @"
Linux build tools are missing.

Start Ubuntu and run:

    sudo apt update
    sudo apt install -y build-essential cmake ninja-build

Then rebuild mmCal.Linux.
"@
}

# ------------------------------------------------------------
# Options
# ------------------------------------------------------------
$native = "OFF"
if ($env:MMCAL_NATIVE -match '^(1|ON|TRUE|YES)$') {
    $native = "ON"
}

$lto = if ($Configuration -eq "Release") { "ON" } else { "OFF" }
$buildDir = "$wslRoot/build/linux/$Configuration"

Write-Host "[mmCal.Linux] Configuration   : $Configuration"
Write-Host "[mmCal.Linux] Source          : $wslRoot"
Write-Host "[mmCal.Linux] Build dir       : $buildDir"
Write-Host "[mmCal.Linux] Native CPU      : $native"
Write-Host "[mmCal.Linux] LTO             : $lto"

# ------------------------------------------------------------
# Configure
# ------------------------------------------------------------
& $wslExe -e cmake `
    -S $wslRoot `
    -B $buildDir `
    "-DCMAKE_BUILD_TYPE=$Configuration" `
    "-DCMAKE_CXX_COMPILER=g++" `
    "-DMMCAL_ENABLE_LTO=$lto" `
    "-DMMCAL_NATIVE_OPTIMIZATION=$native" `
    "-DMMCAL_LINUX_OUT_SUFFIX=ON"

if ($LASTEXITCODE -ne 0) {
    throw "CMake configure failed with exit code $LASTEXITCODE."
}

# ------------------------------------------------------------
# Build
# ------------------------------------------------------------
& $wslExe -e cmake --build $buildDir --target mmCal --parallel
if ($LASTEXITCODE -ne 0) {
    throw "Linux build failed with exit code $LASTEXITCODE."
}

Write-Host "[mmCal.Linux] SUCCESS"
Write-Host "[mmCal.Linux] Output: build\linux\$Configuration\mmCal.out"
