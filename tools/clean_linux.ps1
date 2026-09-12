param(
    [ValidateSet("Debug", "Release")]
    [string]$Configuration = "Release"
)

$ErrorActionPreference = "Stop"

# Keep WSL discovery and invocation consistent with build_linux.ps1.
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

$probe = & $wslExe -e sh -lc "printf MMCAL_WSL_READY" 2>$null
$probeExit = $LASTEXITCODE
if ($probeExit -ne 0 -or (($probe -join "") -notmatch "MMCAL_WSL_READY")) {
    throw "wsl.exe exists, but the default Linux distribution could not be started."
}

$root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$wslRootRaw = & $wslExe -e wslpath -a -- $root 2>$null
if ($LASTEXITCODE -ne 0) {
    throw "wslpath failed while converting project root: $root"
}

$wslRoot = (($wslRootRaw | Select-Object -First 1) -as [string])
if ([string]::IsNullOrWhiteSpace($wslRoot)) {
    throw "wslpath returned an empty project path."
}
$wslRoot = $wslRoot.Trim()

# A Windows drive letter can be convertible by wslpath even when the drive is
# not actually mounted/visible in the active WSL distribution.
& $wslExe -e test -d $wslRoot
if ($LASTEXITCODE -ne 0) {
    throw @"
The project directory is not visible from WSL:

    Windows: $root
    WSL    : $wslRoot

Verify that the drive containing the source tree is mounted in Ubuntu.
For this project, L: would normally be visible below /mnt/l.
"@
}

$buildDir = "$wslRoot/build/linux/$Configuration"
Write-Host "[mmCal.Linux] Cleaning: $buildDir"

& $wslExe -e rm -rf -- $buildDir
if ($LASTEXITCODE -ne 0) {
    throw "Linux clean failed with exit code $LASTEXITCODE."
}

Write-Host "[mmCal.Linux] CLEAN SUCCESS"
