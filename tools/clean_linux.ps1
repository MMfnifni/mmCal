param(
    [ValidateSet("Debug", "Release")]
    [string]$Configuration = "Release"
)

$ErrorActionPreference = "Stop"

if (-not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) {
    throw "wsl.exe was not found."
}

$root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$wslRoot = (& wsl.exe wslpath -a $root | Select-Object -First 1)
if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($wslRoot)) {
    throw "Failed to convert the solution directory to a WSL path."
}
$wslRoot = $wslRoot.Trim()
$buildDir = "$wslRoot/build/linux/$Configuration"

& wsl.exe rm -rf $buildDir
exit $LASTEXITCODE
