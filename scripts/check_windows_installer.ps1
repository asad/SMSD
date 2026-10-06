# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
<# Install the MSI, check the installed console application, then remove it. #>
[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)][string]$ReleaseDir,
    [Parameter(Mandatory = $true)][string]$OutputDir,
    [string]$Python = 'python'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
if (![System.Runtime.InteropServices.RuntimeInformation]::IsOSPlatform([System.Runtime.InteropServices.OSPlatform]::Windows)) {
    throw 'MSI installation checks require native Windows'
}
$output = (Resolve-Path $OutputDir).Path
$release = (Resolve-Path $ReleaseDir).Path
$provenance = Get-Content (Join-Path $output 'installer-provenance.json') -Raw | ConvertFrom-Json
if ($provenance.platform -ne 'windows' -or $provenance.architecture -ne 'amd64') {
    throw 'Expected a Windows AMD64 installer'
}
$msi = Join-Path $output $provenance.installer
if ((Get-FileHash $msi -Algorithm SHA256).Hash.ToLowerInvariant() -ne $provenance.installer_sha256) {
    throw 'MSI checksum differs from its provenance'
}
$install = Join-Path $env:LOCALAPPDATA ('SMSD installer check ' + [guid]::NewGuid().ToString('N') + ' é')
if (Test-Path $install) { throw 'Installation check directory must not exist' }
$msiexec = Join-Path $env:SystemRoot 'System32\msiexec.exe'
$installed = $false
$removed = $false
$verification = Join-Path $output 'installed-image-check.json'
$installExitCode = $null
$removeExitCode = $null
$signature = Get-AuthenticodeSignature $msi
try {
    $arguments = "/i `"$msi`" /qn /norestart INSTALLDIR=`"$install`" /l*v `"$(Join-Path $output 'msi-install.log')`""
    $process = Start-Process $msiexec -ArgumentList $arguments -Wait -PassThru
    $installExitCode = $process.ExitCode
    if ($process.ExitCode -notin @(0, 3010)) { throw "MSI installation failed: $($process.ExitCode)" }
    $installed = $true
    if (!(Test-Path (Join-Path $install 'SMSD.exe'))) { throw 'The installed console launcher is missing' }
    & $Python (Join-Path $PSScriptRoot 'build_java_installer.py') verify-image `
        --image $install --release-dir $release --output-json $verification
    if ($LASTEXITCODE -ne 0) { throw 'Installed MSI application checks failed' }
}
finally {
    if ($installed) {
        $arguments = "/x `"$msi`" /qn /norestart /l*v `"$(Join-Path $output 'msi-remove.log')`""
        $process = Start-Process $msiexec -ArgumentList $arguments -Wait -PassThru
        $removeExitCode = $process.ExitCode
        if ($process.ExitCode -notin @(0, 3010)) { throw "MSI removal failed: $($process.ExitCode)" }
        if (Test-Path (Join-Path $install 'SMSD.exe')) { throw 'MSI removal left the launcher installed' }
        if (Test-Path (Join-Path $install 'app')) { throw 'MSI removal left the application payload installed' }
        if (Test-Path (Join-Path $install 'runtime')) { throw 'MSI removal left the runtime installed' }
        $removed = $true
        if (Test-Path $install) {
            $remaining = @(Get-ChildItem $install -Force)
            if ($remaining.Count -ne 0) { throw 'MSI removal left unexpected installed files' }
            Remove-Item $install
        }
    }
}
if (!$installed -or !$removed) { throw 'MSI installation and removal checks are incomplete' }
$image = Get-Content $verification -Raw | ConvertFrom-Json
if ($image.cli_jar_sha256 -ne $provenance.cli_jar_sha256 -or $image.version -ne $provenance.version) {
    throw 'Installed application differs from the packaged release'
}
$report = [ordered]@{
    version = $provenance.version
    platform = 'windows'
    architecture = 'amd64'
    installer = $provenance.installer
    installer_sha256 = $provenance.installer_sha256
    cli_jar_sha256 = $provenance.cli_jar_sha256
    status = 'passed'
    installation_method = 'native MSI install/run/remove'
    install_verified = $true
    cleanup_verified = $true
    install_exit_code = $installExitCode
    remove_exit_code = $removeExitCode
    authenticode_status = $signature.Status.ToString()
    execution_environment = 'native Windows AMD64'
    installed_image = $image
}
$report | ConvertTo-Json -Depth 20 | Set-Content (Join-Path $output 'installation-check.json') -Encoding utf8
$provenance.package_installation = 'passed: native MSI install/run/remove'
$provenance | ConvertTo-Json -Depth 20 | Set-Content (Join-Path $output 'installer-provenance.json') -Encoding utf8
Write-Host 'MSI installation, standalone chemical searches and removal passed.'
