# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
# SMSD standalone CLI launcher. JAVA_OPTS is split on whitespace.

$base = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)

# Find the fat JAR
$jar = Get-ChildItem "$base/target/smsd-cli-*-jar-with-dependencies.jar" -ErrorAction SilentlyContinue | Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $jar) {
    $jar = Get-ChildItem "$PSScriptRoot/smsd-cli-*-jar-with-dependencies.jar" -ErrorAction SilentlyContinue | Sort-Object LastWriteTime -Descending | Select-Object -First 1
}
if (-not $jar) {
    Write-Error "SMSD JAR not found. Run 'mvn package' first."
    exit 1
}

$javaOptions = @()
if ($env:JAVA_OPTS) { $javaOptions = $env:JAVA_OPTS -split '\s+' }
& java @javaOptions -jar $jar.FullName @args
exit $LASTEXITCODE
