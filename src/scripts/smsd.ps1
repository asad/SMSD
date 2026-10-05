# SMSD — Small Molecule Subgraph Detector
# https://github.com/asad/SMSD
# Apache License 2.0

$base = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)

# Find the fat JAR
$jar = Get-ChildItem "$base/target/smsd-*-jar-with-dependencies.jar" -ErrorAction SilentlyContinue | Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $jar) {
    $jar = Get-ChildItem "$PSScriptRoot/smsd-*-jar-with-dependencies.jar" -ErrorAction SilentlyContinue | Sort-Object LastWriteTime -Descending | Select-Object -First 1
}
if (-not $jar) {
    Write-Error "SMSD JAR not found. Run 'mvn package' first."
    exit 1
}

$javaOptions = @()
if ($env:JAVA_OPTS) { $javaOptions = $env:JAVA_OPTS -split '\s+' }
& java @javaOptions -jar $jar.FullName @args
exit $LASTEXITCODE
