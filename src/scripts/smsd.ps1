# SMSD — Small Molecule Subgraph Detector
# https://github.com/asad/SMSD
# Apache License 2.0

$base = Split-Path -Parent $PSScriptRoot

# Find the fat JAR
$jar = Get-ChildItem "$base/target/smsd-*-jar-with-dependencies.jar" -ErrorAction SilentlyContinue | Select-Object -First 1
if (-not $jar) {
    $jar = Get-ChildItem "$base/smsd-*-jar-with-dependencies.jar" -ErrorAction SilentlyContinue | Select-Object -First 1
}
if (-not $jar) {
    Write-Error "SMSD JAR not found. Run 'mvn package' first."
    exit 1
}

& java $env:JAVA_OPTS -jar $jar.FullName @args
