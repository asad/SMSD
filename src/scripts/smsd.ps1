$base = Split-Path -Parent $PSScriptRoot
$jar  = Join-Path $base "target/${project.build.finalName}-jar-with-dependencies.jar"
& java -jar $jar @args

