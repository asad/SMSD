@echo off
setlocal
set "SMSD_JAR_PATH=%SMSD_JAR%"
if defined SMSD_JAR_PATH goto run
for /f "delims=" %%J in ('dir /b /a-d /o-d "%~dp0..\..\target\smsd-*-jar-with-dependencies.jar" 2^>nul') do (
  set "SMSD_JAR_PATH=%~dp0..\..\target\%%J"
  goto run
)
for /f "delims=" %%J in ('dir /b /a-d /o-d "%~dp0smsd-*-jar-with-dependencies.jar" 2^>nul') do (
  set "SMSD_JAR_PATH=%~dp0%%J"
  goto run
)
echo SMSD JAR not found. Run mvn package first or set SMSD_JAR. 1>&2
exit /b 1
:run
if defined JAVA_HOME (
  "%JAVA_HOME%\bin\java.exe" %JAVA_OPTS% -jar "%SMSD_JAR_PATH%" %*
) else (
  java %JAVA_OPTS% -jar "%SMSD_JAR_PATH%" %*
)
exit /b %ERRORLEVEL%
