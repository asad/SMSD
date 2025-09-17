@echo off
setlocal
set BASE=%~dp0..
set JAR=%BASE%\target\${project.build.finalName}-jar-with-dependencies.jar
java -jar "%JAR%" %*

