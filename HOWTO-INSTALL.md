# How to use this bundle

1) Place **WHITEPAPER.md** and **README.md** at the repository root.
2) Copy `src/test/java/com/bioinception/smsd/SMSDPortedCasesTest.java` into your tree,
   keeping the package `com.bioinception.smsd` (adjust if your API package differs).
3) Build and run tests:
   ```bash
   mvn -U clean test package
   ```
4) Run the CLI (jar-first):
   ```bash
   java -jar target/smsd-cdk-3.0.0-jar-with-dependencies.jar sub  --Q SMI --q "CCN" --T SMI --t "CCCNC" -m --json - --json-pretty
   ```

## Notes
- The tests exercise substructure and MCIS across 100+ cases, including recursive SMARTS.
- Output is plain (no icons). Format output as per `OUTPUT_STYLE_NOTES.md`.
- To apply SPDX/Apache headers automatically, add the `license-maven-plugin` from `POM_SNIPPETS.md` and run:
  ```bash
  mvn license:format
  ```
