# How to Build and Run

## Requirements
- Java 11+ (JDK 17 recommended)
- Maven 3.9+

## Build

```bash
mvn -U clean package
```

This produces `target/smsd-6.10.0-jar-with-dependencies.jar` (fat JAR with all dependencies).

## Run Tests

```bash
mvn clean test
```

## Run the CLI

```bash
java -jar target/smsd-*-jar-with-dependencies.jar \
  --Q SMI --q "CCN" \
  --T SMI --t "CCCNC" \
  -m --json - --json-pretty
```

## Docker

```bash
docker build -t smsd .
docker run --rm smsd --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --json -
```

## Notes
- The test suite exercises substructure and MCS across 420 cases, including recursive SMARTS, adversarial edge cases, and large molecules.
- Cross-platform launcher scripts are generated at `src/scripts/smsd` and `src/scripts/smsd.bat`.
