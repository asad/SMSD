Of course, here is the text formatted as Markdown.

# SMSD (CDK-compat) — Pure Java CLI & API

**A high-performance, pure-Java command-line tool and library for substructure and maximum common substructure (MCS) searching in chemical graphs, built on the Chemistry Development Kit (CDK).**

This tool provides a modern, fast, and dependency-light solution for common cheminformatics tasks. It avoids reliance on external toolkits and prioritizes a simple, powerful command-line interface with structured JSON output.

-----

## Key Features

  * **Substructure Search**: Identify if a query molecule exists as a substructure within a target molecule.
      * Supports SMARTS queries for complex pattern matching.
      * Enumerate all unique mappings.
  * **Maximum Common Substructure (MCS)**: Find the largest shared substructure between two molecules.
      * Tunable algorithms for connected or disconnected MCS.
      * Options for induced vs. non-induced subgraphs.
  * **Rich Input/Output**:
      * Read molecules from SMILES, SMARTS, MOL, SDF, and more.
      * Produce machine-readable JSON output for easy integration into pipelines.
      * Pretty-printing for human-readable inspection.
  * **High Performance**:
      * Based on advanced VF2++ graph matching algorithms.
      * Intelligent pruning and heuristics for rapid searching.
      * Configurable timeouts to prevent runaway calculations.
  * **Modern CLI**:
      * Powered by Picocli for a user-friendly experience.
      * Clear, organized commands (`sub`, `mcs`).
      * Detailed help messages and validation.
  * **Minimal Dependencies**: Optimized to use only essential CDK modules, keeping the footprint small.

-----

## Build

To build the project and create an executable "fat JAR", run:

```bash
mvn -U clean package
```

This will produce `target/smsd-cdk-3.0.0-jar-with-dependencies.jar`.

-----

## Command-Line Usage

The tool is designed to be run as a self-contained JAR.

### General Syntax

```bash
java -jar target/smsd-cdk-*-jar-with-dependencies.jar [COMMAND] [OPTIONS]
```

### Commands

  * `sub`: Perform a substructure search.
  * `mcs`: Find the maximum common substructure.

### Common Options

  * `-q, --query <STRING>`: Query molecule (SMILES, SMARTS, or file path).
  * `-Q, --query-type <TYPE>`: Type of query: `SMI`, `SMARTS`, `MOL`, `SDF`.
  * `-t, --target <STRING>`: Target molecule (SMILES or file path).
  * `-T, --target-type <TYPE>`: Type of target: `SMI`, `MOL`, `SDF`.
  * `--json <FILE>`: Write JSON output to a file. Use `-` for standard output.
  * `--pretty`: Pretty-print the JSON output.
  * `--timeout <ms>`: Set a timeout in milliseconds (default: 10000).

-----

### 1\. Substructure Search (`sub`)

#### Check for existence (simple)

Returns `{"exists":true}` or `{"exists":false}`.

```bash
java -jar target/smsd-cdk-*-jar-with-dependencies.jar sub \
  -q "CCN" -Q SMI \
  -t "CCCNC" -T SMI \
  --json -
```

#### Find all mappings

Use the `--all-mappings` flag to enumerate every unique way the query fits into the target.

```bash
java -jar target/smsd-cdk-*-jar-with-dependencies.jar sub \
  -q "c1ccccc1" -Q SMI \
  -t "c1ccccc1C(=O)O" -T SMI \
  --all-mappings \
  --json - --pretty
```

**Example Output:**

```json
{
  "query": "c1ccccc1",
  "target": "c1ccccc1C(=O)O",
  "mode": "substructure",
  "mappings": [
    {
      "query_atoms": [0, 1, 2, 3, 4, 5],
      "target_atoms": [0, 1, 2, 3, 4, 5],
      "pairs": [[0,0], [1,1], [2,2], [3,3], [4,4], [5,5]]
    }
  ]
}
```

#### Using SMARTS Query

Use the `SMARTS` query type for advanced matching, like finding any primary amine.

```bash
java -jar target/smsd-cdk-*-jar-with-dependencies.jar sub \
  -q "[NH2;!$(N=O)]" -Q SMARTS \
  -t "CCN" -T SMI \
  --json -
```

-----

### 2\. Maximum Common Substructure (`mcs`)

#### Find the MCS between two molecules

```bash
java -jar target/smsd-cdk-*-jar-with-dependencies.jar mcs \
  -q "O=C(C)Oc1ccccc1C(=O)O" -Q SMI \
  -t "CC(=O)Nc1ccc(O)cc1" -T SMI \
  --json - --pretty
```

**Example Output:**

```json
{
  "query": "O=C(C)Oc1ccccc1C(=O)O",
  "target": "CC(=O)Nc1ccc(O)cc1",
  "mode": "mcs",
  "induced": false,
  "connected": true,
  "mcs": {
    "query_atoms": [0, 1, 2, 3, 4, 5, 6, 7],
    "target_atoms": [0, 1, 2, 3, 4, 5, 6, 7],
    "pairs": [[0,2], [1,1], [2,0], [3,3], [4,4], [5,5], [6,6], [7,7]],
    "smiles": "CC(=O)Oc1ccccc1"
  }
}
```

#### MCS Options

  * `--induced`: Require that if a bond exists in the target MCS, it must also exist in the query MCS (stricter).
  * `--connected`: Return only the largest connected component of the MCS.

-----

## Repository Layout

  * `src/main/java/...` — Core Java source code.
  * `src/test/java/...` — JUnit tests.
  * `pom.xml` — Maven build configuration.