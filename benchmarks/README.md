# SMSD benchmarks

Run these tools from the repository root. Use the same molecules, chemistry
options and search budgets when comparing results, and check returned mappings.
The [8.0.0 validation record](../docs/VALIDATION_8.0.0.md) covers the release
checks; the measurement protocol below governs every reported result.

## Prepared-query APIs

Build the source before running these workloads. The three runners read
`search_api_cases.tsv` and
write JSON Lines containing every measurement, parsed graph, atom and selected
bond mapping, fragments, stopping status and proof flags.

| Phase | Work inside the timer |
| --- | --- |
| `prepare` | Construct a prepared query from a parsed graph, or compile SMARTS |
| `prepared` | Search a pre-parsed target with a prepared query |
| `end_to_end` | Parse both inputs, prepare the query and search |
| `batch` | Search a collection of pre-parsed targets sequentially |

The corpus includes MCIS and MCES, connected and disconnected results,
weighted MCIS, molecule and SMARTS substructure, capped enumeration and batch
calls. Chemistry, budgets, limits and batch sizes are explicit in each row.
The `any` profile compares elements and ignores bond order; other profiles
exercise strict chemistry, rings, stereo and tautomer matching.

From the repository root, run the Python and native workloads:

```bash
mkdir -p build/search-benchmarks
python benchmarks/benchmark_search_api.py \
  --output build/search-benchmarks/python.jsonl --warmup 2 --repeats 5
c++ -std=c++17 -O3 -DNDEBUG -Icpp/include \
  benchmarks/benchmark_search_api.cpp -o build/search-benchmarks/native
build/search-benchmarks/native --cases benchmarks/search_api_cases.tsv \
  --output build/search-benchmarks/cpp.jsonl --warmup 2 --repeats 5
```

For Java, compile the runner against the current source-built dependency JAR:

```bash
mvn -f java/pom.xml package -DskipTests
javac --release 8 -cp java/cli/target/smsd-cli-8.0.0-jar-with-dependencies.jar \
  -d build/search-benchmarks benchmarks/BenchmarkSearchAPI.java
java -cp "build/search-benchmarks:java/cli/target/smsd-cli-8.0.0-jar-with-dependencies.jar" \
  BenchmarkSearchAPI --cases benchmarks/search_api_cases.tsv \
  --output build/search-benchmarks/java.jsonl --warmup 2 --repeats 5
```

Use a semicolon instead of a colon in the Java classpath on Windows. With
Java 8, replace `--release 8` with `-source 8 -target 8`. All runners accept
`--case ID` and `--phase PHASE` to select a workload.

Validate the recorded results before analysing timings:

```bash
python benchmarks/search_api_validation.py --cases benchmarks/search_api_cases.tsv \
  --results build/search-benchmarks/cpp.jsonl build/search-benchmarks/java.jsonl \
    build/search-benchmarks/python.jsonl --output build/search-benchmarks/validation.json
```

The independent validator checks returned mappings, selected components and
completed cross-language results. Small eligible graphs also use exhaustive
objective checks. It reports predicates and graph sizes outside its independent
oracle scope explicitly; language regression tests cover those contracts.
Raw embeddings are compared after accounting for parser atom order.

Keep preparation and warm search results separate. Use several fresh processes
in a balanced order, record compiler, runtime, source and power conditions, and
retain every observation. Timing ratios require the same objective and matching
rules with completed proofs. A timeout with a valid incumbent is a separate
quality outcome. Run coverage and builds outside the timing window.

## Python

Install SMSD and RDKit in the same environment. These short runs check the
benchmark tools; they do not measure the full corpus:

```bash
python benchmarks/benchmark_python.py \
  --mcs-only --compare-mode any --timeout-sec 1 --warmup 0 --iters 1 \
  --output build/local-benchmarks/core.tsv
python benchmarks/run_external_benchmarks.py \
  --datasets stress dalke-random dalke-nn tautobase \
  --policies any --timeout-sec 1 --warmup 0 --repeats 1 --limit 1 \
  --output-dir build/local-benchmarks/external-mcs-smoke
python benchmarks/benchmark_bindings.py \
  --iterations 10 --repeats 1 --warmup 0 \
  --output build/local-benchmarks/bindings-smoke.json
```

Remove `--limit 1` for the full external corpus. Use `--resume` with the same
arguments and environment to continue an interrupted run. Results record the
imported package paths and versions.

| Tool | Measures |
|---|---|
| `benchmark_python.py` | The 20 curated MCS pairs in `data/mcs_pairs.tsv` and 20 substructure cases |
| `run_external_benchmarks.py` | Random and neighbouring pairs, tautomer pairs, stress cases and SMARTS |
| `benchmark_bindings.py` | Parsing, RDKit conversion, scalar calls and batch overhead |

The external MCS runner checks mappings and records completion and stopping
status. Stopped searches are excluded from speed ratios.

| Policy | Matching rules |
|---|---|
| `any` | Matching elements, with any bond order |
| `strict` | Exact bonds, charge, ring membership and aromatic atoms; RDKit uses a Python comparator |
| `fmcs` | SMSD strict bond order with flexible aromaticity; RDKit uses `CompareOrder` |

Different MCS edge and aromaticity rules can produce different results. The
mapping checks of the external MCS runner, above, are in `mcs_protocol.py`, and
the comparison limits under [Prepared-query APIs](#prepared-query-apis) apply.
The [data guide](data/README.md) records corpus sources and counts. The derived
Dalke-style pairs are separate from the original Dalke benchmark.

## Java and C++

Build the Java CLI with Maven and a JDK. Java 25 LTS is preferred for the
build, and the JAR runs on Java 8 or later. C++ benchmarks require a C++17
compiler.

```bash
mvn -f java/pom.xml package -DskipTests
python benchmarks/benchmark_leaderboard.py --mode integration
bash benchmarks/benchmark_cpp.sh benchmarks/data/stress_pairs.tsv 1000 1 0
mkdir -p build/local-benchmarks
c++ -std=c++17 -O3 -Icpp/include benchmarks/benchmark_primitives.cpp \
  -o build/local-benchmarks/primitives
build/local-benchmarks/primitives assignment
```

The integration driver measures the Java CLI, including its start-up time,
against RDKit on the curated MCS pairs, and Java substructure matching against
RDKit and CDK on `substructure_pairs.tsv`. The C++ shell benchmark compiles
the SMSD headers and measures the supplied stress pairs. The primitive
benchmark times native assignment and general matching (`assignment`, with an
optional repeat multiplier), NLF3 construction and a complete-subgraph MCS
(`primitives`) and cold-target matcher set-up (`setup`); each selector prints
checksums that must agree between builds. These timings have different
set-up costs from the Python benchmarks.

For opt-in Java corpus tests, see the [build guide](../docs/HOWTO-INSTALL.md).
