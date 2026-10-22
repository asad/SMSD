# SMSD benchmarks

This guide shows how to run the benchmarks and how to report their results.
Run every command from the repository root, after building the current source.
The [data guide](data/README.md) lists the input files, their sources and
licences.

## Prepared-query APIs

The C++, Java and Python runners read `search_api_cases.tsv`, whose rows set
the search, chemistry, budget, limits and batch size of each case. They write
JSON Lines holding every measurement, mapping, stopping status and proof flag.

| Phase | Work inside the timer |
| --- | --- |
| `prepare` | Construct a prepared query from a parsed graph, or compile SMARTS |
| `prepared` | Search a pre-parsed target with a prepared query |
| `end_to_end` | Parse both inputs, prepare the query and search |
| `batch` | Search a collection of pre-parsed targets sequentially |

Python and C++:

```bash
mkdir -p build/search-benchmarks
python benchmarks/benchmark_search_api.py \
  --output build/search-benchmarks/python.jsonl --warmup 2 --repeats 5
c++ -std=c++17 -O3 -DNDEBUG -Icpp/include \
  benchmarks/benchmark_search_api.cpp -o build/search-benchmarks/native
build/search-benchmarks/native --cases benchmarks/search_api_cases.tsv \
  --output build/search-benchmarks/cpp.jsonl --warmup 2 --repeats 5
```

Java, against the source-built dependency JAR:

```bash
mvn -f java/pom.xml package -DskipTests
javac --release 8 -cp java/cli/target/smsd-cli-8.0.1-jar-with-dependencies.jar \
  -d build/search-benchmarks benchmarks/BenchmarkSearchAPI.java
java -cp "build/search-benchmarks:java/cli/target/smsd-cli-8.0.1-jar-with-dependencies.jar" \
  BenchmarkSearchAPI --cases benchmarks/search_api_cases.tsv \
  --output build/search-benchmarks/java.jsonl --warmup 2 --repeats 5
```

On Windows, separate classpath entries with `;`. With a Java 8 compiler,
replace `--release 8` with `-source 8 -target 8`. Every runner accepts
`--case ID` and `--phase PHASE` to select a workload.

Validate the results before analysing any timing; the validator checks every
returned mapping and compares completed results across the languages:

```bash
python benchmarks/search_api_validation.py --cases benchmarks/search_api_cases.tsv \
  --results build/search-benchmarks/cpp.jsonl build/search-benchmarks/java.jsonl \
    build/search-benchmarks/python.jsonl --output build/search-benchmarks/validation.json
```

## Python and RDKit

Install SMSD and RDKit in the same environment. These short runs check the
scripts rather than measure the full inputs:

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

Remove `--limit 1` to run every external pair, and add `--resume` with the
same arguments and environment to continue an interrupted run.

| Script | Measures |
|---|---|
| `benchmark_python.py` | The curated MCS pairs in `data/mcs_pairs.tsv` and substructure cases |
| `run_external_benchmarks.py` | Random and nearest-neighbour pairs, tautomer pairs, stress cases and SMARTS |
| `benchmark_bindings.py` | Parsing, RDKit conversion, scalar calls and batch overhead |

| Policy | Matching rules |
|---|---|
| `any` | Matching elements, with any bond order |
| `strict` | Exact bonds, charge, ring membership and aromatic atoms |
| `fmcs` | Strict bond order with flexible aromaticity; RDKit uses `CompareOrder` |

## Java and C++

```bash
mvn -f java/pom.xml package -DskipTests
python benchmarks/benchmark_leaderboard.py --mode integration
bash benchmarks/benchmark_cpp.sh benchmarks/data/stress_pairs.tsv 1000 1 0
mkdir -p build/local-benchmarks
c++ -std=c++17 -O3 -Icpp/include benchmarks/benchmark_primitives.cpp \
  -o build/local-benchmarks/primitives
build/local-benchmarks/primitives assignment
```

The integration driver times the Java CLI, including its start-up, against
RDKit on the curated MCS pairs, and Java substructure search against RDKit and
CDK on `substructure_pairs.tsv`; its results are not engine rankings.
`benchmark_cpp.sh` builds and times C++ MCS search on the stress pairs.
`benchmark_primitives` takes the selector `assignment`, `primitives` or
`setup` and prints checksums that must agree between builds.

The Java benchmark tests are off by default. To run them:

```bash
mvn -f java/pom.xml -pl core test -Dslow.tests.exclude=nothing \
  '-Dtest=BenchmarkSuiteTest*,ExternalBenchmarkTest*,JavaCdkVsSmsdBenchmarkTest' \
  -Dbenchmark=true -Dsmsd.benchmark=true
```

The system properties `smsd.benchmark.timeoutMillis` (default `1000`),
`smsd.benchmark.rounds` (default `1`) and `smsd.benchmark.warmup` (default `0`)
set the run, and results go to `build/local-benchmarks/java` unless
`smsd.benchmark.outputDir` names another directory.

## Reporting results

- Validate every returned mapping before analysing timings.
- Report each phase separately, and keep preparation apart from warm searches.
- Use several fresh processes in a balanced order and retain every observation.
- Record the compiler, runtime, source revision and power conditions.
- Compare timings only for the same objective, chemistry and matching rules,
  and only between completed searches; report stopped and timed-out searches,
  even with a valid mapping, as separate outcomes.
- Keep builds and coverage runs outside the timing window.
