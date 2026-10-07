# SMSD benchmarks

Run these tools from the repository root. Use the same molecules, chemistry
options and search budgets when comparing results, and check returned mappings.
The [7.2.0 report](RESULTS_7.2.0.md) contains the recorded measurements;
[7.2.2 validation](../docs/VALIDATION_7.2.2.md) covers the current release tests.

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
| `benchmark_python.py` | 20 curated MCS pairs and 20 substructure cases |
| `run_external_benchmarks.py` | Random and neighbouring pairs, tautomer pairs, stress cases and SMARTS |
| `benchmark_bindings.py` | Parsing, RDKit conversion, scalar calls and batch overhead |

The external MCS runner checks mappings and records RDKit cancellation. SMSD's
mapping API does not expose cancellation status. Results near the search budget
are excluded from speed ratios.

| Policy | Matching rules |
|---|---|
| `any` | Matching elements, with any bond order |
| `strict` | Exact bonds, charge, ring membership and aromatic atoms; RDKit uses a Python comparator |
| `fmcs` | SMSD strict bond order with flexible aromaticity; RDKit uses `CompareOrder` |

Different MCS edge and aromaticity rules can produce different results. The
[report](RESULTS_7.2.0.md) explains the validity checks and comparison limits.
The [data guide](data/README.md) records corpus sources and counts. The derived
Dalke-style pairs are separate from the original Dalke benchmark.

## Java and C++

Build the Java CLI with Maven and a JDK. Java 25 LTS is preferred; the 7.2.2 JAR
runs on Java 8 or later. C++ benchmarks require a C++17 compiler.

```bash
mvn -f java/pom.xml package -DskipTests
SMSD_JAR=java/target/smsd-7.2.2-jar-with-dependencies.jar \
  NUM_RUNS=1 TIMEOUT_MS=1000 bash benchmarks/benchmark_java.sh
bash benchmarks/benchmark_cpp.sh benchmarks/data/stress_pairs.tsv 1000 1 0
```

The Java shell benchmark includes CLI startup. The C++ shell benchmark compiles
the SMSD headers and measures the supplied stress pairs. These timings have different
setup costs from the Python benchmarks.

For opt-in Java corpus tests, see the [build guide](../docs/HOWTO-INSTALL.md).
Other `benchmark_*.py`, Java and C++ files provide feature diagnostics or
historical workloads. Existing `results_*` files retain their original scope.
