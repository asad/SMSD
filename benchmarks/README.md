# SMSD local benchmarks

The current measurement protocol and release results are in
[RESULTS_7.2.0.md](RESULTS_7.2.0.md). Performance depends on chemistry options,
corpus, search budget and mapping validity. The maintained comparison records
latency, returned atoms, RDKit cancellation and independent mapping checks
separately. Atom counts alone do not establish equivalent results or optimality.

## 7.2.1 release scope

Version 7.2.1 reorganizes the repository into `java/`, `cpp/` and `python/`
modules and updates release packaging. No new measurements are attributed to
those changes. `RESULTS_7.2.0.md`, its measured numbers, source fingerprints,
input hashes and the `smsd-7.2.0-benchmark-data.tar.gz` archive retain their
original names and scope. Commands below use the current checkout layout;
new 7.2.1 release checks are tracked in
[validation](../docs/VALIDATION_7.2.1.md).

## Maintained measurements

| Entry point | Scope |
|---|---|
| `benchmark_python.py` | 20 curated MCS pairs and 20 substructure cases; in-process native SMSD and RDKit |
| `run_external_benchmarks.py` | All 1,000 random pairs, 1,000 neighbor pairs, 468 Tautobase pairs, 12 stress pairs and 1,400 SMARTS patterns |
| `benchmark_bindings.py` | Parsing, RDKit conversion/cache, graph-property copies, scalar calls, dispatch and single-thread batch overhead |
| `benchmark_cpp.sh` / `benchmark_cpp.cpp` | Actual SMSD C++ headers; six default MCS pairs or a supplied pair TSV |
| `benchmark_assignment_matching.cpp` | Assignment and bipartite-matching primitives |
| `benchmark_mcs_components.cpp` | Connected-mapping filtering |
| `benchmark_substructure_setup.cpp` | Substructure setup and screening |
| `benchmark_search_primitives.cpp` | Search helper primitives |
| `benchmark_tautomer.cpp`, `benchmark_tautomer_zinc.java` | Tautomer feature diagnostics on the curated molecule pool; the historical Java filename does not establish ZINC provenance |
| `benchmark_substructure_java.java` | In-process SMSD/CDK substructure diagnostics on 28 pairs |
| `benchmark_1000_java.java` | In-process Java MCS on generated pairs from the curated molecule pool |
| Python `test_external_benchmarks.py`; Java opt-in benchmark tests | Feature diagnostics; see the report for overlap and comparison limits |

The Dalke-style files are derived from MoleculeNet, rather than the original
ChEMBL-13 experiments. The checked-in neighbor corpus contains self-pairs,
duplicates and neighbors below Tanimoto 0.7. See [data/README.md](data/README.md)
for provenance, counts and the corrected generator. No fixed similarity cutoff
is implied by the word “neighbor.” SMARTS matching on a single ibuprofen target
does not reproduce the full Ehrlich–Rarey experiment or validate entire dialects.

## Run the native comparison

Use an isolated environment containing the intended SMSD wheel and RDKit.
Every JSON result records imported module paths and versions; verify these
before interpreting results. For CPU measurements, build with Metal and CUDA
disabled and use the same compiler configuration for baseline and candidate.

```bash
python -m build --wheel \
  -Ccmake.define.SMSD_BUILD_METAL=OFF \
  -Ccmake.define.SMSD_BUILD_CUDA=OFF
python benchmarks/benchmark_python.py \
  --mcs-only --compare-mode any --timeout-sec 10 --warmup 1 --iters 3 \
  --output build/local-benchmarks/core.tsv
python benchmarks/run_external_benchmarks.py \
  --policies strict fmcs --timeout-sec 1 --warmup 0 --repeats 1 \
  --output-dir build/local-benchmarks/external
python benchmarks/run_external_benchmarks.py \
  --datasets stress dalke-random dalke-nn tautobase \
  --policies any --timeout-sec 1 --warmup 0 --repeats 1 \
  --output-dir build/local-benchmarks/any
python benchmarks/run_external_benchmarks.py \
  --datasets stress dalke-random dalke-nn tautobase \
  --policies fmcs --smsd-strategy auto --timeout-sec 1 --warmup 0 --repeats 1 \
  --output-dir build/local-benchmarks/auto
python benchmarks/benchmark_bindings.py \
  --output build/local-benchmarks/bindings.json
```

The external runner checkpoints each completed pair and each engine call. Use
`--resume` with identical arguments and environment metadata to continue a run.
The defaults are 10 seconds, one warmup and three timed trials; the full release
corpus run explicitly uses one second, no warmup and one trial to measure
bounded-budget quality. These are different protocols and are labelled in the
report. Timings recorded alongside other builds/tests do not support speed claims.

All MCS comparison modes request a connected, non-induced mapping and maximize
atoms. Parsing, graph conversion, warmup and witness validation are outside the
search timer. Engine order alternates by pair and trial.

| Mode | Interpretation |
|---|---|
| `any` | Shared element matching and any bond order; ratio eligibility still requires valid mappings, equal atom counts and uncanceled RDKit results |
| `strict` | Exact bonds, charge, atom/bond ring parity and aromatic atoms; RDKit requires a Python atom comparator, so its timing includes callback overhead and is excluded from speed ratios |
| `fmcs` / `defaults` | SMSD strict bond order with flexible aromaticity versus RDKit `CompareOrder`; aromatic bond semantics differ, so timing is descriptive |

SMSD preserves every query edge between mapped atoms. RDKit FMCS can return a
common bond subgraph that omits some of those edges. The runner seeks a compatible
query-vertex witness among the first 128 embeddings per molecule. “No witness”
is a bounded-search result, not proof that all FMCS embeddings are invalid.
SMSD's returned-mapping API exposes no cancellation status; the report records
its elapsed-budget crossings and keeps cancellation unknown.
Calls within one percent of the requested budget are also excluded from speed
ratios because millisecond budget rounding can finish just below the outer timer.

## Integration and historical entry points

`benchmark_leaderboard.py` aggregates core or Java CLI/CDK outputs;
`benchmark_all.py` measures Java CLI startup, Python and RDKit; `benchmark_java.sh`
and `benchmark_rdkit.py` cover 15 historical pairs. Their timings have different
setup costs or chemistry contracts and do not establish cross-engine rankings.
`benchmark_1000.py` includes Java CLI startup on generated pool pairs.
`benchmark_python_vs_rdkit.py` retains the older ten-pair set and delegates to the
maintained protocol. `benchmark_smsd_v6.py` retains its historical filename and
measures public API, fingerprint and batch features.

`benchmark_present_vs_light_scaling.py` requires a separately built sibling
native extension. Without independent witness validation and matched chemistry,
its historical size/timing ratios are diagnostic only. `compare_benchmark_tsv.py`
compares two saved TSVs; use the associated JSON validity data before making
performance claims. Existing `results_*` files are historical archives, not the
current release evidence. [BENCHMARK_REPORT.md](BENCHMARK_REPORT.md) points to the
current report and does not carry forward older unsupported speed rankings.

Java requires JDK 25 and a built shaded JAR. C++ requires a C++17 compiler; the
standalone shell harness compiles the actual headers with optimization and
neither configures GPU backends nor substitutes a different algorithm.

```bash
mvn -f java/pom.xml package -DskipTests
SMSD_JAR=java/target/smsd-7.2.1-jar-with-dependencies.jar \
  NUM_RUNS=3 TIMEOUT_MS=10000 bash benchmarks/benchmark_java.sh
bash benchmarks/benchmark_cpp.sh
SMSD_BENCHMARK=1 SMSD_BENCHMARK_TIMEOUT_MS=1000 \
  python -m pytest --import-mode=importlib python/tests/test_external_benchmarks.py -v -s
mvn -f java/pom.xml test -Dslow.tests.exclude=nothing -Dbenchmark=true -Dsmsd.benchmark=true \
  '-Dtest=BenchmarkSuiteTest*,ExternalBenchmarkTest*,JavaCdkVsSmsdBenchmarkTest' \
  -Dsmsd.benchmark.timeoutMs=1000 -Dsmsd.benchmark.rounds=1 \
  -Dsmsd.benchmark.warmup=0 -Dsmsd.benchmark.outputDir="$PWD/build/local-benchmarks/java"
```

No hosted runner is required for these commands.
