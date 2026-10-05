# Benchmark input data

These checked-in inputs are repository benchmark corpora. Their presence does
not establish that they reproduce a published paper's original data or protocol.
Current measured results, input hashes, parse failures, matching policies and
quality checks are in [the 7.2.0 benchmark report](../RESULTS_7.2.0.md).

## Checked-in records

Counts below are non-comment records, excluding CSV/TSV headers, before molecule
or SMARTS parsing. A row count is not a count of valid completed comparisons.

| File | Records | Provenance recorded in this repository |
|---|---:|---|
| `stress_pairs.tsv` | 12 | Curated cage, ring, chain and peptide examples; not identified as an original published dataset |
| `chodera_tautobase_subset.txt` | 468 | Tautomer-pair subset labeled Chodera/Tautobase; exact extraction metadata is not retained |
| `tautobase_smirks.txt` | 1,680 | Tautobase-labeled transforms and source-reference fields |
| `dalke_random_pairs.tsv` | 1,000 | Derived pairs from the MoleculeNet collection below, not the original Dalke/Hastings ChEMBL-13 corpus |
| `dalke_nn_pairs.tsv` | 1,000 | Legacy derived neighbor pairs from the same collection; see corpus limitations below |
| `ehrlich_rarey_smarts.txt` | 1,400 | Deduplicated SMARTS; its header records Ehrlich/Rarey attribution, Hamburg source URL and license |
| `chembl_mcs_benchmark.smi` | 5,590 | MoleculeNet-derived pool; the filename does not establish ChEMBL provenance |
| `BBBP.csv`, `bace.csv`, `clintox.csv`, `sider.csv` | — | Source collections used for the derived pool |
| `bzr.sdf` | 163 | Legacy SDF collection; upstream extraction metadata is not recorded here |
| `cdk2.sdf` | 47 | Legacy SDF collection; upstream extraction metadata is not recorded here |

The checked-in neighbor corpus has 53 pairs sharing the same source ID,
39 duplicate unordered source-ID pairs, and 444 rows with recorded Tanimoto
below 0.7 (minimum 0.125). It must not be described as 1,000 distinct,
self-excluded, high-similarity nearest-neighbor pairs. Report results for this
fixed input separately from newly generated pairs.

## Generate derived random and neighbor pairs

```sh
python benchmarks/generate_dalke_pairs.py \
  --input benchmarks/data/chembl_mcs_benchmark.smi \
  --output-dir build/local-benchmarks/generated-pairs \
  --seed 42 --pairs 1000
```

The maintained generator requires RDKit, records its version and the source
SHA-256, and excludes the query index from neighbor selection. It uses Morgan
radius-2/1,024-bit Tanimoto and applies no similarity cutoff. Generated files
are written outside the checked-in data directory by default. They can differ
from the legacy files above and require their own hashes and row counts.

## Run maintained benchmarks

The shared MCS protocol uses atom maximization, a connected mapping, and a
10-second per-call budget. Timing comparisons also record validity, equal
mapping quality and RDKit cancellation; SMSD's mapping API does not expose
cancellation status. Different chemical policies and tautomer feature runs are
reported separately.

```sh
python benchmarks/run_external_benchmarks.py \
  --timeout-sec 10 --warmup 1 --repeats 3 \
  --output-dir build/local-benchmarks

mvn -Dtest=ExternalBenchmarkTest -Dbenchmark=true test
mvn -Dtest=JavaCdkVsSmsdBenchmarkTest -Dbenchmark=true test
```

Java uses the CDK version pinned in `java/pom.xml` (2.13 for this checkout). The
Java CDK comparison measures substructure search rather than MCS. See the
[benchmark guide](../README.md) for the remaining entry points and the
[current report](../RESULTS_7.2.0.md) for coverage and observed limitations.
