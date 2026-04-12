# Random 200 Reaction Backend Comparison

## Sample Profile

| Metric | Value |
|---|---:|
| Sample size | 200 |
| Issue-free reactions | 40.00% |
| Error-free reactions | 100.00% |
| Warning-free reactions | 100.00% |
| Charge-balanced reactions | 94.50% |
| Non-H balanced reactions | 55.50% |
| Median reactant molecules | 2.0 |
| Median product molecules | 1.0 |

## Backend Summary

| Backend | Failures | Coverage | Exact | Mol->Mol | Atom->Atom | Bond-change | Mean abs bond delta | Median ms | P95 ms | Rxn/s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| smsd/smsd | 0 | 100.00% | 73.00% | 83.00% | 94.82% | 67.00% | 3.410 | 198.692 | 7042.495 | 0.52 |
| rdkit/rdkit | 0 | 100.00% | 72.50% | 82.50% | 94.64% | 66.50% | 3.420 | 20.163 | 228.967 | 15.84 |

## Pure Backend Focus

| Backend | Exact | Mol->Mol | Atom->Atom | Bond-change | Median ms | Rxn/s |
|---|---:|---:|---:|---:|---:|---:|
| smsd/smsd | 73.00% | 83.00% | 94.82% | 67.00% | 198.692 | 0.52 |
| rdkit/rdkit | 72.50% | 82.50% | 94.64% | 66.50% | 20.163 | 15.84 |

## Algorithm Distribution

### smsd/smsd

| Algorithm | Count |
|---|---:|
| RINGS_PAIRWISE | 164 |
| MIN_PAIRWISE | 18 |
| MAX_PAIRWISE | 11 |
| COMPONENT_CONSTRAINED_MCS | 3 |
| MIXTURE_PAIRWISE | 2 |
| MULTI_ANCHOR_BEAM | 1 |
| EXPLICIT_PAIRING | 1 |

### rdkit/rdkit

| Algorithm | Count |
|---|---:|
| RINGS_PAIRWISE | 163 |
| MIN_PAIRWISE | 20 |
| MAX_PAIRWISE | 11 |
| COMPONENT_CONSTRAINED_MCS | 2 |
| MIXTURE_PAIRWISE | 2 |
| MULTI_ANCHOR_BEAM | 1 |
| EXPLICIT_PAIRING | 1 |
