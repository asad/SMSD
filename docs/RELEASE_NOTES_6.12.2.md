# SMSD Pro 6.12.2

**Syed Asad Rahman — BioInception PVT LTD**

Performance, quality, and API refinement release.

---

## What's New

### Broader default matching behaviour

The default `ChemOptions` now applies relaxed ring and charge constraints,
making out-of-the-box MCS results applicable to a wider range of workflows
without any configuration. Named profiles are available for stricter settings:

| Profile | Use case |
|---|---|
| `default` | General-purpose, broad applicability |
| `strict` | Ring-constrained, charge-sensitive |
| `pharma` | Drug-discovery, complete-ring preference |
| `reaction` | Atom–atom mapping, relaxed constraints |

### Configurable `mcs()` API

All matching options are now available as direct keyword arguments:

```python
mapping = smsd.mcs(mol1, mol2,
    ring_matches_ring_only=False,
    match_bond_order="loose")
```

### Speed/quality trade-off control

A `max_stage` parameter controls the depth of the MCS search.
Lower values return fast approximate results; higher values
allow the engine to find the provably optimal solution.

### Improved similarity metrics

`overlap_coefficient` and `tanimoto_coefficient` are now separately
available on all MCS results, with corrected formula definitions.

### Java and Python parity

Java and Python APIs now expose the same set of options and return
equivalent results for equivalent inputs across all platforms.

### Security

The Docker runtime image now runs SMSD as a non-privileged user.

---

## Compatibility

- Java 25+, C++17, Python 3.10+
- GPU acceleration: Apple Silicon (Metal), NVIDIA (CUDA Volta+)
- Platforms: macOS (arm64, x86_64), Linux (x86_64, aarch64), Windows (AMD64)

---

## Copyright

Copyright (c) 2018-2026 BioInception PVT LTD
Algorithm copyright (c) 2009-2026 Syed Asad Rahman
Apache License 2.0 — see LICENSE and NOTICE for full terms.
