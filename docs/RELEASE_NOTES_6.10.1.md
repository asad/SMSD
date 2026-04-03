# SMSD 6.10.1 Release Notes

## Summary

SMSD `6.10.1` is a stability and correctness release that hardens the MCS
mapping repair pipeline, eliminates non-deterministic test behaviour, and
fixes CI/CD reliability issues.

## Highlights

- **MCS mapping repair rewritten.** The invalid-mapping recovery path is now
  iterative with bounded termination. Previously, deeply invalid seed
  mappings on large molecules (vancomycin, CoA, paclitaxel) could cause
  unbounded recursion. The new pipeline removes one offending atom per
  round, handles duplicate-target violations correctly, and degrades
  gracefully to substructure-based recovery when repair alone is insufficient.

- **Deterministic test suite.** All timing-dependent assertions have been
  removed from the test corpus. Algorithm correctness is now verified by
  structural invariants (MCS size, mapping validity, substructure membership),
  not wall-clock elapsed time. Tests are reproducible regardless of machine
  speed.

- **CI/CD fixes.** Removed `forkedProcessTimeoutInSeconds` from Maven
  Surefire (was killing the fork JVM prematurely). Switched Python publish
  workflow to manual dispatch only. Fixed GitHub Actions artifact version
  references. Added Node.js 24 opt-in for all workflows.

- **Test corrections.** Fixed `atomWeights` array length for benzene query.
  Relaxed MCS thresholds for drug-pair tests to match validated algorithmic
  output. Adjusted `completeRingsOnly` tests for fused-ring edge cases.

## Availability and Use

SMSD is distributed under the Apache License 2.0. Commercial use, internal
deployment, and redistribution are allowed under that license.

If you redistribute SMSD or derivative works, keep the applicable license and
notice materials with the distribution, including [LICENSE](../LICENSE) and
[NOTICE](../NOTICE), as required.

Trademark use is covered by the [NOTICE](../NOTICE) file.

## Attribution

If you use SMSD Pro in academic work, please cite:

Rahman SA.
SMSD Pro: Coverage-Driven, Tautomer-Aware Maximum Common Substructure Search.
ChemRxiv, 2025.
DOI: 10.26434/chemrxiv.15001534
https://doi.org/10.26434/chemrxiv.15001534/v1

For the original SMSD toolkit, please also cite:

Rahman SA, Bashton M, Holliday GL, Schrader R, Thornton JM.
Small Molecule Subgraph Detector (SMSD) toolkit.
Journal of Cheminformatics, 1:12, 2009.
DOI: 10.1186/1758-2946-1-12
