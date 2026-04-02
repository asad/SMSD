# SMSD 6.9.0 Release Notes

## Summary

SMSD `6.9.0` focuses on chemistry correctness, native I/O coverage, and
clearer benchmark parity across the Java, C++, and Python implementations.

## Highlights

- Native/public MCS handling is more direction-stable on hard asymmetric pairs.
- `ringMatchesRingOnly=true` now enforces symmetric ring/non-ring parity in
  the native core, Python binding, and Java implementation.
- The local benchmark leaderboard now supports explicit comparison modes:
  `defaults`, `strict`, and `fmcs`.
- Native MOL V2000/V3000 handling, metadata round-trip, and patent-style
  `R#` support remain part of the `6.9.0` release.

## Availability and Use

SMSD is distributed under the Apache License 2.0. Commercial use, internal
deployment, and redistribution are allowed under that license.

If you redistribute SMSD or derivative works, keep the applicable license and
notice materials with the distribution, including [LICENSE](../LICENSE) and
[NOTICE](../NOTICE), as required.

Trademark use is covered by the [NOTICE](../NOTICE) file.

## Attribution

If you use SMSD in academic work, please cite:

Rahman SA, Bashton M, Holliday GL, Schrader R, Thornton JM.
Small Molecule Subgraph Detector (SMSD) toolkit.
Journal of Cheminformatics, 1:12, 2009.
DOI: 10.1186/1758-2946-1-12
