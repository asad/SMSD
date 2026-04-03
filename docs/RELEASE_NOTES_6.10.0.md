# SMSD 6.10.0 Release Notes

## Summary

SMSD `6.10.0` improves Python bindings, strengthens native aromaticity
perception, and makes ring-sensitive graph matching more reliable across
native, Java, and Python paths.

## Highlights

- Added explicit aromaticity controls in Python and native APIs.
- Improved aromaticity perception for Kekule and fused-ring systems.
- Added `kekulize()` and `dearomatize()` support for `MolGraph`.
- Harmonized aromaticity behaviour across Java, Python, and native paths using
  a shared parity test corpus.
- Improved MCS timeout handling and self-match behavior.
- Included housekeeping, bug fixes, and internal optimization work.

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
