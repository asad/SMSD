# SMSD Pro 6.10.2 Release Notes

## Summary

SMSD Pro `6.10.2` is a correctness release that fixes the MCS connectivity filter
for non-induced mode and adds regression tests for challenging molecule pairs.

## Highlights

- **MCS connectivity filter corrected.** The connected-component filter now
  enforces reachability through bonds present in *both* query and target
  molecules. Previously, in non-induced mode, a bond existing only in the
  query could incorrectly count toward connectivity, inflating the reported
  MCS size with disconnected atoms. Affected pairs include aspirin vs
  acetaminophen, glucose vs fructose, and other cases where the two molecules
  share atoms but not all intervening bonds.

- **GOLDEN_843 regression tests.** Added regression tests in both Python and
  Java for the GOLDEN_843 molecule pair (a large reaction-style pair that
  previously triggered timeouts). Tests verify both completion within budget
  and minimum MCS size.

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
