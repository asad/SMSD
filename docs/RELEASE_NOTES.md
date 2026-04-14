# SMSD Pro 7.1.1

Bug-fix patch on top of v7.1.0.  No new features, no public API breakage.

## Fixed

- Cross-language ECFP / FCFP fingerprint parity between Java, C++, and
  Python (two long-standing Java drifts at radius ≥ 1).
- Java canonical SMILES writer: bond symbol for aromatic-adjacent
  single bonds and implicit H count inside stereo brackets.
- Python `smsd.canonical_smiles(smi)` / `smsd.to_smiles(smi)` raised
  `TypeError` on string input.  Both now accept `str` or `MolGraph`.
- `MatchResult.overlapCoefficient` returned the wrong similarity
  metric.  Both `MatchResult.overlap` and `MatchResult.overlapCoefficient`
  now return Szymkiewicz-Simpson overlap as documented; the new
  `MatchResult.tanimoto` attribute exposes the Jaccard value.
- Canonical SMILES writer now emits `[nH]` for pyrrole-type aromatic
  nitrogen, so output kekulizes cleanly in downstream readers.
- FP-level `smsd.overlapCoefficient` / `count_overlapCoefficient`
  camelCase aliases now return Simpson overlap as documented.

## Verified

- Python pytest: 603 passed, 6 skipped, 0 failures.
- Java JUnit: 581 tests, 0 failures, 0 errors.

Apache 2.0 — see `NOTICE`.
