"""
SPDX-License-Identifier: Apache-2.0
Copyright (c) 2018-2026 BioInception PVT LTD
Algorithm Copyright (c) 2009-2026 Syed Asad Rahman

Thin Python wrapper over the native SMSD MCS engine.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field



# ---------------------------------------------------------------------------
# Result dataclass (public API — backward compatible)
# ---------------------------------------------------------------------------

@dataclass
class LightMCSResult:
    """Result from the lightweight MCS wrapper."""
    size: int = 0
    mapping: list = field(default_factory=list)
    candidates: list = field(default_factory=list)
    elapsed_ms: float = 0.0


# ---------------------------------------------------------------------------
# Molecule conversion helpers
# ---------------------------------------------------------------------------

def _ensure_molgraph(mol):
    """Convert inputs through the shared parser and RDKit conversion cache."""
    from smsd import _ensure_mol
    return _ensure_mol(mol)


# ---------------------------------------------------------------------------
# Main entry point — delegates to the native engine
# ---------------------------------------------------------------------------

def find_mcs_lightweight(
    mol1, mol2, *,
    timeout: float = 1.0,
    ring_matches_ring: bool = False,
    bond_any: bool = False,
) -> LightMCSResult:
    """Compute a maximum common subgraph between two molecules.

    Accepts SMILES strings, SMSD MolGraph objects, or RDKit Mol objects.
    All matching is delegated to the native engine.

    Args:
        mol1: First molecule (SMILES, MolGraph, or RDKit Mol).
        mol2: Second molecule.
        timeout: Wall-clock timeout in seconds (default 1.0).
        ring_matches_ring: Ring atom only matches ring atom.
        bond_any: If True, any bond matches any bond.

    Returns:
        LightMCSResult with size, mapping, candidates, elapsed_ms.
    """
    t0 = time.monotonic()

    from smsd import _ensure_mol_ex, _auto_translate
    g1, rdkit1 = _ensure_mol_ex(mol1)
    g2, rdkit2 = _ensure_mol_ex(mol2)

    from smsd._smsd import find_mcs_coverage

    mapping = find_mcs_coverage(
        g1, g2,
        ring_match=ring_matches_ring,
        bond_any=bond_any,
        timeout_ms=int(timeout * 1000),
    )

    mapping = _auto_translate(mapping, g1, g2, rdkit1, rdkit2)[0]
    # Historical wrapper returns one-based indices in each input's atom order.
    pairs = [(k + 1, v + 1) for k, v in mapping.items()]
    elapsed = (time.monotonic() - t0) * 1000

    return LightMCSResult(
        size=len(pairs),
        mapping=pairs,
        candidates=[pairs] if pairs else [],
        elapsed_ms=elapsed,
    )
