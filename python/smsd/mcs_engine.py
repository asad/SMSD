"""
SPDX-License-Identifier: Apache-2.0
Copyright (c) 2018-2026 BioInception PVT LTD
Algorithm Copyright (c) 2009-2026 Syed Asad Rahman

Thin Python wrapper over the native SMSD MCS engine.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


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
    """Convert SMILES, MolGraph, or RDKit Mol to raw C++ MolGraph."""
    if isinstance(mol, str):
        import smsd._smsd as _smsd
        return _smsd.parse_smiles(mol)
    # Already a C++ MolGraph
    if hasattr(mol, 'atomic_num') and hasattr(mol, 'bond_order'):
        return mol
    # RDKit Mol
    if HAS_RDKIT and hasattr(mol, 'GetNumAtoms'):
        smi = Chem.MolToSmiles(mol)
        import smsd._smsd as _smsd
        return _smsd.parse_smiles(smi)
    raise TypeError(f"Cannot convert {type(mol)} to MolGraph")


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

    g1 = _ensure_molgraph(mol1)
    g2 = _ensure_molgraph(mol2)

    from smsd._smsd import find_mcs_coverage

    mapping = find_mcs_coverage(
        g1, g2,
        ring_match=ring_matches_ring,
        bond_any=bond_any,
        timeout_ms=int(timeout * 1000),
    )

    # Convert 0-based C++ mapping to 1-based pairs for backward compat
    pairs = [(k + 1, v + 1) for k, v in mapping.items()]
    elapsed = (time.monotonic() - t0) * 1000

    return LightMCSResult(
        size=len(pairs),
        mapping=pairs,
        candidates=[pairs] if pairs else [],
        elapsed_ms=elapsed,
    )
