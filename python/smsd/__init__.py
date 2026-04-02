# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""
SMSD -- Substructure & Maximum Common Substructure search for chemical graphs.

Copyright (c) 2018-2026 BioInception PVT LTD
Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
Licensed under Apache 2.0.

Quick start::

    from smsd import parse_smiles, find_mcs, is_substructure, similarity

    benzene = parse_smiles("c1ccccc1")
    phenol  = parse_smiles("c1ccc(O)cc1")

    assert is_substructure(benzene, phenol)

    mcs = find_mcs(benzene, phenol)
    print(f"MCS size: {len(mcs)}")

    sim = similarity(benzene, phenol)
    print(f"Similarity: {sim:.3f}")
"""

import importlib.machinery
import importlib.util
import sys
from pathlib import Path

__version__ = "6.9.0"
__author__ = "Syed Asad Rahman"


def _load_native_module():
    """Load the package-local native extension when available.

    This avoids accidentally binding the source-tree Python package to an older
    site-packages extension with a different exported API surface.
    """
    pkg_dir = Path(__file__).resolve().parent
    for suffix in importlib.machinery.EXTENSION_SUFFIXES:
        candidate = pkg_dir / f"_smsd{suffix}"
        if not candidate.exists():
            continue
        spec = importlib.util.spec_from_file_location("smsd._smsd", candidate)
        if spec is None or spec.loader is None:
            continue
        module = importlib.util.module_from_spec(spec)
        sys.modules["smsd._smsd"] = module
        spec.loader.exec_module(module)
        return module

    from smsd import _smsd as module  # fallback to regular installed import
    return module


_smsd = _load_native_module()


def parse_smiles(smiles):
    """Parse a SMILES string into a MolGraph.

    The Python API treats empty input as an empty molecule for convenience.
    Strict parser semantics remain unchanged in the native C++ layer.
    """
    if smiles == "":
        return MolGraph()
    return _parse_smiles_raw(smiles)

ChemOptions = _smsd.ChemOptions
McsOptions = _smsd.McsOptions
MolGraph = _smsd.MolGraph
MolGraphBuilder = _smsd.MolGraphBuilder

BondOrderMode = _smsd.BondOrderMode
AromaticityMode = _smsd.AromaticityMode
RingFusionMode = _smsd.RingFusionMode

_parse_smiles_raw = _smsd.parse_smiles
read_mol_block = _smsd.read_mol_block
to_smiles = _smsd.to_smiles
to_smarts = _smsd.to_smarts
write_mol_block = _smsd.write_mol_block
write_mol_block_v3000 = _smsd.write_mol_block_v3000
write_sdf_record = _smsd.write_sdf_record

is_substructure = _smsd.is_substructure
find_substructure = _smsd.find_substructure

find_mcs = _smsd.find_mcs
find_all_mcs = _smsd.find_all_mcs
translate_to_atom_ids = _smsd.translate_to_atom_ids
mcs_to_smiles = _smsd.mcs_to_smiles
find_mcs_smiles = _smsd.find_mcs_smiles

similarity_upper_bound = _smsd.similarity_upper_bound
screen_targets = _smsd.screen_targets

path_fingerprint = _smsd.path_fingerprint
mcs_fingerprint = _smsd.mcs_fingerprint
fingerprint_subset = _smsd.fingerprint_subset
analyze_fp_quality = _smsd.analyze_fp_quality
circular_fingerprint = _smsd.circular_fingerprint
circular_fingerprint_counts = _smsd.circular_fingerprint_counts
topological_torsion = _smsd.topological_torsion
topological_torsion_counts = _smsd.topological_torsion_counts

tanimoto = _smsd.tanimoto
dice = _smsd.dice
cosine = _smsd.cosine
soergel = _smsd.soergel
count_tanimoto = _smsd.count_tanimoto
count_dice = _smsd.count_dice
count_cosine = _smsd.count_cosine

_to_hex_cpp = _smsd.to_hex
_from_hex_cpp = _smsd.from_hex
_to_binary_string_cpp = _smsd.to_binary_string

gpu_is_available = _smsd.gpu_is_available
gpu_device_info = _smsd.gpu_device_info

RSLabel = _smsd.RSLabel
EZLabel = _smsd.EZLabel
CIPDescriptors = _smsd.CIPDescriptors
_assign_rs_raw = _smsd.assign_rs
_assign_ez_raw = _smsd.assign_ez
_assign_cip_raw = _smsd.assign_cip
_assign_cip_from_smiles_raw = _smsd.assign_cip_from_smiles

_compute_sssr_raw = _smsd.compute_sssr
_layout_sssr_raw = _smsd.layout_sssr

Point2D = _smsd.Point2D
_reduce_crossings_raw = _smsd.reduce_crossings
_force_directed_layout_raw = _smsd.force_directed_layout
_stress_majorisation_raw = _smsd.stress_majorisation
_match_template_raw = _smsd.match_template

# Reaction-aware MCS — lazy import (available when compiled with v6.5.0+)
try:
    from smsd._smsd import (
        reaction_aware_mcs as _reaction_aware_mcs_raw,
        map_reaction_aware as _map_reaction_aware_raw,
    )
    _HAS_REACTION_AWARE = True
except ImportError:
    _HAS_REACTION_AWARE = False

# Bond-change scoring + batch constrained MCS (v6.6.0)
try:
    from smsd._smsd import (
        bond_change_score as _bond_change_score_raw,
        batch_mcs_constrained as _batch_mcs_constrained_raw,
    )
    _HAS_BOND_CHANGE = True
except ImportError:
    _HAS_BOND_CHANGE = False

# SMARTS matching — lazy import (available when compiled with smarts_parser.hpp)
try:
    from smsd._smsd import (
        find_mcs_smarts as _find_mcs_smarts_raw,
        smarts_match as _smarts_match_raw,
        smarts_find_all as _smarts_find_all_raw,
    )
    _HAS_SMARTS = True
except ImportError:
    _HAS_SMARTS = False


# ---------------------------------------------------------------------------
# Convenience helpers that wrap the raw C++ bindings with Pythonic defaults
# ---------------------------------------------------------------------------

def similarity(mol1, mol2):
    """Compute RASCAL similarity upper bound between two molecules.

    Accepts either MolGraph objects or SMILES strings.

    Returns:
        float: Tanimoto-like similarity in [0.0, 1.0].
    """
    g1 = _ensure_mol(mol1)
    g2 = _ensure_mol(mol2)
    result = similarity_upper_bound(g1, g2)
    # Clamp to [0.0, 1.0] for safety against floating-point drift
    return max(0.0, min(1.0, result))


def mcs(mol1, mol2, *, tautomer_aware=False, prefer_rare_heteroatoms=False,
        timeout_ms=10000, **kwargs):
    """Find Maximum Common Substructure between two molecules.

    Accepts MolGraph objects, SMILES strings, or RDKit Mol objects.
    When RDKit Mol objects are passed, returned indices are automatically
    translated to RDKit's native atom ordering.

    Args:
        mol1: First molecule (MolGraph, SMILES string, or RDKit Mol).
        mol2: Second molecule (MolGraph, SMILES string, or RDKit Mol).
        tautomer_aware: Use tautomer-aware matching.
        prefer_rare_heteroatoms: Prefer mappings that include rare
            heteroatoms (S, P, Se) even if slightly smaller. Activates
            reaction-aware MCS with heteroatom-weighted scoring.
        timeout_ms: Timeout in milliseconds.
        **kwargs: Additional McsOptions fields (e.g. connected_only=False).

    Returns:
        dict: Mapping from mol1 atom indices to mol2 atom indices.

    .. versionchanged:: 6.5.2
       Now accepts RDKit Mol objects and auto-translates indices.
       Added *prefer_rare_heteroatoms* parameter.
    """
    if prefer_rare_heteroatoms:
        return map_reaction_aware(mol1, mol2, timeout_ms=timeout_ms, **kwargs)
    g1, rdkit1 = _ensure_mol_ex(mol1)
    g2, rdkit2 = _ensure_mol_ex(mol2)
    chem = ChemOptions.tautomer_profile() if tautomer_aware else ChemOptions()
    opts = McsOptions()
    opts.timeout_ms = timeout_ms
    for k, v in kwargs.items():
        setattr(opts, k, v)
    mapping = find_mcs(g1, g2, chem, opts)
    translated, _ = _auto_translate(mapping, g1, g2, rdkit1, rdkit2)
    return translated


def map_reaction_aware(mol1, mol2, *, timeout_ms=10000,
                       bond_change_aware=False, **kwargs):
    """Reaction-aware MCS: find MCS candidates, generate near-MCS variants
    (K-1, K-2), and re-rank by heteroatom coverage, rare-element importance,
    and connectivity.

    Prefers mappings that capture reaction-center heteroatoms (S, P, Se)
    even if they are slightly smaller than the mathematical maximum.

    Accepts MolGraph objects, SMILES strings, or RDKit Mol objects.
    When RDKit Mol objects are passed, returned indices are automatically
    translated to RDKit's native atom ordering.

    Args:
        mol1: First molecule (MolGraph, SMILES string, or RDKit Mol).
        mol2: Second molecule (MolGraph, SMILES string, or RDKit Mol).
        timeout_ms: Timeout in milliseconds.
        bond_change_aware: When True, rank candidates by bond-change
            plausibility instead of heteroatom coverage (v6.6.0).
        **kwargs: Additional McsOptions fields (e.g. near_mcs_delta=3).

    Returns:
        dict: Mapping from mol1 atom indices to mol2 atom indices.

    Raises:
        RuntimeError: If the C++ library was compiled without reaction-aware
            MCS support.

    Example::

        mapping = map_reaction_aware("SAM_smiles", "SAH_smiles")
        # mapping will prefer including the S atom even if MCS is 1 atom smaller

    .. versionadded:: 6.5.0
    .. versionchanged:: 6.5.2
       Now accepts RDKit Mol objects and auto-translates indices.
    """
    if not _HAS_REACTION_AWARE:
        raise RuntimeError(
            "Reaction-aware MCS is not available in this build. "
            "Rebuild with v6.5.0+ C++ sources.")
    g1, rdkit1 = _ensure_mol_ex(mol1)
    g2, rdkit2 = _ensure_mol_ex(mol2)
    chem = ChemOptions()
    opts = McsOptions()
    opts.timeout_ms = timeout_ms
    opts.reaction_aware = True
    if bond_change_aware and hasattr(opts, 'bond_change_aware'):
        opts.bond_change_aware = True
    for k, v in kwargs.items():
        setattr(opts, k, v)
    mapping = _reaction_aware_mcs_raw(g1, g2, chem, opts, timeout_ms)
    translated, _ = _auto_translate(mapping, g1, g2, rdkit1, rdkit2)
    return translated


def find_mcs_progressive(mol1, mol2, *, on_progress=None, tautomer_aware=False,
                         timeout_ms=10000, **kwargs):
    """Find MCS with progressive intermediate reporting via callback.

    The search runs in phases with increasing time budgets. After each phase,
    on_progress is called with the best mapping found so far, its size, and
    the elapsed wall-clock time in milliseconds.

    Args:
        mol1: First molecule (MolGraph or SMILES string).
        mol2: Second molecule (MolGraph or SMILES string).
        on_progress: Callback ``(best_mapping, best_size, elapsed_ms) -> None``.
            Called after each pipeline phase. May be None (no reporting).
        tautomer_aware: Use tautomer-aware matching.
        timeout_ms: Total timeout in milliseconds.
        **kwargs: Additional McsOptions fields (e.g. connected_only=False).

    Returns:
        dict: Final best mapping from mol1 atom indices to mol2 atom indices.

    Example::

        def report(mapping, size, ms):
            print(f"  {ms:>5d} ms  MCS size = {size}")

        mapping = smsd.find_mcs_progressive("c1ccccc1O", "c1ccccc1N",
                                            on_progress=report)
    """
    import time

    g1 = _ensure_mol(mol1)
    g2 = _ensure_mol(mol2)
    chem = ChemOptions.tautomer_profile() if tautomer_aware else ChemOptions()
    opts = McsOptions()
    opts.timeout_ms = timeout_ms
    for k, v in kwargs.items():
        setattr(opts, k, v)

    if on_progress is None:
        return find_mcs(g1, g2, chem, opts)

    phase_fractions = [0.02, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00]
    best = {}
    best_size = 0
    t0 = time.monotonic()
    min_n = min(g1.n, g2.n)

    for frac in phase_fractions:
        phase_opts = McsOptions()
        phase_opts.timeout_ms = int(timeout_ms * frac)
        # Copy relevant fields
        for attr in ("induced", "connected_only", "disconnected_mcs",
                     "maximize_bonds", "min_fragment_size", "max_fragments"):
            if hasattr(opts, attr):
                setattr(phase_opts, attr, getattr(opts, attr))
        phase_opts.extra_seeds = (frac >= 1.0) and getattr(opts, "extra_seeds", True)

        result = find_mcs(g1, g2, chem, phase_opts)
        if len(result) > best_size:
            best = result
            best_size = len(best)

        elapsed_ms = int((time.monotonic() - t0) * 1000)
        on_progress(dict(best), best_size, elapsed_ms)

        if best_size >= min_n:
            break

    return best


def all_mcs(mol1, mol2, *, max_results=10, tautomer_aware=False, timeout_ms=10000, **kwargs):
    """Find multiple distinct MCS mappings of the maximum size.

    Accepts either MolGraph objects or SMILES strings. Returns up to
    ``max_results`` distinct atom-atom mappings that all share the same
    (maximum) MCS size. Useful for scaffold analysis and SAR studies.

    Args:
        mol1: First molecule (MolGraph or SMILES string).
        mol2: Second molecule (MolGraph or SMILES string).
        max_results: Maximum number of distinct mappings to return (default 10).
        tautomer_aware: Use tautomer-aware matching.
        timeout_ms: Timeout in milliseconds.
        **kwargs: Additional McsOptions fields (e.g. connected_only=False).

    Returns:
        list[dict]: List of mappings from mol1 atom indices to mol2 atom indices.
            Each mapping has the same (maximum) size.
    """
    g1 = _ensure_mol(mol1)
    g2 = _ensure_mol(mol2)
    chem = ChemOptions.tautomer_profile() if tautomer_aware else ChemOptions()
    opts = McsOptions()
    opts.timeout_ms = timeout_ms
    for k, v in kwargs.items():
        setattr(opts, k, v)
    return find_all_mcs(g1, g2, chem, opts, max_results)


def mcs_smiles(mol1, mol2, *, tautomer_aware=False, timeout_ms=10000, **kwargs):
    """Find MCS between two molecules and return it as a canonical SMILES string.

    Accepts either MolGraph objects or SMILES strings.

    Args:
        mol1: First molecule (MolGraph or SMILES string).
        mol2: Second molecule (MolGraph or SMILES string).
        tautomer_aware: Use tautomer-aware matching.
        timeout_ms: Timeout in milliseconds.
        **kwargs: Additional McsOptions fields (e.g. connected_only=False).

    Returns:
        str: Canonical SMILES of the MCS, or "" if no common substructure.
    """
    g1 = _ensure_mol(mol1)
    g2 = _ensure_mol(mol2)
    chem = ChemOptions.tautomer_profile() if tautomer_aware else ChemOptions()
    opts = McsOptions()
    opts.timeout_ms = timeout_ms
    for k, v in kwargs.items():
        setattr(opts, k, v)
    return find_mcs_smiles(g1, g2, chem, opts, timeout_ms)


def substructure_search(query, target, *, timeout_ms=10000):
    """Find a substructure mapping of query in target.

    Accepts either MolGraph objects or SMILES strings.

    Returns:
        list: List of (query_atom, target_atom) pairs, or empty list if not
        a substructure.
    """
    q = _ensure_mol(query)
    t = _ensure_mol(target)
    return find_substructure(q, t, ChemOptions(), timeout_ms)


def topological_torsion(mol, fp_size=2048):
    """Compute a topological torsion fingerprint (Nilakantan et al. 1987).

    Enumerates all 4-atom linear paths A-B-C-D in the molecule and hashes
    their atom types (atomicNum, heavyDeg, piElectrons, ring) using FNV-1a
    with canonical ordering.

    Args:
        mol: MolGraph or SMILES string.
        fp_size: Fingerprint size in bits (default 2048).

    Returns:
        list[int]: Sorted list of set bit positions.
    """
    if not isinstance(fp_size, int) or fp_size <= 0:
        raise ValueError(f"fp_size must be a positive integer, got {fp_size!r}")
    from smsd._smsd import topological_torsion as _tt  # type: ignore[import]
    g = _ensure_mol(mol)
    return _tt(g, fp_size)


def fingerprint(mol, *, kind="mcs", path_length=7, fp_size=2048):
    """Compute a molecular fingerprint.

    Args:
        mol: MolGraph or SMILES string.
        kind: "path" for simple path fingerprint, "mcs" for MCS-aware.
        path_length: Maximum bond path length (default 7).
        fp_size: Fingerprint size in bits (default 2048).

    Returns:
        list[int]: Sorted list of set bit positions.
    """
    if not isinstance(fp_size, int) or fp_size <= 0:
        raise ValueError(f"fp_size must be a positive integer, got {fp_size!r}")
    if not isinstance(path_length, int) or path_length <= 0:
        raise ValueError(f"path_length must be a positive integer, got {path_length!r}")
    g = _ensure_mol(mol)
    if kind == "path":
        return path_fingerprint(g, path_length, fp_size)
    elif kind == "mcs":
        return mcs_fingerprint(g, path_length, fp_size)
    else:
        raise ValueError(f"Unknown fingerprint kind: {kind!r}. Use 'path' or 'mcs'.")


def tanimoto(fp1, fp2, fp_size=2048):
    """Compute Tanimoto coefficient between two fingerprints.

    Args:
        fp1: List of set bit positions (from fingerprint()).
        fp2: List of set bit positions (from fingerprint()).
        fp_size: Fingerprint size in bits.

    Returns:
        float: Tanimoto similarity in [0.0, 1.0].
    """
    if fp1 is None or fp2 is None:
        return 0.0
    s1 = set(fp1)
    s2 = set(fp2)
    intersection = len(s1 & s2)
    union = len(s1 | s2)
    if union == 0:
        return 0.0
    return intersection / union


def batch_substructure(query, targets, *, timeout_ms=10000):
    """Check whether query is a substructure of each molecule in targets.

    Automatically dispatches to multi-core OpenMP via the C++ extension.
    Pass pre-parsed MolGraph objects for best performance (avoids per-call
    SMILES parsing).

    Args:
        query:      MolGraph or SMILES string.
        targets:    List of MolGraph or SMILES strings.
        timeout_ms: Per-pair timeout in milliseconds.

    Returns:
        list[bool]: True at index i when query is a substructure of targets[i].
    """
    from smsd._smsd import batch_substructure as _batch_sub  # type: ignore[import]
    q = _ensure_mol(query)
    ts = [_ensure_mol(t) for t in targets]
    return _batch_sub(q, ts, ChemOptions(), timeout_ms)


def batch_mcs(query, targets, *, timeout_ms=10000, **kwargs):
    """Compute MCS between query and every molecule in targets in parallel.

    Dispatches to multi-core OpenMP automatically.

    Args:
        query:      MolGraph or SMILES string.
        targets:    List of MolGraph or SMILES strings.
        timeout_ms: Per-pair timeout in milliseconds.
        **kwargs:   Extra McsOptions fields (e.g. connected_only=False).

    Returns:
        list[dict]: MCS atom mappings (query → target index) for each target.
    """
    from smsd._smsd import batch_mcs as _batch_mcs  # type: ignore[import]
    q = _ensure_mol(query)
    ts = [_ensure_mol(t) for t in targets]
    opts = McsOptions()
    opts.timeout_ms = timeout_ms
    for k, v in kwargs.items():
        setattr(opts, k, v)
    return _batch_mcs(q, ts, ChemOptions(), opts)


def assign_rs(mol):
    """Assign CIP R/S descriptors to all tetrahedral stereocentres.

    Accepts a MolGraph or SMILES string.

    Returns:
        dict: Mapping from atom index (int) to descriptor ('R' or 'S').
            Only stereocentres with resolvable, distinct CIP priorities
            are included. Atoms without chirality annotation are omitted.

    Example::

        result = smsd.assign_rs("N[C@@H](C)C(=O)O")  # L-alanine
        # result: {1: 'S'}
    """
    g = _ensure_mol(mol)
    desc = _assign_cip_raw(g)
    result = {}
    for i, label in enumerate(desc.rs_labels):
        if label == RSLabel.R:
            result[i] = 'R'
        elif label == RSLabel.S:
            result[i] = 'S'
    return result


def assign_ez(mol):
    """Assign CIP E/Z descriptors to all stereogenic double bonds.

    Accepts a MolGraph or SMILES string.

    Returns:
        dict: Mapping from (atom1, atom2) tuple to descriptor ('E' or 'Z').
            Only double bonds with stereo annotations and resolvable
            CIP priorities are included.

    Example::

        result = smsd.assign_ez("C/C=C/C")  # (E)-2-butene
        # result: {(1, 2): 'E'}
    """
    g = _ensure_mol(mol)
    desc = _assign_cip_raw(g)
    result = {}
    for a1, a2, label in desc.ez_bonds:
        if label == EZLabel.E:
            result[(a1, a2)] = 'E'
        elif label == EZLabel.Z:
            result[(a1, a2)] = 'Z'
    return result


def assign_cip(mol):
    """Assign all CIP stereo descriptors (R/S and E/Z) for a molecule.

    Accepts a MolGraph or SMILES string.

    Returns:
        CIPDescriptors: Object with rs_labels (per-atom R/S) and
            ez_bonds (list of (atom1, atom2, EZLabel) tuples).

    Example::

        desc = smsd.assign_cip("N[C@@H](C)C(=O)O")
        print(desc.rs_labels)  # [NONE, S, NONE, ...]
    """
    g = _ensure_mol(mol)
    return _assign_cip_raw(g)


# ---------------------------------------------------------------------------
# Bond-change scoring (v6.6.0)
# ---------------------------------------------------------------------------

def bond_change_score(mapping, mol1, mol2):
    """Score an MCS mapping by bond-change plausibility.

    Lower score = more chemically plausible mapping.  Useful for
    selecting the best MCS candidate in reaction atom-atom mapping.

    Bond-change costs (biochemistry-informed):
      - C-C bond broken/formed: 3.0 (rare in most biochemistry)
      - C-N, C-O bond change:  1.5 (amide hydrolysis, esterification)
      - Heteroatom bonds (S, P): 0.5 (common reaction centres)
      - Bond order change:      1.0

    Args:
        mapping: dict {query_idx: target_idx} from MCS.
        mol1: Query molecule (MolGraph or SMILES string).
        mol2: Target molecule (MolGraph or SMILES string).

    Returns:
        float: Bond-change penalty score (lower is better).

    Example::

        score = smsd.bond_change_score({0: 0, 1: 1}, "CCO", "CC=O")

    .. versionadded:: 6.6.0
    """
    if not _HAS_BOND_CHANGE:
        raise RuntimeError(
            "bond_change_score requires v6.6.0+ C++ sources. "
            "Rebuild with updated sources.")
    g1 = _ensure_mol(mol1)
    g2 = _ensure_mol(mol2)
    return _bond_change_score_raw(mapping, g1, g2)


# ---------------------------------------------------------------------------
# Batch constrained MCS (v6.6.0)
# ---------------------------------------------------------------------------

def batch_mcs_constrained(queries, targets, *, timeout_ms=10000, **kwargs):
    """Find MCS for each query against ALL targets with non-overlapping
    target atom constraints.

    For multi-component reaction mapping: N reactants to M products.
    Each query finds its best MCS across all targets. Target atoms
    used by earlier (larger) queries are excluded for later queries.

    Accepts MolGraph, SMILES strings, or native Mol objects (auto-detected).
    When native Mol objects are passed, returned indices are automatically
    translated to native atom ordering.

    Args:
        queries: list of reactant molecules.
        targets: list of product molecules (often just one combined product).
        timeout_ms: per-pair timeout in milliseconds.
        **kwargs: additional McsOptions fields.

    Returns:
        list[dict]: One mapping per query (same order as input).
            Each dict maps {query_atom_idx: target_atom_idx}.

    Example::

        # Two reactants mapping to one product with non-overlapping atoms
        mappings = smsd.batch_mcs_constrained(
            [aryl_bromide, aniline],  # reactant fragments
            [biaryl_amine],           # combined product
        )
        # mappings[0] = {ArBr atoms → product atoms}
        # mappings[1] = {ArNH2 atoms → product atoms} (non-overlapping with [0])

    .. versionchanged:: 6.6.0
       Redesigned: each query maps against ALL targets (not 1:1 pairing).
       Auto-detects input types and translates indices.
    """
    if not _HAS_BOND_CHANGE:
        raise RuntimeError(
            "batch_mcs_constrained requires v6.6.0+ C++ sources.")

    # Convert inputs, tracking which are native Mol objects
    gq_list, rdkit_q_list = [], []
    for q in queries:
        g, rdkit_mol = _ensure_mol_ex(q)
        gq_list.append(g)
        rdkit_q_list.append(rdkit_mol)

    gt_list, rdkit_t_list = [], []
    for t in targets:
        g, rdkit_mol = _ensure_mol_ex(t)
        gt_list.append(g)
        rdkit_t_list.append(rdkit_mol)

    opts = McsOptions()
    opts.timeout_ms = timeout_ms
    for k, v in kwargs.items():
        setattr(opts, k, v)

    raw_results = _batch_mcs_constrained_raw(gq_list, gt_list, ChemOptions(), opts)

    # Translate indices if native Mol objects were used
    translated_results = []
    for i, mapping in enumerate(raw_results):
        if not mapping:
            translated_results.append(mapping)
            continue
        # Find which target this mapping refers to (by checking which target
        # atoms are present). For now, translate query indices using query's map.
        rdkit_q = rdkit_q_list[i] if i < len(rdkit_q_list) else None
        # Target index translation: find the target that contains the mapped atoms
        # For simplicity, translate using first target's map (common case: 1 target)
        rdkit_t = rdkit_t_list[0] if rdkit_t_list else None
        translated, _ = _auto_translate(mapping, gq_list[i], gt_list[0], rdkit_q, rdkit_t)
        translated_results.append(translated)

    return translated_results


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# SMILES utilities (v6.6.4)
# ---------------------------------------------------------------------------

def strip_atom_maps(smiles):
    """Remove atom map numbers from a SMILES string.

    Correctly handles the ambiguity between atom map ``:N`` inside brackets
    and aromatic ring closure ``:N`` outside brackets.  A naive regex like
    ``re.sub(r':\\d+', '', smi)`` incorrectly strips ring closures such as
    ``:1`` in ``[CH]:[CH]:1``.

    Args:
        smiles: SMILES string with atom map numbers.

    Returns:
        str: SMILES with atom maps removed, ring closures preserved.

    Example::

        clean = smsd.strip_atom_maps("[CH3:1][C:2](=[O:3])[OH:4]")
        # Returns: "[CH3][C](=[O])[OH]"

        # Preserves aromatic ring closures:
        clean = smsd.strip_atom_maps("[Cl:1][C:2]1:[CH:3]:[CH:4]:1")
        # Returns: "[Cl][C]1:[CH]:[CH]:1"  (ring :1 preserved)

    .. versionadded:: 6.6.4
    """
    result = []
    in_bracket = False
    i = 0
    n = len(smiles)
    while i < n:
        c = smiles[i]
        if c == '[':
            in_bracket = True
            result.append(c)
            i += 1
        elif c == ']':
            in_bracket = False
            result.append(c)
            i += 1
        elif c == ':' and in_bracket:
            # Skip :digits inside bracket (atom map number)
            i += 1
            while i < n and smiles[i].isdigit():
                i += 1
        else:
            result.append(c)
            i += 1
    return ''.join(result)


def parse_mapped_smiles(smiles):
    """Parse a SMILES string that contains atom map numbers.

    Strips atom maps first, then parses. This is the safe way to convert
    reaction SMILES fragments to MolGraph for MCS computation.

    Args:
        smiles: SMILES string, may contain atom map numbers like ``:1``.

    Returns:
        MolGraph: Parsed molecular graph.

    Example::

        g = smsd.parse_mapped_smiles("[CH3:1][C:2](=[O:3])[OH:4]")
        print(g.n)  # 4

    .. versionadded:: 6.6.4
    """
    return parse_smiles(strip_atom_maps(smiles))


# Ring layout utilities
# ---------------------------------------------------------------------------

def compute_sssr(mol):
    """Return the Smallest Set of Smallest Rings (minimum cycle basis).

    Rings are sorted ascending by size. No redundant candidates.

    Accepts a MolGraph or SMILES string.

    Returns:
        list[list[int]]: SSSR rings, each ring as a list of atom indices.

    Example::

        rings = smsd.compute_sssr("c1ccc2ccccc2c1")  # naphthalene
        # rings: [[0, 1, 2, 3, 8, 9], [3, 4, 5, 6, 7, 8]]  (two 6-rings)
    """
    g = _ensure_mol(mol)
    return _compute_sssr_raw(g)


def layout_sssr(mol):
    """Return SSSR rings optimized for 2D coordinate generation.

    The rings are ordered for layout purposes:
    - Largest ring system first
    - Fused rings ordered by shared-edge adjacency within each system
    - Within each system, rings sorted by size (smallest first for placement)

    Accepts a MolGraph or SMILES string.

    Returns:
        list[list[int]]: Layout-ordered SSSR rings.

    Example::

        rings = smsd.layout_sssr("c1ccc2ccccc2c1")  # naphthalene
    """
    g = _ensure_mol(mol)
    return _layout_sssr_raw(g)


def reduce_crossings(mol, coords, max_iter=1000):
    """Reduce bond crossings in a 2D layout by optimizing ring orientations.

    Uses simulated annealing to flip ring systems and minimize crossings.

    Accepts a MolGraph or SMILES string for the molecule, and a list of
    [x, y] coordinate pairs for each atom.

    Args:
        mol: MolGraph or SMILES string.
        coords: List of [x, y] pairs for each atom. Modified in place.
        max_iter: Maximum number of simulated annealing iterations (default 1000).

    Returns:
        tuple: (crossings_remaining, updated_coords)

    Example::

        mol = smsd.parse_smiles("c1ccccc1")
        # ... generate initial coords ...
        crossings, new_coords = smsd.reduce_crossings(mol, coords)
    """
    g = _ensure_mol(mol)
    # Ensure coords is a list of [x, y] lists
    coords_list = [[float(c[0]), float(c[1])] for c in coords]
    return _reduce_crossings_raw(g, coords_list, max_iter)


def force_directed_layout(mol, coords, max_iter=500, target_bond_length=1.5):
    """Force-directed 2D layout minimisation.

    Iteratively moves atoms to minimise bond-length stress, non-bonded
    repulsion, and crossing penalties. Suitable for refining an existing
    layout or producing crossing-free coordinates for polycyclic systems.

    Accepts a MolGraph or SMILES string for the molecule, and a list of
    [x, y] coordinate pairs for each atom.

    Args:
        mol: MolGraph or SMILES string.
        coords: List of [x, y] pairs for each atom.
        max_iter: Maximum number of force iterations (default 500).
        target_bond_length: Desired bond length in coordinate units (default 1.5).

    Returns:
        tuple: (final_stress, updated_coords)

    Example::

        mol = smsd.parse_smiles("c1ccc2ccccc2c1")  # naphthalene
        coords = [[i * 1.0, 0.0] for i in range(10)]  # initial linear layout
        stress, new_coords = smsd.force_directed_layout(mol, coords)
    """
    g = _ensure_mol(mol)
    coords_list = [[float(c[0]), float(c[1])] for c in coords]
    return _force_directed_layout_raw(g, coords_list, max_iter, target_bond_length)


def stress_majorisation(mol, coords, max_iter=300, target_bond_length=1.5):
    """Stress majorisation (SMACOF algorithm) for 2D layout.

    Minimises weighted stress: sum_ij w_ij * (d_ij - D_ij)^2
    where D_ij is the graph-distance scaled by target_bond_length.
    Guarantees monotone convergence to a local minimum.

    Accepts a MolGraph or SMILES string for the molecule, and a list of
    [x, y] coordinate pairs for each atom.

    Args:
        mol: MolGraph or SMILES string.
        coords: List of [x, y] pairs for each atom.
        max_iter: Maximum SMACOF iterations (default 300).
        target_bond_length: Desired bond length in coordinate units (default 1.5).

    Returns:
        tuple: (final_stress, updated_coords)

    Example::

        mol = smsd.parse_smiles("c1ccc2ccccc2c1")  # naphthalene
        coords = [[i * 1.0, 0.0] for i in range(10)]  # initial linear layout
        stress, new_coords = smsd.stress_majorisation(mol, coords)
    """
    g = _ensure_mol(mol)
    coords_list = [[float(c[0]), float(c[1])] for c in coords]
    return _stress_majorisation_raw(g, coords_list, max_iter, target_bond_length)


def match_template(mol, target_bond_length=1.5):
    """Check if a molecule matches a known scaffold template.

    Returns pre-computed coordinates scaled to target_bond_length if a
    template match is found, empty list otherwise.

    Supported scaffolds: benzene, naphthalene, indole, purine, steroid,
    morphinan, quinoline, biphenyl, cyclohexane, piperidine.

    Accepts a MolGraph or SMILES string.

    Args:
        mol: MolGraph or SMILES string.
        target_bond_length: Desired bond length for scaling (default 1.5).

    Returns:
        list: List of (x, y) tuples if matched, empty list otherwise.

    Example::

        coords = smsd.match_template("c1ccccc1")  # benzene -> 6 coordinate pairs
        if coords:
            print(f"Template match: {len(coords)} atoms")
    """
    g = _ensure_mol(mol)
    return _match_template_raw(g, target_bond_length)


# ---------------------------------------------------------------------------
# MatchResult -- structured return type for MCS queries
# ---------------------------------------------------------------------------

class MatchResult:
    """Structured result from an MCS computation.

    Attributes:
        mapping (dict): Atom index mapping {mol1_idx: mol2_idx}.
        size (int): Number of matched atoms.
        tanimoto (float): Tanimoto-like similarity = size / min(query_atoms, target_atoms).
        query_atoms (int): Number of atoms in the query molecule.
        target_atoms (int): Number of atoms in the target molecule.
        mcs_smiles (str): Canonical SMILES of the MCS subgraph, or "" if empty.
    """
    __slots__ = ('mapping', 'size', 'overlap', 'tanimoto', 'query_atoms', 'target_atoms', 'mcs_smiles')

    def __init__(self, mapping, query_n, target_n, mcs_smi=""):
        self.mapping = mapping
        self.size = len(mapping)
        self.query_atoms = query_n
        self.target_atoms = target_n
        min_n = min(query_n, target_n)
        self.overlap = self.size / min_n if min_n > 0 else 0.0  # Simpson overlap
        denom = query_n + target_n - self.size
        self.tanimoto = self.size / denom if denom > 0 else 0.0  # standard Tanimoto
        self.mcs_smiles = mcs_smi

    def __len__(self):
        return self.size

    def __repr__(self):
        return f"MatchResult(size={self.size}, tanimoto={self.tanimoto:.3f}, mcs_smiles={self.mcs_smiles!r})"


def mcs_result(mol1, mol2, **kwargs):
    """Find MCS between two molecules and return a MatchResult.

    Like :func:`mcs`, but wraps the raw mapping dict in a
    :class:`MatchResult` that includes size and Tanimoto similarity.

    Accepts MolGraph objects, SMILES strings, or RDKit Mol objects.
    When RDKit Mol objects are passed, returned indices are automatically
    translated to RDKit's native atom ordering.

    Args:
        mol1: First molecule (MolGraph, SMILES string, or RDKit Mol).
        mol2: Second molecule (MolGraph, SMILES string, or RDKit Mol).
        **kwargs: Keyword arguments forwarded to :func:`mcs`
            (e.g. tautomer_aware, timeout_ms, prefer_rare_heteroatoms).

    Returns:
        MatchResult: Structured result with mapping, size, and tanimoto.

    Example::

        import smsd
        r = smsd.mcs_result("c1ccccc1", "c1ccc(O)cc1")
        print(r)          # MatchResult(size=6, tanimoto=1.000)
        print(len(r))     # 6
        print(r.mapping)  # {0: 0, 1: 1, ...}

    .. versionchanged:: 6.5.3
       Fixed: passes original mol objects (not MolGraph) so RDKit
       auto-detection works. Uses _auto_translate tuple to avoid
       redundant MCS re-computation for SMILES extraction.
    """
    g1, rdkit1 = _ensure_mol_ex(mol1)
    g2, rdkit2 = _ensure_mol_ex(mol2)
    # Run MCS on MolGraphs, get both translated and original SMSD mapping
    prefer = kwargs.pop('prefer_rare_heteroatoms', False)
    taut = kwargs.pop('tautomer_aware', False)
    tms = kwargs.pop('timeout_ms', 10000)
    if prefer:
        raw = _reaction_aware_mcs_raw(g1, g2, ChemOptions(), McsOptions(), tms) \
            if _HAS_REACTION_AWARE else find_mcs(g1, g2)
    else:
        chem = ChemOptions.tautomer_profile() if taut else ChemOptions()
        opts = McsOptions()
        opts.timeout_ms = tms
        for k, v in kwargs.items():
            setattr(opts, k, v)
        raw = find_mcs(g1, g2, chem, opts)
    # Translate + keep SMSD mapping for SMILES extraction
    mapping, smsd_mapping = _auto_translate(raw, g1, g2, rdkit1, rdkit2)
    smi = ""
    if smsd_mapping:
        try:
            smi = mcs_to_smiles(g1, smsd_mapping)
        except Exception:
            smi = ""
    return MatchResult(mapping, g1.n, g2.n, smi)


# ---------------------------------------------------------------------------
# mcs_from_smiles -- all-in-one SMILES -> MCS MatchResult
# ---------------------------------------------------------------------------

def mcs_from_smiles(smi1, smi2, **kwargs):
    """All-in-one convenience: parse two SMILES, compute MCS, return MatchResult.

    Combines :func:`parse_smiles`, :func:`find_mcs`, :func:`mcs_to_smiles`
    into a single call that returns a :class:`MatchResult` with mapping, size,
    Tanimoto, **and** ``mcs_smiles``.

    Args:
        smi1 (str): SMILES string for the first molecule.
        smi2 (str): SMILES string for the second molecule.
        **kwargs: Keyword arguments forwarded to :func:`mcs`
            (e.g. tautomer_aware, timeout_ms, connected_only).

    Returns:
        MatchResult: Structured result including ``mcs_smiles``.

    Example::

        r = smsd.mcs_from_smiles("c1ccccc1", "c1ccc(O)cc1")
        print(r.mcs_smiles)  # canonical SMILES of the MCS
        print(r.tanimoto)    # overlap score
    """
    if not isinstance(smi1, str) or not smi1:
        raise ValueError("smi1 must be a non-empty SMILES string")
    if not isinstance(smi2, str) or not smi2:
        raise ValueError("smi2 must be a non-empty SMILES string")
    g1 = parse_smiles(smi1)
    g2 = parse_smiles(smi2)
    mapping = mcs(g1, g2, **kwargs)
    smi = mcs_to_smiles(g1, mapping) if mapping else ""
    return MatchResult(mapping, g1.n, g2.n, smi)


# ---------------------------------------------------------------------------
# fingerprint_from_smiles -- all-in-one SMILES -> ECFP fingerprint
# ---------------------------------------------------------------------------

def fingerprint_from_smiles(smiles, radius=2, fp_size=2048, mode="ecfp"):
    """All-in-one convenience: parse a SMILES string and compute its fingerprint.

    Args:
        smiles (str): SMILES string.
        radius (int): ECFP radius (2 = ECFP4, 3 = ECFP6; default 2).
        fp_size (int): Fingerprint size in bits (default 2048).
        mode (str): ``"ecfp"`` (circular/Morgan), ``"fcfp"`` (pharmacophoric),
                    ``"path"`` (path-based), or ``"mcs"`` (MCS-aware).

    Returns:
        list[int]: Sorted list of set bit positions.

    Example::

        fp = smsd.fingerprint_from_smiles("c1ccccc1")
        print(len(fp))  # number of set bits
    """
    if not isinstance(smiles, str) or not smiles:
        raise ValueError("smiles must be a non-empty SMILES string")
    if not isinstance(fp_size, int) or fp_size <= 0:
        raise ValueError(f"fp_size must be a positive integer, got {fp_size!r}")
    g = parse_smiles(smiles)
    if mode == "ecfp":
        return circular_fingerprint(g, radius, fp_size, mode="ecfp")
    elif mode == "fcfp":
        return circular_fingerprint(g, radius, fp_size, mode="fcfp")
    elif mode == "path":
        return path_fingerprint(g, radius, fp_size)
    elif mode == "mcs":
        return mcs_fingerprint(g, radius, fp_size)
    else:
        raise ValueError(f"Unknown mode: {mode!r}. Use 'ecfp', 'fcfp', 'path', or 'mcs'.")


# ---------------------------------------------------------------------------
# counts_to_array -- sparse count dict -> dense list
# ---------------------------------------------------------------------------

def to_hex(fp, fp_size=2048):
    """Convert a fingerprint (list of set-bit positions) to a hex string.

    Useful for storing fingerprints in databases or transmitting via REST APIs.

    Args:
        fp (list[int]): Sorted list of set bit positions.
        fp_size (int): Total fingerprint size in bits (default 2048).

    Returns:
        str: Hex string representation of the fingerprint.

    Example::

        fp = smsd.path_fingerprint(mol, fp_size=2048)
        hex_str = smsd.to_hex(fp, fp_size=2048)
        fp_back = smsd.from_hex(hex_str)
        assert fp_back == fp
    """
    if not isinstance(fp_size, int) or fp_size <= 0:
        raise ValueError(f"fp_size must be a positive integer, got {fp_size!r}")
    n_words = (fp_size + 63) // 64
    words = [0] * n_words
    if fp:
        for bit in fp:
            if 0 <= bit < fp_size:
                words[bit >> 6] |= 1 << (bit & 63)
    return "".join(f"{w:016x}" for w in words)


def from_hex(hex_str):
    """Convert a hex string back to a fingerprint (list of set-bit positions).

    Inverse of :func:`to_hex`.

    Args:
        hex_str (str): Hex string produced by :func:`to_hex`.

    Returns:
        list[int]: Sorted list of set bit positions.

    Example::

        fp_back = smsd.from_hex(hex_str)
    """
    if not isinstance(hex_str, str) or not hex_str:
        return []
    n_words = (len(hex_str) + 15) // 16
    bits = []
    for i in range(n_words):
        start = i * 16
        end = min(start + 16, len(hex_str))
        word = int(hex_str[start:end], 16)
        for b in range(64):
            if word & (1 << b):
                bits.append(i * 64 + b)
    return bits


def to_binary_string(fp, fp_size=2048):
    """Convert a fingerprint to a binary string of '0' and '1' characters.

    Args:
        fp (list[int]): Sorted list of set bit positions.
        fp_size (int): Total fingerprint size in bits (default 2048).

    Returns:
        str: Binary string of length ``fp_size``.

    Example::

        bits = smsd.to_binary_string(fp, fp_size=2048)
        assert len(bits) == 2048
        assert all(c in ('0', '1') for c in bits)
    """
    if not isinstance(fp_size, int) or fp_size <= 0:
        raise ValueError(f"fp_size must be a positive integer, got {fp_size!r}")
    fp_set = set(fp) if fp else set()
    return "".join("1" if i in fp_set else "0" for i in range(fp_size))


def counts_to_array(counts_dict, fp_size):
    """Convert a sparse count dict ``{bit_position: count}`` into a dense list
    of length ``fp_size``.

    Useful for serialising count-based fingerprints to fixed-length vectors
    compatible with NumPy, database columns, or ML feature arrays.

    Positions outside ``[0, fp_size)`` are silently ignored.

    Args:
        counts_dict (dict): Sparse count map ``{int: int}``.
        fp_size (int): Target array length.

    Returns:
        list[int]: Dense count array of length ``fp_size``.

    Example::

        arr = smsd.counts_to_array({42: 3, 1000: 1}, 2048)
        assert arr[42] == 3
        assert arr[1000] == 1
    """
    if not isinstance(fp_size, int) or fp_size <= 0:
        raise ValueError(f"fp_size must be a positive integer, got {fp_size!r}")
    arr = [0] * fp_size
    if counts_dict is None:
        return arr
    for pos, count in counts_dict.items():
        if 0 <= pos < fp_size:
            arr[pos] = count
    return arr


# ---------------------------------------------------------------------------
# batch_similarity -- sorted pairwise similarity
# ---------------------------------------------------------------------------

def batch_similarity(query, targets, fp_size=2048, radius=2):
    """Compute fingerprint similarity between a query and a list of targets,
    returning results sorted by decreasing Tanimoto score.

    Computes ECFP fingerprints for the query and each target, then returns
    ``(index, tanimoto)`` pairs sorted from most to least similar.

    Args:
        query: Query molecule (MolGraph or SMILES string).
        targets: List of target molecules (MolGraph or SMILES strings).
        fp_size (int): Fingerprint size in bits (default 2048).
        radius (int): ECFP radius (default 2).

    Returns:
        list[tuple[int, float]]: Sorted ``(target_index, tanimoto)`` pairs,
            highest similarity first.

    Example::

        results = smsd.batch_similarity("c1ccccc1", ["c1ccc(O)cc1", "CCCC"])
        for idx, sim in results:
            print(f"Target {idx}: Tanimoto = {sim:.3f}")
    """
    q = _ensure_mol(query)
    q_fp = path_fingerprint(q, radius, fp_size)
    q_set = set(q_fp)

    results = []
    for i, t in enumerate(targets):
        g = _ensure_mol(t)
        t_fp = path_fingerprint(g, radius, fp_size)
        t_set = set(t_fp)
        intersection = len(q_set & t_set)
        union = len(q_set | t_set)
        sim = intersection / union if union > 0 else 0.0
        results.append((i, sim))

    results.sort(key=lambda x: x[1], reverse=True)
    return results


# ---------------------------------------------------------------------------
# R-Group Decomposition
# ---------------------------------------------------------------------------

def decompose_r_groups(core, molecules, *, timeout_ms=10000):
    """Decompose molecules into a core scaffold and R-groups.

    Finds the substructure mapping of core onto each molecule, then
    identifies R-group substituents attached to the core.

    Args:
        core: Core scaffold (MolGraph or SMILES string).
        molecules: List of molecules (MolGraph or SMILES strings).
        timeout_ms: Per-molecule timeout in milliseconds.

    Returns:
        list[dict]: One dict per molecule. Each dict maps:
            - ``"core"``: list of atom indices in the molecule that match the core
            - ``"R1"``, ``"R2"``, ...: lists of atom indices for each R-group
            Empty dict if the core is not a substructure of the molecule.

    Example::

        import smsd
        results = smsd.decompose_r_groups("c1ccccc1", ["c1ccc(O)cc1", "c1ccc(N)cc1"])
        for r in results:
            print(r.get("core"), r.get("R1"))

    Note:
        This is a pure-Python implementation. For CDK-backed R-group
        decomposition with full chemical perception, use the Java API
        ``SearchEngine.decomposeRGroups()`` directly.
    """
    g_core = _ensure_mol(core)
    results = []

    def iter_neighbors(mol, atom_idx):
        for nb in range(mol.n):
            if nb != atom_idx and mol.has_bond(atom_idx, nb):
                yield nb

    for mol in molecules:
        g_mol = _ensure_mol(mol)
        decomposition = {}

        # Find substructure mapping: core -> mol
        mapping_pairs = find_substructure(g_core, g_mol, ChemOptions(), timeout_ms)
        if not mapping_pairs or (isinstance(mapping_pairs, list) and len(mapping_pairs) == 0):
            results.append(decomposition)
            continue

        # Build mapping dict from substructure result
        if isinstance(mapping_pairs, dict):
            mapping = mapping_pairs
        elif isinstance(mapping_pairs, list) and len(mapping_pairs) > 0:
            if isinstance(mapping_pairs[0], tuple):
                mapping = {q: t for q, t in mapping_pairs}
            else:
                results.append(decomposition)
                continue
        else:
            results.append(decomposition)
            continue

        core_atom_indices = set(mapping.values())
        decomposition["core"] = sorted(core_atom_indices)

        # BFS to find connected components of non-core atoms
        r_group_counter = 1
        processed = set()

        for core_idx in sorted(core_atom_indices):
            if core_idx >= g_mol.n:
                continue
            for nb in iter_neighbors(g_mol, core_idx):
                if nb in core_atom_indices or nb in processed:
                    continue
                # BFS from this neighbor to find the full R-group
                r_group = []
                queue = [nb]
                visited = {nb}
                while queue:
                    current = queue.pop(0)
                    if current in core_atom_indices:
                        continue
                    r_group.append(current)
                    processed.add(current)
                    for next_nb in iter_neighbors(g_mol, current):
                        if next_nb not in visited and next_nb not in core_atom_indices:
                            visited.add(next_nb)
                            queue.append(next_nb)
                if r_group:
                    decomposition[f"R{r_group_counter}"] = sorted(r_group)
                    r_group_counter += 1

        results.append(decomposition)

    return results


# ---------------------------------------------------------------------------
# SMARTS-based MCS and matching
# ---------------------------------------------------------------------------

def find_mcs_smarts(smarts, target, *, max_matches=1000):
    """Find the largest SMARTS substructure match in a target molecule.

    Uses the C++ SMARTS matcher to enumerate all embeddings of the SMARTS
    pattern in the target, then returns the mapping with the most atoms.

    Args:
        smarts: SMARTS query string (e.g. "[CX3](=O)[NX3]").
        target: Target molecule (MolGraph or SMILES string).
        max_matches: Maximum embeddings to enumerate.

    Returns:
        dict: Mapping from SMARTS atom index to target atom index.
            Empty dict if no match or SMARTS support not compiled.

    Example::

        import smsd
        mapping = smsd.find_mcs_smarts("[CX3](=O)[NX3]", "CC(=O)Nc1ccc(O)cc1")
        print(f"Matched {len(mapping)} atoms")
    """
    if not _HAS_SMARTS:
        raise NotImplementedError(
            "SMARTS matching not available in this build. "
            "Rebuild with smarts_parser.hpp included."
        )
    g = _ensure_mol(target)
    return _find_mcs_smarts_raw(smarts, g, max_matches)


def smarts_match(smarts, target):
    """Check if a SMARTS pattern matches a target molecule.

    Args:
        smarts: SMARTS query string.
        target: Target molecule (MolGraph or SMILES string).

    Returns:
        bool: True if the pattern matches.

    Example::

        import smsd
        assert smsd.smarts_match("[#6]~[#7]", "CCN")
        assert not smsd.smarts_match("[#6]~[#7]", "CCC")
    """
    if not _HAS_SMARTS:
        raise NotImplementedError(
            "SMARTS matching not available in this build."
        )
    g = _ensure_mol(target)
    return _smarts_match_raw(smarts, g)


def smarts_find_all(smarts, target, *, max_matches=1000):
    """Find all SMARTS substructure matches in a target molecule.

    Args:
        smarts: SMARTS query string.
        target: Target molecule (MolGraph or SMILES string).
        max_matches: Maximum number of embeddings to return.

    Returns:
        list[dict]: List of mappings (SMARTS atom idx -> target atom idx).

    Example::

        import smsd
        matches = smsd.smarts_find_all("[#6]", "CCCC")
        assert len(matches) == 4  # 4 carbon atoms
    """
    if not _HAS_SMARTS:
        raise NotImplementedError(
            "SMARTS matching not available in this build."
        )
    g = _ensure_mol(target)
    return _smarts_find_all_raw(smarts, g, max_matches)


def _ensure_mol(obj):
    """Convert a SMILES string to MolGraph, or return MolGraph as-is."""
    if obj is None:
        raise TypeError("Expected MolGraph or SMILES string, got None")
    if isinstance(obj, str):
        if not obj:
            raise ValueError("SMILES string must not be empty")
        return parse_smiles(obj)
    if isinstance(obj, MolGraph):
        return obj
    # Auto-detect RDKit Mol via duck-typing (avoids hard import)
    if hasattr(obj, 'GetNumAtoms') and hasattr(obj, 'GetAtomWithIdx'):
        return from_rdkit(obj)
    raise TypeError(
        f"Expected MolGraph, SMILES string, or RDKit Mol, got {type(obj).__name__}"
    )


def _ensure_mol_ex(obj):
    """Like _ensure_mol, but also returns the original RDKit Mol (if any).

    Returns:
        tuple: (MolGraph, rdkit_mol_or_None)
    """
    if obj is None:
        raise TypeError("Expected MolGraph, SMILES string, or RDKit Mol, got None")
    if isinstance(obj, str):
        if not obj:
            raise ValueError("SMILES string must not be empty")
        return parse_smiles(obj), None
    if isinstance(obj, MolGraph):
        return obj, None
    # Auto-detect RDKit Mol via duck-typing
    if hasattr(obj, 'GetNumAtoms') and hasattr(obj, 'GetAtomWithIdx'):
        return from_rdkit(obj), obj
    raise TypeError(
        f"Expected MolGraph, SMILES string, or RDKit Mol, got {type(obj).__name__}"
    )


def _auto_translate(mapping, g1, g2, rdkit1, rdkit2):
    """Translate SMSD mapping to RDKit indices if either input was RDKit Mol.

    Returns a tuple: (translated_mapping, original_smsd_mapping).
    The original SMSD mapping is needed for mcs_to_smiles().

    .. versionchanged:: 6.5.3
       Returns tuple; validates even with mixed input types.
    """
    if rdkit1 is None and rdkit2 is None:
        return mapping, mapping  # no translation, both are the same
    if not mapping:
        return mapping, mapping
    smsd_mapping = dict(mapping)  # preserve original before translation
    translated = translate_mapping(mapping, g1, g2)
    # Element validation — check whichever RDKit Mol is available
    validated = {}
    for q, t in translated.items():
        ok = True
        if rdkit1 is not None:
            if q >= rdkit1.GetNumAtoms():
                ok = False
            elif rdkit2 is not None and t < rdkit2.GetNumAtoms():
                if rdkit1.GetAtomWithIdx(q).GetAtomicNum() != \
                   rdkit2.GetAtomWithIdx(t).GetAtomicNum():
                    ok = False
        if rdkit2 is not None and t >= rdkit2.GetNumAtoms():
            ok = False
        if ok:
            validated[q] = t
    return validated, smsd_mapping


# ---------------------------------------------------------------------------
# RDKit interoperability helpers (RDKit is an optional dependency)
# ---------------------------------------------------------------------------

class _RdkitResult:
    """Wrapper holding a MolGraph + its SMSD→RDKit index mapping."""
    __slots__ = ('graph', 'index_map')
    def __init__(self, graph, index_map):
        self.graph = graph
        self.index_map = index_map  # list: index_map[smsd_idx] = rdkit_idx


# Cache keyed by id(rdkit_mol) -> _RdkitResult.
# Per-Mol-object cache avoids the problem of different Mol objects with
# the same canonical SMILES having different atom orderings.
# NOTE: These caches are NOT thread-safe. In multi-threaded applications,
# guard from_rdkit() and clear_cache() calls with an external lock.
import threading as _threading
_rdkit_cache_lock = _threading.Lock()
_rdkit_cache = {}
_mol_graph_cache = {}  # SMILES -> MolGraph (for parse_smiles caching only)


def _compute_index_map(rdkit_mol, canonical_smi):
    """Compute deterministic SMSD→RDKit atom index mapping.

    Uses RDKit's _smilesAtomOutputOrder property which records the exact
    order atoms appear in the canonical SMILES string. This is the ground
    truth for SMILES parse order = SMSD's internal atom order.

    The mapping: index_map[smsd_idx] = rdkit_idx
    """
    from rdkit import Chem

    n = rdkit_mol.GetNumAtoms()

    # RDKit records which original atom appears at each position in the
    # canonical SMILES via the _smilesAtomOutputOrder property.
    # We must call MolToSmiles first to populate this property.
    _ = Chem.MolToSmiles(rdkit_mol)  # ensures property is set
    try:
        output_order = list(rdkit_mol.GetPropsAsDict()
                            .get('_smilesAtomOutputOrder', []))
    except Exception:
        output_order = []

    if output_order and len(output_order) >= n:
        # output_order[smiles_position] = rdkit_atom_idx
        # SMSD parses the canonical SMILES, so smsd_idx = smiles_position
        # But output_order may include implicit H indices beyond n — filter them
        index_map = []
        for rdkit_idx in output_order:
            if rdkit_idx < n:
                index_map.append(rdkit_idx)
        # Pad if needed (shouldn't happen for valid molecules)
        while len(index_map) < n:
            index_map.append(len(index_map))
        return index_map
    else:
        # Fallback: use GetSubstructMatch on re-parsed canonical SMILES
        mol_from_smi = Chem.MolFromSmiles(canonical_smi)
        if mol_from_smi is not None:
            match = rdkit_mol.GetSubstructMatch(mol_from_smi)
            if match and len(match) == n:
                return list(match)
        # Last resort: identity
        return list(range(n))


def from_rdkit(mol, use_cache=True):
    """Convert an RDKit Mol to SMSD MolGraph via SMILES round-trip.

    Returns a MolGraph. The SMSD→RDKit index mapping is stored internally
    and used by :func:`translate_mapping` and :func:`mcs_rdkit_native`.

    Each distinct RDKit Mol object gets its own mapping (no cross-Mol
    cache collisions). The cache is keyed by ``id(mol)`` and automatically
    evicted when the cache exceeds 10K entries.

    Example::

        from rdkit import Chem
        import smsd

        mol = Chem.MolFromSmiles("c1ccccc1")
        g = smsd.from_rdkit(mol)
    """
    try:
        from rdkit import Chem
    except ImportError:
        raise ImportError(
            "RDKit is required for from_rdkit(). "
            "Install with: pip install rdkit"
        )

    if mol is None:
        raise TypeError("mol must not be None")
    if mol.GetNumAtoms() == 0:
        raise ValueError("mol has 0 atoms; cannot convert empty molecule")

    # Check if molecule is only implicit/explicit H (no heavy atoms)
    heavy = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() != 1)
    if heavy == 0:
        raise ValueError("mol contains only hydrogen atoms; no heavy atoms to convert")

    mol_id = id(mol)
    if use_cache:
        with _rdkit_cache_lock:
            cached = _rdkit_cache.get(mol_id)
            if cached is not None:
                cached_mol, cached_result = cached
                if cached_mol is mol:
                    return cached_result.graph

    smi = Chem.MolToSmiles(mol)
    if not smi:
        raise ValueError("RDKit returned empty SMILES for the molecule")

    # Each RDKit Mol gets its OWN MolGraph instance to prevent index map
    # collisions.  Two Mols with the same canonical SMILES but different
    # atom orderings must have separate index maps, and get_index_map()
    # uses id(MolGraph) for reverse lookup.  Sharing a MolGraph across
    # Mols would cause the second Mol's index map to overwrite the first.
    g = parse_smiles(smi)

    # Compute deterministic index mapping for THIS specific Mol object
    index_map = _compute_index_map(mol, smi)
    result = _RdkitResult(g, index_map)

    with _rdkit_cache_lock:
        _graph_to_result[id(g)] = result  # safe: unique g per Mol

        if use_cache:
            if len(_rdkit_cache) > 10000:
                _rdkit_cache.clear()
                _graph_to_result.clear()
            _rdkit_cache[mol_id] = (mol, result)

    return g


_graph_to_result = {}  # id(MolGraph) → _RdkitResult (O(1) reverse lookup)


def get_index_map(g):
    """Return the SMSD→RDKit atom index mapping for a MolGraph from from_rdkit().

    Returns a list where ``result[smsd_idx] = rdkit_idx``, or None if the
    MolGraph was not created by from_rdkit().
    """
    r = _graph_to_result.get(id(g))
    return r.index_map if r is not None else None


def translate_mapping(mapping, g_query=None, g_target=None):
    """Translate an MCS mapping from SMSD indices to RDKit indices.

    Args:
        mapping: dict of {smsd_query_idx: smsd_target_idx} from find_mcs()
        g_query: the query MolGraph (from from_rdkit())
        g_target: the target MolGraph (from from_rdkit())

    Returns:
        dict of {rdkit_query_idx: rdkit_target_idx}
    """
    if mapping is None or len(mapping) == 0:
        return {}
    if not isinstance(mapping, dict):
        raise TypeError(f"mapping must be a dict, got {type(mapping).__name__}")

    q_map = get_index_map(g_query) if g_query is not None else None
    t_map = get_index_map(g_target) if g_target is not None else None

    result = {}
    for smsd_q, smsd_t in mapping.items():
        if not isinstance(smsd_q, int) or not isinstance(smsd_t, int):
            raise TypeError(f"mapping keys/values must be int, got ({type(smsd_q).__name__}, {type(smsd_t).__name__})")
        if smsd_q < 0 or smsd_t < 0:
            continue  # skip negative indices silently
        # Translate with bounds checking — skip (not crash) on out-of-bounds
        # which can happen with molecules containing explicit H atoms.
        if q_map is not None:
            if smsd_q >= len(q_map):
                continue  # SMSD index out of range for this index map
            rdkit_q = q_map[smsd_q]
        else:
            rdkit_q = smsd_q
        if t_map is not None:
            if smsd_t >= len(t_map):
                continue
            rdkit_t = t_map[smsd_t]
        else:
            rdkit_t = smsd_t
        result[rdkit_q] = rdkit_t
    return result


def mcs_rdkit_native(mol1, mol2, **kwargs):
    """Find MCS between two RDKit Mols, returning RDKit-compatible atom indices.

    Uses SMSD's 11-level MCS funnel for the search, then re-matches the MCS
    SMILES against both RDKit molecules using RDKit's own GetSubstructMatch
    to get correct native indices. This guarantees element-correct mappings.

    Returns a dict where keys are atom indices in mol1 and values are atom
    indices in mol2 — matching RDKit's atom ordering.

    Example::

        from rdkit import Chem
        import smsd

        m1 = Chem.MolFromSmiles("c1ccccc1")
        m2 = Chem.MolFromSmiles("c1ccc(O)cc1")
        mapping = smsd.mcs_rdkit_native(m1, m2)
        for q, t in mapping.items():
            assert m1.GetAtomWithIdx(q).GetAtomicNum() == m2.GetAtomWithIdx(t).GetAtomicNum()
    """
    if mol1 is None or mol2 is None:
        raise ValueError("mol1 and mol2 must not be None")
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        raise ImportError("RDKit is required for mcs_rdkit_native()")

    g1 = from_rdkit(mol1)
    g2 = from_rdkit(mol2)

    # Step 1: Find MCS using SMSD's full 11-level pipeline
    smsd_mapping = find_mcs(g1, g2, **kwargs)
    if not smsd_mapping:
        return {}

    # Step 2: Extract MCS as SMILES via SMSD
    try:
        mcs_smi = mcs_to_smiles(g1, smsd_mapping)
    except Exception:
        mcs_smi = ""

    if not mcs_smi:
        # Fallback: try index-based translation
        rdkit_mapping = translate_mapping(smsd_mapping, g1, g2)
        n1, n2 = mol1.GetNumAtoms(), mol2.GetNumAtoms()
        validated = {}
        for q, t in rdkit_mapping.items():
            if q < n1 and t < n2:
                if mol1.GetAtomWithIdx(q).GetAtomicNum() == mol2.GetAtomWithIdx(t).GetAtomicNum():
                    validated[q] = t
        return validated

    # Step 3: Parse MCS SMILES in RDKit and match against both molecules
    mcs_mol = Chem.MolFromSmarts(mcs_smi)
    if mcs_mol is None:
        # Try as SMILES if SMARTS parse fails
        mcs_mol = Chem.MolFromSmiles(mcs_smi)
    if mcs_mol is None:
        return {}

    match1 = mol1.GetSubstructMatch(mcs_mol)
    match2 = mol2.GetSubstructMatch(mcs_mol)

    if not match1 or not match2:
        # SMARTS/SMILES match failed — fall back to index translation
        rdkit_mapping = translate_mapping(smsd_mapping, g1, g2)
        n1, n2 = mol1.GetNumAtoms(), mol2.GetNumAtoms()
        validated = {}
        for q, t in rdkit_mapping.items():
            if q < n1 and t < n2:
                if mol1.GetAtomWithIdx(q).GetAtomicNum() == mol2.GetAtomWithIdx(t).GetAtomicNum():
                    validated[q] = t
        return validated

    # Step 4: Build RDKit→RDKit mapping via the MCS bridge
    # match1[mcs_idx] = mol1_atom_idx, match2[mcs_idx] = mol2_atom_idx
    result = {}
    for mcs_idx in range(min(len(match1), len(match2))):
        q_idx = match1[mcs_idx]
        t_idx = match2[mcs_idx]
        # Element validation (should always pass with proper SMARTS match)
        if mol1.GetAtomWithIdx(q_idx).GetAtomicNum() == mol2.GetAtomWithIdx(t_idx).GetAtomicNum():
            result[q_idx] = t_idx
    return result


def clear_cache():
    """Clear all from_rdkit() caches. Thread-safe."""
    with _rdkit_cache_lock:
        _rdkit_cache.clear()
        _mol_graph_cache.clear()
        _graph_to_result.clear()


def batch_mcs_rdkit(query_mol, target_mols, **kwargs):
    """Batch MCS with correct RDKit atom indices.

    Computes MCS between a single query and multiple targets using SMSD's
    MCS engine, returning atom-index mappings in RDKit's native atom ordering
    with full element validation.

    Args:
        query_mol: RDKit Mol object for the query.
        target_mols: List of RDKit Mol objects for the targets.
        **kwargs: Keyword arguments forwarded to find_mcs()
            (e.g. tautomer_aware, timeout_ms).

    Returns:
        list[dict]: One mapping per target. Each dict maps
            query RDKit atom index -> target RDKit atom index.
            Mismatched elements are excluded.

    Example::

        from rdkit import Chem
        import smsd

        query = Chem.MolFromSmiles("c1ccccc1")
        targets = [Chem.MolFromSmiles(s) for s in ["c1ccc(O)cc1", "c1ccncc1"]]
        results = smsd.batch_mcs_rdkit(query, targets)
        for mapping in results:
            for q, t in mapping.items():
                assert query.GetAtomWithIdx(q).GetAtomicNum() == targets[0].GetAtomWithIdx(t).GetAtomicNum()
    """
    g_query = from_rdkit(query_mol)
    results = []
    for t_mol in target_mols:
        g_target = from_rdkit(t_mol)
        smsd_mapping = find_mcs(g_query, g_target, **kwargs)
        rdkit_mapping = translate_mapping(smsd_mapping, g_query, g_target)
        # Element validation
        n1, n2 = query_mol.GetNumAtoms(), t_mol.GetNumAtoms()
        validated = {}
        for q, t in rdkit_mapping.items():
            if q < n1 and t < n2:
                if query_mol.GetAtomWithIdx(q).GetAtomicNum() == t_mol.GetAtomWithIdx(t).GetAtomicNum():
                    validated[q] = t
        results.append(validated)
    return results


def mcs_rdkit(mol1, mol2, **kwargs):
    """Find Maximum Common Substructure between two RDKit Mol objects.

    Accepts any keyword arguments supported by :func:`mcs` (e.g.
    ``tautomer_aware``, ``timeout_ms``).

    Example::

        from rdkit import Chem
        import smsd

        m1 = Chem.MolFromSmiles("c1ccccc1")
        m2 = Chem.MolFromSmiles("c1ccc(O)cc1")
        result = smsd.mcs_rdkit(m1, m2)
    """
    g1 = from_rdkit(mol1)
    g2 = from_rdkit(mol2)
    return mcs(g1, g2, **kwargs)


def substructure_rdkit(query_mol, target_mol, **kwargs):
    """Check substructure relationship between two RDKit Mol objects.

    Accepts any keyword arguments supported by :func:`substructure_search`
    (e.g. ``timeout_ms``).

    Example::

        from rdkit import Chem
        import smsd

        benzene = Chem.MolFromSmiles("c1ccccc1")
        phenol  = Chem.MolFromSmiles("c1ccc(O)cc1")
        result = smsd.substructure_rdkit(benzene, phenol)
    """
    gq = from_rdkit(query_mol)
    gt = from_rdkit(target_mol)
    return substructure_search(gq, gt, **kwargs)


# =========================================================================
# Export & Depiction — plug into RDKit, matplotlib, or Jupyter
# =========================================================================


def to_rdkit(mol_or_smiles):
    """Convert SMSD MolGraph (or SMILES string) to an RDKit Mol object.

    Example::

        import smsd
        g = smsd.parse_smiles("c1ccccc1")
        rdkit_mol = smsd.to_rdkit(g)

        # Now use RDKit for drawing, fingerprints, descriptors, etc.
        from rdkit.Chem import Draw
        Draw.MolToImage(rdkit_mol)
    """
    try:
        from rdkit import Chem
    except ImportError:
        raise ImportError(
            "RDKit is required for to_rdkit(). "
            "Install with: pip install rdkit"
        )
    if isinstance(mol_or_smiles, str):
        smi = mol_or_smiles
    elif hasattr(mol_or_smiles, "__class__") and "MolGraph" in type(mol_or_smiles).__name__:
        smi = to_smiles(mol_or_smiles)
    else:
        raise TypeError(f"Expected MolGraph or SMILES string, got {type(mol_or_smiles)}")
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smi}")
    return mol


def depict_mcs(mol1, mol2, *, tautomer_aware=False, timeout_ms=10000,
               size=(600, 300), highlight_color=(0.5, 0.8, 1.0)):
    """Find MCS and return an RDKit image with matched atoms highlighted.

    Works in Jupyter notebooks (auto-displays) and scripts (returns PIL Image).

    Example::

        import smsd
        img = smsd.depict_mcs("c1ccccc1", "c1ccc(O)cc1")
        img.save("mcs.png")  # or just display in Jupyter
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, AllChem
    except ImportError:
        raise ImportError(
            "RDKit is required for depict_mcs(). "
            "Install with: pip install rdkit"
        )

    # Parse inputs
    g1 = _ensure_mol(mol1)
    g2 = _ensure_mol(mol2)

    # Find MCS
    mapping = mcs(g1, g2, tautomer_aware=tautomer_aware, timeout_ms=timeout_ms)

    # Convert to RDKit mols for drawing
    smi1 = to_smiles(g1) if not isinstance(mol1, str) else mol1
    smi2 = to_smiles(g2) if not isinstance(mol2, str) else mol2
    rdmol1 = Chem.MolFromSmiles(smi1)
    rdmol2 = Chem.MolFromSmiles(smi2)
    if rdmol1 is None or rdmol2 is None:
        raise ValueError("Could not convert molecules for depiction")

    AllChem.Compute2DCoords(rdmol1)
    AllChem.Compute2DCoords(rdmol2)

    # Highlight matched atoms
    q_atoms = list(mapping.keys()) if mapping else []
    t_atoms = list(mapping.values()) if mapping else []

    img = Draw.MolsToGridImage(
        [rdmol1, rdmol2],
        molsPerRow=2,
        subImgSize=(size[0] // 2, size[1]),
        highlightAtomLists=[q_atoms, t_atoms],
        legends=[f"Query ({len(q_atoms)} matched)", f"Target ({len(t_atoms)} matched)"],
    )
    return img


def depict_substructure(query, target, *, timeout_ms=10000,
                        size=(600, 300)):
    """Find substructure and return an RDKit image with matched atoms highlighted.

    Example::

        import smsd
        img = smsd.depict_substructure("c1ccccc1", "c1ccc(O)cc1")
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, AllChem
    except ImportError:
        raise ImportError(
            "RDKit is required for depict_substructure(). "
            "Install with: pip install rdkit"
        )

    gq = _ensure_mol(query)
    gt = _ensure_mol(target)

    result = substructure_search(gq, gt, timeout_ms=timeout_ms)
    mapping = result if isinstance(result, dict) else {}

    smi_q = to_smiles(gq) if not isinstance(query, str) else query
    smi_t = to_smiles(gt) if not isinstance(target, str) else target
    rdmol_q = Chem.MolFromSmiles(smi_q)
    rdmol_t = Chem.MolFromSmiles(smi_t)
    if rdmol_q is None or rdmol_t is None:
        raise ValueError("Could not convert molecules for depiction")

    AllChem.Compute2DCoords(rdmol_q)
    AllChem.Compute2DCoords(rdmol_t)

    q_atoms = list(mapping.keys()) if mapping else []
    t_atoms = list(mapping.values()) if mapping else []

    img = Draw.MolsToGridImage(
        [rdmol_q, rdmol_t],
        molsPerRow=2,
        subImgSize=(size[0] // 2, size[1]),
        highlightAtomLists=[q_atoms, t_atoms],
        legends=[f"Query", f"Target ({len(t_atoms)} matched)"],
    )
    return img


def to_svg(mol_or_smiles, *, size=(300, 200)):
    """Generate SVG depiction of a molecule using RDKit.

    Example::

        svg = smsd.to_svg("c1ccccc1")
        with open("benzene.svg", "w") as f:
            f.write(svg)
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, AllChem, rdMolDraw2D
    except ImportError:
        raise ImportError(
            "RDKit is required for to_svg(). "
            "Install with: pip install rdkit"
        )

    if isinstance(mol_or_smiles, str):
        mol = Chem.MolFromSmiles(mol_or_smiles)
    elif hasattr(mol_or_smiles, "__class__") and "MolGraph" in type(mol_or_smiles).__name__:
        mol = Chem.MolFromSmiles(to_smiles(mol_or_smiles))
    else:
        mol = mol_or_smiles  # assume RDKit mol

    if mol is None:
        raise ValueError("Could not parse molecule for SVG depiction")

    AllChem.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def export_sdf(molecules, filename):
    """Export molecules to an SDF file using the native SMSD writer.

    Example::

        mols = [smsd.parse_smiles(s) for s in ["CCO", "c1ccccc1", "CC(=O)O"]]
        smsd.export_sdf(mols, "output.sdf")
    """
    path = Path(filename)
    records = [write_sdf_record(_ensure_mol(mol)) for mol in molecules]
    path.write_text("".join(records), encoding="utf-8")


def read_molfile(filename):
    """Read a MOL/SDF record from disk into a MolGraph."""
    return read_mol_block(Path(filename).read_text(encoding="utf-8"))


def write_molfile(mol, filename, *, v3000=False, sdf=False):
    """Write a MolGraph to disk as MOL V2000, MOL V3000, or a single SDF record."""
    g = _ensure_mol(mol)
    if sdf:
        text = write_sdf_record(g)
    elif v3000:
        text = write_mol_block_v3000(g)
    else:
        text = write_mol_block(g)
    Path(filename).write_text(text, encoding="utf-8")


def from_openbabel(mol):
    """Convert an Open Babel molecule into an SMSD MolGraph."""
    try:
        from openbabel import openbabel as ob  # type: ignore[import]
        from openbabel import pybel  # type: ignore[import]
    except ImportError:
        try:
            import pybel  # type: ignore[import]
            ob = None
        except ImportError as exc:  # pragma: no cover - optional dependency
            raise ImportError(
                "Open Babel is required for from_openbabel(). "
                "Install openbabel-wheel or the system Open Babel bindings."
            ) from exc

    if hasattr(mol, "write"):
        smi = mol.write("smi").strip().splitlines()[0].split()[0]
        return parse_smiles(smi)

    if ob is not None and hasattr(mol, "NumAtoms"):
        conv = ob.OBConversion()
        conv.SetOutFormat("smi")
        smi = conv.WriteString(mol).strip().splitlines()[0].split()[0]
        return parse_smiles(smi)

    raise TypeError(
        f"Expected pybel.Molecule or openbabel.OBMol, got {type(mol).__name__}"
    )


def to_openbabel(mol):
    """Convert a MolGraph, SMILES string, or RDKit molecule into a pybel Molecule."""
    try:
        from openbabel import pybel  # type: ignore[import]
    except ImportError:
        try:
            import pybel  # type: ignore[import]
        except ImportError as exc:  # pragma: no cover - optional dependency
            raise ImportError(
                "Open Babel is required for to_openbabel(). "
                "Install openbabel-wheel or the system Open Babel bindings."
            ) from exc

    g = _ensure_mol(mol)
    return pybel.readstring("smi", to_smiles(g))


__all__ = [
    # Types
    "ChemOptions",
    "McsOptions",
    "MolGraph",
    "MolGraphBuilder",
    # Enums
    "BondOrderMode",
    "AromaticityMode",
    "RingFusionMode",
    # SMILES / SMARTS
    "parse_smiles",
    "read_mol_block",
    "read_molfile",
    "parse_mapped_smiles",
    "strip_atom_maps",
    "to_smiles",
    "to_smarts",
    "write_mol_block",
    "write_mol_block_v3000",
    "write_sdf_record",
    "write_molfile",
    # Substructure
    "is_substructure",
    "find_substructure",
    "substructure_search",
    # MCS
    "find_mcs",
    "find_all_mcs",
    "mcs",
    "all_mcs",
    # Index translation
    "translate_to_atom_ids",
    "find_mcs_progressive",
    # Reaction-aware MCS
    "map_reaction_aware",
    # Bond-change scoring (v6.6.0)
    "bond_change_score",
    # Batch constrained MCS (v6.6.0)
    "batch_mcs_constrained",
    # MCS SMILES extraction
    "mcs_to_smiles",
    "find_mcs_smiles",
    "mcs_smiles",
    # Similarity
    "similarity",
    "similarity_upper_bound",
    "screen_targets",
    # Fingerprints
    "fingerprint",
    "path_fingerprint",
    "mcs_fingerprint",
    "fingerprint_subset",
    "analyze_fp_quality",
    "tanimoto",
    "topological_torsion",
    # Structured result
    "MatchResult",
    "mcs_result",
    # Convenience SMILES-in APIs
    "mcs_from_smiles",
    "fingerprint_from_smiles",
    "counts_to_array",
    "to_hex",
    "from_hex",
    "to_binary_string",
    "batch_similarity",
    # R-Group decomposition
    "decompose_r_groups",
    # Parallel batch
    "batch_substructure",
    "batch_mcs",
    "batch_mcs_rdkit",
    # GPU / compute backend
    "gpu_is_available",
    "gpu_device_info",
    # RDKit interoperability
    "from_rdkit",
    "from_openbabel",
    "get_index_map",
    "translate_mapping",
    "mcs_rdkit_native",
    "clear_cache",
    "to_openbabel",
    "to_rdkit",
    "mcs_rdkit",
    "substructure_rdkit",
    # Export & depiction (RDKit-powered)
    "depict_mcs",
    "depict_substructure",
    "to_svg",
    "export_sdf",
    # SMARTS matching
    "find_mcs_smarts",
    "smarts_match",
    "smarts_find_all",
    # CIP stereo descriptors
    "RSLabel",
    "EZLabel",
    "CIPDescriptors",
    "assign_rs",
    "assign_ez",
    "assign_cip",
    # Ring perception
    "compute_sssr",
    "layout_sssr",
]
