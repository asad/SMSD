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

    import smsd

    # Substructure check — accepts SMILES strings directly
    assert smsd.is_substructure("c1ccccc1", "c1ccc(O)cc1")

    # MCS — returns {query_atom: target_atom} mapping
    mapping = smsd.find_mcs("c1ccc(O)cc1", "c1ccc(N)cc1")
    print(f"MCS size: {len(mapping)}")

    # Multiple MCS mappings for SAR analysis
    mappings = smsd.find_mcs("c1ccc(O)cc1", "c1ccc(N)cc1", max_results=10)

    # Substructure with atom mapping
    mapping = smsd.find_substructure("c1ccccc1", "c1ccc(O)cc1")
    if mapping:
        print("Benzene found in phenol")
"""

import importlib.machinery
import importlib.util
import sys
from pathlib import Path

__version__ = "7.2.1"
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
MCSOptions = _smsd.MCSOptions
MolGraph = _smsd.MolGraph
MolGraphBuilder = _smsd.MolGraphBuilder

BondOrderMode = _smsd.BondOrderMode
AromaticityMode = _smsd.AromaticityMode
AromaticityModel = _smsd.AromaticityModel
RingFusionMode = _smsd.RingFusionMode
MatcherEngine = _smsd.MatcherEngine
Solvent = _smsd.Solvent

_parse_smiles_raw = _smsd.parse_smiles
read_mol_block = _smsd.read_mol_block
_to_smiles_native = _smsd.to_smiles
to_smarts = _smsd.to_smarts


def to_smiles(mol):
    """Write a canonical SMILES string for a molecule.

    Accepts either a :class:`MolGraph` or a SMILES string (which is
    parsed to a :class:`MolGraph` first).  The output is the
    canonical SMILES defined by SMSD's Morgan-based ordering, and is
    byte-identical across the Java, C++, and Python engines.

    .. versionchanged:: 7.1.1
       Now accepts a SMILES string in addition to :class:`MolGraph`.
       Previously raised ``TypeError: incompatible function arguments``
       because the wrapper was aliased directly to a pybind11 overload
       that only accepted :class:`MolGraph`.
    """
    if isinstance(mol, str):
        mol = _parse_smiles_raw(mol)
    return _to_smiles_native(mol)


def canonical_smiles(mol):
    """Stereo-aware canonical SMILES, accepting a SMILES string or :class:`MolGraph`.

    Thin alias of :func:`to_smiles`.  See :func:`to_smiles` for details.

    .. versionchanged:: 7.1.1
       Now accepts a SMILES string in addition to :class:`MolGraph`.
    """
    return to_smiles(mol)


def canonical_hash(mol):
    """Return the canonical hash of a molecule as a Python int.

    The hash encodes connectivity, bond order, aromaticity, ring membership,
    formal charge, mass number, and stereo (tetrahedral chirality and E/Z
    double-bond configuration). Configurable stereo inclusion is controlled
    by the stereo data stored in the MolGraph itself.

    Hashes are compact screening keys. Collisions and bounded canonical
    labeling mean they are not proofs of molecular identity or difference;
    use a validated graph mapping when identity matters.

    Example::

        >>> import smsd
        >>> r = smsd.parse_smiles("N[C@@H](C)C(=O)O")   # L-alanine
        >>> s = smsd.parse_smiles("N[C@H](C)C(=O)O")    # D-alanine
        >>> smsd.canonical_hash(r) == smsd.canonical_hash(s)
        False
    """
    return mol.get_canonical_hash()
write_mol_block = _smsd.write_mol_block
write_mol_block_v3000 = _smsd.write_mol_block_v3000
write_sdf_record = _smsd.write_sdf_record

# Low-level C++ bindings (internal — use find_mcs / find_substructure instead)
_native_is_substructure = _smsd.is_substructure
_native_find_substructure = _smsd.find_substructure
_native_find_mcs = _smsd.find_mcs
_native_find_all_mcs = _smsd.find_all_mcs
canonicalize_mapping = _smsd.canonicalize_mapping
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

dice = _smsd.dice
cosine = _smsd.cosine
soergel = _smsd.soergel
count_dice = _smsd.count_dice
count_cosine = _smsd.count_cosine
# v7.1.1: `smsd.overlapCoefficient` and `smsd.count_overlapCoefficient`
# are now redirected to the correct Szymkiewicz-Simpson implementation
# in `overlap_coefficient` / `count_overlap_coefficient` (defined below).
# The underlying C++ binding with the same name actually computed
# Tanimoto, matching the `MatchResult.overlapCoefficient` semantic bug
# in v7.1.0.  The rebound aliases are declared after the pure-Python
# helpers so they pick up the corrected implementation.

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

# -----------------------------------------------------------------------
# Extended API — v6.11.0
# -----------------------------------------------------------------------

# Multi-molecule MCS
find_nmcs = _smsd.find_nmcs
find_scaffold_mcs = _smsd.find_scaffold_mcs

# Mapping validation & maximality
validate_mapping = _smsd.validate_mapping
is_mapping_maximal = _smsd.is_mapping_maximal

# R-group decomposition
decompose_rgroups = _smsd.decompose_rgroups

# Graph utilities
extract_subgraph = _smsd.extract_subgraph
murcko_scaffold = _smsd.murcko_scaffold
count_components = _smsd.count_components
split_components = _smsd.split_components
same_canonical_graph = _smsd.same_canonical_graph

# All substructure matches
find_all_substructures = _smsd.find_all_substructures

# File I/O
read_mol_file = _smsd.read_mol_file
read_sdf = _smsd.read_sdf
write_sdf = _smsd.write_sdf

# Chemistry utilities
classify_pharmacophore = _smsd.classify_pharmacophore
implicit_h = _smsd.implicit_h

# Layout utilities
count_crossings = _smsd.count_crossings
is_degenerate_layout = _smsd.is_degenerate_layout

# Fingerprint utilities
counts_to_array = _smsd.counts_to_array
prewarm_graph = _smsd.prewarm_graph

# Layout engine v6.11.0
Point3D = _smsd.Point3D
generate_coords_2d = _smsd.generate_coords_2d
generate_coords_3d = _smsd.generate_coords_3d
layout_quality = _smsd.layout_quality
resolve_overlaps = _smsd.resolve_overlaps

# Coordinate transforms
translate_2d = _smsd.translate_2d
rotate_2d = _smsd.rotate_2d
scale_2d = _smsd.scale_2d
mirror_x = _smsd.mirror_x
mirror_y = _smsd.mirror_y
center_2d = _smsd.center_2d
align_2d = _smsd.align_2d
bounding_box_2d = _smsd.bounding_box_2d
normalise_bond_length = _smsd.normalise_bond_length
canonical_orientation = _smsd.canonical_orientation

# ── Depiction Engine v6.11.0 ──────────────────────────────────────────────
# SVG rendering with configurable bond and atom-label styling.
DepictOptions = _smsd.DepictOptions


def depict_svg(mol_or_smiles, opts=None, **kwargs):
    """Render a molecule as an SVG string.

    Parameters
    ----------
    mol_or_smiles : MolGraph or str
        Molecule to render (accepts MolGraph, SMILES string, or RDKit Mol).
    opts : DepictOptions, optional
        Full options object. If None, uses the default depiction settings.
    **kwargs
        Shorthand overrides applied to opts (e.g., bond_length=40, width=600).

    Returns
    -------
    str
        SVG string suitable for file output or Jupyter display.
    """
    if opts is None:
        opts = DepictOptions()
    for k, v in kwargs.items():
        setattr(opts, k, v)
    if isinstance(mol_or_smiles, str):
        return _smsd.depict_smiles(mol_or_smiles, opts)
    mol = _ensure_mol(mol_or_smiles) if not isinstance(mol_or_smiles, _smsd.MolGraph) else mol_or_smiles
    return _smsd.depict(mol, opts)


def depict_mapping(mol_or_smiles, mapping, opts=None, **kwargs):
    """Render a molecule with MCS/substructure atoms highlighted and numbered.

    Parameters
    ----------
    mol_or_smiles : MolGraph or str
        Molecule to render.
    mapping : dict
        Atom index -> partner atom index (from MCS or substructure result).
    opts : DepictOptions, optional
    **kwargs
        Shorthand overrides (e.g., show_atom_indices=True).

    Returns
    -------
    str
        SVG string with matched atoms highlighted in green.
    """
    if opts is None:
        opts = DepictOptions()
    for k, v in kwargs.items():
        setattr(opts, k, v)
    if isinstance(mol_or_smiles, str):
        mol_or_smiles = parse_smiles(mol_or_smiles)
    mol = _ensure_mol(mol_or_smiles) if not isinstance(mol_or_smiles, _smsd.MolGraph) else mol_or_smiles
    return _smsd.depict_with_mapping(mol, mapping, opts)


def depict_pair(mol1, mol2, mapping, opts=None, **kwargs):
    """Render two molecules side-by-side with MCS mapping highlighted.

    Parameters
    ----------
    mol1, mol2 : MolGraph or str
        Molecules to compare.
    mapping : dict
        mol1 atom -> mol2 atom mapping (from find_mcs or find_all_mcs).
    opts : DepictOptions, optional
    **kwargs
        Shorthand overrides (e.g., width=800, height=400).

    Returns
    -------
    str
        SVG string showing both molecules with matched regions highlighted.
    """
    if opts is None:
        opts = DepictOptions()
    for k, v in kwargs.items():
        setattr(opts, k, v)
    if isinstance(mol1, str):
        mol1 = parse_smiles(mol1)
    if isinstance(mol2, str):
        mol2 = parse_smiles(mol2)
    m1 = _ensure_mol(mol1) if not isinstance(mol1, _smsd.MolGraph) else mol1
    m2 = _ensure_mol(mol2) if not isinstance(mol2, _smsd.MolGraph) else mol2
    return _smsd.depict_pair(m1, m2, mapping, opts)


def save_svg(svg_str, path):
    """Write an SVG string to a file.

    Parameters
    ----------
    svg_str : str
        SVG content (from depict_svg, depict_pair, etc.).
    path : str
        Output file path (e.g., 'molecule.svg').
    """
    with open(path, 'w') as f:
        f.write(svg_str)


# Bond-change scoring + batch constrained MCS (v6.6.0)
try:
    from smsd._smsd import (
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
    try:
        from smsd._smsd import (
            SmartsQuery,
            compile_smarts as _compile_smarts_raw,
        )
    except ImportError:
        SmartsQuery = None
        _compile_smarts_raw = None
    _HAS_SMARTS = True
except ImportError:
    SmartsQuery = None
    _compile_smarts_raw = None
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


def find_mcs(mol1, mol2, *,
             max_results=1,
             # Chemistry matching flags
             ring_matches_ring_only=False,
             complete_rings_only=False,
             match_bond_order="strict",
             match_atom_type=True,
             match_formal_charge=False,
             match_isotope=False,
             use_chirality=False,
             use_bond_stereo=False,
             tautomer_aware=False,
             # MCS search options
             timeout_ms=10000,
             connected_only=True,
             induced=False,
             maximize_bonds=False,
             max_stage=5,
             # Strategy
             strategy="auto",
             prefer_rare_heteroatoms=False,
             # Advanced
             **kwargs):
    """Find Maximum Common Substructure between two molecules.

    Unified MCS entry point — accepts SMILES strings, MolGraph objects, or
    RDKit Mol objects. Uses the fast lightweight engine by default with
    automatic fallback to the full native C++ pipeline.

    Args:
        mol1: First molecule (MolGraph, SMILES string, or RDKit Mol).
        mol2: Second molecule.
        max_results: Number of distinct MCS mappings to return.
            When 1 (default), returns a single ``dict``.
            When > 1, returns a ``list[dict]`` of up to *max_results*
            mappings that all share the maximum MCS size.
        ring_matches_ring_only: Ring atom only matches ring atom (default False).
        complete_rings_only: MCS must include full rings (default False).
        match_bond_order: "strict" (exact), "loose" (ignore bond order),
                          "any" (ignore bond order). Default "strict".
        match_atom_type: Match element type (default True).
        match_formal_charge: Match formal charge (default False).
        match_isotope: Match mass number (default False).
        use_chirality: Consider R/S stereochemistry (default False).
        use_bond_stereo: Consider E/Z stereochemistry (default False).
        tautomer_aware: Use tautomer-aware matching (default False).
        timeout_ms: Timeout in milliseconds (default 10000).
        connected_only: Connected MCS only (default True).
        induced: Induced subgraph MCS (default False).
        maximize_bonds: Edge MCS / MCES (default False).
        max_stage: Internal effort knob 0-5 (default 5 — highest effort).
            Advanced tuning only; intended for callers willing to trade
            MCS size for speed.
        strategy: "auto", "lightweight", or "native". "auto" uses lightweight
            for a single result with default chemistry/search settings, plus
            timeout, ring matching, and "strict"/"any" bond order. Other requested
            options use native. "lightweight" raises ValueError for unsupported
            options, including "loose" bond order and advanced kwargs.
        prefer_rare_heteroatoms: Accepted for API stability; currently a no-op.

    Returns:
        dict when max_results=1, list[dict] when max_results > 1.
        Each dict maps mol1 atom indices to mol2 atom indices.

    Examples::

        # Single best MCS (default)
        mapping = smsd.find_mcs("c1ccc(O)cc1", "c1ccc(N)cc1")

        # Multiple MCS mappings for SAR analysis
        mappings = smsd.find_mcs(mol1, mol2, max_results=10)

        # Strict matching
        mapping = smsd.find_mcs(mol1, mol2,
                                ring_matches_ring_only=True,
                                match_formal_charge=True)

        # Minimum-effort mode (faster, possibly smaller MCS)
        mapping = smsd.find_mcs(mol1, mol2, max_stage=1)
    """
    # prefer_rare_heteroatoms is accepted for API stability but currently
    # has no effect; retained so existing callers do not need to be edited.
    del prefer_rare_heteroatoms

    if strategy not in ("auto", "lightweight", "native"):
        raise ValueError("strategy must be 'auto', 'lightweight', or 'native'")
    if not isinstance(match_bond_order, str) or match_bond_order.lower() not in ("strict", "loose", "any"):
        raise ValueError("match_bond_order must be 'strict', 'loose', or 'any'")
    bond_mode = match_bond_order.lower()
    bond_any = bond_mode == "any"

    # The coverage wrapper accepts only timeout, ring matching, and a boolean
    # bond-order mode. Never let its early return bypass requested constraints.
    lightweight_unsupported = [
        name for name, unsupported in (
            ("max_results", max_results != 1),
            ("complete_rings_only", complete_rings_only),
            ("match_bond_order", bond_mode == "loose"),
            ("match_atom_type", not match_atom_type),
            ("match_formal_charge", match_formal_charge),
            ("match_isotope", match_isotope),
            ("use_chirality", use_chirality),
            ("use_bond_stereo", use_bond_stereo),
            ("tautomer_aware", tautomer_aware),
            ("connected_only", not connected_only),
            ("induced", induced),
            ("maximize_bonds", maximize_bonds),
            ("max_stage", max_stage != 5),
        ) if unsupported
    ]
    lightweight_unsupported.extend(sorted(kwargs))
    if strategy == "lightweight" and lightweight_unsupported:
        raise ValueError(
            "strategy='lightweight' does not support requested options: "
            + ", ".join(lightweight_unsupported)
            + "; use strategy='native' or 'auto'"
        )

    # Build ChemOptions / MCSOptions for native path
    def _build_opts():
        chem = ChemOptions.tautomer_profile() if tautomer_aware else ChemOptions()
        chem.match_atom_type = match_atom_type
        chem.match_formal_charge = match_formal_charge
        chem.match_isotope = match_isotope
        chem.use_chirality = use_chirality
        chem.use_bond_stereo = use_bond_stereo
        chem.ring_matches_ring_only = ring_matches_ring_only
        chem.complete_rings_only = complete_rings_only
        bond_mode_map = {"strict": BondOrderMode.STRICT, "loose": BondOrderMode.LOOSE, "any": BondOrderMode.ANY}
        chem.match_bond_order = bond_mode_map[bond_mode]
        opts = MCSOptions()
        opts.timeout_ms = timeout_ms
        opts.connected_only = connected_only
        opts.induced = induced
        opts.maximize_bonds = maximize_bonds
        opts.max_stage = max_stage
        for k, v in kwargs.items():
            setattr(opts, k, v)
        return chem, opts

    if not isinstance(max_results, int) or isinstance(max_results, bool) or max_results < 1:
        raise ValueError("max_results must be a positive integer")

    # --- Multi-result path (max_results > 1) → native find_all_mcs ---
    if max_results > 1:
        g1, rdkit1 = _ensure_mol_ex(mol1)
        g2, rdkit2 = _ensure_mol_ex(mol2)
        chem, opts = _build_opts()
        _remap_query_weights(opts, g1, rdkit1)
        return [_auto_translate(mapping, g1, g2, rdkit1, rdkit2)[0]
                for mapping in _native_find_all_mcs(g1, g2, chem, opts, max_results)]

    # --- Single-result path (max_results == 1) ---
    # Coverage search for the supported default settings
    import time
    started = time.monotonic()
    light_mapping = {}
    if strategy in ("auto", "lightweight") and not lightweight_unsupported:
        try:
            result = _find_mcs_light(
                mol1, mol2,
                timeout=timeout_ms / 1000.0,
                ring_matches_ring=ring_matches_ring_only,
                bond_any=bond_any,
            )
            if result.size > 0:
                light_mapping = {a - 1: b - 1 for a, b in result.mapping}
                if strategy == "lightweight":
                    return light_mapping
                # In auto: for dot-disconnected (salts), also try native
                smi1 = mol1 if isinstance(mol1, str) else ""
                smi2 = mol2 if isinstance(mol2, str) else ""
                if '.' not in smi1 and '.' not in smi2:
                    return light_mapping
        except Exception:
            if strategy == "lightweight":
                return {}

    # Native C++ pipeline — for advanced options, salts, or explicit "native".
    g1, rdkit1 = _ensure_mol_ex(mol1)
    g2, rdkit2 = _ensure_mol_ex(mol2)
    chem, opts = _build_opts()
    _remap_query_weights(opts, g1, rdkit1)
    if strategy == "auto" and not lightweight_unsupported and timeout_ms > 0:
        remaining = timeout_ms - int((time.monotonic() - started) * 1000)
        if remaining <= 0:
            return light_mapping
        opts.timeout_ms = remaining
    mapping = _native_find_mcs(g1, g2, chem, opts)
    translated, _ = _auto_translate(mapping, g1, g2, rdkit1, rdkit2)
    if strategy == "auto" and len(light_mapping) > len(translated):
        return light_mapping
    return translated




def find_mcs_progressive(mol1, mol2, *, on_progress=None, tautomer_aware=False,
                         timeout_ms=10000, **kwargs):
    """Find MCS with progressive intermediate reporting via callback.

    The search uses one deadline across its native pipeline. The current
    implementation calls on_progress once with the final mapping, its size,
    and elapsed wall-clock time in milliseconds.

    Args:
        mol1: First molecule (MolGraph or SMILES string).
        mol2: Second molecule (MolGraph or SMILES string).
        on_progress: Callback ``(best_mapping, best_size, elapsed_ms) -> None``.
            Called with the final mapping. May be None (no reporting).
        tautomer_aware: Use tautomer-aware matching.
        timeout_ms: Total timeout in milliseconds.
        **kwargs: Additional MCSOptions fields (e.g. connected_only=False).

    Returns:
        dict: Final best mapping from mol1 atom indices to mol2 atom indices.

    Example::

        def report(mapping, size, ms):
            print(f"  {ms:>5d} ms  MCS size = {size}")

        mapping = smsd.find_mcs_progressive("c1ccccc1O", "c1ccccc1N",
                                            on_progress=report)
    """
    import time

    g1, rdkit1 = _ensure_mol_ex(mol1)
    g2, rdkit2 = _ensure_mol_ex(mol2)
    chem = ChemOptions.tautomer_profile() if tautomer_aware else ChemOptions()
    chem, opts = _search_options(timeout_ms, dict(kwargs, tautomer_aware=tautomer_aware))
    _remap_query_weights(opts, g1, rdkit1)

    if on_progress is None:
        result = _native_find_mcs(g1, g2, chem, opts)
        return _auto_translate(result, g1, g2, rdkit1, rdkit2)[0]

    last = None
    t0 = time.monotonic()

    def report(mapping, size, elapsed_ms):
        nonlocal last
        translated = _auto_translate(mapping, g1, g2, rdkit1, rdkit2)[0]
        last = dict(translated)
        on_progress(translated, len(translated), elapsed_ms)

    result = _smsd.find_mcs_progressive(g1, g2, chem, opts, report)
    translated = _auto_translate(result, g1, g2, rdkit1, rdkit2)[0]
    if last != translated:
        on_progress(translated, len(translated), int((time.monotonic() - t0) * 1000))
    return translated


def mcs_smiles(mol1, mol2, *, tautomer_aware=False, timeout_ms=10000, **kwargs):
    """Find MCS between two molecules and return it as a canonical SMILES string.

    Accepts either MolGraph objects or SMILES strings.

    Args:
        mol1: First molecule (MolGraph or SMILES string).
        mol2: Second molecule (MolGraph or SMILES string).
        tautomer_aware: Use tautomer-aware matching.
        timeout_ms: Timeout in milliseconds.
        **kwargs: Additional MCSOptions fields (e.g. connected_only=False).

    Returns:
        str: Canonical SMILES of the MCS, or "" if no common substructure.
    """
    g1 = _ensure_mol(mol1)
    g2 = _ensure_mol(mol2)
    chem = ChemOptions.tautomer_profile() if tautomer_aware else ChemOptions()
    opts = MCSOptions()
    opts.timeout_ms = timeout_ms
    for k, v in kwargs.items():
        setattr(opts, k, v)
    return find_mcs_smiles(g1, g2, chem, opts, timeout_ms)


def find_substructure(query, target, *, max_results=1, timeout_ms=10000):
    """Find substructure mapping(s) of *query* in *target*.

    Unified substructure entry point — accepts SMILES strings, MolGraph
    objects, or RDKit Mol objects.

    Args:
        query: Query molecule (MolGraph, SMILES string, or RDKit Mol).
        target: Target molecule.
        max_results: Number of distinct mappings to return.
            When 1 (default), returns a single ``dict`` (empty if no match).
            When > 1, returns a ``list[dict]`` of up to *max_results*
            distinct embeddings.
        timeout_ms: Timeout in milliseconds (default 10000).

    Returns:
        dict when max_results=1, list[dict] when max_results > 1.
        Each dict maps query atom indices to target atom indices.

    Examples::

        # Single mapping (boolean-like: truthy if match found)
        mapping = smsd.find_substructure("c1ccccc1", "c1ccc(O)cc1")
        if mapping:
            print("Benzene is a substructure of phenol")

        # All distinct embeddings
        mappings = smsd.find_substructure("c1ccccc1", target, max_results=100)
    """
    if not isinstance(max_results, int) or isinstance(max_results, bool) or max_results < 1:
        raise ValueError("max_results must be a positive integer")
    q, rdkit_q = _ensure_mol_ex(query)
    t, rdkit_t = _ensure_mol_ex(target)
    if max_results > 1:
        all_maps = find_all_substructures(q, t, ChemOptions(), timeout_ms)
        return [_auto_translate(dict(mapping), q, t, rdkit_q, rdkit_t)[0]
                for mapping in all_maps[:max_results]]
    mapping = dict(_native_find_substructure(q, t, ChemOptions(), timeout_ms))
    return _auto_translate(mapping, q, t, rdkit_q, rdkit_t)[0]


def is_substructure(query, target, *, timeout_ms=10000):
    """Check whether *query* is a substructure of *target*.

    Convenience boolean wrapper around :func:`find_substructure`.

    Args:
        query: Query molecule (MolGraph, SMILES string, or RDKit Mol).
        target: Target molecule.
        timeout_ms: Timeout in milliseconds (default 10000).

    Returns:
        bool: True if query is a substructure of target.
    """
    q = _ensure_mol(query)
    t = _ensure_mol(target)
    return _native_is_substructure(q, t, ChemOptions(), timeout_ms)




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
    elif kind in ("ecfp", "circular"):
        return circular_fingerprint(g, radius=2, fp_size=fp_size)
    elif kind == "fcfp":
        return circular_fingerprint(g, radius=2, fp_size=fp_size, mode="fcfp")
    elif kind == "torsion":
        return topological_torsion(g, fp_size=fp_size)
    else:
        raise ValueError(
            f"Unknown fingerprint kind: {kind!r}. "
            f"Use 'path', 'mcs', 'ecfp', 'circular', 'fcfp', or 'torsion'."
        )


def _fp_to_positions(fp):
    """Convert fingerprint to a set of bit positions.

    Accepts:
    - list[int]: set bit positions (returned as-is)
    - list[tuple[int, int]]: sparse counts [(pos, count), ...] → extract positions
    - None: returns empty set
    """
    if fp is None:
        return set()
    if fp and isinstance(fp[0], tuple):
        return {pos for pos, _count in fp}
    return set(fp)


def overlap_coefficient(fp1, fp2, fp_size=2048):
    """Overlap coefficient (Szymkiewicz-Simpson): |A∩B| / min(|A|, |B|).

    Measures how much the smaller fingerprint is contained in the larger.
    Returns 1.0 when A is a subset of B (or vice versa).
    Useful for substructure-like similarity where size difference is expected.

    Accepts binary fingerprints (list of set bit positions) or sparse count
    fingerprints (list of (position, count) tuples).

    The empty-vs-empty corner case returns 1.0 (trivially identical
    empty sets), matching the convention used elsewhere in the library.
    If exactly one fingerprint is empty, the result is 0.0.

    Args:
        fp1: Fingerprint as list of set bit positions or sparse counts.
        fp2: Fingerprint as list of set bit positions or sparse counts.
        fp_size: Fingerprint size in bits (used when fp is bit positions).

    Returns:
        float: Overlap coefficient in [0.0, 1.0].
    """
    if fp1 is None or fp2 is None:
        return 0.0
    s1 = _fp_to_positions(fp1)
    s2 = _fp_to_positions(fp2)
    if not s1 and not s2:
        return 1.0   # two empty sets are trivially identical
    min_size = min(len(s1), len(s2))
    if min_size == 0:
        return 0.0   # exactly one is empty → no overlap
    intersection = len(s1 & s2)
    return intersection / min_size


def tanimoto_coefficient(fp1, fp2, fp_size=2048):
    """Tanimoto coefficient (Jaccard index): |A∩B| / |A∪B|.

    The standard fingerprint similarity metric in cheminformatics.
    Returns 1.0 when A and B are identical (including the empty-vs-empty
    corner case, which is conventionally 1.0 — two empty sets are
    trivially identical).  Symmetric and bounded in [0.0, 1.0].

    Accepts binary fingerprints (list of set bit positions) or sparse
    count fingerprints (list of (position, count) tuples).

    Args:
        fp1: Fingerprint as list of set bit positions or sparse counts.
        fp2: Fingerprint as list of set bit positions or sparse counts.
        fp_size: Fingerprint size in bits.

    Returns:
        float: Tanimoto similarity in [0.0, 1.0].
    """
    if fp1 is None or fp2 is None:
        return 0.0
    s1 = _fp_to_positions(fp1)
    s2 = _fp_to_positions(fp2)
    if not s1 and not s2:
        return 1.0   # two empty sets are trivially identical
    intersection = len(s1 & s2)
    union = len(s1 | s2)
    if union == 0:
        return 0.0
    return intersection / union


# v7.1.1: redirect the camelCase aliases to the Szymkiewicz-Simpson
# implementation.  In v7.1.0 they were bound directly to a C++ lambda
# that silently computed Tanimoto / Jaccard, so any caller using
# `smsd.overlapCoefficient(fp1, fp2)` was getting a different value
# from the documented semantic "overlap coefficient = |A∩B| / min(|A|,|B|)".
# Callers who relied on that wrong value should switch to
# `smsd.tanimoto_coefficient(fp1, fp2)`.
overlapCoefficient = overlap_coefficient


def count_overlap_coefficient(fp1, fp2):
    """Count-vector Szymkiewicz-Simpson overlap: sum(min(a,b)) / min(sum(a), sum(b)).

    Generalisation of the binary overlap coefficient for integer
    count fingerprints (e.g. ECFP counts, topological torsion counts).

    .. versionadded:: 7.1.1
       v7.1.0 had a `count_overlapCoefficient` binding that silently
       computed count-based Tanimoto (`sum(min(a,b)) / sum(max(a,b))`).
       Call :func:`count_tanimoto_coefficient` for the Tanimoto value.
    """
    def _to_counts(fp):
        counts = {}
        for item in fp:
            if isinstance(item, (tuple, list)) and len(item) == 2:
                counts[item[0]] = counts.get(item[0], 0) + int(item[1])
            else:
                counts[int(item)] = counts.get(int(item), 0) + 1
        return counts
    a = _to_counts(fp1)
    b = _to_counts(fp2)
    inter = sum(min(a.get(k, 0), b[k]) for k in b)
    s1 = sum(a.values())
    s2 = sum(b.values())
    mn = min(s1, s2)
    if mn == 0:
        return 0.0
    return inter / mn


def count_tanimoto_coefficient(fp1, fp2):
    """Count-vector Tanimoto (generalised Jaccard):
    sum(min(a,b)) / sum(max(a,b)).

    The standard count-based fingerprint similarity metric.
    """
    def _to_counts(fp):
        counts = {}
        for item in fp:
            if isinstance(item, (tuple, list)) and len(item) == 2:
                counts[item[0]] = counts.get(item[0], 0) + int(item[1])
            else:
                counts[int(item)] = counts.get(int(item), 0) + 1
        return counts
    a = _to_counts(fp1)
    b = _to_counts(fp2)
    keys = set(a) | set(b)
    inter = sum(min(a.get(k, 0), b.get(k, 0)) for k in keys)
    union = sum(max(a.get(k, 0), b.get(k, 0)) for k in keys)
    if union == 0:
        return 0.0
    return inter / union


# v7.1.1: redirect the camelCase count alias.  Same semantic fix
# as the binary overlapCoefficient alias above.
count_overlapCoefficient = count_overlap_coefficient




def _remap_query_weights(opts, graph, rdkit_mol):
    """Convert weights from original RDKit atom order into native graph order."""
    weights = opts.atom_weights
    if rdkit_mol is not None and weights:
        if len(weights) < rdkit_mol.GetNumAtoms():
            raise ValueError("atom_weights size is smaller than query atom count")
        index_map = get_index_map(graph)
        if index_map is None or len(index_map) != len(graph):
            raise ValueError("Cannot translate query atom weights")
        opts.atom_weights = [weights[index] for index in index_map]


def _search_options(timeout_ms, options):
    """Build native chemistry/search options without dropping caller settings."""
    options = dict(options)
    tautomer = options.pop("tautomer_aware", False)
    chem = ChemOptions.tautomer_profile() if tautomer else ChemOptions()
    opts = MCSOptions()
    opts.timeout_ms = timeout_ms
    for key, value in options.items():
        if key == "match_bond_order" and isinstance(value, str):
            modes = {"strict": BondOrderMode.STRICT, "loose": BondOrderMode.LOOSE,
                     "any": BondOrderMode.ANY}
            if value.lower() not in modes:
                raise ValueError("match_bond_order must be 'strict', 'loose', or 'any'")
            value = modes[value.lower()]
        if hasattr(opts, key):
            setattr(opts, key, value)
        elif hasattr(chem, key):
            setattr(chem, key, value)
        else:
            raise AttributeError(f"Unknown search option: {key}")
    return chem, opts


def _batch_inputs(query, targets):
    q, rdkit_q = _ensure_mol_ex(query)
    converted = [_ensure_mol_ex(target) for target in targets]
    return q, rdkit_q, [g for g, _ in converted], [mol for _, mol in converted]


def batch_substructure(query, targets, *, timeout_ms=10000, num_threads=0):
    """Check one query against many targets using a per-pair time budget.

    Inputs may be SMILES, MolGraph or RDKit Mol objects. ``num_threads=0``
    uses the available OpenMP workers; ``1`` runs sequentially.
    """
    q, _, ts, _ = _batch_inputs(query, targets)
    return _smsd.batch_substructure(q, ts, ChemOptions(), num_threads, timeout_ms)


def batch_find_substructure(query, targets, *, num_threads=0, timeout_ms=10000):
    """Return substructure pairs in each input molecule's atom ordering.

    The output contains one list of ``(query_atom, target_atom)`` pairs per
    target. An empty list means no match. Inputs may be SMILES, MolGraph or
    RDKit Mol objects. Each pair receives ``timeout_ms`` milliseconds.
    """
    q, rdkit_q, ts, rdkit_ts = _batch_inputs(query, targets)
    results = _smsd.batch_find_substructure(q, ts, ChemOptions(), num_threads, timeout_ms)
    return [list(_auto_translate(dict(mapping), q, t, rdkit_q, rdkit_t)[0].items())
            for mapping, t, rdkit_t in zip(results, ts, rdkit_ts)]


def batch_mcs(query, targets, *, timeout_ms=10000, num_threads=0, **kwargs):
    """Compute native MCS mappings for one query against many targets.

    Inputs may be SMILES, MolGraph or RDKit Mol objects. Mappings use each
    input's atom ordering. Chemistry and MCSOptions keywords are preserved;
    ``timeout_ms`` applies to each pair. Pre-parsed graphs avoid repeat parsing.
    """
    q, rdkit_q, ts, rdkit_ts = _batch_inputs(query, targets)
    chem, opts = _search_options(timeout_ms, kwargs)
    _remap_query_weights(opts, q, rdkit_q)
    results = _smsd.batch_mcs(q, ts, chem, opts, num_threads)
    return [_auto_translate(mapping, q, t, rdkit_q, rdkit_t)[0]
            for mapping, t, rdkit_t in zip(results, ts, rdkit_ts)]


def batch_mcs_size(query, targets, *, timeout_ms=10000, num_threads=0, **kwargs):
    """Compute native MCS atom counts in target order without Python mappings.

    Accepts the same chemistry/search keywords as :func:`batch_mcs`.
    ``timeout_ms`` applies to each pair; ``num_threads=0`` uses OpenMP defaults.
    """
    q, rdkit_q, ts, _ = _batch_inputs(query, targets)
    chem, opts = _search_options(timeout_ms, kwargs)
    _remap_query_weights(opts, q, rdkit_q)
    return _smsd.batch_mcs_size(q, ts, chem, opts, num_threads)


def screen_and_match(query, targets, threshold, *, timeout_ms=10000, num_threads=0, **kwargs):
    """RASCAL pre-screen followed by exact MCS on screened-in targets.

    Args:
        query: MolGraph or SMILES string.
        targets: List of MolGraph or SMILES strings.
        threshold: Minimum RASCAL upper-bound score to keep a target.
        timeout_ms: Per-pair timeout in milliseconds for the exact MCS phase.
        num_threads: OpenMP worker count. `0` uses all available processors.
        **kwargs: Extra MCSOptions fields.

    Returns:
        list[tuple[int, dict]]: `(target_index, mapping)` pairs for hits.
    """
    from smsd._smsd import screen_and_match as _screen_and_match  # type: ignore[import]
    q = _ensure_mol(query)
    ts = [_ensure_mol(t) for t in targets]
    opts = MCSOptions()
    opts.timeout_ms = timeout_ms
    for k, v in kwargs.items():
        setattr(opts, k, v)
    return _screen_and_match(q, ts, threshold, ChemOptions(), opts, num_threads)


def screen_and_mcs_size(query, targets, threshold, *, timeout_ms=10000, num_threads=0, **kwargs):
    """RASCAL pre-screen followed by exact MCS size on screened-in targets.

    Args:
        query: MolGraph or SMILES string.
        targets: List of MolGraph or SMILES strings.
        threshold: Minimum RASCAL upper-bound score to keep a target.
        timeout_ms: Per-pair timeout in milliseconds for the exact MCS phase.
        num_threads: OpenMP worker count. `0` uses all available processors.
        **kwargs: Extra MCSOptions fields.

    Returns:
        list[tuple[int, int]]: `(target_index, mcs_size)` pairs for hits.
    """
    from smsd._smsd import screen_and_mcs_size as _screen_and_mcs_size  # type: ignore[import]
    q = _ensure_mol(query)
    ts = [_ensure_mol(t) for t in targets]
    opts = MCSOptions()
    opts.timeout_ms = timeout_ms
    for k, v in kwargs.items():
        setattr(opts, k, v)
    return _screen_and_mcs_size(q, ts, threshold, ChemOptions(), opts, num_threads)


class TargetCorpus:
    """Pre-warmed, fingerprinted target collection for repeated queries.

    Amortises the cost of graph invariant computation and fingerprint
    generation across many queries.  Call :meth:`prewarm` once after
    loading all targets, then run :meth:`substructure`, :meth:`find_substructure`,
    :meth:`mcs_size`, or :meth:`screen` as many times as needed.

    Example::

        corpus = smsd.TargetCorpus.from_smiles(["CCO", "c1ccccc1", "CC(=O)O"])
        corpus.prewarm()
        hits = corpus.substructure("CC")
        maps = corpus.find_substructure("CC")
    """

    def __init__(self):
        from smsd._smsd import TargetCorpus as _TargetCorpus  # type: ignore[import]
        self._impl = _TargetCorpus()

    @classmethod
    def from_smiles(cls, smiles_list):
        """Create a TargetCorpus from a list of SMILES strings.

        Args:
            smiles_list: Iterable of SMILES strings or MolGraph objects.

        Returns:
            TargetCorpus: A new corpus with all targets added (not yet prewarmed).
        """
        corpus = cls()
        mols = [_ensure_mol(s) for s in smiles_list]
        corpus._impl.add_targets(mols)
        return corpus

    def add_target(self, mol):
        """Add a single target molecule (MolGraph or SMILES string)."""
        self._impl.add_target(_ensure_mol(mol))

    def add_targets(self, mols):
        """Add multiple target molecules (MolGraph or SMILES strings)."""
        self._impl.add_targets([_ensure_mol(m) for m in mols])

    def prewarm(self, *, num_threads=0):
        """Pre-compute graph invariants and ECFP fingerprints for all targets.

        Args:
            num_threads: OpenMP worker count. ``0`` uses all available processors.
        """
        self._impl.prewarm(num_threads)

    def __len__(self):
        return len(self._impl)

    @property
    def is_prewarmed(self):
        """True if :meth:`prewarm` has been called since last target addition."""
        return self._impl.is_prewarmed

    def substructure(self, query, *, num_threads=0):
        """Check if query is a substructure of each target.

        Args:
            query: MolGraph or SMILES string.
            num_threads: OpenMP worker count. ``0`` uses all available processors.

        Returns:
            list[bool]: True at index i when query is substructure of target i.
        """
        return self._impl.substructure(_ensure_mol(query), ChemOptions(), num_threads)

    def find_substructure(self, query, *, num_threads=0):
        """Find substructure atom-atom mappings for query against each target.

        Args:
            query: MolGraph or SMILES string.
            num_threads: OpenMP worker count. ``0`` uses all available processors.

        Returns:
            list[list[tuple[int, int]]]: Atom mappings for each target.
                Empty inner list means no match.
        """
        return self._impl.find_substructure(
            _ensure_mol(query), ChemOptions(), num_threads)

    def mcs_size(self, query, *, num_threads=0, **kwargs):
        """Compute MCS atom counts for query against each target.

        Args:
            query: MolGraph or SMILES string.
            num_threads: OpenMP worker count. ``0`` uses all available processors.
            **kwargs: Extra MCSOptions fields.

        Returns:
            list[int]: MCS sizes in target order.
        """
        opts = MCSOptions()
        for k, v in kwargs.items():
            setattr(opts, k, v)
        return self._impl.mcs_size(
            _ensure_mol(query), ChemOptions(), opts, num_threads)

    def screen(self, query, threshold, *, num_threads=0):
        """Tanimoto fingerprint screening.

        Requires :meth:`prewarm` to have been called (computes fingerprints).

        Args:
            query: MolGraph or SMILES string.
            threshold: Minimum Tanimoto similarity to include a target.
            num_threads: OpenMP worker count. ``0`` uses all available processors.

        Returns:
            list[int]: Indices of targets above the threshold.
        """
        return self._impl.screen(_ensure_mol(query), threshold, num_threads)


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
# Batch constrained MCS (v6.6.0)
# ---------------------------------------------------------------------------

def batch_mcs_constrained(queries, targets, *, timeout_ms=10000,
                          return_target_indices=False, **kwargs):
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
        return_target_indices: Return ``(target_index, mapping)`` pairs when True.
            The target index is -1 for an empty result. Use this for multiple targets.
        **kwargs: additional chemistry/search options.

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

    chem, opts = _search_options(timeout_ms, kwargs)
    if opts.atom_weights and any(mol is not None for mol in rdkit_q_list):
        if len(gq_list) != 1:
            raise ValueError("Weighted constrained batches with RDKit inputs require one query; "
                             "convert multiple queries to MolGraph and supply weights in graph order")
        _remap_query_weights(opts, gq_list[0], rdkit_q_list[0])
    raw_results, selected_targets = _smsd.batch_mcs_constrained_with_targets(
        gq_list, gt_list, chem, opts)
    translated = []
    for i, (mapping, target_index) in enumerate(zip(raw_results, selected_targets)):
        if target_index >= 0:
            mapping = _auto_translate(mapping, gq_list[i], gt_list[target_index],
                                      rdkit_q_list[i], rdkit_t_list[target_index])[0]
        translated.append((target_index, mapping) if return_target_indices else mapping)
    return translated


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
        overlap (float): Szymkiewicz-Simpson overlap = size / min(query_atoms, target_atoms).
            Returns 1.0 when the smaller molecule is fully contained in the larger.
        overlapCoefficient (float): Alias of ``overlap`` (Szymkiewicz-Simpson).
            Kept for backward compatibility with v7.1.0 code that read this field.
        tanimoto (float): Jaccard / Tanimoto = size / (query_atoms + target_atoms - size).
            Symmetric and penalised by size mismatch.
        query_atoms (int): Number of atoms in the query molecule.
        target_atoms (int): Number of atoms in the target molecule.
        mcs_smiles (str): Canonical SMILES of the MCS subgraph, or "" if empty.

    .. versionchanged:: 7.1.1
       Fixed a long-standing semantic bug: ``overlapCoefficient`` was
       documented as Szymkiewicz-Simpson but silently computed Tanimoto.
       Both ``overlap`` and ``overlapCoefficient`` now return Simpson
       overlap as documented.  The Jaccard / Tanimoto value is available
       via the new ``tanimoto`` attribute.
    """
    __slots__ = ('mapping', 'size', 'overlap', 'overlapCoefficient', 'tanimoto',
                 'query_atoms', 'target_atoms', 'mcs_smiles')

    def __init__(self, mapping, query_n, target_n, mcs_smi=""):
        self.mapping = mapping
        self.size = len(mapping)
        self.query_atoms = query_n
        self.target_atoms = target_n
        min_n = min(query_n, target_n)
        # Szymkiewicz-Simpson overlap = |A∩B| / min(|A|,|B|); equals 1.0 on containment.
        self.overlap = self.size / min_n if min_n > 0 else 0.0
        # Back-compat: the documented semantics of `overlapCoefficient` is
        # Simpson overlap, matching the class docstring and the usage example.
        self.overlapCoefficient = self.overlap
        # Jaccard / Tanimoto = |A∩B| / |A∪B|; penalises size mismatch.
        denom = query_n + target_n - self.size
        self.tanimoto = self.size / denom if denom > 0 else 0.0
        self.mcs_smiles = mcs_smi

    def __len__(self):
        return self.size

    def __repr__(self):
        return (f"MatchResult(size={self.size}, overlap={self.overlap:.3f}, "
                f"tanimoto={self.tanimoto:.3f}, mcs_smiles={self.mcs_smiles!r})")


def mcs_result(mol1, mol2, **kwargs):
    """Find MCS between two molecules and return a MatchResult.

    Like :func:`mcs`, but wraps the raw mapping dict in a
    :class:`MatchResult` that includes size and overlap coefficient similarity.

    Accepts MolGraph objects, SMILES strings, or RDKit Mol objects.
    When RDKit Mol objects are passed, returned indices are automatically
    translated to RDKit's native atom ordering.

    Args:
        mol1: First molecule (MolGraph, SMILES string, or RDKit Mol).
        mol2: Second molecule (MolGraph, SMILES string, or RDKit Mol).
        **kwargs: Keyword arguments forwarded to :func:`mcs`
            (e.g. tautomer_aware, timeout_ms, prefer_rare_heteroatoms).

    Returns:
        MatchResult: Structured result with mapping, size, and overlapCoefficient.

    Example::

        import smsd
        r = smsd.mcs_result("c1ccccc1", "c1ccc(O)cc1")
        print(r)          # MatchResult(size=6, overlapCoefficient=1.000)
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
    kwargs.pop('prefer_rare_heteroatoms', False)  # flag kept but scorer removed
    taut = kwargs.pop('tautomer_aware', False)
    tms = kwargs.pop('timeout_ms', 10000)
    chem = ChemOptions.tautomer_profile() if taut else ChemOptions()
    opts = MCSOptions()
    opts.timeout_ms = tms
    for k, v in kwargs.items():
        setattr(opts, k, v)
    raw = _native_find_mcs(g1, g2, chem, opts)
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
        print(r.overlapCoefficient)    # overlap score
    """
    if not isinstance(smi1, str) or not smi1:
        raise ValueError("smi1 must be a non-empty SMILES string")
    if not isinstance(smi2, str) or not smi2:
        raise ValueError("smi2 must be a non-empty SMILES string")
    g1 = parse_smiles(smi1)
    g2 = parse_smiles(smi2)
    mapping = find_mcs(g1, g2, **kwargs)
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
    # Accept both dict {pos: count} and list-of-tuples [(pos, count), ...]
    items = counts_dict.items() if isinstance(counts_dict, dict) else counts_dict
    for pos, count in items:
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
    ``(index, overlapCoefficient)`` pairs sorted from most to least similar.

    Args:
        query: Query molecule (MolGraph or SMILES string).
        targets: List of target molecules (MolGraph or SMILES strings).
        fp_size (int): Fingerprint size in bits (default 2048).
        radius (int): ECFP radius (default 2).

    Returns:
        list[tuple[int, float]]: Sorted ``(target_index, overlapCoefficient)`` pairs,
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
        mapping_pairs = _native_find_substructure(g_core, g_mol, ChemOptions(), timeout_ms)
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

def compile_smarts(smarts, *, max_recursion_depth=20):
    """Parse SMARTS once and return a reusable compiled query object."""
    if not _HAS_SMARTS:
        raise NotImplementedError(
            "SMARTS matching not available in this build. "
            "Rebuild with smarts_parser.hpp included."
        )
    if SmartsQuery is not None and isinstance(smarts, SmartsQuery):
        return smarts
    if not isinstance(smarts, str):
        raise TypeError(
            f"Expected SMARTS string or SmartsQuery, got {type(smarts).__name__}"
        )
    if _compile_smarts_raw is None:
        raise RuntimeError(
            "compile_smarts requires a rebuilt extension with the 6.11.0 SMARTS bindings."
        )
    return _compile_smarts_raw(smarts, max_recursion_depth)


def prewarm(mol):
    """Compute and cache common lazy invariants on a MolGraph once."""
    g = _ensure_mol(mol)
    g.prewarm()
    return g


def perceive_aromaticity(mol, *, model=AromaticityModel.DAYLIGHT_LIKE):
    """Recompute ring membership and aromaticity on a MolGraph."""
    g = _ensure_mol(mol)
    g.perceive_aromaticity(model)
    return g


def kekulize(mol):
    """Convert an aromatic MolGraph into an explicit Kekule form."""
    g = _ensure_mol(mol)
    g.kekulize()
    return g


def dearomatize(mol):
    """Alias of kekulize() for explicit aromatic-flag removal."""
    g = _ensure_mol(mol)
    g.dearomatize()
    return g


def _ensure_smarts_query(smarts):
    if not _HAS_SMARTS:
        raise NotImplementedError(
            "SMARTS matching not available in this build. "
            "Rebuild with smarts_parser.hpp included."
        )
    return compile_smarts(smarts)

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
    g = _ensure_mol(target)
    if isinstance(smarts, str) and _compile_smarts_raw is None:
        return _find_mcs_smarts_raw(smarts, g, max_matches)
    query = _ensure_smarts_query(smarts)
    matches = query.find_all(g, max_matches=max_matches)
    if not matches:
        return {}
    return max(matches, key=len)


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
    g = _ensure_mol(target)
    if isinstance(smarts, str) and _compile_smarts_raw is None:
        return _smarts_match_raw(smarts, g)
    query = _ensure_smarts_query(smarts)
    return query.matches(g)


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
    g = _ensure_mol(target)
    if isinstance(smarts, str) and _compile_smarts_raw is None:
        return _smarts_find_all_raw(smarts, g, max_matches)
    query = _ensure_smarts_query(smarts)
    return query.find_all(g, max_matches=max_matches)


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
    translated = translate_mapping(mapping, g1 if rdkit1 is not None else None,
                                   g2 if rdkit2 is not None else None)
    # Translation preserves the native chemistry policy, including atom-type
    # relaxation explicitly requested by the caller.
    validated = {}
    for q, t in translated.items():
        ok = True
        if rdkit1 is not None:
            if q < 0 or q >= rdkit1.GetNumAtoms():
                ok = False
        if rdkit2 is not None and (t < 0 or t >= rdkit2.GetNumAtoms()):
            ok = False
        if ok:
            validated[q] = t
    return validated, smsd_mapping


# ---------------------------------------------------------------------------
# RDKit interoperability helpers (RDKit is an optional dependency)
# ---------------------------------------------------------------------------

class _RdkitResult:
    """Wrapper holding a MolGraph + its SMSD→RDKit index mapping."""
    __slots__ = ('graph', 'index_map', 'signature')
    def __init__(self, graph, index_map, signature):
        self.graph = graph
        self.index_map = index_map  # list: index_map[smsd_idx] = rdkit_idx
        self.signature = signature


# Cache each live molecule independently so reordered molecules cannot share
# atom-index metadata. Weak keys release entries when their inputs are discarded.
import threading as _threading
import weakref as _weakref
_rdkit_cache_lock = _threading.Lock()
_rdkit_cache = _weakref.WeakKeyDictionary()
_mol_graph_cache = {}  # SMILES -> MolGraph (for parse_smiles caching only)


def _compute_index_map(rdkit_mol, canonical_smi):
    """Compute deterministic SMSD→RDKit atom index mapping.

    Uses RDKit's _smilesAtomOutputOrder property which records the exact
    order atoms appear in the canonical SMILES string. This is the ground
    truth for SMILES parse order = SMSD's internal atom order.

    The mapping: index_map[smsd_idx] = rdkit_idx

    NOTE: The caller (from_rdkit) must have already called MolToSmiles()
    on this molecule, which populates _smilesAtomOutputOrder as a side effect.
    """
    from rdkit import Chem

    n = rdkit_mol.GetNumAtoms()

    # _smilesAtomOutputOrder is populated by the MolToSmiles() call in from_rdkit().
    try:
        output_order = list(rdkit_mol.GetPropsAsDict(includePrivate=True, includeComputed=True)
                            .get('_smilesAtomOutputOrder', []))
    except Exception:
        output_order = []

    if len(output_order) == n and sorted(output_order) == list(range(n)):
        # output_order[smiles_position] = rdkit_atom_idx
        # SMSD parses the canonical SMILES, so smsd_idx = smiles_position
        return output_order
    else:
        # Fallback: use GetSubstructMatch on re-parsed canonical SMILES
        mol_from_smi = Chem.MolFromSmiles(canonical_smi)
        if mol_from_smi is not None:
            match = rdkit_mol.GetSubstructMatch(mol_from_smi, useChirality=True)
            if match and len(match) == n:
                return list(match)
        raise ValueError("Cannot determine the SMSD-to-RDKit atom-index mapping")


def from_rdkit(mol, use_cache=True):
    """Convert an RDKit Mol to SMSD MolGraph via SMILES round-trip.

    Returns a MolGraph. The SMSD→RDKit index mapping is stored internally
    and used by :func:`translate_mapping` and :func:`mcs_rdkit_native`.

    Each distinct RDKit Mol object gets its own mapping (no cross-Mol
    cache collisions). Weak cache keys release discarded molecules. In-place
    molecule changes invalidate their cached conversion. Index metadata remains
    available while a converted graph is alive, including after clear_cache().

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
    if not isinstance(mol, Chem.Mol):
        raise TypeError("mol must be an RDKit Mol")
    if mol.GetNumAtoms() == 0:
        raise ValueError("mol has 0 atoms; cannot convert empty molecule")

    # Check if molecule is only implicit/explicit H (no heavy atoms)
    heavy = mol.GetNumHeavyAtoms()
    if heavy == 0:
        raise ValueError("mol contains only hydrogen atoms; no heavy atoms to convert")

    signature = mol.ToBinary()
    if use_cache:
        with _rdkit_cache_lock:
            cached = _rdkit_cache.get(mol)
            if cached is not None and cached.signature == signature:
                return cached.graph

    smi = Chem.MolToSmiles(mol)
    if not smi:
        raise ValueError("RDKit returned empty SMILES for the molecule")

    # Each RDKit Mol gets its OWN MolGraph instance to prevent index map
    # collisions.  Two Mols with the same canonical SMILES but different
    # atom orderings must have separate index maps, and get_index_map()
    # associates metadata with the live graph.  Sharing a MolGraph across
    # Mols would cause the second Mol's index map to overwrite the first.
    g = parse_smiles(smi)

    # Compute deterministic index mapping for THIS specific Mol object
    index_map = _compute_index_map(mol, smi)
    result = _RdkitResult(g, index_map, mol.ToBinary())

    with _rdkit_cache_lock:
        _graph_to_result[g] = tuple(index_map)

        if use_cache:
            if len(_rdkit_cache) >= 10000:
                _rdkit_cache.clear()
            _rdkit_cache[mol] = result

    return g


_graph_to_result = _weakref.WeakKeyDictionary()  # live graph -> immutable index map


def get_index_map(g):
    """Return the SMSD→RDKit atom index mapping for a MolGraph from from_rdkit().

    Returns a list where ``result[smsd_idx] = rdkit_idx``, or None if the
    MolGraph was not created by from_rdkit().
    """
    with _rdkit_cache_lock:
        index_map = _graph_to_result.get(g)
        return list(index_map) if index_map is not None else None


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
    """Run native MCS and return indices in the original RDKit molecules.

    Accepts :func:`find_mcs` options. Index translation preserves the selected
    mapping, including explicitly relaxed atom or bond matching.
    """
    from_rdkit(mol1)
    from_rdkit(mol2)
    kwargs["strategy"] = "native"
    return find_mcs(mol1, mol2, **kwargs)


def clear_cache():
    """Clear reusable conversions without invalidating live graphs' index maps."""
    with _rdkit_cache_lock:
        _rdkit_cache.clear()
        _mol_graph_cache.clear()


def batch_mcs_rdkit(query_mol, target_mols, **kwargs):
    """Run native batch MCS using original RDKit atom indices.

    Accepts the chemistry/search options and thread count of :func:`batch_mcs`.
    """
    target_mols = list(target_mols)
    from_rdkit(query_mol)
    for mol in target_mols:
        from_rdkit(mol)
    return batch_mcs(query_mol, target_mols, **kwargs)


def mcs_rdkit(mol1, mol2, **kwargs):
    """Find MCS using :func:`find_mcs` options and original RDKit indices."""
    from_rdkit(mol1)
    from_rdkit(mol2)
    return find_mcs(mol1, mol2, **kwargs)


def substructure_rdkit(query_mol, target_mol, **kwargs):
    """Find substructure mappings using original RDKit atom indices."""
    from_rdkit(query_mol)
    from_rdkit(target_mol)
    return find_substructure(query_mol, target_mol, **kwargs)


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
    mapping = find_mcs(g1, g2, tautomer_aware=tautomer_aware, timeout_ms=timeout_ms)

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

    result = find_substructure(gq, gt, timeout_ms=timeout_ms)
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
    "MCSOptions",
    "MolGraph",
    "MolGraphBuilder",
    # Enums
    "BondOrderMode",
    "AromaticityMode",
    "AromaticityModel",
    "RingFusionMode",
    "MatcherEngine",
    "Solvent",
    # SMILES / SMARTS
    "parse_smiles",
    "read_mol_block",
    "read_molfile",
    "parse_mapped_smiles",
    "strip_atom_maps",
    "to_smiles",
    "to_smarts",
    "compile_smarts",
    "write_mol_block",
    "write_mol_block_v3000",
    "write_sdf_record",
    "write_molfile",
    # Substructure
    "find_substructure",
    "is_substructure",
    # MCS
    "find_mcs",
    "canonicalize_mapping",
    "prewarm",
    "perceive_aromaticity",
    "kekulize",
    "dearomatize",
    # Index translation
    "translate_to_atom_ids",
    "find_mcs_progressive",
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
    "overlap_coefficient",
    "tanimoto_coefficient",
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
    "batch_find_substructure",
    "batch_mcs",
    "batch_mcs_size",
    "screen_and_match",
    "screen_and_mcs_size",
    "batch_mcs_rdkit",
    # Target corpus
    "TargetCorpus",
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
    # Lightweight MCS engine (clique solver + Python pipeline)
    "mcs_engine",
]

# Lightweight Python MCS engine — uses clique solver C++ backend
# with RDKit chemistry for reaction mapping workflows.
from smsd import mcs_engine
try:
    _find_mcs_light = mcs_engine.find_mcs_lightweight
except Exception:
    _find_mcs_light = None
