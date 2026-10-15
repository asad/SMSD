# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Chemical graph I/O, fingerprints, MCS and substructure searches.

    import smsd

    result = smsd.Substructure.find("c1ccccc1").match("c1ccc(O)cc1")
    assert result.found

    result = smsd.MCS.find("c1ccc(O)cc1").match("c1ccc(N)cc1")
    if result.found:
        print(result.best.atom_mapping, result.optimal)
"""

import importlib.machinery
import importlib.util
import sys
from pathlib import Path

__version__ = "7.2.2"
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
MolGraph = _smsd.MolGraph
MolGraphBuilder = _smsd.MolGraphBuilder

BondOrderMode = _smsd.BondOrderMode
AromaticityMode = _smsd.AromaticityMode
AromaticityModel = _smsd.AromaticityModel
RingFusionMode = _smsd.RingFusionMode
Solvent = _smsd.Solvent

_parse_smiles_raw = _smsd.parse_smiles
read_mol_block = _smsd.read_mol_block
_to_smiles_native = _smsd.to_smiles
to_smarts = _smsd.to_smarts


def to_smiles(mol):
    """Write a canonical SMILES string for a molecule.

    Accepts either a :class:`MolGraph` or a SMILES string. Canonical
    ordering is defined by the writer and may differ between languages.

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
    labelling mean they are not proofs of molecular identity or difference;
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

# Graph index translation
translate_to_atom_ids = _smsd.translate_to_atom_ids

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

# Graph utilities
extract_subgraph = _smsd.extract_subgraph
murcko_scaffold = _smsd.murcko_scaffold
count_components = _smsd.count_components
split_components = _smsd.split_components
same_canonical_graph = _smsd.same_canonical_graph

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
        mol1 atom -> mol2 atom mapping from a Match.
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


class TargetCorpus:
    """Snapshot targets for repeated query.match_many(corpus) and screening."""

    def __init__(self):
        self._impl = _smsd.TargetCorpus()
        self._indices = []
        self._lock = _threading.RLock()

    @classmethod
    def from_smiles(cls, targets):
        """Create a corpus from molecule objects or SMILES strings."""
        corpus = cls()
        corpus.add_targets(targets)
        return corpus

    def add_target(self, target):
        """Capture a molecule and its current input atom numbering."""
        self.add_targets([target])

    def add_targets(self, targets):
        """Capture targets together; later input changes do not alter them."""
        converted = [_ensure_mol_ex(target) for target in targets]
        indices = [tuple(get_index_map(graph)) if source is not None else None
                   for graph, source in converted]
        with self._lock:
            self._impl.add_targets([graph for graph, _ in converted])
            self._indices.extend(indices)

    def _snapshot(self):
        with self._lock:
            return list(zip(self._impl.targets(), self._indices))

    def prewarm(self, *, num_threads=0):
        """Compute graph invariants and ECFP fingerprints once."""
        self._impl.prewarm(num_threads)

    def __len__(self):
        return len(self._impl)

    @property
    def is_prewarmed(self):
        return self._impl.is_prewarmed

    def screen(self, query, threshold, *, num_threads=0):
        """Return target indices meeting the ECFP Tanimoto threshold."""
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
    """Return SSSR rings optimised for 2D coordinate generation.

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
    """Reduce bond crossings in a 2D layout by optimising ring orientations.

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


# ---------------------------------------------------------------------------
# mcs_from_smiles -- all-in-one SMILES -> MCS MatchResult
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# fingerprint_from_smiles -- all-in-one SMILES -> ECFP fingerprint
# ---------------------------------------------------------------------------

def fingerprint_from_smiles(smiles, radius=2, fp_size=2048, mode="ecfp", *, path_length=7):
    """All-in-one convenience: parse a SMILES string and compute its fingerprint.

    Args:
        smiles (str): SMILES string.
        radius (int): ECFP radius (2 = ECFP4, 3 = ECFP6; default 2).
        fp_size (int): Fingerprint size in bits (default 2048).
        mode (str): ``"ecfp"`` (circular/Morgan), ``"fcfp"`` (pharmacophoric),
                    ``"path"`` (path-based), or ``"mcs"`` (MCS-aware).
        path_length (int): Maximum bond path length for path/MCS modes (default 7).

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
        return path_fingerprint(g, path_length, fp_size)
    elif mode == "mcs":
        return mcs_fingerprint(g, path_length, fp_size)
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
    q_fp = circular_fingerprint(q, radius, fp_size, mode="ecfp")
    q_set = set(q_fp)

    results = []
    for i, t in enumerate(targets):
        g = _ensure_mol(t)
        t_fp = circular_fingerprint(g, radius, fp_size, mode="ecfp")
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


# ---------------------------------------------------------------------------
# SMARTS-based MCS and matching
# ---------------------------------------------------------------------------


def prewarm(mol):
    """Compute and cache common lazy invariants on a MolGraph once."""
    g = _ensure_mol(mol)
    g.prewarm()
    return g


def perceive_aromaticity(mol, *, model=AromaticityModel.DAYLIGHT_LIKE):
    """Recompute ring membership and aromaticity on a MolGraph."""
    g = _ensure_mol(mol) if isinstance(mol, (str, MolGraph)) else from_rdkit(mol, use_cache=False)
    g.perceive_aromaticity(model)
    return g


def kekulize(mol):
    """Convert an aromatic MolGraph into an explicit Kekule form."""
    g = _ensure_mol(mol) if isinstance(mol, (str, MolGraph)) else from_rdkit(mol, use_cache=False)
    g.kekulize()
    return g


def dearomatize(mol):
    """Alias of kekulize() for explicit aromatic-flag removal."""
    g = _ensure_mol(mol) if isinstance(mol, (str, MolGraph)) else from_rdkit(mol, use_cache=False)
    g.dearomatize()
    return g


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
    and used by :func:`translate_mapping` and query results.

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
        mapping: dict of {smsd_query_idx: smsd_target_idx} from a Match
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


def clear_cache():
    """Clear reusable conversions without invalidating live graphs' index maps."""
    with _rdkit_cache_lock:
        _rdkit_cache.clear()
        _mol_graph_cache.clear()


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


def _molecule_for_depiction(graph):
    from rdkit import Chem
    return Chem.MolFromMolBlock(write_mol_block(graph), removeHs=False)


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
    options = MCS.Options()
    options.chemistry.tautomer_aware = tautomer_aware
    result = MCS.find(g1).with_options(options).with_timeout(timeout_ms).match(g2)
    mapping = result.best.atom_mapping if result.found else {}

    # Convert to RDKit mols for drawing
    rdmol1 = _molecule_for_depiction(g1)
    rdmol2 = _molecule_for_depiction(g2)
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

    result = Substructure.find(gq).with_timeout(timeout_ms).match(gt)
    mapping = result.best.atom_mapping if result.found else {}

    rdmol_q = _molecule_for_depiction(gq)
    rdmol_t = _molecule_for_depiction(gt)
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


from .search import MCS, Substructure, Match, SearchResult, ConsensusResult

SearchStatus = _smsd.SearchStatus
CancellationToken = _smsd.CancellationToken

__all__ = ['AromaticityMode', 'AromaticityModel', 'BondOrderMode', 'CIPDescriptors', 'CancellationToken', 'ConsensusResult', 'ChemOptions', 'DepictOptions', 'EZLabel', 'MCS', 'Match', 'MolGraph', 'MolGraphBuilder', 'Point2D', 'Point3D', 'RSLabel', 'RingFusionMode', 'SearchResult', 'SearchStatus', 'Solvent', 'Substructure', 'TargetCorpus', 'align_2d', 'analyze_fp_quality', 'assign_cip', 'assign_ez', 'assign_rs', 'batch_similarity', 'bounding_box_2d', 'canonical_hash', 'canonical_orientation', 'canonical_smiles', 'center_2d', 'circular_fingerprint', 'circular_fingerprint_counts', 'classify_pharmacophore', 'clear_cache', 'compute_sssr', 'cosine', 'count_components', 'count_cosine', 'count_crossings', 'count_dice', 'count_overlapCoefficient', 'count_overlap_coefficient', 'count_tanimoto_coefficient', 'counts_to_array', 'dearomatize', 'depict_mapping', 'depict_mcs', 'depict_pair', 'depict_substructure', 'depict_svg', 'dice', 'export_sdf', 'extract_subgraph', 'fingerprint', 'fingerprint_from_smiles', 'fingerprint_subset', 'force_directed_layout', 'from_hex', 'from_openbabel', 'from_rdkit', 'generate_coords_2d', 'generate_coords_3d', 'get_index_map', 'gpu_device_info', 'gpu_is_available', 'implicit_h', 'is_degenerate_layout', 'kekulize', 'layout_quality', 'layout_sssr', 'match_template', 'mcs_fingerprint', 'mirror_x', 'mirror_y', 'murcko_scaffold', 'normalise_bond_length', 'overlapCoefficient', 'overlap_coefficient', 'parse_mapped_smiles', 'parse_smiles', 'path_fingerprint', 'perceive_aromaticity', 'prewarm', 'prewarm_graph', 'read_mol_block', 'read_mol_file', 'read_molfile', 'read_sdf', 'reduce_crossings', 'resolve_overlaps', 'rotate_2d', 'same_canonical_graph', 'save_svg', 'scale_2d', 'screen_targets', 'similarity', 'similarity_upper_bound', 'soergel', 'split_components', 'stress_majorisation', 'strip_atom_maps', 'tanimoto_coefficient', 'to_binary_string', 'to_hex', 'to_openbabel', 'to_rdkit', 'to_smarts', 'to_smiles', 'to_svg', 'topological_torsion', 'topological_torsion_counts', 'translate_2d', 'translate_mapping', 'translate_to_atom_ids', 'write_mol_block', 'write_mol_block_v3000', 'write_molfile', 'write_sdf', 'write_sdf_record']
