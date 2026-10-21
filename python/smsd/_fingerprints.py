# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Fingerprints, similarity coefficients and fingerprint encodings."""

import re

from . import _smsd
from ._molecule import _checked_int, _deprecated, _ensure_mol, parse_smiles

__all__ = [
    "analyze_fp_quality", "circular_fingerprint", "circular_fingerprint_counts", "cosine",
    "count_cosine", "count_dice", "count_overlapCoefficient", "count_overlap_coefficient",
    "count_tanimoto_coefficient", "counts_to_array", "dice", "fingerprint", "fingerprint_from_smiles",
    "fingerprint_subset", "from_hex", "mcs_fingerprint", "overlapCoefficient", "overlap_coefficient",
    "path_fingerprint", "soergel", "tanimoto_coefficient", "to_binary_string", "to_hex",
    "topological_torsion", "topological_torsion_counts",
]

# The fingerprint natives are bound directly: each runs in well under a microsecond, and a
# Python frame in front of them would cost a measurable share of that.
path_fingerprint = _smsd.path_fingerprint
mcs_fingerprint = _smsd.mcs_fingerprint
fingerprint_subset = _smsd.fingerprint_subset
analyze_fp_quality = _smsd.analyze_fp_quality
circular_fingerprint = _smsd.circular_fingerprint
circular_fingerprint_counts = _smsd.circular_fingerprint_counts
topological_torsion_counts = _smsd.topological_torsion_counts

dice = _smsd.dice
cosine = _smsd.cosine
soergel = _smsd.soergel

_HEX_DIGITS = re.compile("[0-9A-Fa-f]+")


def _positions(fp):
    """The bit positions of a fingerprint given as positions or as (position, count) pairs; None is empty."""
    if fp is None:
        return []
    if not isinstance(fp, (list, tuple)):
        fp = list(fp)
    if fp and isinstance(fp[0], (tuple, list)):
        return [pair[0] for pair in fp]
    return fp


def _counts(fp):
    """Counts per position of a fingerprint given as (position, count) pairs or as positions."""
    counts = {}
    for item in fp:
        if isinstance(item, (tuple, list)) and len(item) == 2:
            counts[item[0]] = counts.get(item[0], 0) + int(item[1])
        else:
            counts[int(item)] = counts.get(int(item), 0) + 1
    return counts


def topological_torsion(mol, fp_size=2048):
    """Compute a topological torsion fingerprint (Nilakantan et al. 1987).

    Enumerates all 4-atom linear paths A-B-C-D in the molecule and hashes
    their atom types (atomic number, heavy degree, pi electrons, ring
    membership) in canonical order.

    Args:
        mol: MolGraph or SMILES string.
        fp_size: Fingerprint size in bits (default 2048).

    Returns:
        list[int]: Sorted list of set bit positions.
    """
    _checked_int(fp_size, "fp_size", 1)
    return _smsd.topological_torsion(_ensure_mol(mol), fp_size)


def _circular(mode):
    return lambda g, radius, path_length, fp_size: circular_fingerprint(g, radius, fp_size, mode)


# Every fingerprint kind, by the name used in fingerprint() and fingerprint_from_smiles().
_KINDS = {
    "path": lambda g, radius, path_length, fp_size: path_fingerprint(g, path_length, fp_size),
    "mcs": lambda g, radius, path_length, fp_size: mcs_fingerprint(g, path_length, fp_size),
    "ecfp": _circular("ecfp"),
    "circular": _circular("ecfp"),
    "fcfp": _circular("fcfp"),
    "torsion": lambda g, radius, path_length, fp_size: _smsd.topological_torsion(g, fp_size),
}


def fingerprint(mol, *, kind="mcs", path_length=7, fp_size=2048, radius=2, options=None):
    """Compute a molecular fingerprint.

    Args:
        mol: MolGraph or SMILES string.
        kind: ``"path"``, ``"mcs"`` (MCS-aware paths), ``"ecfp"`` or ``"circular"``,
            ``"fcfp"`` (pharmacophoric) or ``"torsion"``.
        path_length: Maximum bond path length for path and MCS kinds (default 7).
        fp_size: Fingerprint size in bits (default 2048).
        radius: Radius of the circular kinds (default 2, which is ECFP4).
        options: A :class:`FingerprintOptions` specification. When given, it defines the
            fingerprint entirely and the other keyword arguments are not read.

    Returns:
        list[int]: Sorted list of set bit positions, or the ``(position, count)`` pairs of
        the non-zero bins when ``options.counts`` is set.

    Example::

        import smsd

        chiral = smsd.FingerprintOptions(family="ecfp", radius=2, size=2048, cip_atoms=True)
        assert smsd.fingerprint("C[C@H](N)C(=O)O", options=chiral) != smsd.fingerprint(
            "C[C@@H](N)C(=O)O", options=chiral)
    """
    if options is not None:
        return _smsd.fingerprint_options(_ensure_mol(mol), options)
    _checked_int(fp_size, "fp_size", 1)
    _checked_int(path_length, "path_length", 1)
    compute = _KINDS.get(kind)
    if compute is None:
        raise ValueError(
            f"Unknown fingerprint kind: {kind!r}. "
            f"Use 'path', 'mcs', 'ecfp', 'circular', 'fcfp', or 'torsion'."
        )
    return compute(_ensure_mol(mol), radius, path_length, fp_size)


def fingerprint_from_smiles(smiles, radius=2, fp_size=2048, mode="ecfp", *, path_length=7, options=None):
    """Parse a SMILES string and compute its fingerprint in one call.

    Args:
        smiles (str): SMILES string.
        radius (int): ECFP radius (2 = ECFP4, 3 = ECFP6; default 2).
        fp_size (int): Fingerprint size in bits (default 2048).
        mode (str): ``"ecfp"`` or ``"circular"``, ``"fcfp"`` (pharmacophoric),
            ``"path"``, ``"mcs"`` (MCS-aware paths) or ``"torsion"``.
        path_length (int): Maximum bond path length for path and MCS modes (default 7).
        options: A :class:`FingerprintOptions` specification. When given, it defines the
            fingerprint entirely and the other arguments are not read.

    Returns:
        list[int]: Sorted list of set bit positions, or the ``(position, count)`` pairs of
        the non-zero bins when ``options.counts`` is set.

    Example::

        import smsd

        fp = smsd.fingerprint_from_smiles("c1ccccc1")
        print(len(fp))  # number of set bits
        keys = smsd.fingerprint_from_smiles("CC(=O)Oc1ccccc1C(=O)O", options=smsd.FingerprintOptions.maccs())
        assert 162 in keys  # MACCS key 163, a six-membered ring, is position 162
    """
    if not isinstance(smiles, str) or not smiles:
        raise ValueError("smiles must be a non-empty SMILES string")
    if options is not None:
        return _smsd.fingerprint_options(parse_smiles(smiles), options)
    _checked_int(fp_size, "fp_size", 1)
    compute = _KINDS.get(mode)
    if compute is None:
        raise ValueError(f"Unknown mode: {mode!r}. Use 'ecfp', 'fcfp', 'path', 'mcs', 'circular', or 'torsion'.")
    return compute(parse_smiles(smiles), radius, path_length, fp_size)


def overlap_coefficient(fp1, fp2, fp_size=2048):
    """Overlap coefficient (Szymkiewicz-Simpson): |A∩B| / min(|A|, |B|).

    Measures how far the smaller fingerprint is contained in the larger one,
    so it suits substructure-like comparisons where a size difference is
    expected; it is 1.0 when either fingerprint is a subset of the other.

    Accepts binary fingerprints (lists of set bit positions) or sparse count
    fingerprints (lists of (position, count) pairs).

    Two empty fingerprints return 1.0 (identical empty sets); if exactly one
    fingerprint is empty, the result is 0.0. The native metric is shared with
    C++ and Java.

    Args:
        fp1: Fingerprint as a list of set bit positions or sparse counts.
        fp2: Fingerprint as a list of set bit positions or sparse counts.
        fp_size: Fingerprint size in bits (not read: every position counts
            once whatever its value).

    Returns:
        float: Overlap coefficient in [0.0, 1.0].
    """
    if fp1 is None or fp2 is None:
        return 0.0
    return _smsd._set_overlap(_positions(fp1), _positions(fp2))


def tanimoto_coefficient(fp1, fp2, fp_size=2048):
    """Tanimoto coefficient (Jaccard index): |A∩B| / |A∪B|.

    The standard fingerprint similarity metric in cheminformatics. Returns
    1.0 when A and B are identical, including two empty fingerprints
    (identical empty sets). Symmetric and bounded in [0.0, 1.0]. The native
    metric is shared with C++ and Java.

    Accepts binary fingerprints (lists of set bit positions) or sparse count
    fingerprints (lists of (position, count) pairs).

    Args:
        fp1: Fingerprint as a list of set bit positions or sparse counts.
        fp2: Fingerprint as a list of set bit positions or sparse counts.
        fp_size: Fingerprint size in bits (not read: every position counts
            once whatever its value).

    Returns:
        float: Tanimoto similarity in [0.0, 1.0].
    """
    if fp1 is None or fp2 is None:
        return 0.0
    return _smsd._set_tanimoto(_positions(fp1), _positions(fp2))


def overlapCoefficient(fp1, fp2, fp_size=2048):
    """Deprecated alias of :func:`overlap_coefficient`.

    Use :func:`overlap_coefficient`.
    """
    _deprecated("overlapCoefficient", "overlap_coefficient")
    return overlap_coefficient(fp1, fp2, fp_size)


def count_overlap_coefficient(fp1, fp2):
    """Count-vector Szymkiewicz-Simpson overlap: sum(min(a,b)) / min(sum(a), sum(b)).

    Generalisation of the binary overlap coefficient for integer count
    fingerprints (e.g. ECFP counts, topological torsion counts), given as
    (position, count) pairs.
    """
    a = _counts(fp1)
    b = _counts(fp2)
    inter = sum(min(a.get(k, 0), b[k]) for k in b)
    mn = min(sum(a.values()), sum(b.values()))
    if mn == 0:
        return 0.0
    return inter / mn


def count_tanimoto_coefficient(fp1, fp2):
    """Count-vector Tanimoto (generalised Jaccard): sum(min(a,b)) / sum(max(a,b)).

    The standard count-based fingerprint similarity metric, for count
    fingerprints given as (position, count) pairs.
    """
    a = _counts(fp1)
    b = _counts(fp2)
    keys = set(a) | set(b)
    inter = sum(min(a.get(k, 0), b.get(k, 0)) for k in keys)
    union = sum(max(a.get(k, 0), b.get(k, 0)) for k in keys)
    if union == 0:
        return 0.0
    return inter / union


def count_overlapCoefficient(fp1, fp2):
    """Deprecated alias of :func:`count_overlap_coefficient`.

    Use :func:`count_overlap_coefficient`.
    """
    _deprecated("count_overlapCoefficient", "count_overlap_coefficient")
    return count_overlap_coefficient(fp1, fp2)


def _sparse(fp):
    """Whether a count fingerprint is given as (position, count) pairs or a dict rather than densely."""
    return isinstance(fp, dict) or (len(fp) > 0 and isinstance(fp[0], (tuple, list)))


def _dense_pair(fp1, fp2):
    """Two count vectors of equal length; sparse pairs or dicts are expanded, dense lists pass through."""
    if not (_sparse(fp1) or _sparse(fp2)):
        return fp1, fp2
    a = dict(fp1) if isinstance(fp1, dict) else _counts(fp1)
    b = dict(fp2) if isinstance(fp2, dict) else _counts(fp2)
    size = max([*a, *b], default=-1) + 1
    if size == 0:
        return [], []
    return _smsd.counts_to_array(a, size), _smsd.counts_to_array(b, size)


def count_dice(fp1, fp2):
    """Count-vector Dice: 2*sum(min(a,b)) / (sum(a) + sum(b)).

    Accepts dense count vectors of equal length, (position, count) pairs or
    ``{position: count}`` dictionaries.
    """
    return _smsd.count_dice(*_dense_pair(fp1, fp2))


def count_cosine(fp1, fp2):
    """Count-vector cosine: dot(a,b) / (|a| * |b|).

    Accepts dense count vectors of equal length, (position, count) pairs or
    ``{position: count}`` dictionaries.
    """
    return _smsd.count_cosine(*_dense_pair(fp1, fp2))


def to_hex(fp, fp_size=2048):
    """Convert a fingerprint (list of set-bit positions) to a hex string.

    Useful for storing fingerprints in databases or sending them through REST APIs.

    Args:
        fp (list[int]): Sorted list of set bit positions.
        fp_size (int): Total fingerprint size in bits (default 2048).

    Returns:
        str: Hex string representation of the fingerprint.

    Example::

        import smsd

        fp = smsd.fingerprint_from_smiles("CCO", mode="path", fp_size=2048)
        hex_str = smsd.to_hex(fp, fp_size=2048)
        fp_back = smsd.from_hex(hex_str)
        assert fp_back == fp
    """
    _checked_int(fp_size, "fp_size", 1)
    return _smsd.to_hex(_positions(fp), fp_size)


def from_hex(hex_str):
    """Convert a hex string back to a fingerprint (list of set-bit positions).

    Inverse of :func:`to_hex`.

    Args:
        hex_str (str): Hex string produced by :func:`to_hex`.

    Returns:
        list[int]: Sorted list of set bit positions.

    Raises:
        ValueError: the string holds a character that is not a hexadecimal digit.

    Example::

        import smsd

        fp = smsd.fingerprint_from_smiles("CCO")
        hex_str = smsd.to_hex(fp)
        fp_back = smsd.from_hex(hex_str)
        assert fp_back == fp
    """
    if not isinstance(hex_str, str) or not hex_str:
        return []
    if not _HEX_DIGITS.fullmatch(hex_str):
        raise ValueError("hex_str must contain only hexadecimal digits")
    return _smsd.from_hex(hex_str)


def to_binary_string(fp, fp_size=2048):
    """Convert a fingerprint to a binary string of '0' and '1' characters.

    Args:
        fp (list[int]): Sorted list of set bit positions.
        fp_size (int): Total fingerprint size in bits (default 2048).

    Returns:
        str: Binary string of length ``fp_size``.

    Example::

        import smsd

        fp = smsd.fingerprint_from_smiles("CCO")
        bits = smsd.to_binary_string(fp, fp_size=2048)
        assert len(bits) == 2048
        assert all(c in ('0', '1') for c in bits)
    """
    _checked_int(fp_size, "fp_size", 1)
    return _smsd.to_binary_string(_positions(fp), fp_size)


def counts_to_array(counts_dict, fp_size):
    """Convert a sparse count dict ``{bit_position: count}`` into a dense list
    of length ``fp_size``.

    Useful for serialising count-based fingerprints to fixed-length vectors
    for NumPy, database columns or machine-learning feature arrays.

    Positions outside ``[0, fp_size)`` are ignored.

    Args:
        counts_dict (dict): Sparse count map ``{int: int}``, or a list of
            ``(position, count)`` pairs.
        fp_size (int): Target array length.

    Returns:
        list[int]: Dense count array of length ``fp_size``.

    Example::

        import smsd

        arr = smsd.counts_to_array({42: 3, 1000: 1}, 2048)
        assert arr[42] == 3
        assert arr[1000] == 1
    """
    _checked_int(fp_size, "fp_size", 1)
    if counts_dict is None:
        counts_dict = {}
    elif not isinstance(counts_dict, dict):
        counts_dict = dict(counts_dict)
    return _smsd.counts_to_array(counts_dict, fp_size)
