# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""SMILES, SMARTS, MOL and SDF text, canonical forms, graph utilities and input conversion."""

import os
import warnings

from . import _smsd

__all__ = [
    "canonical_hash", "canonical_smiles", "count_components", "export_sdf", "extract_subgraph",
    "murcko_scaffold", "parse_mapped_smiles", "parse_smiles", "prewarm", "prewarm_graph",
    "read_mol_block", "read_mol_file", "read_molfile", "read_sdf", "same_canonical_graph",
    "split_components", "strip_atom_maps", "to_smarts", "to_smiles", "translate_to_atom_ids",
    "write_mol_block", "write_mol_block_v3000", "write_molfile", "write_sdf", "write_sdf_record",
]

MolGraph = _smsd.MolGraph

same_canonical_graph = _smsd.same_canonical_graph
translate_to_atom_ids = _smsd.translate_to_atom_ids
read_mol_block = _smsd.read_mol_block


def _deprecated(name, replacement):
    """Warn, at the caller of ``smsd.name()``, that ``smsd.replacement()`` supersedes it."""
    warnings.warn(f"smsd.{name}() is deprecated; use smsd.{replacement}()",
                  DeprecationWarning, stacklevel=3)


def _checked_int(value, name, minimum):
    """``value`` when it is an int of at least ``minimum``; a bool is not an integer here."""
    kind = "positive" if minimum > 0 else "non-negative"
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValueError(f"{name} must be a {kind} integer, got {value!r}")
    if value < minimum:
        raise ValueError(f"{name} must be {kind}")
    return value


def _ensure_mol_ex(obj):
    """A MolGraph for a MolGraph, SMILES string or RDKit Mol, with the RDKit Mol when there is one.

    Returns:
        tuple: (MolGraph, rdkit_mol_or_None)
    """
    if isinstance(obj, MolGraph):
        return obj, None
    if isinstance(obj, str):
        if not obj:
            raise ValueError("SMILES string must not be empty")
        return parse_smiles(obj), None
    # An RDKit Mol is recognised by its interface, so RDKit is imported only when it is used.
    if hasattr(obj, "GetNumAtoms") and hasattr(obj, "GetAtomWithIdx"):
        from ._interop import from_rdkit
        return from_rdkit(obj), obj
    raise TypeError(f"Expected MolGraph, SMILES string, or RDKit Mol, got {type(obj).__name__}")


def _ensure_mol(obj):
    """A MolGraph for a MolGraph, SMILES string or RDKit Mol."""
    return obj if type(obj) is MolGraph else _ensure_mol_ex(obj)[0]


def _molecule_first(native):
    """``native`` accepting a MolGraph, SMILES string or RDKit Mol as its first argument."""
    def call(mol, *args, **kwargs):
        return native(mol if type(mol) is MolGraph else _ensure_mol_ex(mol)[0], *args, **kwargs)
    call.__name__ = call.__qualname__ = native.__name__
    call.__doc__ = native.__doc__
    return call


# Graph utilities and text writers; each accepts a MolGraph, SMILES string or RDKit Mol.
to_smarts = _molecule_first(_smsd.to_smarts)
extract_subgraph = _molecule_first(_smsd.extract_subgraph)
murcko_scaffold = _molecule_first(_smsd.murcko_scaffold)
count_components = _molecule_first(_smsd.count_components)
split_components = _molecule_first(_smsd.split_components)
write_mol_block = _molecule_first(_smsd.write_mol_block)
write_mol_block_v3000 = _molecule_first(_smsd.write_mol_block_v3000)
write_sdf_record = _molecule_first(_smsd.write_sdf_record)


def parse_smiles(smiles, *, lenient=False):
    """Parse a SMILES string into a MolGraph.

    An empty string gives an empty molecule. With ``lenient`` the parser recovers
    from malformed text that would otherwise raise ``ValueError``.

    Example::

        >>> import smsd
        >>> graph = smsd.parse_smiles("CCO")
        >>> assert graph.n == 3 and graph.atomic_num == [6, 6, 8]
    """
    if smiles == "":
        return MolGraph()
    if lenient:
        options = _smsd.ParseOptions()
        options.lenient = True
        return _smsd.parse_smiles_lenient(smiles, options)
    return _smsd.parse_smiles(smiles)


def to_smiles(mol):
    """Canonical SMILES of a MolGraph, SMILES string or RDKit Mol.

    The canonical order is defined by the writer and may differ between languages.

    Example::

        >>> import smsd
        >>> graph = smsd.parse_smiles("OCC")
        >>> restored = smsd.parse_smiles(smsd.to_smiles(graph))
        >>> assert smsd.same_canonical_graph(graph, restored)
    """
    return _smsd.to_smiles(_ensure_mol(mol))


def canonical_smiles(mol):
    """Deprecated alias of :func:`to_smiles`.

    Use :func:`to_smiles`, which returns the same text.

    Example::

        >>> import smsd, warnings
        >>> with warnings.catch_warnings():
        ...     warnings.simplefilter("ignore", DeprecationWarning)
        ...     text = smsd.canonical_smiles("OCC")
        >>> assert text == smsd.to_smiles("OCC")
    """
    _deprecated("canonical_smiles", "to_smiles")
    return to_smiles(mol)


def canonical_hash(mol):
    """The canonical hash of a molecule as a Python int.

    The hash covers connectivity, bond order, aromaticity, ring membership,
    formal charge, mass number and stereo (tetrahedral chirality and E/Z
    double-bond configuration), as recorded in the graph.

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
    return _ensure_mol(mol).get_canonical_hash()


def strip_atom_maps(smiles):
    """Remove atom map numbers from a SMILES string.

    Distinguishes an atom map ``:N`` inside brackets from an aromatic ring
    closure ``:N`` outside them, which a plain pattern replacement such as
    ``re.sub(r':\\d+', '', smi)`` would also remove.

    Args:
        smiles: SMILES string with atom map numbers.

    Returns:
        str: SMILES with atom maps removed and ring closures preserved.

    Example::

        import smsd

        clean = smsd.strip_atom_maps("[CH3:1][C:2](=[O:3])[OH:4]")
        # Returns: "[CH3][C](=[O])[OH]"

        # Preserves aromatic ring closures:
        clean = smsd.strip_atom_maps("[Cl:1][C:2]1:[CH:3]:[CH:4]:1")
        # Returns: "[Cl][C]1:[CH]:[CH]:1"  (ring :1 preserved)
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
            # Skip the digits of an atom map number.
            i += 1
            while i < n and smiles[i].isdigit():
                i += 1
        else:
            result.append(c)
            i += 1
    return ''.join(result)


def parse_mapped_smiles(smiles):
    """Parse a SMILES string that contains atom map numbers.

    Removes the atom maps first, then parses; this is the safe way to turn
    reaction SMILES fragments into graphs for MCS searches.

    Args:
        smiles: SMILES string, which may contain atom map numbers such as ``:1``.

    Returns:
        MolGraph: The parsed molecular graph.

    Example::

        import smsd

        g = smsd.parse_mapped_smiles("[CH3:1][C:2](=[O:3])[OH:4]")
        print(g.n)  # 4
    """
    return parse_smiles(strip_atom_maps(smiles))


def prewarm(mol):
    """Compute and cache the common lazy invariants of a molecule once; returns its graph.

    Example::

        >>> import smsd
        >>> graph = smsd.parse_smiles("CCO")
        >>> assert smsd.prewarm(graph) is graph
        >>> assert smsd.MCS.find(graph).match("CCO").best.atom_count == 3
    """
    g = _ensure_mol(mol)
    g.prewarm()
    return g


def prewarm_graph(mol):
    """Deprecated alias of :func:`prewarm` that returns None.

    Use :func:`prewarm` or :meth:`MolGraph.prewarm`.

    Example::

        >>> import smsd, warnings
        >>> graph = smsd.parse_smiles("CCO")
        >>> with warnings.catch_warnings():
        ...     warnings.simplefilter("ignore", DeprecationWarning)
        ...     assert smsd.prewarm_graph(graph) is None
        >>> assert graph.n == 3
    """
    _deprecated("prewarm_graph", "prewarm")
    _ensure_mol(mol).prewarm()


def read_mol_file(filename):
    """Read a MOL file (V2000 or V3000) into a MolGraph; one that cannot be opened or read raises ValueError.

    Example::

        >>> import smsd, warnings
        >>> from pathlib import Path
        >>> from tempfile import TemporaryDirectory
        >>> with TemporaryDirectory() as directory, warnings.catch_warnings():
        ...     warnings.simplefilter("ignore", DeprecationWarning)
        ...     destination = Path(directory) / "molecule.mol"
        ...     _ = destination.write_text(smsd.write_mol_block("CCO"), encoding="utf-8")
        ...     assert smsd.read_mol_file(destination).n == 3
    """
    return _smsd.read_mol_file(os.fspath(filename))


def read_sdf(filename):
    """Read every record of an SDF file at once; malformed records become empty graphs.

    Example::

        >>> import smsd
        >>> from pathlib import Path
        >>> from tempfile import TemporaryDirectory
        >>> with TemporaryDirectory() as directory:
        ...     destination = Path(directory) / "molecules.sdf"
        ...     smsd.write_sdf(["CCO", "c1ccccc1"], destination)
        ...     assert [graph.n for graph in smsd.read_sdf(destination)] == [3, 6]
    """
    return _smsd.read_sdf(os.fspath(filename))


def read_molfile(filename):
    """Deprecated alias of :func:`read_mol_file` that reads the file as UTF-8 text.

    Use :func:`read_mol_file`.

    Example::

        >>> import smsd, warnings
        >>> from pathlib import Path
        >>> from tempfile import TemporaryDirectory
        >>> with TemporaryDirectory() as directory, warnings.catch_warnings():
        ...     warnings.simplefilter("ignore", DeprecationWarning)
        ...     destination = Path(directory) / "molecule.mol"
        ...     _ = destination.write_text(smsd.write_mol_block("CCO"), encoding="utf-8")
        ...     assert smsd.read_molfile(destination).n == 3
    """
    _deprecated("read_molfile", "read_mol_file")
    from pathlib import Path
    return read_mol_block(Path(filename).read_text(encoding="utf-8"))


def write_sdf(molecules, filename):
    """Write molecules as SDF records, each a V2000 MOL block followed by ``$$$$``.

    Example::

        >>> import smsd, warnings
        >>> from pathlib import Path
        >>> from tempfile import TemporaryDirectory
        >>> with TemporaryDirectory() as directory, warnings.catch_warnings():
        ...     warnings.simplefilter("ignore", DeprecationWarning)
        ...     destination = Path(directory) / "molecules.sdf"
        ...     smsd.write_sdf(["CCO", "CCN"], destination)
        ...     assert [graph.n for graph in smsd.read_sdf(destination)] == [3, 3]
    """
    _smsd.write_sdf([_ensure_mol(mol) for mol in molecules], os.fspath(filename))


def export_sdf(molecules, filename):
    """Deprecated alias of :func:`write_sdf`.

    Use :func:`write_sdf`, which writes the same records.

    Example::

        >>> import smsd, warnings
        >>> from pathlib import Path
        >>> from tempfile import TemporaryDirectory
        >>> with TemporaryDirectory() as directory, warnings.catch_warnings():
        ...     warnings.simplefilter("ignore", DeprecationWarning)
        ...     destination = Path(directory) / "molecules.sdf"
        ...     smsd.export_sdf(["CCO", "CCN"], destination)
        ...     assert [graph.n for graph in smsd.read_sdf(destination)] == [3, 3]
    """
    _deprecated("export_sdf", "write_sdf")
    from pathlib import Path
    Path(filename).write_text("".join(write_sdf_record(mol) for mol in molecules), encoding="utf-8")


def write_molfile(mol, filename, *, v3000=False, sdf=False):
    """Write a molecule to disk as MOL V2000, MOL V3000 or a single SDF record.

    Example::

        >>> import smsd
        >>> from pathlib import Path
        >>> from tempfile import TemporaryDirectory
        >>> with TemporaryDirectory() as directory:
        ...     destination = Path(directory) / "molecule.mol"
        ...     smsd.write_molfile("CCO", destination, v3000=True)
        ...     assert smsd.read_mol_file(destination).n == 3
    """
    if sdf:
        text = write_sdf_record(mol)
    elif v3000:
        text = write_mol_block_v3000(mol)
    else:
        text = write_mol_block(mol)
    from pathlib import Path
    Path(filename).write_text(text, encoding="utf-8")
