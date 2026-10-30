# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""RDKit and Open Babel conversion with atom-index bookkeeping.

Both toolkits are optional dependencies, imported only when a conversion
needs them.
"""

import threading
import weakref

from ._molecule import MolGraph, _ensure_mol, parse_smiles, to_smiles

__all__ = [
    "clear_cache", "from_openbabel", "from_rdkit", "get_index_map", "to_openbabel", "to_rdkit",
    "translate_mapping",
]


def _require_rdkit(function, *, drawing=False):
    """Return the RDKit chemistry module, naming the calling function if it is missing.

    ``drawing`` also requires the coordinate and drawing modules.
    """
    try:
        from rdkit import Chem
        if drawing:
            from rdkit.Chem import AllChem, Draw  # noqa: F401
            from rdkit.Chem.Draw import rdMolDraw2D  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            f"RDKit is required for {function}(). "
            "Install with: pip install rdkit"
        ) from exc
    return Chem


# Each live molecule is cached on its own, so reordered molecules cannot share atom-index
# metadata. Weak keys release entries when their inputs are discarded.
_rdkit_cache_lock = threading.Lock()
_rdkit_cache = weakref.WeakKeyDictionary()  # live RDKit Mol -> (graph, signature)
_index_maps = weakref.WeakKeyDictionary()   # live graph -> tuple: SMSD index -> RDKit index


def _compute_index_map(rdkit_mol, canonical_smi):
    """Compute the deterministic SMSD-to-RDKit atom index mapping.

    RDKit records the order in which atoms appear in a canonical SMILES string
    in the ``_smilesAtomOutputOrder`` property, and SMSD parses that string in
    the same order, so ``index_map[smsd_idx] = rdkit_idx``.

    The caller must already have called ``MolToSmiles()`` on the molecule,
    which sets the property.
    """
    from rdkit import Chem

    n = rdkit_mol.GetNumAtoms()

    try:
        output_order = list(rdkit_mol.GetPropsAsDict(includePrivate=True, includeComputed=True)
                            .get('_smilesAtomOutputOrder', []))
    except Exception:
        output_order = []

    if len(output_order) == n and sorted(output_order) == list(range(n)):
        # output_order[smiles_position] = rdkit_atom_idx, and smsd_idx = smiles_position.
        return output_order
    # Otherwise match the molecule against the canonical SMILES parsed afresh.
    mol_from_smi = Chem.MolFromSmiles(canonical_smi)
    if mol_from_smi is not None:
        match = rdkit_mol.GetSubstructMatch(mol_from_smi, useChirality=True)
        if match and len(match) == n:
            return list(match)
    raise ValueError("Cannot determine the SMSD-to-RDKit atom-index mapping")


def from_rdkit(mol, use_cache=True):
    """Convert an RDKit Mol to an SMSD MolGraph through a SMILES round trip.

    The SMSD-to-RDKit index mapping is kept with the graph and used by
    :func:`translate_mapping` and by query results.

    Each distinct RDKit Mol object gets its own graph and mapping. Weak cache
    keys release discarded molecules, and in-place changes to a molecule
    invalidate its cached conversion. Index metadata stays available while a
    converted graph is alive, including after :func:`clear_cache`. Changing a
    returned graph in place, for example with :func:`kekulize`, drops it from
    the cache, so the next conversion of the same molecule starts afresh.

    Example::

        from rdkit import Chem
        import smsd

        mol = Chem.MolFromSmiles("c1ccccc1")
        g = smsd.from_rdkit(mol)
    """
    Chem = _require_rdkit("from_rdkit")

    if mol is None:
        raise TypeError("mol must not be None")
    if not isinstance(mol, Chem.Mol):
        raise TypeError("mol must be an RDKit Mol")
    if mol.GetNumAtoms() == 0:
        raise ValueError("mol has 0 atoms; cannot convert empty molecule")
    if mol.GetNumHeavyAtoms() == 0:
        raise ValueError("mol contains only hydrogen atoms; no heavy atoms to convert")

    if use_cache:
        signature = mol.ToBinary()
        with _rdkit_cache_lock:
            cached = _rdkit_cache.get(mol)
            if cached is not None and cached[1] == signature:
                return cached[0]

    smi = Chem.MolToSmiles(mol)
    if not smi:
        raise ValueError("RDKit returned empty SMILES for the molecule")

    # Each RDKit Mol gets its own graph: two molecules with the same canonical SMILES
    # but different atom orders need separate index maps, which live with the graph.
    graph = parse_smiles(smi)
    indices = tuple(_compute_index_map(mol, smi))
    # The signature is taken after MolToSmiles(), which may change the pickled properties.
    signature = mol.ToBinary() if use_cache else None

    with _rdkit_cache_lock:
        _index_maps[graph] = indices
        if use_cache:
            if len(_rdkit_cache) >= 10000:
                _rdkit_cache.clear()
            _rdkit_cache[mol] = (graph, signature)

    return graph


def _index_tuple(graph):
    """The stored SMSD-to-RDKit index tuple of a graph from :func:`from_rdkit`, or None."""
    with _rdkit_cache_lock:
        return _index_maps.get(graph)


def _forget(graph):
    """Drop the cached conversions that return ``graph``, after it has been changed in place."""
    with _rdkit_cache_lock:
        if graph not in _index_maps:
            return
        for mol in [mol for mol, entry in _rdkit_cache.items() if entry[0] is graph]:
            del _rdkit_cache[mol]


def get_index_map(g):
    """The SMSD-to-RDKit atom index mapping of a graph made by :func:`from_rdkit`.

    Returns a list where ``result[smsd_idx] = rdkit_idx``, or None when the
    graph did not come from :func:`from_rdkit`.

    Example::

        >>> import smsd
        >>> from rdkit import Chem
        >>> graph = smsd.from_rdkit(Chem.MolFromSmiles("OCC"))
        >>> assert sorted(smsd.get_index_map(graph)) == [0, 1, 2]
        >>> assert smsd.get_index_map(smsd.parse_smiles("CCO")) is None
    """
    indices = _index_tuple(g)
    return list(indices) if indices is not None else None


def translate_mapping(mapping, g_query=None, g_target=None):
    """Translate an MCS mapping from SMSD indices to RDKit indices.

    Args:
        mapping: dict of {smsd_query_idx: smsd_target_idx} from a Match
        g_query: the query MolGraph (from from_rdkit())
        g_target: the target MolGraph (from from_rdkit())

    Returns:
        dict of {rdkit_query_idx: rdkit_target_idx}

    Raises:
        IndexError: an index is negative or outside its molecule, as in
            search results.

    Example::

        >>> import smsd
        >>> from rdkit import Chem
        >>> graph = smsd.from_rdkit(Chem.MolFromSmiles("OCC"))
        >>> indices = smsd.get_index_map(graph)
        >>> assert smsd.translate_mapping({0: 1}, graph, graph) == {indices[0]: indices[1]}
    """
    if mapping is None or len(mapping) == 0:
        return {}
    if not isinstance(mapping, dict):
        raise TypeError(f"mapping must be a dict, got {type(mapping).__name__}")

    def numbering(graph):
        if graph is None:
            return None, None
        indices = _index_tuple(graph)
        return indices, (len(indices) if indices is not None else getattr(graph, "n", None))

    q_map, q_size = numbering(g_query)
    t_map, t_size = numbering(g_target)

    def translate(index, index_map, size, role):
        if index < 0 or (size is not None and index >= size):
            raise IndexError(f"{role} atom index out of range")
        return index_map[index] if index_map is not None else index

    result = {}
    for smsd_q, smsd_t in mapping.items():
        if not isinstance(smsd_q, int) or not isinstance(smsd_t, int):
            raise TypeError(f"mapping keys/values must be int, got ({type(smsd_q).__name__}, {type(smsd_t).__name__})")
        result[translate(smsd_q, q_map, q_size, "Query")] = translate(smsd_t, t_map, t_size, "Target")
    return result


def clear_cache():
    """Clear reusable conversions without invalidating live graphs' index maps.

    Example::

        >>> import smsd
        >>> from rdkit import Chem
        >>> graph = smsd.from_rdkit(Chem.MolFromSmiles("OCC"))
        >>> indices = smsd.get_index_map(graph)
        >>> smsd.clear_cache()
        >>> assert smsd.get_index_map(graph) == indices
    """
    with _rdkit_cache_lock:
        _rdkit_cache.clear()


def to_rdkit(mol_or_smiles):
    """Convert an SMSD MolGraph or SMILES string to an RDKit Mol.

    Example::

        import smsd
        g = smsd.parse_smiles("c1ccccc1")
        rdkit_mol = smsd.to_rdkit(g)

        # RDKit can then draw it or compute descriptors.
        from rdkit.Chem import Draw
        Draw.MolToImage(rdkit_mol)
    """
    Chem = _require_rdkit("to_rdkit")
    if isinstance(mol_or_smiles, str):
        smi = mol_or_smiles
    elif isinstance(mol_or_smiles, MolGraph):
        smi = to_smiles(mol_or_smiles)
    else:
        raise TypeError(f"Expected MolGraph or SMILES string, got {type(mol_or_smiles)}")
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smi}")
    return mol


def from_openbabel(mol):
    """Convert an Open Babel molecule into an SMSD MolGraph.

    Example::

        >>> import smsd
        >>> from openbabel import pybel
        >>> graph = smsd.from_openbabel(pybel.readstring("smi", "CCO"))
        >>> assert graph.n == 3
    """
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
    """Convert a MolGraph, SMILES string or RDKit molecule into a pybel Molecule.

    Example::

        >>> import smsd
        >>> molecule = smsd.to_openbabel("CCO")
        >>> assert molecule.OBMol.NumAtoms() == 3
    """
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

    return pybel.readstring("smi", to_smiles(_ensure_mol(mol)))
