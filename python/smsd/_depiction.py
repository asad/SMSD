# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Native SVG depiction, and RDKit images of search results."""

import warnings

from . import _smsd
from ._interop import _require_rdkit
from ._molecule import MolGraph, _deprecated, _ensure_mol, parse_smiles, to_smiles, write_mol_block
from .search import MCS, Substructure

__all__ = [
    "DepictOptions", "depict_mapping", "depict_mcs", "depict_pair", "depict_substructure", "depict_svg",
    "save_svg", "to_svg",
]

# SVG rendering with configurable bond and atom-label styling.
DepictOptions = _smsd.DepictOptions

_DEFAULT_HIGHLIGHT = (0.5, 0.8, 1.0)


def _options(opts, overrides):
    """Default options when none are given, with keyword overrides applied."""
    if opts is None:
        opts = DepictOptions()
    for k, v in overrides.items():
        setattr(opts, k, v)
    return opts


def _graph(mol_or_smiles):
    """A MolGraph for depiction; an empty SMILES string draws an empty molecule."""
    if isinstance(mol_or_smiles, str):
        return parse_smiles(mol_or_smiles)
    return _ensure_mol(mol_or_smiles)


def depict_svg(mol_or_smiles, opts=None, **kwargs):
    """Render a molecule as an SVG string.

    Args:
        mol_or_smiles: MolGraph, SMILES string or RDKit Mol to render.
        opts: DepictOptions; the default settings when None.
        **kwargs: Overrides applied to ``opts``, such as ``bond_length=40`` or ``width=600``.

    Returns:
        str: SVG text for a file or a notebook.
    """
    opts = _options(opts, kwargs)
    return _smsd.depict(_ensure_mol(mol_or_smiles), opts)


def depict_mapping(mol_or_smiles, mapping, opts=None, **kwargs):
    """Render a molecule with MCS or substructure atoms highlighted and numbered.

    Args:
        mol_or_smiles: MolGraph, SMILES string or RDKit Mol to render.
        mapping: Atom index to partner atom index, from an MCS or substructure result.
        opts: DepictOptions; the default settings when None.
        **kwargs: Overrides applied to ``opts``, such as ``show_atom_indices=True``.

    Returns:
        str: SVG text with the matched atoms highlighted in green.
    """
    opts = _options(opts, kwargs)
    return _smsd.depict_with_mapping(_graph(mol_or_smiles), mapping, opts)


def depict_pair(mol1, mol2, mapping, opts=None, **kwargs):
    """Render two molecules side by side with an MCS mapping highlighted.

    Args:
        mol1: First molecule, as a MolGraph, SMILES string or RDKit Mol.
        mol2: Second molecule, in the same forms.
        mapping: Atoms of ``mol1`` to atoms of ``mol2``, from a Match.
        opts: DepictOptions; the default settings when None.
        **kwargs: Overrides applied to ``opts``, such as ``width=800`` or ``height=400``.

    Returns:
        str: SVG text showing both molecules with the matched regions highlighted.
    """
    opts = _options(opts, kwargs)
    return _smsd.depict_pair(_graph(mol1), _graph(mol2), mapping, opts)


def save_svg(svg_str, path):
    """Write SVG text to a file as UTF-8.

    Args:
        svg_str: SVG text, from :func:`depict_svg`, :func:`depict_pair` or another renderer.
        path: Output file path, such as ``'molecule.svg'``.
    """
    from pathlib import Path
    Path(path).write_text(svg_str, encoding="utf-8")


def _molecule_for_depiction(graph):
    from rdkit import Chem
    return Chem.MolFromMolBlock(write_mol_block(graph), removeHs=False)


def _highlighted_pair(first, second, mapping, size, query_legend):
    """Draw two graphs side by side, highlighting the mapped atoms of each.

    ``query_legend`` may hold a ``{}`` field for the number of matched query atoms.
    """
    from rdkit.Chem import AllChem, Draw

    rdmol1 = _molecule_for_depiction(first)
    rdmol2 = _molecule_for_depiction(second)
    if rdmol1 is None or rdmol2 is None:
        raise ValueError("Could not convert molecules for depiction")

    AllChem.Compute2DCoords(rdmol1)
    AllChem.Compute2DCoords(rdmol2)

    q_atoms = list(mapping.keys()) if mapping else []
    t_atoms = list(mapping.values()) if mapping else []

    return Draw.MolsToGridImage(
        [rdmol1, rdmol2],
        molsPerRow=2,
        subImgSize=(size[0] // 2, size[1]),
        highlightAtomLists=[q_atoms, t_atoms],
        legends=[query_legend.format(len(q_atoms)), f"Target ({len(t_atoms)} matched)"],
    )


def depict_mcs(mol1, mol2, *, tautomer_aware=False, timeout_ms=10000,
               size=(600, 300), highlight_color=_DEFAULT_HIGHLIGHT):
    """Find the MCS and return an RDKit image with the matched atoms highlighted.

    Works in notebooks (displayed inline) and scripts (a PIL image).

    ``highlight_color`` is deprecated: it has no effect, and RDKit draws the
    highlights in its default colour.

    Example::

        import smsd
        img = smsd.depict_mcs("c1ccccc1", "c1ccc(O)cc1")
        img.save("mcs.png")  # or just display in a notebook
    """
    if highlight_color != _DEFAULT_HIGHLIGHT:
        warnings.warn("the highlight_color argument of smsd.depict_mcs() has no effect and is deprecated",
                      DeprecationWarning, stacklevel=2)
    _require_rdkit("depict_mcs", drawing=True)
    g1 = _ensure_mol(mol1)
    g2 = _ensure_mol(mol2)

    options = MCS.Options()
    options.chemistry.tautomer_aware = tautomer_aware
    result = MCS.find(g1).with_options(options).with_timeout(timeout_ms).match(g2)
    mapping = result.best.atom_mapping if result.found else {}

    return _highlighted_pair(g1, g2, mapping, size, "Query ({} matched)")


def depict_substructure(query, target, *, timeout_ms=10000,
                        size=(600, 300)):
    """Find a substructure and return an RDKit image with the matched atoms highlighted.

    Example::

        import smsd
        img = smsd.depict_substructure("c1ccccc1", "c1ccc(O)cc1")
    """
    _require_rdkit("depict_substructure", drawing=True)
    gq = _ensure_mol(query)
    gt = _ensure_mol(target)

    result = Substructure.find(gq).with_timeout(timeout_ms).match(gt)
    mapping = result.best.atom_mapping if result.found else {}

    return _highlighted_pair(gq, gt, mapping, size, "Query")


def to_svg(mol_or_smiles, *, size=(300, 200)):
    """Deprecated RDKit rendering of a molecule as SVG.

    Use :func:`depict_svg`, which draws natively and needs no RDKit.

    Example::

        import smsd

        svg = smsd.to_svg("c1ccccc1")
        with open("benzene.svg", "w") as f:
            f.write(svg)
    """
    _deprecated("to_svg", "depict_svg")
    Chem = _require_rdkit("to_svg", drawing=True)
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D

    if isinstance(mol_or_smiles, str):
        mol = Chem.MolFromSmiles(mol_or_smiles)
    elif isinstance(mol_or_smiles, MolGraph):
        mol = Chem.MolFromSmiles(to_smiles(mol_or_smiles))
    else:
        mol = mol_or_smiles  # an RDKit Mol

    if mol is None:
        raise ValueError("Could not parse molecule for SVG depiction")

    AllChem.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()
