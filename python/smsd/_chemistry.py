# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Stereo descriptors, ring perception, aromaticity and atom chemistry."""

from . import _smsd
from ._interop import _forget, from_rdkit
from ._molecule import MolGraph, _deprecated, _ensure_mol, _molecule_first

__all__ = [
    "CIPDescriptors", "EZLabel", "RSLabel", "assign_cip", "assign_ez", "assign_rs",
    "classify_pharmacophore", "compute_sssr", "dearomatize", "implicit_h", "kekulize",
    "perceive_aromaticity",
]

RSLabel = _smsd.RSLabel
EZLabel = _smsd.EZLabel
CIPDescriptors = _smsd.CIPDescriptors

classify_pharmacophore = _smsd.classify_pharmacophore
implicit_h = _smsd.implicit_h

# Each accepts a MolGraph, SMILES string or RDKit Mol.
assign_cip = _molecule_first(_smsd.assign_cip)
compute_sssr = _molecule_first(_smsd.compute_sssr)


def assign_rs(mol):
    """Assign CIP R/S descriptors to all tetrahedral stereocentres.

    Accepts a MolGraph or SMILES string.

    Returns:
        dict: Mapping from atom index (int) to descriptor ('R' or 'S').
            Only stereocentres with resolvable, distinct CIP priorities
            are included. Atoms without chirality annotation are omitted.

    Example::

        import smsd

        result = smsd.assign_rs("N[C@@H](C)C(=O)O")  # L-alanine
        # result: {1: 'S'}
    """
    desc = assign_cip(mol)
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

        import smsd

        result = smsd.assign_ez("C/C=C/C")  # (E)-2-butene
        # result: {(1, 2): 'E'}
    """
    desc = assign_cip(mol)
    result = {}
    for a1, a2, label in desc.ez_bonds:
        if label == EZLabel.E:
            result[(a1, a2)] = 'E'
        elif label == EZLabel.Z:
            result[(a1, a2)] = 'Z'
    return result


def _own_graph(mol):
    """A graph that may be changed in place; foreign molecules are converted afresh."""
    if isinstance(mol, MolGraph):
        # A changed graph does not represent the RDKit molecule it came from.
        _forget(mol)
        return mol
    if isinstance(mol, str):
        return _ensure_mol(mol)
    return from_rdkit(mol, use_cache=False)


def perceive_aromaticity(mol, *, model=_smsd.AromaticityModel.DAYLIGHT_LIKE):
    """Recompute ring membership and aromaticity on a MolGraph."""
    g = _own_graph(mol)
    g.perceive_aromaticity(model)
    return g


def kekulize(mol):
    """Convert an aromatic MolGraph into an explicit Kekule form."""
    g = _own_graph(mol)
    g.kekulize()
    return g


def dearomatize(mol):
    """Deprecated alias of :func:`kekulize`.

    Use :func:`kekulize`, which performs the same conversion.
    """
    _deprecated("dearomatize", "kekulize")
    g = _own_graph(mol)
    g.dearomatize()
    return g
