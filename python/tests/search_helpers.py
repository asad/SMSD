# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
"""Small search configurations shared by chemistry regression cases."""
import smsd


def configured_mcs(query, *, timeout_ms=-1, chemistry=None, fragments=False,
                   objective=None, atom_weights=(), min_fragment_atoms=1,
                   max_fragments=2147483647, **chemistry_flags):
    options = smsd.MCS.Options()
    options.timeout_ms = timeout_ms
    if chemistry is not None:
        options.chemistry = chemistry
    for name, value in chemistry_flags.items():
        if name == "match_bond_order" and isinstance(value, str):
            value = getattr(smsd.BondOrderMode, value.upper())
        if not hasattr(options.chemistry, name):
            raise ValueError(f"Unknown chemistry option: {name}")
        setattr(options.chemistry, name, value)
    options.objective = objective or smsd.MCS.Objective.MCIS
    options.connectivity = smsd.MCS.Connectivity.FRAGMENTS if fragments else smsd.MCS.Connectivity.CONNECTED
    options.atom_weights = list(atom_weights)
    options.min_fragment_atoms = min_fragment_atoms
    options.max_fragments = max_fragments
    return smsd.MCS.find(query).with_options(options)


def atom_map(result):
    return dict(result.best.atom_mapping) if result.found else {}


def selected_atoms(graph, atoms):
    """Build an induced selected-edge fragment for writer regression cases."""
    atom_mapping = {atom: index for index, atom in enumerate(atoms)}
    bonds = tuple((first, second, atom_mapping[first], atom_mapping[second])
                  for first, second, _ in graph.bonds()
                  if first in atom_mapping and second in atom_mapping)
    return smsd.Match(atom_mapping, bonds, ()).to_smiles(graph)
