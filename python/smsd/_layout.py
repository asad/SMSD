# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Coordinate generation, layout refinement and coordinate transforms.

Coordinates are sequences of (x, y) pairs: lists, tuples or arrays. Functions
that take a molecule accept a MolGraph, SMILES string or RDKit Mol.
"""

from . import _smsd
from ._molecule import _molecule_first

__all__ = [
    "Point2D", "Point3D", "align_2d", "bounding_box_2d", "canonical_orientation", "center_2d",
    "count_crossings", "force_directed_layout", "generate_coords_2d", "generate_coords_3d",
    "is_degenerate_layout", "layout_quality", "layout_sssr", "match_template", "mirror_x", "mirror_y",
    "normalise_bond_length", "reduce_crossings", "resolve_overlaps", "rotate_2d", "scale_2d",
    "stress_majorisation", "translate_2d",
]

Point2D = _smsd.Point2D
Point3D = _smsd.Point3D

# Generation, refinement and quality measures of a molecule's layout.
generate_coords_2d = _molecule_first(_smsd.generate_coords_2d)
generate_coords_3d = _molecule_first(_smsd.generate_coords_3d)
layout_sssr = _molecule_first(_smsd.layout_sssr)
match_template = _molecule_first(_smsd.match_template)
reduce_crossings = _molecule_first(_smsd.reduce_crossings)
force_directed_layout = _molecule_first(_smsd.force_directed_layout)
stress_majorisation = _molecule_first(_smsd.stress_majorisation)
layout_quality = _molecule_first(_smsd.layout_quality)
count_crossings = _molecule_first(_smsd.count_crossings)
normalise_bond_length = _molecule_first(_smsd.normalise_bond_length)
canonical_orientation = _molecule_first(_smsd.canonical_orientation)

# Transforms of coordinates alone.
resolve_overlaps = _smsd.resolve_overlaps
is_degenerate_layout = _smsd.is_degenerate_layout
translate_2d = _smsd.translate_2d
rotate_2d = _smsd.rotate_2d
scale_2d = _smsd.scale_2d
mirror_x = _smsd.mirror_x
mirror_y = _smsd.mirror_y
center_2d = _smsd.center_2d
align_2d = _smsd.align_2d
bounding_box_2d = _smsd.bounding_box_2d
