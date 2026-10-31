# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""SVG highlights of an existing MCS result, without another search."""

from html import escape

from ._depiction import depict_pair
from .search import SearchResult, _convert, _native_correspondence

__all__ = ["mcs_svg"]


def mcs_svg(query, target, result, *, opts=None, **kwargs):
    """Render the two molecules with the best result's atom correspondence.

    The native SVG renderer draws mapped atoms and every bond joining them.
    A selected-edge result that omits such bonds raises ``ValueError`` rather
    than highlighting an edge outside the common graph. Use the same molecules
    and atom numbering as the search. No optional drawing toolkit is required.

    The SVG title and ``data-search-status``, ``data-optimal`` and
    ``data-complete`` attributes retain proof status, including stopped searches.
    A found mapping need not be maximum: inspect ``result.optimal`` separately.
    Depiction options are those of :func:`smsd.depict_pair`.

    Example::

        >>> import smsd
        >>> from smsd.highlight import mcs_svg
        >>> result = smsd.MCS.find("CCO").match("CCN")
        >>> svg = mcs_svg("CCO", "CCN", result)
        >>> assert '<svg' in svg and 'data-search-status="COMPLETED"' in svg
    """
    if not isinstance(result, SearchResult):
        raise TypeError("result must be a SearchResult")
    query_graph, query_indices = _convert(query)
    target_graph, target_indices = _convert(target)
    mapping = {}
    if result.best is not None:
        match = _native_correspondence(result.best, query_indices, target_indices)
        mapping = match.atom_mapping
        query_bonds = {tuple(sorted((first, second))) for first, second, _, _ in match.bond_mapping}
        target_bonds = {tuple(sorted((first, second))) for _, _, first, second in match.bond_mapping}
        _check_selected_bonds(query_graph, set(mapping), query_bonds)
        _check_selected_bonds(target_graph, set(mapping.values()), target_bonds)
    svg = depict_pair(query_graph, target_graph, mapping, opts=opts, **kwargs)
    status = result.status.name
    optimal, complete = str(result.optimal).lower(), str(result.complete).lower()
    attributes = (f'data-search-status="{escape(status, quote=True)}" '
                  f'data-optimal="{optimal}" data-complete="{complete}" ')
    svg = svg.replace("<svg ", "<svg " + attributes, 1)
    title = (f"MCS: {len(mapping)} mapped atoms; status {status}; "
             f"optimal {optimal}; complete {complete}")
    start = svg.index(">", svg.index("<svg ")) + 1
    return svg[:start] + "\n<title>" + escape(title) + "</title>" + svg[start:]


def _check_selected_bonds(graph, mapped_atoms, selected_bonds):
    """The renderer cannot distinguish a selected edge from another mapped-atom edge."""
    for atom in mapped_atoms:
        if not 0 <= atom < graph.n:
            raise ValueError("The result's atom mapping does not belong to these molecules")
        for neighbour in graph.neighbors[atom]:
            if neighbour in mapped_atoms and tuple(sorted((atom, neighbour))) not in selected_bonds:
                raise ValueError("The SVG pair renderer cannot depict this selected-edge MCS; "
                                 "render result.best.to_graph(query) to show its selected bonds")
