# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""SVG highlights of an existing MCS or substructure result, without another search."""

from html import escape

from ._depiction import depict_mapping, depict_pair
from .search import SearchResult, Substructure, _convert, _native_correspondence

__all__ = ["mcs_svg", "substructure_svg"]


def mcs_svg(query, target, result, *, highlight="gold", opts=None, **kwargs):
    """Render the two molecules with the best result's atom correspondence.

    The native SVG renderer draws mapped atoms and each bond between them that
    the mapping carries onto a bond of the other molecule. A selected-edge
    result that omits such a bond raises ``ValueError`` rather than highlighting
    an edge outside the common graph. Use the same molecules and atom numbering
    as the search. No optional drawing toolkit is required.

    The SVG title and ``data-search-status``, ``data-optimal`` and
    ``data-complete`` attributes retain proof status, including stopped searches.
    A found mapping need not be maximum: inspect ``result.optimal`` separately.
    ``highlight`` is ``"gold"`` or ``"mint"``; depiction options are those of
    :func:`smsd.depict_pair`.

    Example::

        >>> import smsd
        >>> from smsd.highlight import mcs_svg
        >>> result = smsd.MCS.find("CCO").match("CCN")
        >>> svg = mcs_svg("CCO", "CCN", result)
        >>> assert '<svg' in svg and 'data-search-status="COMPLETED"' in svg
    """
    query_graph, target_graph, mapping = _mapped_graphs(
        query, target, result, "The SVG pair renderer cannot depict this selected-edge MCS; "
        "render result.best.to_graph(query) to show its selected bonds", paired=True)
    svg = depict_pair(query_graph, target_graph, mapping, opts=opts, highlight=highlight, **kwargs)
    return _with_status(svg, result, f"MCS: {len(mapping)} mapped atoms")


def substructure_svg(query, target, result, *, highlight="mint", opts=None, **kwargs):
    """Render the target with the best hit's atoms and the bonds the hit carries.

    The native SVG renderer draws hit atoms and every bond joining them. A hit
    whose target atoms share a bond that the query lacks raises ``ValueError``
    rather than highlighting that bond. Pass the query molecule or SMARTS pattern
    and the target as searched, with the same atom numbering. No optional
    drawing toolkit is required.

    The SVG title and ``data-search-status``, ``data-optimal`` and
    ``data-complete`` attributes retain search status, including stopped
    searches, which may have missed a hit. ``highlight`` is ``"mint"`` or
    ``"gold"``; depiction options are those of :func:`smsd.depict_mapping`.

    Example::

        >>> import smsd
        >>> from smsd.highlight import substructure_svg
        >>> result = smsd.Substructure.find("C=O").match("CC(=O)O")
        >>> svg = substructure_svg("C=O", "CC(=O)O", result)
        >>> assert '<svg' in svg and 'data-search-status="COMPLETED"' in svg
    """
    _, target_graph, mapping = _mapped_graphs(
        query, target, result, "The SVG renderer cannot depict this hit: its target atoms "
        "share a bond that the query lacks", _substructure_query)
    hit = {target_atom: query_atom for query_atom, target_atom in mapping.items()}
    svg = depict_mapping(target_graph, hit, opts=opts, highlight=highlight, **kwargs)
    return _with_status(svg, result, f"Substructure: {len(hit)} hit atoms")


def _substructure_query(query):
    """The query's graph and numbering, or no graph and the written numbering for a SMARTS
    pattern that does not parse as SMILES."""
    try:
        return _convert(query)
    except ValueError:
        if not isinstance(query, str):
            raise
        Substructure.smarts(query)
        return None, None


def _mapped_graphs(query, target, result, undrawable, convert_query=_convert, paired=False):
    """Both graphs and the best atom mapping in native indices, raising ``undrawable`` unless
    the drawing highlights exactly the bonds the result carries: those joining mapped atoms
    or, when ``paired``, those the mapping carries onto a bond of the other molecule."""
    if not isinstance(result, SearchResult):
        raise TypeError("result must be a SearchResult")
    query_graph, query_indices = convert_query(query)
    target_graph, target_indices = _convert(target)
    mapping = {}
    if result.best is not None:
        match = _native_correspondence(result.best, query_indices, target_indices)
        mapping = match.atom_mapping
        query_bonds = {tuple(sorted((first, second))) for first, second, _, _ in match.bond_mapping}
        target_bonds = {tuple(sorted((first, second))) for _, _, first, second in match.bond_mapping}
        partners = {second: first for first, second in mapping.items()}
        if query_graph is not None:
            _check_selected_bonds(query_graph, mapping, query_bonds, undrawable,
                                  target_graph if paired else None)
        _check_selected_bonds(target_graph, partners, target_bonds, undrawable,
                              query_graph if paired else None)
    return query_graph, target_graph, mapping


def _with_status(svg, result, summary):
    """``svg`` with the search status as root attributes and in a title after ``summary``."""
    status = result.status.name
    optimal, complete = str(result.optimal).lower(), str(result.complete).lower()
    attributes = (f'data-search-status="{escape(status, quote=True)}" '
                  f'data-optimal="{optimal}" data-complete="{complete}" ')
    svg = svg.replace("<svg ", "<svg " + attributes, 1)
    title = f"{summary}; status {status}; optimal {optimal}; complete {complete}"
    start = svg.index(">", svg.index("<svg ")) + 1
    return svg[:start] + "\n<title>" + escape(title) + "</title>" + svg[start:]


def _check_selected_bonds(graph, mapping, selected_bonds, undrawable, other=None):
    """Raise ``undrawable`` if the renderer would highlight a bond outside ``selected_bonds``:
    any bond joining mapped atoms or, given ``other``, one the mapping carries onto its bonds."""
    for atom, partner in mapping.items():
        if not (0 <= atom < graph.n and (other is None or 0 <= partner < other.n)):
            raise ValueError("The result's atom mapping does not belong to these molecules")
    for atom, partner in mapping.items():
        for neighbour in graph.neighbors[atom]:
            if (neighbour in mapping and tuple(sorted((atom, neighbour))) not in selected_bonds
                    and (other is None or other.bond_order(partner, mapping[neighbour]) != 0)):
                raise ValueError(undrawable)
