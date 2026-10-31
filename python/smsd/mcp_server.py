# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Bounded molecular tools served by the official MCP SDK over standard input/output.

Importing this module neither loads the SDK nor starts a server. Molecular
arguments are text, never filenames or URLs to fetch.
"""
import json
import math
import re
import time
import xml.etree.ElementTree as ET

from . import ChemOptions, FingerprintOptions, __version__
from ._depiction import depict_pair
from ._fingerprints import (fingerprint, tanimoto_coefficient, count_tanimoto_coefficient,
                            dice, count_dice, cosine, count_cosine)
from ._molecule import _checked_int, parse_smiles, strip_atom_maps, to_smiles
from .search import MCS, Substructure

__all__ = ["create_server", "find_mcs", "find_substructure", "molecular_fingerprint",
           "search_smarts", "molecular_similarity", "screen_molecules",
           "standardise_molecule", "depict_mcs", "main"]

_LIMITS = {"text_characters": 16384, "atoms_per_molecule": 1024, "targets": 64,
           "batch_text_characters": 65536, "batch_atoms": 4096, "timeout_ms": 30000,
           "node_limit": 1000000, "matches": 32, "fingerprint_bits": 65536,
           "fingerprint_radius": 8, "path_length": 8, "path_walks": 2000000,
           "svg_bytes": 1048576, "response_bytes": 2097152, "depiction_atoms": 256,
           "standardise_atoms": 128}
_KINDS = ("mcs", "ecfp", "circular", "fcfp", "path", "torsion")
_METRICS = ("tanimoto", "dice", "cosine")


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    if len(value) > _LIMITS["text_characters"]:
        raise ValueError(f"{name} exceeds the MCP text limit of 16384 characters")
    return value


def _switch(value, name):
    if not isinstance(value, bool):
        raise ValueError(f"{name} must be True or False")
    return value


def _integer(value, name, minimum, maximum):
    _checked_int(value, name, minimum)
    if value > maximum:
        raise ValueError(f"{name} must be at most {maximum}")
    return value


def _molecule(text, name="smiles", maximum=1024):
    graph = parse_smiles(_text(text, name))
    if graph.n > maximum:
        raise ValueError(f"MCP molecular input supports at most {maximum} atoms")
    return graph


def _budget(options, timeout_ms, node_limit):
    _integer(timeout_ms, "timeout_ms", 0, _LIMITS["timeout_ms"])
    _integer(node_limit, "node_limit", 1, _LIMITS["node_limit"])
    options.timeout_ms, options.node_limit = timeout_ms, node_limit


def _chemistry(match_stereo, complete_rings_only, ring_matches_ring_only):
    chemistry = ChemOptions()
    chemistry.use_chirality = _switch(match_stereo, "match_stereo")
    chemistry.use_bond_stereo = match_stereo
    chemistry.complete_rings_only = _switch(complete_rings_only, "complete_rings_only")
    chemistry.ring_matches_ring_only = _switch(ring_matches_ring_only, "ring_matches_ring_only")
    return chemistry


def _option_value(value):
    return value.name if hasattr(value, "name") else value


def _response(output, started, native_ms=None):
    output["limits"] = dict(_LIMITS)
    output["timings"] = {"wall_elapsed_ms": (time.perf_counter_ns() - started) / 1000000,
                         "native_search_ms": native_ms,
                         "wall_scope": "tool_entry_to_pre_serialisation"}
    if len(json.dumps(output, allow_nan=False).encode("utf-8")) > _LIMITS["response_bytes"]:
        raise ValueError("MCP response exceeds the byte limit; reduce the input or result limit")
    return output


def _serialise(result, options):
    chemistry = options.chemistry
    return {"engine_version": __version__, "found": result.found, "status": result.status.name,
        "optimal": result.optimal, "complete": result.complete, "score": result.score,
        "nodes_visited": result.nodes_visited, "elapsed_ms": result.elapsed_ms,
        "proof": {"scope": "native_search", "optimal": result.optimal, "complete": result.complete},
        "budget": {"timeout_ms": options.timeout_ms, "node_limit": options.node_limit},
        "chemistry": {name: _option_value(getattr(chemistry, name)) for name in (
            "match_atom_type", "match_formal_charge", "match_isotope", "match_bond_order",
            "use_chirality", "use_bond_stereo", "complete_rings_only", "ring_matches_ring_only", "tautomer_aware")},
        "matches": [{"atom_mapping": [[query_atom, target_atom] for query_atom, target_atom
                                      in sorted(match.atom_mapping.items())],
                     "bond_mapping": [list(bond) for bond in match.bond_mapping],
                     "fragments": [list(fragment) for fragment in match.fragments],
                     "atom_count": match.atom_count, "bond_count": match.bond_count}
                    for match in result.matches]}


def _operation(scope):
    return {"engine_version": __version__, "status": "COMPLETED", "complete": True,
            "proof": {"scope": scope, "optimal": None, "search_proof_applicable": False}}


def _mcs(query_smiles, target_smiles, objective, connected, timeout_ms, node_limit,
         match_stereo, complete_rings_only, ring_matches_ring_only, maximum=1024):
    if objective not in ("mcis", "mces"):
        raise ValueError("objective must be 'mcis' or 'mces'")
    options = MCS.Options()
    _budget(options, timeout_ms, node_limit)
    options.objective = MCS.Objective.MCIS if objective == "mcis" else MCS.Objective.MCES
    options.connectivity = MCS.Connectivity.CONNECTED if _switch(connected, "connected") else MCS.Connectivity.FRAGMENTS
    options.chemistry = _chemistry(match_stereo, complete_rings_only, ring_matches_ring_only)
    query = _molecule(query_smiles, "query_smiles", maximum)
    target = _molecule(target_smiles, "target_smiles", maximum)
    return MCS.find(query).with_options(options).match(target), options, query, target


def find_mcs(query_smiles: str, target_smiles: str, objective: str = "mcis",
             connected: bool = True, timeout_ms: int = 1000, node_limit: int = 100000,
             match_stereo: bool = False, complete_rings_only: bool = False,
             ring_matches_ring_only: bool = False) -> dict:
    """Find one common graph and retain native proof and stopping flags.

    MCIS maximises atoms; MCES maximises selected bonds. Search limits apply
    after parsing. An unfinished match need not be maximum.

    Example::

        result = find_mcs("CCO", "CCN")
        assert result["found"] and result["matches"][0]["atom_count"] == 2
    """
    started = time.perf_counter_ns()
    result, options, _, _ = _mcs(query_smiles, target_smiles, objective, connected, timeout_ms,
                                node_limit, match_stereo, complete_rings_only, ring_matches_ring_only)
    output = _serialise(result, options)
    output.update(objective=objective, connected=connected)
    return _response(output, started, result.elapsed_ms)


def _substructure(query, target_smiles, query_format, induced, timeout_ms, node_limit,
                  match_stereo, complete_rings_only, ring_matches_ring_only, limit=None):
    if query_format not in ("smarts", "smiles"):
        raise ValueError("query_format must be 'smarts' or 'smiles'")
    options = Substructure.Options()
    _budget(options, timeout_ms, node_limit)
    options.induced = _switch(induced, "induced")
    options.chemistry = _chemistry(match_stereo, complete_rings_only, ring_matches_ring_only)
    query = _text(query, "query")
    prepared = Substructure.smarts(query) if query_format == "smarts" else Substructure.find(_molecule(query, "query"))
    prepared = prepared.with_options(options)
    target = _molecule(target_smiles, "target_smiles")
    result = prepared.match(target) if limit is None else prepared.match_all(target, limit=limit)
    return result, options


def find_substructure(query: str, target_smiles: str, query_format: str = "smarts",
                      induced: bool = False, timeout_ms: int = 1000, node_limit: int = 100000,
                      match_stereo: bool = False, complete_rings_only: bool = False,
                      ring_matches_ring_only: bool = False) -> dict:
    """Find one SMARTS or molecule embedding, retaining native stopping status.

    An unfinished result with found=False does not prove absence. Budgets apply
    to native search after parsing or SMARTS compilation.

    Example::

        assert find_substructure("[OX2H]", "CCO")["found"]
    """
    started = time.perf_counter_ns()
    result, options = _substructure(query, target_smiles, query_format, induced, timeout_ms,
        node_limit, match_stereo, complete_rings_only, ring_matches_ring_only)
    output = _serialise(result, options)
    output.update(query_format=query_format, induced=induced)
    return _response(output, started, result.elapsed_ms)


def search_smarts(query_smarts: str, target_smiles: str, max_matches: int = 16,
                  timeout_ms: int = 1000, node_limit: int = 100000,
                  match_stereo: bool = False, complete_rings_only: bool = False,
                  ring_matches_ring_only: bool = False) -> dict:
    """Enumerate bounded SMARTS embeddings without turning incomplete results into absence.

    Example::

        result = search_smarts("[#8]", "OCCO")
        assert len(result["matches"]) == 2
    """
    started = time.perf_counter_ns()
    _integer(max_matches, "max_matches", 1, _LIMITS["matches"])
    result, options = _substructure(query_smarts, target_smiles, "smarts", False, timeout_ms,
        node_limit, match_stereo, complete_rings_only, ring_matches_ring_only, max_matches)
    output = _serialise(result, options)
    output.update(query_format="smarts", max_matches=max_matches, truncated=not result.complete)
    return _response(output, started, result.elapsed_ms)


def _fingerprint_parameters(kind, fp_size, radius, path_length, counts):
    _checked_int(fp_size, "fp_size", 1)
    _checked_int(radius, "radius", 0)
    _checked_int(path_length, "path_length", 1)
    if fp_size > _LIMITS["fingerprint_bits"] or radius > 8 or path_length > 8:
        raise ValueError("MCP fingerprint limits are fp_size <= 65536, radius <= 8 and path_length <= 8")
    _switch(counts, "counts")
    if kind not in _KINDS:
        raise ValueError("Unknown fingerprint kind")
    if counts and kind == "mcs":
        raise ValueError("The MCS fingerprint has no count mode")
    return FingerprintOptions(family=kind, size=fp_size, radius=radius,
                              path_length=path_length, counts=True) if counts else None


def _walk_work(graph, kind, path_length):
    if kind not in ("mcs", "path", "torsion"):
        return 0
    neighbors = graph.neighbors
    walks = {(atom, neighbor): 1 for atom, adjacent in enumerate(neighbors) for neighbor in adjacent}
    total = graph.n
    for _ in range(3 if kind == "torsion" else path_length):
        total += sum(walks.values())
        if total > _LIMITS["path_walks"]:
            raise ValueError("MCP fingerprint path work bound exceeds 2000000 walks")
        outgoing = [sum(walks[(atom, neighbor)] for neighbor in adjacent) for atom, adjacent in enumerate(neighbors)]
        walks = {(atom, neighbor): outgoing[neighbor] - walks[(neighbor, atom)] for atom, neighbor in walks}
    return total


def molecular_fingerprint(smiles: str, kind: str = "ecfp", fp_size: int = 2048,
                          radius: int = 2, path_length: int = 7, counts: bool = False) -> dict:
    """Compute native set-bit positions or sparse position/count pairs, not an MCS proof.

    Example::

        assert molecular_fingerprint("CCO", fp_size=64)["values"]
    """
    started = time.perf_counter_ns()
    options = _fingerprint_parameters(kind, fp_size, radius, path_length, counts)
    graph = _molecule(smiles)
    _walk_work(graph, kind, path_length)
    values = fingerprint(graph, kind=kind, fp_size=fp_size, radius=radius, path_length=path_length, options=options)
    output = _operation("fingerprint_generation")
    output.update(kind=kind, fp_size=fp_size, radius=radius, path_length=path_length, counts=counts,
                  values=[list(pair) for pair in values] if counts else list(values))
    return _response(output, started)


def _fingerprints(smiles_values, kind, fp_size, radius, path_length, counts):
    options = _fingerprint_parameters(kind, fp_size, radius, path_length, counts)
    graphs = [_molecule(value) for value in smiles_values]
    if sum(graph.n for graph in graphs) > _LIMITS["batch_atoms"]:
        raise ValueError("MCP batch exceeds 4096 atoms")
    if sum(_walk_work(graph, kind, path_length) for graph in graphs) > _LIMITS["path_walks"]:
        raise ValueError("MCP batch fingerprint path work bound exceeds 2000000 walks")
    return [fingerprint(graph, kind=kind, fp_size=fp_size, radius=radius,
                        path_length=path_length, options=options) for graph in graphs]


def _similarity(first, second, metric, counts, fp_size):
    functions = {"tanimoto": (tanimoto_coefficient, count_tanimoto_coefficient),
                 "dice": (dice, count_dice), "cosine": (cosine, count_cosine)}
    if metric not in functions:
        raise ValueError("metric must be 'tanimoto', 'dice' or 'cosine'")
    function = functions[metric][int(counts)]
    return float(function(first, second) if counts else function(first, second, fp_size=fp_size))


def molecular_similarity(query_smiles: str, target_smiles: str, metric: str = "tanimoto",
                          kind: str = "ecfp", fp_size: int = 2048, radius: int = 2,
                          path_length: int = 7, counts: bool = False) -> dict:
    """Compare native fingerprints; the coefficient is not an MCS or biological-accuracy proof.

    Example::

        assert molecular_similarity("CCO", "CCO")["score"] == 1.0
    """
    started = time.perf_counter_ns()
    values = _fingerprints([query_smiles, target_smiles], kind, fp_size, radius, path_length, counts)
    output = _operation("fingerprint_similarity")
    output.update(score=_similarity(*values, metric, counts, fp_size), metric=metric, kind=kind,
                  fp_size=fp_size, radius=radius, path_length=path_length, counts=counts)
    return _response(output, started)


def screen_molecules(query_smiles: str, targets_smiles: list[str], threshold: float = 0.7,
                     max_hits: int = 32, offset: int = 0, kind: str = "ecfp",
                     fp_size: int = 2048, radius: int = 2, path_length: int = 7,
                     counts: bool = False) -> dict:
    """Screen a bounded in-memory list by Tanimoto, with stable ranked ties and explicit pagination.

    Example::

        result = screen_molecules("CCO", ["CCO", "CCN", "CCO"], threshold=1.0)
        assert [hit["target_index"] for hit in result["hits"]] == [0, 2]
    """
    started = time.perf_counter_ns()
    if not isinstance(targets_smiles, list) or len(targets_smiles) > _LIMITS["targets"]:
        raise ValueError("targets_smiles must be a list of at most 64 molecular strings")
    texts = [_text(query_smiles, "query_smiles")] + [_text(value, "target_smiles") for value in targets_smiles]
    if sum(map(len, texts)) > _LIMITS["batch_text_characters"]:
        raise ValueError("MCP batch exceeds 65536 text characters")
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) or not 0 <= threshold <= 1 or not math.isfinite(threshold):
        raise ValueError("threshold must be a finite number between 0 and 1")
    _integer(max_hits, "max_hits", 1, _LIMITS["matches"])
    _integer(offset, "offset", 0, _LIMITS["targets"])
    values = _fingerprints(texts, kind, fp_size, radius, path_length, counts)
    hits = [{"target_index": index, "score": _similarity(values[0], value, "tanimoto", counts, fp_size)}
            for index, value in enumerate(values[1:])]
    hits = sorted((hit for hit in hits if hit["score"] >= threshold), key=lambda hit: -hit["score"])
    page = hits[offset:offset + max_hits]
    next_offset = offset + len(page) if offset + len(page) < len(hits) else None
    output = _operation("fingerprint_screening")
    output.update(hits=page, total_hits=len(hits), total_targets=len(targets_smiles), threshold=threshold,
                  max_hits=max_hits, offset=offset, next_offset=next_offset, truncated=len(page) < len(hits),
                  ordering="score_descending_then_input_index", metric="tanimoto", kind=kind,
                  fp_size=fp_size, radius=radius, path_length=path_length, counts=counts)
    return _response(output, started)


def standardise_molecule(smiles: str, remove_atom_maps: bool = False, cleanup: bool = False,
                         largest_fragment: bool = False, neutralise: bool = False,
                         canonical_tautomer: bool = False) -> dict:
    """Canonicalise representation; chemical transformations require explicit options and RDKit.

    The native default does not neutralise, remove fragments or select a tautomer.
    A canonical tautomer is a representation choice, not a stability prediction.

    Example::

        result = standardise_molecule("[Na+].[Cl-]")
        assert result["backend"] == "smsd" and not result["chemical_transformations"]
    """
    started = time.perf_counter_ns()
    flags = {name: _switch(value, name) for name, value in (
        ("remove_atom_maps", remove_atom_maps), ("cleanup", cleanup), ("largest_fragment", largest_fragment),
        ("neutralise", neutralise), ("canonical_tautomer", canonical_tautomer))}
    text = _text(smiles, "smiles")
    if remove_atom_maps:
        text = strip_atom_maps(text)
    graph = _molecule(text)
    if graph.n > _LIMITS["standardise_atoms"]:
        raise ValueError(f"MCP standardisation supports at most {_LIMITS['standardise_atoms']} atoms")
    chemical = any(flags[name] for name in ("cleanup", "largest_fragment", "neutralise", "canonical_tautomer"))
    output = _operation("representation_standardisation")
    operations = ["remove_atom_maps"] if remove_atom_maps else []
    if chemical:
        try:
            from rdkit import Chem, __version__ as backend_version
            from rdkit.Chem.MolStandardize import rdMolStandardize
        except ImportError as error:
            raise ImportError("Chemical standardisation requires the optional smsd[rdkit] dependency") from error
        molecule = Chem.MolFromSmiles(text)
        if molecule is None:
            raise ValueError("RDKit could not parse the molecular input")
        if cleanup:
            parameters = rdMolStandardize.CleanupParameters()
            parameters.maxRestarts = 32
            molecule = rdMolStandardize.Cleanup(molecule, parameters)
            operations.append("cleanup")
        if largest_fragment:
            molecule = rdMolStandardize.LargestFragmentChooser().choose(molecule)
            operations.append("largest_fragment")
        if neutralise:
            molecule = rdMolStandardize.Uncharger().uncharge(molecule)
            operations.append("neutralise")
        if canonical_tautomer:
            enumerator = rdMolStandardize.TautomerEnumerator()
            enumerator.SetMaxTautomers(128)
            enumerator.SetMaxTransforms(256)
            tautomers = enumerator.Enumerate(molecule)
            molecule = enumerator.PickCanonical(tautomers)
            complete = tautomers.status == rdMolStandardize.TautomerEnumeratorStatus.Completed
            output.update(status="COMPLETED" if complete else "LIMIT_REACHED", complete=complete,
                          tautomer_status=str(tautomers.status), tautomer_limits={"tautomers": 128, "transforms": 256})
            operations.append("canonical_tautomer" if complete else "canonical_from_partial_tautomers")
        normalised = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        backend = "rdkit"
    else:
        normalised, backend, backend_version = to_smiles(graph), "smsd", __version__
    operations.append("canonical_smiles")
    output.update(smiles=normalised, backend=backend, backend_version=backend_version,
                  requested=flags, operations=operations, chemical_transformations=chemical)
    return _response(output, started)


def _safe_svg(svg):
    if not isinstance(svg, str) or len(svg.encode("utf-8")) > _LIMITS["svg_bytes"]:
        raise ValueError("MCP SVG exceeds the byte limit")
    if re.search(r"<!\s*(DOCTYPE|ENTITY)", svg, re.IGNORECASE):
        raise ValueError("SVG document declarations are not supported")
    root = ET.fromstring(svg)
    allowed = {"svg", "g", "defs", "style", "path", "line", "polyline", "polygon", "circle",
               "ellipse", "rect", "text", "tspan", "title", "desc", "clipPath"}
    for element in root.iter():
        tag = element.tag.rsplit("}", 1)[-1]
        if (element.tag.startswith("{") and not element.tag.startswith("{http://www.w3.org/2000/svg}")) or tag not in allowed or (element is root and tag != "svg"):
            raise ValueError("Unsupported element in generated SVG")
        if tag == "style" and re.search(r"@import|url\s*\(|expression\s*\(", element.text or "", re.IGNORECASE):
            raise ValueError("External content in generated SVG")
        for name, value in element.attrib.items():
            local = name.rsplit("}", 1)[-1].lower()
            if (name.startswith("{") and name != "{http://www.w3.org/XML/1998/namespace}space") or local.startswith("on") or local in ("href", "src") or re.search(
                    r"javascript:|data:|https?:|file:|@import|expression\s*\(|url\s*\(", value, re.IGNORECASE):
                raise ValueError("External or executable content in generated SVG")
    return svg


def depict_mcs(query_smiles: str, target_smiles: str, objective: str = "mcis",
                connected: bool = True, timeout_ms: int = 1000, node_limit: int = 100000,
                match_stereo: bool = False, complete_rings_only: bool = False,
                ring_matches_ring_only: bool = False, width: int = 800, height: int = 360) -> dict:
    """Return a native mapping-highlight SVG and the unchanged search proof/status.

    No SVG is returned when no match was found; a stopped match stays unproved.
    Input bonds between highlighted atoms are drawn too; for MCES the JSON
    bond_mapping, rather than the drawing, defines the selected common edges.

    Example::

        result = depict_mcs("CCO", "CCN")
        assert result["found"] and "<svg" in result["svg"]
    """
    started = time.perf_counter_ns()
    _integer(width, "width", 200, 1600)
    _integer(height, "height", 150, 1200)
    result, options, query, target = _mcs(query_smiles, target_smiles, objective, connected, timeout_ms,
        node_limit, match_stereo, complete_rings_only, ring_matches_ring_only, _LIMITS["depiction_atoms"])
    output = _serialise(result, options)
    svg = None if result.best is None else _safe_svg(depict_pair(query, target, dict(result.best.atom_mapping),
                                                                width=width, height=height))
    output.update(svg=svg, mime_type="image/svg+xml", objective=objective, connected=connected,
                  width=width, height=height, highlight_scope="mapped_atoms_and_input_bonds_between_them")
    return _response(output, started, result.elapsed_ms)


def create_server():
    """Create the official SDK stdio server without starting it.

    Example::

        server = create_server()  # requires smsd[mcp]
    """
    from importlib import import_module
    from .tool_definitions import TOOL_REGISTRY, _mcp_callable, _validate_arguments
    try:
        module = import_module("mcp.server")
    except ModuleNotFoundError as error:
        if error.name != "mcp":
            raise
        raise ImportError("The MCP server requires the official SDK; install smsd[mcp]") from error
    server_type = getattr(module, "MCPServer", None)
    if server_type is None:
        from mcp.server.fastmcp import FastMCP as server_type
        from mcp.server.fastmcp.exceptions import ToolError
    else:
        from mcp.server.mcpserver.exceptions import ToolError

    class _MolecularServer(server_type):
        async def list_tools(self):
            definitions = await super().list_tools()
            return [definition.model_copy(update={
                "inputSchema" if "inputSchema" in type(definition).model_fields else "input_schema":
                    TOOL_REGISTRY[definition.name].input_schema}) for definition in definitions]

        async def call_tool(self, name, arguments, *args, **kwargs):
            # Check the original request before SDK argument normalisation.
            try:
                _validate_arguments(name, arguments)
            except ValueError as error:
                raise ToolError(str(error)) from error
            return await super().call_tool(name, arguments, *args, **kwargs)

    server = _MolecularServer("SMSD", instructions=(
        "Bounded molecular-text tools. Inspect native status, optimal and complete independently. "
        "An unfinished search proves neither a maximum nor absence. Fingerprint scores are not "
        "MCS or biological accuracy. Default standardisation changes representation only. "
        "Tools do not read arbitrary files or contact network/model services."))
    for tool in TOOL_REGISTRY.values():
        server.tool(description=tool.description)(_mcp_callable(tool))
    return server


def main():
    """Serve over standard input/output, with diagnostics on standard error.

    Example::

        # Run from a shell: python -m smsd.mcp_server
    """
    create_server().run(transport="stdio")


if __name__ == "__main__":
    main()
