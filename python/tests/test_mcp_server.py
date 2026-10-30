# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""MCP tools run native chemistry and preserve completion and proof status."""

import asyncio
import json
import os
import sys
import xml.etree.ElementTree as ET

import pytest
import smsd
from smsd import mcp_server


def test_mcs_tool_returns_the_native_mapping_and_proof():
    expected = smsd.MCS.find("CCO").match("CCN")
    output = mcp_server.find_mcs("CCO", "CCN")
    assert output["engine_version"] == smsd.__version__
    assert output["status"] == expected.status.name
    assert output["optimal"] == expected.optimal
    assert output["complete"] == expected.complete
    assert output["matches"][0]["atom_mapping"] == [list(pair) for pair in sorted(expected.best.atom_mapping.items())]
    assert output["matches"][0]["bond_mapping"] == [list(bond) for bond in expected.best.bond_mapping]
    assert json.loads(json.dumps(output)) == output


@pytest.mark.parametrize("search,arguments", [
    (mcp_server.find_mcs, ("CCO", "CCCO")),
    (mcp_server.find_substructure, ("CO", "CCCO")),
])
def test_zero_time_budget_stays_unfinished(search, arguments):
    output = search(*arguments, timeout_ms=0, node_limit=1)
    assert output["status"] == "TIMED_OUT"
    assert not output["optimal"] and not output["complete"]
    assert output["budget"] == {"timeout_ms": 0, "node_limit": 1}


def test_node_budget_matches_the_native_stopping_status():
    options = smsd.MCS.Options()
    options.timeout_ms = 1000
    options.node_limit = 1
    expected = smsd.MCS.find("CCOc1ccc(CC)cc1").with_options(options).match("CCNc1ccc(OC)cc1")
    output = mcp_server.find_mcs("CCOc1ccc(CC)cc1", "CCNc1ccc(OC)cc1", node_limit=1)
    assert output["status"] == expected.status.name == "LIMIT_REACHED"
    assert output["optimal"] == expected.optimal
    assert output["complete"] == expected.complete
    assert output["nodes_visited"] == expected.nodes_visited


def test_substructure_tool_distinguishes_smarts_and_smiles():
    assert mcp_server.find_substructure("[OX2H]", "CCO")["found"]
    output = mcp_server.find_substructure("CO", "CCN", query_format="smiles")
    assert not output["found"] and output["complete"]
    assert output["status"] == "COMPLETED"
    assert output["query_format"] == "smiles"


def test_fingerprint_tool_uses_native_bits_and_counts():
    binary = mcp_server.molecular_fingerprint("CCCCCC", fp_size=64)
    assert binary["values"] == smsd.fingerprint("CCCCCC", kind="ecfp", fp_size=64)
    counted = mcp_server.molecular_fingerprint("CCCCCC", fp_size=64, counts=True)
    options = smsd.FingerprintOptions(family="ecfp", size=64, counts=True)
    assert counted["values"] == [list(pair) for pair in smsd.fingerprint("CCCCCC", options=options)]
    assert any(count > 1 for _, count in counted["values"])


@pytest.mark.parametrize("options", [{"timeout_ms": -1}, {"node_limit": 0},
                                    {"node_limit": True}, {"objective": "other"}])
def test_mcs_tool_rejects_unbounded_or_invalid_options(options):
    with pytest.raises(ValueError):
        mcp_server.find_mcs("CC", "CCC", **options)


@pytest.mark.parametrize("options", [{"fp_size": 1073741824, "counts": True},
                                    {"radius": 10000}, {"path_length": 10000}])
def test_fingerprint_tool_rejects_excessive_allocation_or_work_options(options):
    with pytest.raises(ValueError, match="MCP fingerprint limits"):
        mcp_server.molecular_fingerprint("CCO", **options)


def test_fingerprint_tool_rejects_excessive_input_and_atom_inventory():
    with pytest.raises(ValueError, match="text limit"):
        mcp_server.molecular_fingerprint("C" * 16385)
    with pytest.raises(ValueError, match="1024 atoms"):
        mcp_server.molecular_fingerprint("C" * 1025)


def test_fingerprint_tool_bounds_path_expansion():
    labels = [[] for _ in range(8)]
    ring = 1
    for first in range(8):
        for second in range(first + 2, 8):
            label = str(ring) if ring < 10 else "%" + str(ring)
            labels[first].append(label)
            labels[second].append(label)
            ring += 1
    densely_connected = "".join("C" + "".join(atom_labels) for atom_labels in labels)
    with pytest.raises(ValueError, match="path work bound"):
        mcp_server.molecular_fingerprint(densely_connected, kind="path", path_length=8)


def test_official_sdk_serves_real_tools_over_stdio():
    pytest.importorskip("mcp")
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    script = ("import smsd; "
              f"smsd.__path__.insert(0, {os.path.dirname(mcp_server.__file__)!r}); "
              "from smsd.mcp_server import main; main()")

    async def exercise():
        parameters = StdioServerParameters(command=sys.executable, args=["-c", script], env=dict(os.environ))
        async with stdio_client(parameters) as (reader, writer):
            async with ClientSession(reader, writer) as client:
                await client.initialize()
                tools = await client.list_tools()
                assert {tool.name for tool in tools.tools} == {
                    "find_mcs", "find_substructure", "molecular_fingerprint", "search_smarts",
                    "molecular_similarity", "screen_molecules", "standardise_molecule", "depict_mcs"}
                for tool in tools.tools:
                    schema = getattr(tool, "inputSchema", None) or tool.input_schema
                    assert schema["additionalProperties"] is False
                response = await client.call_tool("find_mcs", {"query_smiles": "CCO",
                                                                "target_smiles": "CCN"})
                output = json.loads(response.content[0].text)
                assert output["found"] and output["optimal"]
                assert output["matches"][0]["atom_count"] == 2
                stopped = await client.call_tool("find_substructure", {
                    "query": "[OX2H]", "target_smiles": "CCO", "timeout_ms": 0})
                output = json.loads(stopped.content[0].text)
                assert output["status"] == "TIMED_OUT" and not output["complete"]
                fingerprint = await client.call_tool("molecular_fingerprint", {"smiles": "CCO", "fp_size": 64})
                output = json.loads(fingerprint.content[0].text)
                assert output["values"] == smsd.fingerprint("CCO", kind="ecfp", fp_size=64)
                matches = await client.call_tool("search_smarts", {
                    "query_smarts": "[#8]", "target_smiles": "OCCO"})
                output = json.loads(matches.content[0].text)
                assert {tuple(match["atom_mapping"][0]) for match in output["matches"]} == {(0, 0), (0, 3)}
                invalid = await client.call_tool("find_mcs", {
                    "query_smiles": "CC", "target_smiles": "CC", "node_limit": True})
                assert getattr(invalid, "isError", getattr(invalid, "is_error", False))
                misspelt = await client.call_tool("find_mcs", {
                    "query_smiles": "CCO", "target_smiles": "CCN", "timeout_mss": 0})
                assert getattr(misspelt, "isError", getattr(misspelt, "is_error", False))
                assert "Unknown tool argument" in misspelt.content[0].text

    asyncio.run(exercise())


def test_smarts_enumeration_keeps_mappings_and_enumeration_limit():
    complete = mcp_server.search_smarts("[#8]", "OCCO")
    assert {tuple(match["atom_mapping"][0]) for match in complete["matches"]} == {(0, 0), (0, 3)}
    assert complete["complete"] and not complete["truncated"]
    partial = mcp_server.search_smarts("[#8]", "OCCO", max_matches=1)
    assert len(partial["matches"]) == 1 and partial["truncated"]
    assert not partial["complete"]
    stopped = mcp_server.search_smarts("[#8]", "OCCO", timeout_ms=0)
    assert stopped["status"] == "TIMED_OUT" and not stopped["complete"]


@pytest.mark.parametrize("metric", ["tanimoto", "dice", "cosine"])
@pytest.mark.parametrize("counts", [False, True])
def test_similarity_identical_and_disjoint_element_fingerprints(metric, counts):
    identical = mcp_server.molecular_similarity("CCO", "CCO", metric=metric, counts=counts)
    assert identical["score"] == pytest.approx(1.0)
    assert identical["proof"]["optimal"] is None
    assert not identical["proof"]["search_proof_applicable"]
    disjoint = mcp_server.molecular_similarity("C", "[He]", metric=metric, fp_size=65536, counts=counts)
    assert disjoint["score"] == 0.0


def test_screening_keeps_input_ties_and_explicit_pages():
    first = mcp_server.screen_molecules("CCO", ["CCO", "CCN", "CCO"], threshold=1, max_hits=1)
    assert first["hits"] == [{"target_index": 0, "score": 1.0}]
    assert first["total_hits"] == 2 and first["total_targets"] == 3
    assert first["next_offset"] == 1 and first["truncated"]
    second = mcp_server.screen_molecules("CCO", ["CCO", "CCN", "CCO"], threshold=1, max_hits=1, offset=1)
    assert second["hits"] == [{"target_index": 2, "score": 1.0}]
    assert second["next_offset"] is None
    empty = mcp_server.screen_molecules("CCO", [])
    assert empty["hits"] == [] and empty["total_hits"] == 0 and not empty["truncated"]


def test_count_screening_uses_the_count_tanimoto_coefficient():
    targets = ["CCCC", "CC", "CCCCCC", "CCCC"]
    options = smsd.FingerprintOptions(family="ecfp", size=128, counts=True)
    query = smsd.fingerprint("CCCC", options=options)
    expected = sorted(enumerate(smsd.count_tanimoto_coefficient(query, smsd.fingerprint(target, options=options))
                                for target in targets), key=lambda pair: -pair[1])
    output = mcp_server.screen_molecules("CCCC", targets, threshold=0, counts=True, fp_size=128)
    assert [(hit["target_index"], hit["score"]) for hit in output["hits"]] == expected


@pytest.mark.parametrize("threshold", [-0.1, 1.1, float("nan"), float("inf"), True, "0.7"])
def test_screening_rejects_invalid_thresholds(threshold):
    with pytest.raises(ValueError):
        mcp_server.screen_molecules("CC", ["CC"], threshold=threshold)


def test_batch_limits_apply_to_all_targets_before_pagination():
    with pytest.raises(ValueError, match="at most 64"):
        mcp_server.screen_molecules("C", ["C"] * 65)
    with pytest.raises(ValueError, match="4096 atoms"):
        mcp_server.screen_molecules("C", ["C" * 1024] * 4, max_hits=1)
    with pytest.raises(ValueError, match="65536 text"):
        mcp_server.screen_molecules("C", ["[CH4:123456789]" * 1024] * 5, max_hits=1)


def test_native_standardisation_preserves_chemistry_and_fragments():
    smiles = "[13CH3:7][C@H](O)F.[Na+]"
    output = mcp_server.standardise_molecule(smiles)
    graph = smsd.parse_smiles(output["smiles"])
    assert graph.n == 5 and sorted(graph.formal_charge) == [0, 0, 0, 0, 1]
    assert sorted(graph.mass_number) == [0, 0, 0, 0, 13]
    assert 7 in graph.atom_class and len(smsd.assign_rs(graph)) == 1
    assert smsd.count_components(graph) == 2
    assert sorted(smsd.assign_rs(graph).values()) == sorted(smsd.assign_rs(smsd.parse_smiles(smiles)).values())
    assert output["backend"] == "smsd" and not output["chemical_transformations"]
    assert output["operations"] == ["canonical_smiles"]
    without_maps = mcp_server.standardise_molecule(smiles, remove_atom_maps=True)
    assert not any(smsd.parse_smiles(without_maps["smiles"]).atom_class)


@pytest.mark.parametrize("options", [{}, {"remove_atom_maps": True}, {"neutralise": True},
                                     {"canonical_tautomer": True}])
def test_standardisation_applies_one_atom_limit_on_every_route(options):
    limit = mcp_server._LIMITS["standardise_atoms"]
    assert limit == 128
    with pytest.raises(ValueError, match=f"at most {limit} atoms"):
        mcp_server.standardise_molecule("C" * (limit + 1), **options)
    if options.keys() & {"neutralise", "canonical_tautomer"}:
        pytest.importorskip("rdkit")
    output = mcp_server.standardise_molecule("C" * limit, **options)
    assert smsd.parse_smiles(output["smiles"]).n == limit
    assert output["limits"]["standardise_atoms"] == limit


def test_optional_chemical_standardisation_is_requested_and_attributed():
    pytest.importorskip("rdkit")
    output = mcp_server.standardise_molecule("CC(=O)[O-].[Na+]", largest_fragment=True, neutralise=True)
    graph = smsd.parse_smiles(output["smiles"])
    assert graph.n == 4 and all(charge == 0 for charge in graph.formal_charge)
    assert output["backend"] == "rdkit" and output["backend_version"]
    assert output["chemical_transformations"]
    assert output["operations"] == ["largest_fragment", "neutralise", "canonical_smiles"]
    tautomer = mcp_server.standardise_molecule("CC(=O)C", canonical_tautomer=True)
    assert tautomer["tautomer_status"] and tautomer["tautomer_limits"] == {"tautomers": 128, "transforms": 256}
    assert tautomer["complete"] == (tautomer["status"] == "COMPLETED")


def test_mcs_depiction_returns_safe_svg_and_original_search_proof():
    output = mcp_server.depict_mcs("CCO", "CCN")
    root = ET.fromstring(output["svg"])
    assert root.tag == "{http://www.w3.org/2000/svg}svg"
    assert output["matches"][0]["atom_count"] == 2
    assert output["mime_type"] == "image/svg+xml" and output["optimal"]
    assert output["timings"]["wall_elapsed_ms"] >= 0
    assert output["timings"]["native_search_ms"] == output["elapsed_ms"]
    stopped = mcp_server.depict_mcs("CCO", "CCN", timeout_ms=0)
    assert stopped["svg"] is None and stopped["status"] == "TIMED_OUT" and not stopped["optimal"]


@pytest.mark.parametrize("svg", [
    '<svg xmlns="http://www.w3.org/2000/svg"><script>alert(1)</script></svg>',
    '<svg xmlns="http://www.w3.org/2000/svg" onload="alert(1)"/>',
    '<svg xmlns="http://www.w3.org/2000/svg"><image href="https://example.org/a"/></svg>',
    '<svg xmlns="http://www.w3.org/2000/svg"><style>@import "evil.css";</style></svg>',
    '<!DOCTYPE svg [<!ENTITY x "hello">]><svg>&x;</svg>',
    '<svg xmlns="http://www.w3.org/2000/svg"><text xmlns="http://www.w3.org/1999/xhtml">x</text></svg>',
])
def test_generated_svg_guard_rejects_executable_or_external_content(svg):
    with pytest.raises(ValueError):
        mcp_server._safe_svg(svg)


@pytest.mark.parametrize("function,arguments,options", [
    (mcp_server.find_mcs, ("CC", "CC"), {"timeout_ms": 30001}),
    (mcp_server.find_substructure, ("C", "CC"), {"node_limit": 1000001}),
    (mcp_server.search_smarts, ("C", "CC"), {"max_matches": 33}),
    (mcp_server.depict_mcs, ("CC", "CC"), {"width": 1601}),
    (mcp_server.standardise_molecule, ("CC",), {"neutralise": "false"}),
])
def test_tool_specific_limits_reject_unbounded_or_coerced_arguments(function, arguments, options):
    with pytest.raises(ValueError):
        function(*arguments, **options)
