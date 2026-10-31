# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Generated tool definitions and optional adapters keep one molecular contract."""
import inspect
import json

import pytest

from smsd import mcp_server
from smsd.tool_definitions import (TOOL_REGISTRY, call_tool, get_tool_definitions,
                                    langchain_tools, llamaindex_tools)


def test_registry_and_all_exported_schemas_keep_function_arguments_and_defaults():
    assert set(TOOL_REGISTRY) == {"find_mcs", "find_substructure", "molecular_fingerprint", "search_smarts",
                                 "molecular_similarity", "screen_molecules", "standardise_molecule", "depict_mcs"}
    for name, tool in TOOL_REGISTRY.items():
        assert tool.function is getattr(mcp_server, name)
        schema = tool.input_schema
        assert not schema["additionalProperties"]
        assert set(schema["properties"]) == set(inspect.signature(tool.function).parameters)
        for argument, parameter in inspect.signature(tool.function).parameters.items():
            if parameter.default is inspect.Parameter.empty:
                assert argument in schema["required"]
            else:
                assert argument not in schema["required"]
                assert schema["properties"][argument]["default"] == parameter.default
    for provider in ("mcp", "input-schema", "function", "chat-function"):
        definitions = get_tool_definitions(provider)
        assert len(definitions) == 8 and json.loads(json.dumps(definitions)) == definitions
        for definition, tool in zip(definitions, TOOL_REGISTRY.values()):
            if provider == "mcp":
                assert definition["inputSchema"] == tool.input_schema
            elif provider == "input-schema":
                assert definition["input_schema"] == tool.input_schema
            else:
                function = definition if provider == "function" else definition["function"]
                assert function["parameters"] == tool.input_schema and function["strict"] is False
                assert definition["type"] == "function"


def test_exported_schema_mutation_does_not_change_the_registry():
    exported = get_tool_definitions("input-schema")
    exported[0]["input_schema"]["properties"]["node_limit"]["maximum"] = 10**30
    assert TOOL_REGISTRY["find_mcs"].input_schema["properties"]["node_limit"]["maximum"] == 1000000
    with pytest.raises(TypeError):
        TOOL_REGISTRY["other"] = TOOL_REGISTRY["find_mcs"]


def test_local_dispatch_runs_native_chemistry_with_default_options():
    result = call_tool("find_mcs", {"query_smiles": "CCO", "target_smiles": "CCN"})
    assert result["found"] and result["matches"][0]["atom_count"] == 2
    assert result["budget"] == {"timeout_ms": 1000, "node_limit": 100000}
    result = call_tool("screen_molecules", {"query_smiles": "CCO", "targets_smiles": ["CCO", "CCN", "CCO"],
                                             "threshold": 1, "max_hits": 1})
    assert result["hits"] == [{"target_index": 0, "score": 1.0}] and result["next_offset"] == 1


@pytest.mark.parametrize("name,arguments", [
    ("__import__", {}), ("find_mcs", {}),
    ("find_mcs", {"query_smiles": "CC", "target_smiles": "CC", "unknown": 1}),
    ("find_mcs", {"query_smiles": "CC", "target_smiles": "CC", "node_limit": True}),
    ("find_mcs", {"query_smiles": "CC", "target_smiles": "CC", "timeout_ms": "1000"}),
    ("find_mcs", {"query_smiles": "CC", "target_smiles": "CC", "timeout_ms": 30001}),
    ("find_mcs", {"query_smiles": "CC", "target_smiles": "CC", "objective": "other"}),
    ("molecular_fingerprint", {"smiles": "C" * 16385}),
    ("screen_molecules", {"query_smiles": "CC", "targets_smiles": [1]}),
    ("screen_molecules", {"query_smiles": "CC", "targets_smiles": [], "threshold": float("nan")}),
])
def test_dispatch_rejects_unknown_coerced_or_unbounded_json(name, arguments):
    with pytest.raises(ValueError):
        call_tool(name, arguments)


def test_refusals_name_the_argument_and_list_what_is_allowed():
    pair = {"query_smiles": "CCO", "target_smiles": "CCN"}
    arguments = ", ".join(TOOL_REGISTRY["find_mcs"].input_schema["properties"])
    cases = (
        ("find_mcs", {**pair, "objective": "MCS"},
         "objective is not an allowed value; allowed values: 'mcis', 'mces'"),
        ("molecular_similarity", {**pair, "metric": "jaccard"},
         "metric is not an allowed value; allowed values: 'tanimoto', 'dice', 'cosine'"),
        ("find_mcs", {**pair, "timeout": 500},
         f"Unknown tool argument 'timeout' for find_mcs; allowed arguments: {arguments}"),
        ("find_mcs", {**pair, "timeout": 500, "nodes": 1},
         f"Unknown tool arguments 'timeout', 'nodes' for find_mcs; allowed arguments: {arguments}"),
        ("find_mcs", {"query_smiles": "CCO"},
         "Missing required tool argument 'target_smiles' for find_mcs; "
         "required arguments: query_smiles, target_smiles"),
        ("find_mcs", {},
         "Missing required tool arguments 'query_smiles', 'target_smiles' for find_mcs; "
         "required arguments: query_smiles, target_smiles"),
        ("find_mcs", {**pair, "timeout_ms": 30001},
         "timeout_ms is outside its supported range; allowed range: 0 to 30000"),
        ("find_mc", pair, "Unknown molecular tool; available tools: " + ", ".join(TOOL_REGISTRY)),
    )
    for name, call, message in cases:
        with pytest.raises(ValueError) as refusal:
            call_tool(name, call)
        assert str(refusal.value) == message


def test_refusals_echo_a_bounded_part_of_unknown_argument_names():
    pair = {"query_smiles": "CCO", "target_smiles": "CCN"}
    with pytest.raises(ValueError) as refusal:
        call_tool("find_mcs", {**pair, "x" * 1000: 1})
    assert str(refusal.value).startswith("Unknown tool argument '" + "x" * 63 + " for find_mcs;")
    with pytest.raises(ValueError) as refusal:
        call_tool("find_mcs", {**pair, **{f"extra{index}": index for index in range(20)}})
    assert str(refusal.value).startswith(
        "Unknown tool arguments " + ", ".join(repr(f"extra{index}") for index in range(8)) + ", ... for find_mcs;")


@pytest.mark.parametrize("name", sorted(TOOL_REGISTRY))
def test_each_fixed_choice_and_unknown_argument_refusal_lists_the_choices(name):
    schema = TOOL_REGISTRY[name].input_schema
    required = {key: ["C"] if schema["properties"][key]["type"] == "array" else "C" for key in schema["required"]}
    for argument, choices in ((key, value["enum"]) for key, value in schema["properties"].items() if "enum" in value):
        with pytest.raises(ValueError) as refusal:
            call_tool(name, {**required, argument: "other"})
        assert str(refusal.value) == (f"{argument} is not an allowed value; allowed values: "
                                      + ", ".join(repr(choice) for choice in choices))
    with pytest.raises(ValueError) as refusal:
        call_tool(name, {**required, "other": 1})
    assert str(refusal.value).endswith("allowed arguments: " + ", ".join(schema["properties"]))
    assert all(key in str(refusal.value) for key in schema["required"])


def test_unknown_export_provider_is_rejected():
    with pytest.raises(ValueError):
        get_tool_definitions("other")


def test_langchain_tools_invoke_real_search_and_enforce_the_same_budget():
    pytest.importorskip("langchain_core")
    from pydantic import ValidationError
    tools = {tool.name: tool for tool in langchain_tools()}
    output = tools["molecular_similarity"].invoke({"query_smiles": "CCO", "target_smiles": "CCO"})
    assert output["score"] == 1.0
    stopped = tools["find_mcs"].invoke({"query_smiles": "CCO", "target_smiles": "CCN", "timeout_ms": 0})
    assert stopped["status"] == "TIMED_OUT" and not stopped["optimal"]
    with pytest.raises(ValidationError):
        tools["find_mcs"].invoke({"query_smiles": "CC", "target_smiles": "CC", "node_limit": True})


def test_llamaindex_tools_call_real_search_without_a_model():
    pytest.importorskip("llama_index.core")
    tools = {tool.metadata.name: tool for tool in llamaindex_tools()}
    result = tools["find_mcs"].call(query_smiles="CCO", target_smiles="CCN").raw_output
    assert result["found"] and result["matches"][0]["atom_count"] == 2
    result = tools["screen_molecules"].call(query_smiles="CCO", targets_smiles=["CCO"], threshold=1).raw_output
    assert result["hits"] == [{"target_index": 0, "score": 1.0}]
