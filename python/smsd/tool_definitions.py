# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Local function definitions and optional adapters for the bounded molecular tools.

Exporting definitions does not contact a model provider. The registry contains
the same callables used by the MCP server; optional frameworks are imported only
when their adapter is requested.
"""
import inspect
import math
from dataclasses import dataclass
from types import MappingProxyType
from typing import Annotated, Literal, get_args, get_origin

from . import mcp_server as _tools

__all__ = ["ToolDefinition", "TOOL_REGISTRY", "get_tool_definitions", "call_tool",
           "langchain_tools", "llamaindex_tools"]


def _parameter_schema(name, annotation):
    if annotation is str:
        schema = {"type": "string"}
    elif annotation is int:
        schema = {"type": "integer"}
    elif annotation is float:
        schema = {"type": "number"}
    elif annotation is bool:
        schema = {"type": "boolean"}
    elif get_origin(annotation) is list and get_args(annotation) == (str,):
        schema = {"type": "array", "items": {"type": "string", "minLength": 1,
                                                      "maxLength": _tools._LIMITS["text_characters"]},
                  "maxItems": _tools._LIMITS["targets"]}
    else:
        raise TypeError(f"Unsupported tool argument annotation: {name}")
    if name in ("query_smiles", "target_smiles", "query", "query_smarts", "smiles"):
        schema.update(minLength=1, maxLength=_tools._LIMITS["text_characters"])
    enums = {"objective": ("mcis", "mces"), "query_format": ("smiles", "smarts"),
             "kind": _tools._KINDS, "metric": _tools._METRICS}
    if name in enums:
        schema["enum"] = list(enums[name])
    ranges = {"timeout_ms": (0, _tools._LIMITS["timeout_ms"]),
              "node_limit": (1, _tools._LIMITS["node_limit"]),
              "max_matches": (1, _tools._LIMITS["matches"]),
              "max_hits": (1, _tools._LIMITS["matches"]),
              "offset": (0, _tools._LIMITS["targets"]), "threshold": (0, 1),
              "fp_size": (1, _tools._LIMITS["fingerprint_bits"]),
              "radius": (0, _tools._LIMITS["fingerprint_radius"]),
              "path_length": (1, _tools._LIMITS["path_length"]),
              "width": (200, 1600), "height": (150, 1200)}
    if name in ranges:
        schema.update(zip(("minimum", "maximum"), ranges[name]))
    return schema


@dataclass(frozen=True)
class ToolDefinition:
    """A registered molecular callable with a fresh JSON schema on each access.

    Example::

        definition = TOOL_REGISTRY["find_mcs"]
        assert "query_smiles" in definition.input_schema["required"]
    """
    function: object

    @property
    def name(self):
        """The public callable name.

        Example::

            assert TOOL_REGISTRY["find_mcs"].name == "find_mcs"
        """
        return self.function.__name__

    @property
    def description(self):
        """The callable's opening documentation paragraph.

        Example::

            assert TOOL_REGISTRY["find_mcs"].description
        """
        return inspect.getdoc(self.function).split("\n\n", 1)[0].replace("\n", " ")

    @property
    def input_schema(self):
        """Object schema with finite single-argument bounds and Python defaults.

        Example::

            assert TOOL_REGISTRY["find_mcs"].input_schema["additionalProperties"] is False
        """
        properties, required = {}, []
        for name, parameter in inspect.signature(self.function).parameters.items():
            schema = _parameter_schema(name, parameter.annotation)
            if parameter.default is inspect.Parameter.empty:
                required.append(name)
            else:
                schema["default"] = parameter.default
            properties[name] = schema
        return {"type": "object", "properties": properties, "required": required,
                "additionalProperties": False}


TOOL_REGISTRY = MappingProxyType({function.__name__: ToolDefinition(function) for function in (
    _tools.find_mcs, _tools.find_substructure, _tools.molecular_fingerprint,
    _tools.search_smarts, _tools.molecular_similarity, _tools.screen_molecules,
    _tools.standardise_molecule, _tools.depict_mcs)})


def get_tool_definitions(provider: str = "mcp") -> list[dict]:
    """Export definitions in MCP, input-schema, function or chat-function format.

    ``function`` uses flat definitions; ``chat-function`` nests each function.
    ``strict=False`` preserves optional Python defaults. The schemas do not use
    strict mode. No client, credential or network is used.

    Example::

        definitions = get_tool_definitions("function")
        assert definitions[0]["type"] == "function" and not definitions[0]["strict"]
    """
    if provider not in ("mcp", "input-schema", "function", "chat-function"):
        raise ValueError("provider must be 'mcp', 'input-schema', 'function' or 'chat-function'")
    definitions = []
    for tool in TOOL_REGISTRY.values():
        definition = {"name": tool.name, "description": tool.description}
        if provider in ("mcp", "input-schema"):
            definition["inputSchema" if provider == "mcp" else "input_schema"] = tool.input_schema
        else:
            definition.update(parameters=tool.input_schema, strict=False)
            definition = ({"type": "function", **definition} if provider == "function"
                          else {"type": "function", "function": definition})
        definitions.append(definition)
    return definitions


def _validate(value, schema, name):
    kind = schema["type"]
    valid = {"string": lambda: isinstance(value, str),
             "boolean": lambda: isinstance(value, bool),
             "integer": lambda: isinstance(value, int) and not isinstance(value, bool),
             "number": lambda: isinstance(value, (float, int)) and not isinstance(value, bool)
                                      and (isinstance(value, int) or math.isfinite(value)),
             "array": lambda: isinstance(value, list)}[kind]()
    if not valid:
        raise ValueError(f"{name} must have JSON type {kind}")
    if "enum" in schema and value not in schema["enum"]:
        raise ValueError(f"{name} is not an allowed value")
    if "minimum" in schema and not schema["minimum"] <= value <= schema["maximum"]:
        raise ValueError(f"{name} is outside its supported range")
    if kind == "string" and not schema.get("minLength", 0) <= len(value) <= schema.get("maxLength", len(value)):
        raise ValueError(f"{name} is outside its supported text length")
    if kind == "array":
        if len(value) > schema["maxItems"]:
            raise ValueError(f"{name} contains too many items")
        for index, item in enumerate(value):
            _validate(item, schema["items"], f"{name}[{index}]")


def call_tool(name: str, arguments: dict) -> dict:
    """Invoke a registered tool with checked JSON arguments and ordinary defaults.

    This is a local dispatcher, not a model request. Unknown tools, extra fields
    and JSON type coercions are rejected before molecular work begins.

    Example::

        result = call_tool("molecular_similarity", {"query_smiles": "CCO", "target_smiles": "CCO"})
        assert result["score"] == 1.0
    """
    tool = _validate_arguments(name, arguments)
    return tool.function(**arguments)


def _validate_arguments(name, arguments):
    if not isinstance(name, str) or name not in TOOL_REGISTRY:
        raise ValueError("Unknown molecular tool")
    if not isinstance(arguments, dict):
        raise ValueError("arguments must be a JSON object")
    tool = TOOL_REGISTRY[name]
    schema = tool.input_schema
    if any(key not in schema["properties"] for key in arguments):
        raise ValueError("Unknown tool argument")
    if any(key not in arguments for key in schema["required"]):
        raise ValueError("Missing required tool argument")
    for key, value in arguments.items():
        _validate(value, schema["properties"][key], key)
    return tool


def _argument_annotation(parameter, schema):
    from pydantic import Field
    annotation = parameter.annotation
    if "enum" in schema:
        return Literal[tuple(schema["enum"])]
    if schema["type"] == "array":
        annotation = list[Annotated[str, Field(strict=True, min_length=schema["items"]["minLength"],
                                               max_length=schema["items"]["maxLength"])]]
    constraints = {target: schema[source] for source, target in (
        ("minimum", "ge"), ("maximum", "le"), ("minLength", "min_length"),
        ("maxLength", "max_length"), ("maxItems", "max_length")) if source in schema}
    return Annotated[annotation, Field(strict=True, **constraints)]


def _argument_model(tool):
    from pydantic import ConfigDict, create_model
    fields = {}
    for name, parameter in inspect.signature(tool.function).parameters.items():
        annotation = _argument_annotation(parameter, tool.input_schema["properties"][name])
        default = ... if parameter.default is inspect.Parameter.empty else parameter.default
        fields[name] = (annotation, default)
    return create_model(tool.name + "Arguments", __config__=ConfigDict(extra="forbid", strict=True), **fields)


def _mcp_callable(tool):
    from functools import wraps

    @wraps(tool.function)
    def invoke(**arguments):
        return call_tool(tool.name, arguments)

    signature = inspect.signature(tool.function)
    parameters = [parameter.replace(annotation=_argument_annotation(parameter,
        tool.input_schema["properties"][name])) for name, parameter in signature.parameters.items()]
    invoke.__signature__ = signature.replace(parameters=parameters)
    invoke.__annotations__ = {parameter.name: parameter.annotation for parameter in parameters}
    invoke.__annotations__["return"] = dict
    return invoke


def langchain_tools() -> list:
    """Make LangChain StructuredTools with lazy optional imports and the same bounds.

    Requires ``langchain-core`` with Pydantic 2 support. The caller chooses and
    operates any model; these tools do not create a client.

    Example::

        tools = langchain_tools()  # requires langchain-core
        result = tools[0].invoke({"query_smiles": "CCO", "target_smiles": "CCN"})
        assert result["found"]
    """
    try:
        from langchain_core.tools import StructuredTool
    except ImportError as error:
        raise ImportError("The LangChain adapter requires langchain-core") from error
    return [StructuredTool.from_function(func=tool.function, name=tool.name,
                description=tool.description, args_schema=_argument_model(tool))
            for tool in TOOL_REGISTRY.values()]


def llamaindex_tools() -> list:
    """Make LlamaIndex FunctionTools with lazy optional imports and the same bounds.

    Requires ``llama-index-core`` with Pydantic 2 support. No index, model or
    external service is created by this adapter.

    Example::

        tools = llamaindex_tools()  # requires llama-index-core
        result = tools[0].call(query_smiles="CCO", target_smiles="CCN").raw_output
        assert result["found"]
    """
    try:
        from llama_index.core.tools import FunctionTool
    except ImportError as error:
        raise ImportError("The LlamaIndex adapter requires llama-index-core") from error
    return [FunctionTool.from_defaults(fn=tool.function, name=tool.name,
                description=tool.description, fn_schema=_argument_model(tool))
            for tool in TOOL_REGISTRY.values()]
