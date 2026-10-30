# Molecular tools and MCP

SMSD provides eight local molecular tools through the official Model Context
Protocol Python SDK. The same registry exports function definitions for model
clients and optional LangChain and LlamaIndex adapters. Exporting a definition
or calling a tool does not make a model request, need an API key or contact a
network service.

Install SMSD 8.1.0 with the optional MCP SDK, then start the server:

```sh
python -m pip install 'smsd[mcp]==8.1.0'
smsd-mcp
# Equivalent: python -m smsd.mcp_server
```

From a source checkout, install with `python -m pip install '.[mcp]'`. Releases
before 8.1.0 do not include these tools.

Configure an MCP client to start `smsd-mcp` with no arguments. The server has
no HTTP listener. Molecular arguments are SMILES or SMARTS text, never paths,
URLs, Python expressions or commands. Client configuration, logging and any
model connection belong to the client.

See [Connecting AI assistants](AI_ASSISTANTS.md) for client set-up.

| Tool | Operation |
| --- | --- |
| `find_mcs` | One common graph: atom-maximising MCIS or bond-maximising MCES; connected or fragmented. |
| `find_substructure` | One molecular or SMARTS embedding, with optional induced matching. |
| `search_smarts` | Bounded SMARTS embedding enumeration. |
| `molecular_fingerprint` | Set-bit positions or sparse position/count pairs. |
| `molecular_similarity` | Tanimoto, Dice or cosine on native fingerprints. |
| `screen_molecules` | Tanimoto screening of a small in-memory target list, with stable ties and explicit pages. |
| `standardise_molecule` | Native representation canonicalisation; explicitly requested optional chemical transformations. |
| `depict_mcs` | Native mapping-highlight SVG together with the MCS result and proof flags. |

## Search results and limits

Every result is JSON-compatible and includes the engine version, limits and
wall elapsed time. Search results additionally retain the native `status`,
`optimal`, `complete`, `score`, `nodes_visited`, `elapsed_ms`, chemistry options,
budgets and mappings. `proof` describes these engine flags; it is not an
independent chemical validation certificate. Atom indices are zero-based in
the parsed inputs. `wall_elapsed_ms` ends before the JSON size/serialisation
check and SDK/transport handling; `wall_scope` records that boundary.

Inspect proof and stopping status separately from `found`. An unfinished match
need not be maximum; an unfinished `found: false` does not prove absence.
`complete` also describes whether requested enumeration finished. A result
limit can leave enumeration incomplete even when an optimum is known.
Fingerprint and representation operations mark search proof as inapplicable.
A fingerprint coefficient is not an MCS proof or a biological-accuracy score.

Search defaults are 1,000 ms and 100,000 nodes per call. The permitted maxima
are 30,000 ms and 1,000,000 nodes; zero time requests an immediate stop. A zero
node cap is rejected rather than selecting an unlimited search. These budgets
cover native search, after parsing and SMARTS compilation. Wall time also
includes preparation and any depiction, but excludes JSON serialisation and
SDK/transport handling; it is not a hard wall
deadline for the whole Python call. Use process/container limits for that.

Molecular text is limited to 16,384 characters and parsed molecules to 1,024
atoms. SMARTS has the text limit; its compiled predicate size is not a molecule
atom count. Enumeration and returned screening pages have at most 32 matches.
Fingerprint limits are 65,536 positions, radius 8 and path length 8, with at
most 2,000,000 estimated path walks. Search chemistry switches expose stereo,
complete-ring and ring-to-ring restrictions; unsupported options are not
silently inferred from natural language.

Screening accepts at most 64 targets, 65,536 combined text characters and 4,096
combined query/target atoms. It scores all supplied targets, sorts by descending
score and then input index, and returns `total_hits`, `offset`, `next_offset`
and `truncated`. A page is incomplete when it does not contain the whole hit
list; `next_offset: null` means there is no later page. This small-list tool
does not promise a GPU path or an exact structure search.

```python
from smsd.tool_definitions import call_tool

result = call_tool("find_mcs", {
    "query_smiles": "CCO", "target_smiles": "CCN", "timeout_ms": 1000,
})
assert result["found"] and result["matches"][0]["atom_count"] == 2

hits = call_tool("screen_molecules", {
    "query_smiles": "CCO", "targets_smiles": ["CCO", "CCN", "CCO"],
    "threshold": 1.0, "max_hits": 1,
})
assert hits["hits"][0]["target_index"] == 0 and hits["next_offset"] == 1
```

## Standardisation and depiction

The default `standardise_molecule` parses and writes native canonical SMILES.
It preserves the represented charge, fragments, isotope, atom maps and stereo;
it does not neutralise, remove salts or choose a tautomer. `remove_atom_maps`
is a separate explicit representation option. Whatever options are set,
`standardise_molecule` accepts at most 128 atoms, and each result's
`limits` reports this as `standardise_atoms`.

The optional `smsd[rdkit]` dependency is required only when `cleanup`,
`largest_fragment`, `neutralise` or `canonical_tautomer` is explicitly true.
The result records the backend/version and operation sequence; there is no
automatic fallback. These operations change chemical graphs and can change
charge, fragments or stereo. Canonical tautomer selection is a representation
choice, not a prediction of stability. Tautomer enumeration is limited to 128
tautomers and 256 transformations. If that enumeration stops, the result reports
`LIMIT_REACHED`, `complete: false` and a canonical choice from the partial set.

`depict_mcs` accepts at most 256 atoms per input and returns at most 1 MiB of
SVG, with dimensions bounded to 200–1,600 by 150–1,200. Its native SVG is
checked for executable elements, external references and document entities.
The drawing highlights mapped atoms and input bonds between them. For MCES,
the returned `bond_mapping` is the selected-edge witness; the drawing can also
show other input bonds between highlighted atoms. No match gives `svg: null`.
A stopped search remains stopped when depicted. Each result is limited to
2 MiB as JSON text, measured before the MCP SDK encodes it.

## Function definitions and framework adapters

```python
from smsd.tool_definitions import get_tool_definitions, langchain_tools, llamaindex_tools

input_schema_tools = get_tool_definitions("input-schema")
function_tools = get_tool_definitions("function")
chat_function_tools = get_tool_definitions("chat-function")
mcp_tools = get_tool_definitions("mcp")
```

Every definition carries the tool's `name` and `description`. `mcp` gives the
argument schema as `inputSchema`, and `input-schema` gives it as `input_schema`.
`function` gives flat definitions with `type: "function"` and `parameters`;
`chat-function` nests `name`, `description` and `parameters` in a `function`
object beside `type: "function"`. Both function formats put `strict: false`
beside `parameters`, because omitted optional arguments keep their Python
defaults; the schemas do not use strict mode. All exports reject
additional fields and include finite single-argument bounds; combined work/atom
limits are also checked by the functions. Both `call_tool` and the MCP server
check the original JSON object before SDK normalisation, rejecting unknown
fields and type coercions. The server advertises the registry's
`additionalProperties: false` schema.

Install `langchain-core` or `llama-index-core` separately when needed; the
adapters require their Pydantic 2 interfaces. They load lazily and create only
local tools, not clients, models or indexes:

```python
from smsd.tool_definitions import langchain_tools, llamaindex_tools

tools = langchain_tools()  # optional langchain-core dependency
assert tools[0].invoke({"query_smiles": "CCO", "target_smiles": "CCN"})["found"]

tools = llamaindex_tools()  # optional llama-index-core dependency
assert tools[0].call(query_smiles="CCO", target_smiles="CCN").raw_output["found"]
```

See the [container recipe](../packaging/mcp/README.md) and the
[Hugging Face demonstration source](../examples/huggingface/README.md). Neither
deploys a service. The official protocol reference is the
[MCP Python SDK](https://github.com/modelcontextprotocol/python-sdk).
