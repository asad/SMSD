# Chemistry tools for AI assistants

SMSD Pro 8.2.0 offers eight molecular tools that an AI assistant can call:
maximum common substructure (MCS), substructure and SMARTS matching,
fingerprints, similarity, small-list screening, standardisation and MCS
depiction. The tools run on your own machine. They are served to an assistant
by a local Model Context Protocol (MCP) server, `smsd-mcp`, and the same tools
are available from Python as a registry with JSON schemas, so that your own
code can offer them to a model or call them with no model at all. This chapter
shows which route suits which task, what each tool takes and returns, how to
start the server and connect a client, how to read a tool result, and the
limits that keep each call bounded.

## Choose a route

| Task | Use | Needs | Model involved |
| --- | --- | --- | --- |
| Give a desktop or terminal assistant the tools | `smsd-mcp`, started by the assistant's MCP client | The `mcp` extra | Yours, through the client |
| Offer the tools to a model API from your own program | `get_tool_definitions(provider)` and `call_tool(name, arguments)` | The core package | Yours, through your code |
| Run the tools from a script, a test or a workflow | `call_tool(name, arguments)` | The core package | None |
| Search a library, enumerate every match or set detailed chemistry | The Python API: [Substructure search](SUBSTRUCTURE.md), [MCS](MCS.md), [Screening](SCREENING.md) | The core package | None |

The tools are deliberately narrow: each call takes SMILES or SMARTS text,
runs one bounded operation and returns a JSON-compatible dictionary. They suit
questions about a few molecules at a time. For a library of thousands of
molecules, a saved index or a GPU, use the Python API directly; the tools
expose a subset of its options and cap every input.

## The tools

| Tool | Required inputs | Main optional inputs | Returns |
| --- | --- | --- | --- |
| `find_mcs` | `query_smiles`, `target_smiles` | `objective` (`mcis` or `mces`), `connected`, `timeout_ms`, `node_limit`, chemistry switches | One common graph as `matches`, with `status`, `optimal`, `complete`, `score` and the budget used |
| `find_substructure` | `query`, `target_smiles` | `query_format` (`smarts` or `smiles`), `induced`, `timeout_ms`, `node_limit`, chemistry switches | One embedding of the query in the target, with the same proof fields |
| `search_smarts` | `query_smarts`, `target_smiles` | `max_matches` (1 to 32), `timeout_ms`, `node_limit`, chemistry switches | Every embedding up to `max_matches`, and `truncated` when the list stopped at the cap |
| `molecular_fingerprint` | `smiles` | `kind`, `fp_size`, `radius`, `path_length`, `counts` | `values`: set-bit positions, or position and count pairs |
| `molecular_similarity` | `query_smiles`, `target_smiles` | `metric` (`tanimoto`, `dice` or `cosine`) and the fingerprint inputs | `score`, the coefficient of the two fingerprints |
| `screen_molecules` | `query_smiles`, `targets_smiles` (a list) | `threshold`, `max_hits`, `offset` and the fingerprint inputs | One page of `hits`, each a `target_index` and Tanimoto `score`, with `total_hits` and `next_offset` |
| `standardise_molecule` | `smiles` | `remove_atom_maps`; `cleanup`, `largest_fragment`, `neutralise` and `canonical_tautomer`, which need RDKit | Canonical `smiles`, the `backend` used and the `operations` applied |
| `depict_mcs` | `query_smiles`, `target_smiles` | The `find_mcs` inputs, `width`, `height` | The `find_mcs` result plus `svg`, a drawing of the pair with the common atoms highlighted |

The chemistry switches are `match_stereo`, `complete_rings_only` and
`ring_matches_ring_only`, all off by default. Every other chemistry option
keeps its default, described in
[Getting started](GETTING_STARTED.md#default-chemistry), and each search
result lists the options in force under `chemistry`. The fingerprint `kind` is
one of `ecfp`, `fcfp`, `circular`, `path`, `torsion` or `mcs`; the
[Options](#options) section gives every default.

## Start the server

Install SMSD Pro with the `mcp` extra, which adds the MCP Python SDK. The
server and the tool definitions need SMSD 8.1.0 or later and Python 3.12 or
later; the refusal messages that list allowed values need 8.2.0:

```bash
python -m pip install "smsd[mcp]==8.2.0"
smsd-mcp
# Equivalent: python -m smsd.mcp_server
```

`smsd-mcp` takes no arguments. It talks to its client over standard input and
output, writes diagnostics to standard error and opens no network port, so it
is normally started by the client rather than from a shell; run by hand, it
waits silently for protocol messages. Without the `mcp` extra it stops with
`ImportError: The MCP server requires the official SDK; install smsd[mcp]`.

To check an installation without an assistant, run
[`mcp_stdio_client.py`](../../python/examples/mcp_stdio_client.py) from the
examples. It starts the server with the same interpreter, lists the eight
tools and runs one `find_mcs` call through the protocol.

## Connect an MCP client

An MCP client needs two things: the command that starts the server and its
arguments, which here are none. Many clients take a JSON entry of this form:

```json
{
  "mcpServers": {
    "smsd": {
      "command": "/absolute/path/to/smsd-mcp",
      "args": []
    }
  }
}
```

Give the absolute path that `which smsd-mcp` (on Windows,
`where.exe smsd-mcp`) prints, because a desktop client may not search your
shell's PATH; the path also fixes which Python environment, and so which SMSD
version, the server uses. The server needs no API key and no environment
settings. Once connected, the client lists the eight tools with their
schemas, and the server also sends the client a short note on reading
results: inspect `status`, `optimal` and `complete`, and do not read an
unfinished search as a maximum or as proof of absence.

The configuration file, its location and the command that adds a server
differ from client to client.
[Connecting AI assistants](../AI_ASSISTANTS.md#connect-the-mcp-server) gives
the exact steps for several desktop and terminal assistants and for an open
model run on your own machine. Every result carries `engine_version`; use it
to confirm the SMSD version, because some clients show the SDK's version for
the server instead.

## Call the tools from Python

`smsd.tool_definitions.call_tool(name, arguments)` runs a tool locally, with
the arguments as a dictionary of JSON values, and returns the same dictionary
that the server sends to a client. It needs neither the `mcp` extra nor a
model. The program below runs an MCS, two single substructure searches and a
SMARTS enumeration:

```python
from smsd.tool_definitions import call_tool

# MCS of N-methylacetamide and propanamide: atom-maximising, connected.
mcs = call_tool("find_mcs", {"query_smiles": "CC(=O)NC",
                             "target_smiles": "CCC(=O)N"})
print("find_mcs:", mcs["found"], mcs["status"], mcs["optimal"], mcs["complete"])
print("  atoms:", mcs["matches"][0]["atom_count"],
      "mapping:", mcs["matches"][0]["atom_mapping"])

# One embedding of a SMARTS pattern: a hydroxyl on an aromatic carbon.
hit = call_tool("find_substructure", {"query": "[OX2H]c",
                                      "target_smiles": "CC(C)(C)c1ccc(O)cc1"})
print("find_substructure:", hit["found"], hit["status"],
      hit["matches"][0]["atom_mapping"])

# The same call with a molecule as the query.
phenol = call_tool("find_substructure", {"query": "Oc1ccccc1",
                                         "target_smiles": "CC(C)(C)c1ccc(O)cc1",
                                         "query_format": "smiles"})
print("molecule query:", phenol["found"], len(phenol["matches"][0]["atom_mapping"]))

# Every embedding, up to max_matches: the three hydroxyls of glycerol.
every = call_tool("search_smarts", {"query_smarts": "[OX2H]",
                                    "target_smiles": "OCC(O)CO"})
print("search_smarts:", [m["atom_mapping"] for m in every["matches"]],
      "complete:", every["complete"], "truncated:", every["truncated"])

assert mcs["optimal"] and mcs["matches"][0]["atom_mapping"] == [[0, 1], [1, 2], [2, 3], [3, 4]]
assert hit["matches"][0]["atom_mapping"] == [[0, 8], [1, 7]]
assert phenol["found"] and len(phenol["matches"][0]["atom_mapping"]) == 7
assert len(every["matches"]) == 3 and every["complete"] and not every["truncated"]
```

The program prints:

```text
find_mcs: True COMPLETED True True
  atoms: 4 mapping: [[0, 1], [1, 2], [2, 3], [3, 4]]
find_substructure: True COMPLETED [[0, 8], [1, 7]]
molecule query: True 7
search_smarts: [[[0, 0]], [[0, 3]], [[0, 5]]] complete: True truncated: False
```

Mappings are lists of `[query_atom, target_atom]` pairs, sorted by query atom,
with zero-based atom numbers in the order of the input SMILES, as described in
[Atom numbering](GETTING_STARTED.md#atom-numbering). `find_substructure`
expects SMARTS unless `query_format` is `"smiles"`. Each entry of `matches`
also holds `bond_mapping`, `fragments`, `atom_count` and `bond_count`.

The tool functions can also be imported from `smsd.mcp_server` and called with
keyword arguments; importing that module neither loads the SDK nor starts a
server. Use `call_tool` for arguments that come from a model or a file: it
checks the whole argument object against the schema first, as the server does.

## Read a tool result

A search result keeps the native fields of a `SearchResult`: `found`,
`status`, `optimal`, `complete`, `score`, `nodes_visited` and `elapsed_ms`,
with the same meanings as in
[Reading a result](GETTING_STARTED.md#reading-a-result). It adds:

- `engine_version`, the SMSD version that ran the call;
- `proof`, which repeats `optimal` and `complete` with the scope
  `native_search`; it describes the search, and is not an independent
  chemical validation;
- `budget`, the `timeout_ms` and `node_limit` in force, and `chemistry`, the
  chemistry options;
- `limits`, the fixed input and output limits of the server;
- `timings`, with `wall_elapsed_ms` from tool entry to just before the result
  is serialised, and `native_search_ms`.

The rule for an assistant is the rule for any search: a returned match is
always a valid match, but `found: false` proves absence only when `status` is
`COMPLETED` and `optimal` is true, and an MCS is proved maximum only when
`optimal` is true. The program below shows a search stopped by its time limit,
an enumeration stopped by its result cap, and the messages with which
`call_tool` refuses bad arguments:

```python
from smsd.tool_definitions import call_tool

anthracene = "c1ccc2cc3ccccc3cc2c1"

# A time limit of 0 stops the search at once: found=False proves nothing.
stopped = call_tool("find_substructure", {"query": "c1ccccc1",
                                          "target_smiles": anthracene,
                                          "timeout_ms": 0})
print("stopped:", stopped["found"], stopped["status"],
      stopped["optimal"], stopped["complete"])

# A result cap: the matches returned are real, but the list is not the whole.
capped = call_tool("search_smarts", {"query_smarts": "[OX2H]",
                                     "target_smiles": "OCC(O)CO",
                                     "max_matches": 2})
print("capped:", len(capped["matches"]), capped["status"],
      "complete:", capped["complete"], "truncated:", capped["truncated"])


def proves_absence(result):
    return (not result["found"] and result["status"] == "COMPLETED"
            and result["optimal"])


absent = call_tool("find_substructure", {"query": "[N]", "target_smiles": anthracene})
print("absence proved:", proves_absence(absent), "| after a stop:",
      proves_absence(stopped))

# Bad arguments are refused before any molecular work; the message says what is accepted.
for arguments in ({"query_smiles": "CCO"},
                  {"query_smiles": "CCO", "target_smiles": "CCN", "timeout_ms": 60000},
                  {"query_smiles": "CCO", "target_smiles": "CCN", "timeout_ms": "500"},
                  {"query_smiles": "CCO", "target_smiles": "CCN", "path": "a.sdf"},
                  {"query_smiles": "C1CC", "target_smiles": "CCN"}):
    try:
        call_tool("find_mcs", arguments)
    except ValueError as error:
        print("refused:", error)

assert stopped["status"] == "TIMED_OUT" and not stopped["complete"]
assert capped["status"] == "LIMIT_REACHED" and capped["truncated"]
assert proves_absence(absent) and not proves_absence(stopped)
```

The program prints:

```text
stopped: False TIMED_OUT False False
capped: 2 LIMIT_REACHED complete: False truncated: True
absence proved: True | after a stop: False
refused: Missing required tool argument 'target_smiles' for find_mcs; required arguments: query_smiles, target_smiles
refused: timeout_ms is outside its supported range; allowed range: 0 to 30000
refused: timeout_ms must have JSON type integer
refused: Unknown tool argument 'path' for find_mcs; allowed arguments: query_smiles, target_smiles, objective, connected, timeout_ms, node_limit, match_stereo, complete_rings_only, ring_matches_ring_only
refused: Unclosed ring(s): 1
```

`search_smarts` reports a capped list as `LIMIT_REACHED` with
`complete: false` and `truncated: true`; the matches it does return are valid,
and `optimal` stays true because a match was found. Arguments are refused
before any molecular work: unknown tools, unknown or missing arguments, JSON
type changes such as a string for an integer, and values out of range each
raise `ValueError`. The message names the argument and, where a list helps,
what is accepted: the available tools, the tool's allowed or required
arguments, the allowed values or the allowed range. A SMILES or SMARTS parse
error also raises `ValueError`. Return that message to the model as the tool's
answer, so that it can correct the call.

Fingerprint, similarity, screening and standardisation results are not
searches. They report `status: "COMPLETED"`, except a canonical tautomer
chosen from a stopped enumeration, described below, and a `proof` whose scope
names the operation and whose `search_proof_applicable` is false. A similarity
score says how alike two fingerprints are; it is not an MCS size or a measure
of biological activity.

## Fingerprints, similarity and screening

`molecular_fingerprint` returns the fingerprint as positions, which keeps the
result small and readable. `molecular_similarity` compares two molecules with
one coefficient, and `screen_molecules` scores a short list against one query
and returns the hits at or above `threshold` one page at a time:

```python
from smsd.tool_definitions import call_tool

# ECFP4 set bits in a 64-position fingerprint, then position/count pairs.
bits = call_tool("molecular_fingerprint", {"smiles": "OCCO", "fp_size": 64})
counts = call_tool("molecular_fingerprint", {"smiles": "OCCO", "fp_size": 64,
                                             "counts": True})
print("bits:", bits["values"])
print("counts:", counts["values"])
print("proof:", bits["proof"])

# Phenol against 2-methylphenol with each coefficient.
pair = {"query_smiles": "Oc1ccccc1", "target_smiles": "Cc1ccccc1O"}
for metric in ("tanimoto", "dice", "cosine"):
    score = call_tool("molecular_similarity", {**pair, "metric": metric})["score"]
    print(f"{metric}: {score:.3f}")

# A small list, two hits per page; next_offset leads to the next page.
library = ["Oc1ccccc1", "Cc1ccccc1O", "Nc1ccccc1", "CCO", "Oc1ccccc1", "c1ccccc1"]
request = {"query_smiles": "Oc1ccccc1", "targets_smiles": library,
           "threshold": 0.3, "max_hits": 2}
offset, pages = 0, []
while offset is not None:
    page = call_tool("screen_molecules", {**request, "offset": offset})
    pages.append([(hit["target_index"], round(hit["score"], 3)) for hit in page["hits"]])
    offset = page["next_offset"]
print("total hits:", page["total_hits"], "of", page["total_targets"])
for number, hits in enumerate(pages, 1):
    print(f"page {number}:", hits)

assert bits["values"] == [4, 5, 6, 27, 40, 62]
assert counts["values"] == [[4, 2], [5, 2], [6, 2], [27, 2], [40, 2], [62, 2]]
assert bits["proof"]["search_proof_applicable"] is False
assert pages == [[(0, 1.0), (4, 1.0)], [(1, 0.346)]]
```

The program prints:

```text
bits: [4, 5, 6, 27, 40, 62]
counts: [[4, 2], [5, 2], [6, 2], [27, 2], [40, 2], [62, 2]]
proof: {'scope': 'fingerprint_generation', 'optimal': None, 'search_proof_applicable': False}
tanimoto: 0.346
dice: 0.514
cosine: 0.525
total hits: 3 of 6
page 1: [(0, 1.0), (4, 1.0)]
page 2: [(1, 0.346)]
```

`screen_molecules` scores every target, keeps those whose Tanimoto score is at
least `threshold`, sorts them by descending score and then by input index, and
returns at most `max_hits` from `offset`. `next_offset` is the offset of the
next page, or `None` (`null` in JSON) on the last page, and `truncated` is
true when the page does not hold the whole hit list. The screen always uses
Tanimoto; `metric` belongs to `molecular_similarity` only. It runs on the CPU
and is not an exact structure search. For a larger library, use the screens
in [Screening a library](SCREENING.md), which are not capped in this way.

## Standardisation and depiction

`standardise_molecule` by default parses the molecule and writes canonical
SMILES. It changes the representation only: charges, fragments, isotopes,
stereochemistry and atom maps are kept, and `remove_atom_maps` is a separate
switch. `depict_mcs` runs the same search as `find_mcs` and adds an SVG
drawing of the pair:

```python
from smsd.tool_definitions import call_tool

# The default changes the representation only: charges and fragments stay.
salt = call_tool("standardise_molecule", {"smiles": "[Na+].[O-]C(=O)c1ccccc1"})
print("smiles:", salt["smiles"], "| backend:", salt["backend"])
print("operations:", salt["operations"],
      "| chemical transformations:", salt["chemical_transformations"])

mapped = call_tool("standardise_molecule", {"smiles": "[CH3:1][OH:2]",
                                            "remove_atom_maps": True})
print("atom maps removed:", mapped["smiles"], mapped["operations"])

# The MCS of phenol and aniline, with the search result and a drawing.
drawing = call_tool("depict_mcs", {"query_smiles": "Oc1ccccc1",
                                   "target_smiles": "Nc1ccccc1",
                                   "width": 600, "height": 300})
print("depict_mcs:", drawing["status"], drawing["optimal"],
      drawing["matches"][0]["atom_count"], "atoms")
print("svg:", "<svg" in drawing["svg"], drawing["mime_type"],
      drawing["width"], "x", drawing["height"])

# No common graph, no drawing.
none = call_tool("depict_mcs", {"query_smiles": "O", "target_smiles": "N"})
print("no match:", none["found"], none["svg"])

assert not salt["chemical_transformations"] and "." in salt["smiles"]
assert mapped["operations"] == ["remove_atom_maps", "canonical_smiles"]
assert drawing["optimal"] and drawing["matches"][0]["atom_count"] == 6
assert none["svg"] is None
```

The program prints:

```text
smiles: C(c1ccccc1)([O-])=O.[Na+] | backend: smsd
operations: ['canonical_smiles'] | chemical transformations: False
atom maps removed: CO ['remove_atom_maps', 'canonical_smiles']
depict_mcs: COMPLETED True 6 atoms
svg: True image/svg+xml 600 x 300
no match: False None
```

`cleanup`, `largest_fragment`, `neutralise` and `canonical_tautomer` change
the chemical graph, so each must be asked for explicitly, and they need the
`rdkit` extra; without it the call raises `ImportError`, and there is no
silent fallback. With any of them set, `backend` is `rdkit`, `operations`
lists each step in order and `chemical_transformations` is true. A canonical
tautomer is a choice of representation, not a prediction of the most stable
form. Tautomer enumeration stops at 128 tautomers or 256 transformations; when
it stops, the result reports `LIMIT_REACHED` and `complete: false`, and the
choice is made from the tautomers found.

The drawing highlights the mapped atoms and each bond between them that the
mapping carries onto a bond of the other molecule. For an MCES search, the
selected common bonds are those in `bond_mapping`; the drawing can also
highlight such bonds outside it. A search with no
match returns `svg: None`, and a stopped search keeps its stopped status when
drawn. The SVG is checked before it is returned and contains no scripts,
external references or document declarations.

## Tool definitions for a model API

`get_tool_definitions(provider)` returns the eight tools as a list of
dictionaries in one of four shapes. Pass the list to a model API as its tools,
then run each tool call the model makes with `call_tool` and send the result
back. No client is created and nothing is sent anywhere:

```python
import json

from smsd.tool_definitions import TOOL_REGISTRY, get_tool_definitions

for provider in ("mcp", "input-schema", "function", "chat-function"):
    first = get_tool_definitions(provider)[0]
    print(f"{provider}:", sorted(first))

schema = TOOL_REGISTRY["find_mcs"].input_schema
print("required:", schema["required"])
print("timeout_ms:", json.dumps(schema["properties"]["timeout_ms"]))
print("objective:", json.dumps(schema["properties"]["objective"]))
print("additionalProperties:", schema["additionalProperties"])
print("tools:", len(TOOL_REGISTRY), "| first:", next(iter(TOOL_REGISTRY)))

assert len(get_tool_definitions()) == 8
assert get_tool_definitions("chat-function")[0]["function"]["strict"] is False
```

The program prints:

```text
mcp: ['description', 'inputSchema', 'name']
input-schema: ['description', 'input_schema', 'name']
function: ['description', 'name', 'parameters', 'strict', 'type']
chat-function: ['function', 'type']
required: ['query_smiles', 'target_smiles']
timeout_ms: {"type": "integer", "minimum": 0, "maximum": 30000, "default": 1000}
objective: {"type": "string", "enum": ["mcis", "mces"], "default": "mcis"}
additionalProperties: False
tools: 8 | first: find_mcs
```

| `provider` | Shape | Typical use |
| --- | --- | --- |
| `mcp` (default) | `name`, `description`, `inputSchema` | MCP tool lists |
| `input-schema` | `name`, `description`, `input_schema` | APIs that take an input schema per tool |
| `function` | flat: `type: "function"`, `name`, `description`, `parameters`, `strict` | APIs that take flat function definitions |
| `chat-function` | nested: `type: "function"`, `function: {name, description, parameters, strict}` | Chat-completion APIs and most local model runtimes |

Every schema sets `additionalProperties: false`, marks the required
arguments, gives each optional argument its default, and bounds every number,
text length and list length. The function shapes set `strict: false` so that
a model may leave out optional arguments and the defaults apply. A client that
drops schema keywords cannot widen the limits, because `call_tool` and the
server check every argument themselves. The schemas hold single-argument
bounds only; combined limits, such as the atoms of a whole screen, are checked
when the tool runs.

Which shape each vendor's API expects, the changes some need, and complete
loops for several hosted and local models are in
[Use the tool definitions from Python](../AI_ASSISTANTS.md#use-the-tool-definitions-from-python).
Optional LangChain and LlamaIndex adapters, `langchain_tools()` and
`llamaindex_tools()`, wrap the same callables with the same bounds; see
[Function definitions and framework adapters](../MCP_SERVER.md#function-definitions-and-framework-adapters).

## A tool loop with fixed replies

A tool loop sends the conversation and the tool definitions to a model, runs
each tool call the model asks for, appends the results and asks again, until
the model answers without a tool call. The loop below has that structure, but
a small class returns prepared replies in place of a model, so the program
runs offline and gives the same output every time. This is also a convenient
way to test your own loop before connecting a model. The replies use the
`chat-function` shape: the first asks for one valid call and one with an
unknown argument, the second corrects the call after reading the refusal, and
the third answers:

```python
import json

from smsd.tool_definitions import call_tool, get_tool_definitions


class FixedReplies:
    """Stands in for a model: returns the next prepared reply on each request."""

    def __init__(self, replies):
        self.replies, self.requests = list(replies), []

    def create(self, messages, tools):
        self.requests.append({"messages": list(messages), "tools": tools})
        return self.replies.pop(0)


def tool_call(call_id, name, arguments):
    return {"id": call_id, "type": "function",
            "function": {"name": name, "arguments": json.dumps(arguments)}}


def run_tool(name, arguments):
    """Run one tool locally; a refusal goes back to the model as text."""
    try:
        return json.dumps(call_tool(name, json.loads(arguments)))
    except (ValueError, ImportError) as error:
        return json.dumps({"error": str(error)})


pair = {"query_smiles": "c1ccccc1O", "target_smiles": "c1ccccc1N"}
model = FixedReplies([
    # Turn 1: one valid call, and one with an argument the tool does not take.
    {"role": "assistant", "content": None, "tool_calls": [
        tool_call("call_1", "find_mcs", pair),
        tool_call("call_2", "find_mcs", {**pair, "timeout": 500})]},
    # Turn 2: the corrected call, after reading the refusal.
    {"role": "assistant", "content": None, "tool_calls": [
        tool_call("call_3", "find_mcs", {**pair, "timeout_ms": 500})]},
    # Turn 3: the answer, which reports the proof flags.
    {"role": "assistant", "content": "The MCS has 6 atoms; COMPLETED, "
     "optimal and complete, so it is proved maximum.", "tool_calls": None},
])

tools = get_tool_definitions("chat-function")
messages = [{"role": "user", "content": "Find the MCS of phenol and aniline."}]
while True:
    reply = model.create(messages=messages, tools=tools)
    messages.append(reply)
    if not reply["tool_calls"]:
        break
    for call in reply["tool_calls"]:
        text = run_tool(call["function"]["name"], call["function"]["arguments"])
        messages.append({"role": "tool", "tool_call_id": call["id"], "content": text})
        result = json.loads(text)
        print(call["id"], "->", result.get("error") or
              f"status={result['status']} optimal={result['optimal']} "
              f"atoms={result['matches'][0]['atom_count']}")

print("model requests:", len(model.requests))
print("answer:", messages[-1]["content"])

results = [json.loads(m["content"]) for m in messages if m["role"] == "tool"]
assert results[0]["optimal"] and results[0]["matches"][0]["atom_count"] == 6
assert results[1]["error"].startswith("Unknown tool argument 'timeout' for find_mcs")
assert results[2]["budget"]["timeout_ms"] == 500
assert model.requests[0]["tools"][0]["function"]["name"] == "find_mcs"
```

The program prints:

```text
call_1 -> status=COMPLETED optimal=True atoms=6
call_2 -> Unknown tool argument 'timeout' for find_mcs; allowed arguments: query_smiles, target_smiles, objective, connected, timeout_ms, node_limit, match_stereo, complete_rings_only, ring_matches_ring_only
call_3 -> status=COMPLETED optimal=True atoms=6
model requests: 3
answer: The MCS has 6 atoms; COMPLETED, optimal and complete, so it is proved maximum.
```

Two points carry over to a loop with a real model. Send every result back in
full, including `status`, `optimal` and `complete`, so that the model can say
whether its answer is proved; and catch `ValueError` and `ImportError` from
`call_tool` and return the message as the tool result, rather than ending the
loop, so that the model can correct its call. To replace the fixed replies
with a model, replace `model.create` with the client call of your model API;
[Connecting AI assistants](../AI_ASSISTANTS.md#use-the-tool-definitions-from-python)
shows that change for each API. For a fixed sequence of calls with no model at
all, see [`agent_workflow.py`](../../python/examples/agent_workflow.py), which
screens a list, runs an MCS on each hit and draws the best one.

## Limits and safety

The limits are fixed in the server and apply equally to `call_tool`. Each
result reports them under `limits`.

| Limit | Value | Applies to |
| --- | --- | --- |
| Text argument | 16,384 characters | Every SMILES or SMARTS argument |
| Atoms per molecule | 1,024 | Searches and fingerprints |
| Atoms per molecule, depiction | 256 | `depict_mcs` |
| Atoms per molecule, standardisation | 128 | `standardise_molecule` |
| Search time, `timeout_ms` | Default 1,000; at most 30,000; 0 stops at once | `find_mcs`, `find_substructure`, `search_smarts`, `depict_mcs` |
| Search nodes, `node_limit` | Default 100,000; from 1 to 1,000,000 | The same four tools |
| Matches or hits returned | At most 32 | `search_smarts`, a `screen_molecules` page |
| Screen targets | At most 64 | `screen_molecules` |
| Combined text | 65,536 characters | The query and targets of `screen_molecules` |
| Combined atoms | 4,096 | `molecular_similarity`, `screen_molecules` |
| Fingerprint | `fp_size` up to 65,536; `radius` up to 8; `path_length` up to 8; at most 2,000,000 estimated path walks | The fingerprint tools |
| Drawing | 200 to 1,600 by 150 to 1,200 pixels; SVG at most 1 MiB | `depict_mcs` |
| Result | 2 MiB as JSON text | Every tool |

Points to keep in mind when deploying the tools:

- Local only. The tools read no files, fetch no URLs and contact no network
  or model service. Molecular arguments are SMILES or SMARTS text, never
  paths, URLs, expressions or commands. The server has no network listener.
- Hosted assistants. A hosted assistant sends your prompts, the tool
  arguments and the results, including any structures, to its provider. Use
  a model run on your own machine when structures must not leave it.
- Time budgets cover the search, not the whole call. They start after the
  input is parsed and the SMARTS compiled; preparation, depiction and
  serialisation come on top, and standardisation has no time budget. For a
  hard deadline, apply process or container limits. The
  [container recipe](../../packaging/mcp/README.md) runs the server with no
  network access, a read-only file system and memory and CPU caps.
- No zero node limit. A `node_limit` of 0 is refused rather than read as
  unlimited.
- Natural-language requests are not mapped to options. Only the switches in
  the schema exist; an assistant cannot ask for an option the tools lack.
- Approval prompts. Keep your client's tool approval prompts on unless you
  trust every tool it can call.

## Options

### Search tools

`find_mcs`, `depict_mcs`, `find_substructure` and `search_smarts` share these
arguments where they apply:

| Argument | Type | Default | Meaning |
| --- | --- | --- | --- |
| `objective` | string | `"mcis"` | `mcis` maximises common atoms; `mces` maximises common bonds (`find_mcs`, `depict_mcs`) |
| `connected` | boolean | `true` | A single connected common graph, or fragments when false (`find_mcs`, `depict_mcs`) |
| `query_format` | string | `"smarts"` | How `query` is read: `smarts` or `smiles` (`find_substructure`) |
| `induced` | boolean | `false` | Induced matching: unbonded query atoms must map to unbonded target atoms (`find_substructure`) |
| `max_matches` | integer | `16` | Embeddings to return, 1 to 32 (`search_smarts`) |
| `timeout_ms` | integer | `1000` | Search time limit, 0 to 30,000 ms |
| `node_limit` | integer | `100000` | Search node limit, 1 to 1,000,000 |
| `match_stereo` | boolean | `false` | Match tetrahedral and double-bond stereochemistry |
| `complete_rings_only` | boolean | `false` | Matched ring atoms and bonds lie in fully matched rings |
| `ring_matches_ring_only` | boolean | `false` | Ring atoms and bonds match only ring atoms and bonds |
| `width`, `height` | integer | `800`, `360` | Drawing size, 200 to 1,600 and 150 to 1,200 (`depict_mcs`) |

### Fingerprint tools

`molecular_fingerprint`, `molecular_similarity` and `screen_molecules` share
the fingerprint arguments:

| Argument | Type | Default | Meaning |
| --- | --- | --- | --- |
| `kind` | string | `"ecfp"` | `ecfp` (also `circular`), `fcfp`, `path`, `torsion` or `mcs` |
| `fp_size` | integer | `2048` | Fingerprint width, 1 to 65,536 |
| `radius` | integer | `2` | Radius of `ecfp`, `circular` and `fcfp`, 0 to 8 |
| `path_length` | integer | `7` | Longest path, in bonds, of `path` and `mcs`, 1 to 8 |
| `counts` | boolean | `false` | Position and count pairs instead of set bits; not available for `mcs` |
| `metric` | string | `"tanimoto"` | `tanimoto`, `dice` or `cosine` (`molecular_similarity`) |
| `threshold` | number | `0.7` | Lowest Tanimoto score kept, 0 to 1, inclusive (`screen_molecules`) |
| `max_hits` | integer | `32` | Hits per page, 1 to 32 (`screen_molecules`) |
| `offset` | integer | `0` | First hit of the page, 0 to 64 (`screen_molecules`) |

[Fingerprints and similarity](FINGERPRINTS.md) describes the fingerprint
families and coefficients.

### standardise_molecule

| Argument | Type | Default | Meaning |
| --- | --- | --- | --- |
| `remove_atom_maps` | boolean | `false` | Remove atom-map numbers before writing |
| `cleanup` | boolean | `false` | RDKit clean-up of the structure |
| `largest_fragment` | boolean | `false` | Keep the largest fragment |
| `neutralise` | boolean | `false` | Neutralise charges where possible |
| `canonical_tautomer` | boolean | `false` | Write the canonical tautomer |

The last four need the `rdkit` extra.

## Notes

- Each call is independent. The tools keep no state between calls, so a
  query is prepared again on every call; for many searches with one query,
  use the Python API.
- `screen_molecules` fingerprints the query and every target on each call,
  and so does `molecular_similarity` for its pair.
- `elapsed_ms` and the `timings` vary from run to run; do not compare them
  across machines or ask an assistant to.
- `find_mcs` returns one common graph. For every optimal mapping, constraints
  or the common graph of a series, use the Python API in
  [Maximum common substructure](MCS.md).
- A capped enumeration or a capped screen page is incomplete, not wrong; ask
  for the next page or a higher cap, up to the limits above.
- An empty or whitespace-only argument is refused, and so is a `threshold`
  that is not a finite number from 0 to 1.

## See also

- [Molecular tools and MCP](../MCP_SERVER.md): the reference for the tools,
  their results and limits.
- [Connecting AI assistants](../AI_ASSISTANTS.md): set-up for specific
  assistants, model APIs and local model runtimes, with complete tool loops.
- [Command line and containers](COMMAND_LINE.md#molecular-tools-server): the
  server among the other command-line tools.
- [API reference, molecular tools](../API_REFERENCE.md#molecular-tools): the
  signature of every tool function.
- [Getting started](GETTING_STARTED.md#reading-a-result): what a search
  result proves.
- The [MCP stdio client](../../python/examples/mcp_stdio_client.py) and
  [agent workflow](../../python/examples/agent_workflow.py) examples.
