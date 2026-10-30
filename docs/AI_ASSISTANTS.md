# Connecting AI assistants

SMSD lets an AI assistant run molecular comparisons on your machine: maximum
common substructure, substructure and SMARTS matching, fingerprints,
similarity, screening, standardisation and MCS depiction. Each search result
reports its `status` and whether it is `optimal` and `complete`, so the
assistant can tell a proven answer from a search that stopped at its limit.

You can use a hosted assistant or an open model that runs on your own machine.
Connect SMSD through its local MCP server, `smsd-mcp`, or pass its tool
definitions to a model API from your own Python code. SMSD itself sends nothing
over the network. A hosted assistant sends your prompts, tool arguments and
results, including any structures, to its model provider.

The eight tools and their results are described in
[Molecular tools and MCP](MCP_SERVER.md).

## Install

```sh
python -m pip install "smsd[mcp]>=8.1.0"
which smsd-mcp    # Windows: where.exe smsd-mcp
```

The MCP server and the tool definitions need SMSD 8.1.0 or later and Python
3.12 or later. The `[mcp]` extra adds the MCP SDK that `smsd-mcp` needs. The
Python loops below need `smsd` without the extra, plus the client package that
each section installs.

`smsd-mcp` (or `python -m smsd.mcp_server`) starts the server over standard
input and output. It takes no arguments, opens no network port and writes
diagnostics to standard error. Assistants start it themselves when they need
the tools. Use the path that `which` or `where.exe` prints wherever a client
needs an absolute path.

To check the server without an assistant, run
[`python/examples/mcp_stdio_client.py`](../python/examples/mcp_stdio_client.py);
it lists the tools and runs one search over standard input and output.

## Connect the MCP server

Clients are listed in alphabetical order. Each block follows the syntax in that
product's documentation as of 29 October 2026. Each result carries
`engine_version`; use it to confirm the SMSD version, because some clients show
the MCP SDK's version for the server.

### Claude Code

```sh
claude mcp add --transport stdio smsd -- smsd-mcp
claude mcp get smsd
```

The server is added to the current project for you only. Add `--scope user`
before `--` to use it in all your projects, or `--scope project` to share it
with your team through `.mcp.json`. Inside a session, `/mcp` shows its status.
Reference: [Claude Code MCP](https://code.claude.com/docs/en/mcp).

### Claude Desktop

Add the server to `claude_desktop_config.json`:

- macOS: `~/Library/Application Support/Claude/claude_desktop_config.json`
- Windows: `%APPDATA%\Claude\claude_desktop_config.json`

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

Replace the command with the absolute path printed by `which smsd-mcp` or
`where.exe smsd-mcp`, because Claude Desktop may not find commands on your
shell's PATH. In JSON, write Windows paths with doubled backslashes or forward
slashes. Quit Claude Desktop fully and restart it. The server then appears
under Connectors in the chat `+` menu. Server output is logged to
`mcp-server-smsd.log` in `~/Library/Logs/Claude` (macOS) or
`%APPDATA%\Claude\logs` (Windows). Reference:
[Connect to local MCP servers](https://modelcontextprotocol.io/docs/develop/connect-local-servers).

### Codex

```sh
codex mcp add smsd -- smsd-mcp
codex mcp list
```

The command writes an `[mcp_servers.smsd]` table to `~/.codex/config.toml`.
The Codex CLI, the Codex IDE extension and the ChatGPT desktop app share this
file, so the server appears in all three. To add it from the IDE extension
instead, open the gear menu, choose MCP servers, then Add server; enter `smsd`,
choose STDIO, enter `smsd-mcp` as the command, save and select Restart
extension. In the terminal interface, `/mcp` lists active servers.

Codex waits a short time for optional servers when it builds the first tool
list. To wait for each server's `startup_timeout_sec` instead, set
`mcp_optional_startup_grace_ms = 0` at the top level of `config.toml`. To make
start-up fail when the SMSD server cannot start, add `required = true` to the
`[mcp_servers.smsd]` table. References:
[Codex MCP](https://learn.chatgpt.com/docs/extend/mcp?surface=cli),
[configuration reference](https://learn.chatgpt.com/docs/config-file/config-reference).

### Gemini CLI and Antigravity CLI

```sh
gemini mcp add --scope user smsd smsd-mcp
gemini mcp list
```

Without `--scope user`, the server is added to `.gemini/settings.json` in the
current folder. Gemini CLI connects stdio servers only in trusted folders; run
`gemini trust` in a folder to trust it. Google moved some account types from
Gemini CLI to Antigravity CLI (`agy`) on 18 June 2026;
[Google's announcement](https://developers.googleblog.com/an-important-update-transitioning-gemini-cli-to-antigravity-cli/)
lists which. Antigravity CLI reads servers from
`~/.gemini/config/mcp_config.json`, or you can add them with its `/mcp`
manager:

```json
{
  "mcpServers": {
    "smsd": {
      "command": "smsd-mcp",
      "args": []
    }
  }
}
```

References: [Gemini CLI MCP](https://geminicli.com/docs/tools/mcp-server/),
[Antigravity MCP](https://antigravity.google/docs/mcp).

### Goose with a local model

Goose is an open-source agent that can use a model served by Ollama. Run
`goose configure`, choose Ollama as the provider, confirm the API host
`http://localhost:11434` and choose a model with tool support, such as `qwen3`
([models with tools](https://ollama.com/search?c=tools)). On machines with less
than 24 GiB of VRAM, Ollama's default context is 4k, which is small for eight
tool definitions; start Ollama with `OLLAMA_CONTEXT_LENGTH=8192 ollama serve`.
Then add SMSD under `extensions:` in `~/.config/goose/config.yaml`
(Windows: `%APPDATA%\Block\goose\config\config.yaml`):

```yaml
extensions:
  smsd:
    type: stdio
    name: smsd
    enabled: true
    cmd: smsd-mcp
    args: []
    env_keys: []
    envs: {}
    timeout: 300
```

References: [Goose configuration files](https://goose-docs.ai/docs/guides/config-files),
[Ollama context length](https://docs.ollama.com/context-length). Other local
MCP clients are listed under [Open models](#open-models).

## Use the tool definitions from Python

`smsd.tool_definitions.get_tool_definitions(provider)` returns the eight tools
in the shape a model API expects. Pass one of the four values below as
`provider`; the default is `mcp`. `call_tool(name, arguments)` runs a tool
locally and returns its JSON-compatible result. It rejects unknown tools,
unknown or missing arguments, wrong JSON types and out-of-range values with a
`ValueError`; the loops below return that message to the model so that it can
correct the call.
[`python/examples/agent_workflow.py`](../python/examples/agent_workflow.py) runs
the same tools in a fixed sequence (screening, MCS and depiction) without a
model.

| Value | Shape | Use with |
| --- | --- | --- |
| `mcp` | `name`, `description`, `inputSchema` | MCP tool lists |
| `input-schema` | `name`, `description`, `input_schema` | Claude Messages API |
| `function` | flat: `type: "function"`, `name`, `description`, `parameters` | OpenAI Responses API; Gemini API after the change shown below |
| `chat-function` | nested: `type: "function"`, `function: {name, description, parameters}` | Chat Completions API, DeepSeek API, llama.cpp, LM Studio, Ollama, vLLM |

The `function` and `chat-function` values set `strict: false`, so the model can
leave out optional arguments and SMSD applies its defaults. Each loop below
reads any API key it needs from the environment and sends every tool result
back to the model in full, including `status`, `optimal` and `complete`. The
APIs are listed in alphabetical order of vendor.

### Claude Messages API

```sh
python -m pip install anthropic
export ANTHROPIC_API_KEY=...    # your own key
```

```python
import json

import anthropic
from smsd.tool_definitions import call_tool, get_tool_definitions


def run_tool(name, arguments):
    """Run one SMSD tool locally and return (text, is_error)."""
    try:
        return json.dumps(call_tool(name, arguments)), False
    except (ValueError, ImportError) as error:
        return str(error), True


client = anthropic.Anthropic()  # reads ANTHROPIC_API_KEY
tools = get_tool_definitions("input-schema")
messages = [{"role": "user", "content": "Find the maximum common substructure of "
             "c1ccccc1O and c1ccccc1N. Is the result proven optimal?"}]

while True:
    response = client.messages.create(
        model="claude-opus-5-5", max_tokens=4096, tools=tools, messages=messages)
    messages.append({"role": "assistant", "content": response.content})
    if response.stop_reason != "tool_use":
        break
    results = []
    for block in response.content:
        if block.type == "tool_use":
            text, is_error = run_tool(block.name, block.input)
            results.append({"type": "tool_result", "tool_use_id": block.id,
                            "content": text, "is_error": is_error})
    messages.append({"role": "user", "content": results})

print("".join(block.text for block in response.content if block.type == "text"))
```

`claude-opus-5-5` was a current model on 29 October 2026; the
[models overview](https://platform.claude.com/docs/en/models/overview) lists
current IDs. Reference:
[Claude tool use](https://platform.claude.com/docs/en/agents-and-tools/tool-use/overview).

### Gemini API

```sh
python -m pip install google-genai
export GEMINI_API_KEY=...    # your own key
```

```python
from google import genai
from google.genai import types
from smsd.tool_definitions import call_tool, get_tool_definitions


def run_tool(name, arguments):
    """Run one SMSD tool locally; errors go back to the model."""
    try:
        return {"result": call_tool(name, arguments)}
    except (ValueError, ImportError) as error:
        return {"error": str(error)}


# Gemini's parameters schema has no additionalProperties field, so drop it.
# call_tool still rejects unknown arguments.
declarations = [
    {"name": tool["name"], "description": tool["description"],
     "parameters": {key: value for key, value in tool["parameters"].items()
                    if key != "additionalProperties"}}
    for tool in get_tool_definitions("function")
]
client = genai.Client()  # reads GEMINI_API_KEY or GOOGLE_API_KEY
config = types.GenerateContentConfig(
    tools=[types.Tool(function_declarations=declarations)], store=False)
contents = [types.Content(role="user", parts=[types.Part(
    text="Find the maximum common substructure of c1ccccc1O and c1ccccc1N. "
         "Is the result proven optimal?")])]

while True:
    response = client.models.generate_content(
        model="gemini-3.8-flash", contents=contents, config=config)
    if not response.function_calls:
        break
    contents.append(response.candidates[0].content)
    parts = [types.Part(function_response=types.FunctionResponse(
        id=call.id, name=call.name, response=run_tool(call.name, call.args or {})))
        for call in response.function_calls]
    contents.append(types.Content(role="user", parts=parts))

print(response.text)
```

Gemini accepts `default` without applying it and may not enforce every limit
in the schema; `call_tool` checks every argument itself. `store=False` opts the
request out of storage, which is on by default for new paid projects from
26 October 2026; it needs google-genai 2.33.0 or later. `gemini-3.8-flash` was
a current model on 29 October 2026; see
[Gemini models](https://ai.google.dev/gemini-api/docs/models). Google
recommends its Interactions API for new projects
([Generate Content](https://ai.google.dev/gemini-api/docs/generate-content/get-started));
its function tools use the flat `type`, `name`, `description` and `parameters`
shape of the `function` value. References:
[function calling with generateContent](https://ai.google.dev/gemini-api/docs/generate-content/function-calling),
[function calling with Interactions](https://ai.google.dev/gemini-api/docs/function-calling).

### OpenAI Responses API

```sh
python -m pip install openai
export OPENAI_API_KEY=...    # your own key
```

```python
import json

from openai import OpenAI
from smsd.tool_definitions import call_tool, get_tool_definitions


def run_tool(name, arguments):
    """Run one SMSD tool locally; errors go back to the model as text."""
    try:
        return json.dumps(call_tool(name, json.loads(arguments)))
    except (ValueError, ImportError) as error:
        return f"Error: {error}"


client = OpenAI()  # reads OPENAI_API_KEY
tools = get_tool_definitions("function")
input_list = [{"role": "user", "content": "Find the maximum common substructure of "
               "c1ccccc1O and c1ccccc1N. Is the result proven optimal?"}]

while True:
    response = client.responses.create(model="gpt-6-astra", tools=tools, input=input_list)
    input_list += response.output
    calls = [item for item in response.output if item.type == "function_call"]
    if not calls:
        break
    for item in calls:
        input_list.append({"type": "function_call_output", "call_id": item.call_id,
                           "output": run_tool(item.name, item.arguments)})

print(response.output_text)
```

`gpt-6-astra` was a current model on 29 October 2026; see
[OpenAI models](https://developers.openai.com/api/docs/models). OpenAI
recommends the Responses API for new projects
([migration guide](https://developers.openai.com/api/docs/guides/migrate-to-responses)).
Reference:
[Function calling](https://developers.openai.com/api/docs/guides/function-calling).

### Chat Completions API

```python
import json

from openai import OpenAI
from smsd.tool_definitions import call_tool, get_tool_definitions


def run_tool(name, arguments):
    """Run one SMSD tool locally; errors go back to the model as text."""
    try:
        return json.dumps(call_tool(name, json.loads(arguments)))
    except (ValueError, ImportError) as error:
        return f"Error: {error}"


client = OpenAI()  # reads OPENAI_API_KEY
tools = get_tool_definitions("chat-function")
messages = [{"role": "user", "content": "Find the maximum common substructure of "
             "c1ccccc1O and c1ccccc1N. Is the result proven optimal?"}]

while True:
    response = client.chat.completions.create(
        model="gpt-5.6-terra", reasoning_effort="none", messages=messages, tools=tools)
    message = response.choices[0].message
    messages.append(message)
    if not message.tool_calls:
        break
    for tool_call in message.tool_calls:
        messages.append({"role": "tool", "tool_call_id": tool_call.id,
                         "content": run_tool(tool_call.function.name,
                                             tool_call.function.arguments)})

print(message.content)
```

On 29 October 2026, OpenAI's
[migration guide](https://developers.openai.com/api/docs/guides/migrate-to-responses)
stated that, from GPT-5.4 onwards, Chat Completions supports tool calls only
with `reasoning_effort="none"`. `gpt-5.6-terra` was available on that date and
supports this setting. GPT-6 Astra, GPT-6.1 Sol and GPT-6.1 Terra need the
Responses API for tool calls; see
[OpenAI models](https://developers.openai.com/api/docs/models). The same loop
works with other Chat Completions services; see [Open models](#open-models).

## Open models

### One loop for DeepSeek and local runtimes

The hosted DeepSeek API and the local runtimes below accept tools in the
`chat-function` shape through Chat Completions-compatible endpoints. This loop
works with each of them; change only `BASE_URL` and `MODEL`. Runtimes are
listed in alphabetical order.

```sh
python -m pip install openai
export DEEPSEEK_API_KEY=...    # DeepSeek only; local runtimes need no key
```

```python
import json
import os
from urllib.parse import urlparse

from openai import OpenAI
from smsd.tool_definitions import call_tool, get_tool_definitions

# Change only these two values; the table below lists each runtime.
BASE_URL = "https://api.deepseek.com"
MODEL = "deepseek-flash"


def run_tool(name, arguments):
    """Run one SMSD tool locally; errors go back to the model as text."""
    try:
        return json.dumps(call_tool(name, json.loads(arguments)))
    except (ValueError, ImportError) as error:
        return f"Error: {error}"


# Local runtimes need no key, but the client needs a value.
# The DeepSeek key is sent only to the DeepSeek API.
if urlparse(BASE_URL).hostname == "api.deepseek.com":
    api_key = os.environ["DEEPSEEK_API_KEY"]
else:
    api_key = "local"
client = OpenAI(base_url=BASE_URL, api_key=api_key)
tools = get_tool_definitions("chat-function")
messages = [{"role": "user", "content": "Find the maximum common substructure of "
             "c1ccccc1O and c1ccccc1N. Is the result proven optimal?"}]

while True:
    response = client.chat.completions.create(model=MODEL, messages=messages, tools=tools)
    message = response.choices[0].message
    messages.append(message)
    if not message.tool_calls:
        break
    for tool_call in message.tool_calls:
        messages.append({"role": "tool", "tool_call_id": tool_call.id,
                         "content": run_tool(tool_call.function.name,
                                             tool_call.function.arguments)})

print(message.content)
```

| Runtime | Start command | `BASE_URL` | `MODEL` |
| --- | --- | --- | --- |
| DeepSeek API | none: set `DEEPSEEK_API_KEY` | `https://api.deepseek.com` | `deepseek-flash` or `deepseek-v4-pro` |
| llama.cpp | `llama-server --jinja -hf bartowski/Qwen2.5-7B-Instruct-GGUF:Q4_K_M --port 9931` | `http://localhost:9931/v1` | the model name listed at `http://localhost:9931/v1/models` |
| LM Studio | `lms server start` | `http://localhost:1234/v1` | a model identifier listed at `http://localhost:1234/v1/models` |
| Ollama | `ollama pull qwen3`, then `ollama serve` | `http://localhost:11434/v1/` | `qwen3` |
| vLLM | `vllm serve Qwen/Qwen3-8B --host 127.0.0.1 --enable-auto-tool-choice --tool-call-parser hermes --reasoning-parser deepseek_r1` | `http://localhost:8000/v1` | `Qwen/Qwen3-8B` |

Notes on each runtime, as of 29 October 2026:

- DeepSeek API: thinking mode is on by default and supports tool calls. The
  loop appends each assistant message whole, which passes back its
  `reasoning_content` as DeepSeek requires. Current model IDs:
  [Models and pricing](https://api-docs.deepseek.com/quick_start/pricing).
  DeepSeek also publishes open weights, such as `DeepSeek-V4.1-Flash-1017`
  (763B parameters;
  [model card](https://huggingface.co/deepseek-ai/DeepSeek-V4.1-Flash-1017)).
- llama.cpp: native tool-call formats include Hermes 2 and 3, Llama 3.1 to
  3.3, Mistral Nemo and Qwen 2.5. DeepSeek R1 is listed as work in progress.
  Other models use a generic format. The default port changed from 8080 to
  9931, so the commands on this page set it. Reference:
  [function calling](https://github.com/ggml-org/llama.cpp/blob/master/docs/function-calling.md).
- LM Studio: Qwen, Llama 3.1, Llama 3.2 and Mistral models have native tool
  use; other models use LM Studio's default tool format. Reference:
  [LM Studio tool use](https://lmstudio.ai/docs/developer/openai-compat/tools).
- Ollama: models with tool support in its library include `gemma4`, `gpt-oss`,
  `llama3.1`, `mistral-small3.2`, `muse-glimmer`, `qwen3` and `qwen3.8`. Some
  entries, such as `deepseek-v4.1-flash`, have only cloud tags; choose a local
  tag to keep data on your machine. Ollama sets its default context from the
  available VRAM: 4k below 24 GiB, 32k from 24 to 48 GiB and 256k from 48 GiB.
  A 4k context is small for eight tool definitions; on such machines start the
  server with `OLLAMA_CONTEXT_LENGTH=8192 ollama serve`. Ollama does not
  support `tool_choice`, and the loop does not send it. References:
  [tools listing](https://ollama.com/search?c=tools),
  [context length](https://docs.ollama.com/context-length).
- vLLM: tool calls need `--enable-auto-tool-choice` and a `--tool-call-parser`
  that matches the model. The command above follows
  [Qwen's vLLM instructions](https://qwen.readthedocs.io/en/latest/framework/function_call.html);
  vLLM's own
  [reasoning page](https://docs.vllm.ai/en/latest/features/reasoning_outputs.html)
  gives `--reasoning-parser qwen3` for the Qwen3 series. vLLM's
  [tool-calling page](https://docs.vllm.ai/en/latest/features/tool_calling.html)
  gives the parser for other families, including DeepSeek, Llama and Mistral.
  For Gemma 4, the
  [Gemma 4 recipe](https://docs.vllm.ai/projects/recipes/en/stable/Google/Gemma4.html)
  gives `--enable-auto-tool-choice --tool-call-parser gemma4 --reasoning-parser gemma4`, together with a chat template.

### MCP with a local model

[Goose](#goose-with-a-local-model) connects `smsd-mcp` to a model served by
Ollama. llama.cpp's `llama-server` can also start stdio MCP servers itself.
Save this as `mcp.json`:

```json
{
  "mcpServers": {
    "smsd": { "command": "smsd-mcp", "args": [], "timeout_ms": 60000 }
  }
}
```

```sh
llama-server --jinja -hf bartowski/Qwen2.5-7B-Instruct-GGUF:Q4_K_M --port 9931 --mcp-servers-config mcp.json
```

llama.cpp stops each tool call after 30 s unless the entry sets `timeout_ms`.
The value above leaves room for an SMSD search at its 30,000 ms maximum. The
SMSD tools then appear as `smsd_find_mcs` and so on in the built-in web
interface at `http://127.0.0.1:9931`. This option is experimental. The server
starts `smsd-mcp` with its own privileges, and the llama.cpp documentation
advises against enabling it in untrusted environments. Reference:
[llama.cpp server](https://github.com/ggml-org/llama.cpp/blob/master/tools/server/README.md).

LM Studio (version 0.3.17 or later) reads MCP servers from its `mcp.json` file:
on the Program tab, choose Install > Edit mcp.json. The file is
`~/.lmstudio/mcp.json` on macOS and Linux and
`%USERPROFILE%/.lmstudio/mcp.json` on Windows
([LM Studio MCP](https://lmstudio.ai/docs/app/mcp),
[LM Studio 0.3.17](https://lmstudio.ai/blog/lmstudio-v0.3.17)). Jan adds MCP
servers under Settings > MCP Servers ([Jan MCP](https://docs.jan.ai/desktop/mcp)).
In either, enter the absolute path to `smsd-mcp` as the command.

## Keeping this current

The commands, settings and model IDs on this page follow each vendor's
documentation as of 29 October 2026. Assistants and model lists change
often, so check the vendor page before you rely on a command or model ID. The
links follow the order of this page:

- Claude Code: [MCP](https://code.claude.com/docs/en/mcp)
- Claude Desktop: [local MCP servers](https://support.claude.com/en/articles/10949351-getting-started-with-local-mcp-servers-on-claude-desktop),
  [configuration file](https://modelcontextprotocol.io/docs/develop/connect-local-servers)
- Codex: [MCP in the CLI](https://learn.chatgpt.com/docs/extend/mcp?surface=cli),
  [MCP in the IDE extension](https://learn.chatgpt.com/docs/extend/mcp?surface=ide),
  [configuration reference](https://learn.chatgpt.com/docs/config-file/config-reference)
- Gemini CLI: [MCP servers](https://geminicli.com/docs/tools/mcp-server/);
  Antigravity CLI: [MCP](https://antigravity.google/docs/mcp),
  [move from Gemini CLI](https://developers.googleblog.com/an-important-update-transitioning-gemini-cli-to-antigravity-cli/)
- Goose: [configuration files](https://goose-docs.ai/docs/guides/config-files)
- Claude API: [tool use](https://platform.claude.com/docs/en/agents-and-tools/tool-use/overview),
  [models](https://platform.claude.com/docs/en/models/overview)
- Gemini API: [function calling](https://ai.google.dev/gemini-api/docs/generate-content/function-calling),
  [Interactions function calling](https://ai.google.dev/gemini-api/docs/function-calling),
  [models](https://ai.google.dev/gemini-api/docs/models)
- OpenAI API: [function calling](https://developers.openai.com/api/docs/guides/function-calling),
  [migration to Responses](https://developers.openai.com/api/docs/guides/migrate-to-responses),
  [models](https://developers.openai.com/api/docs/models)
- DeepSeek API: [tool calls](https://api-docs.deepseek.com/guides/tool_calls),
  [models](https://api-docs.deepseek.com/quick_start/pricing)
- llama.cpp: [server](https://github.com/ggml-org/llama.cpp/blob/master/tools/server/README.md),
  [function calling](https://github.com/ggml-org/llama.cpp/blob/master/docs/function-calling.md)
- LM Studio: [tool use](https://lmstudio.ai/docs/developer/openai-compat/tools),
  [MCP](https://lmstudio.ai/docs/app/mcp),
  [0.3.17 release notes](https://lmstudio.ai/blog/lmstudio-v0.3.17)
- Ollama: [compatibility endpoint](https://docs.ollama.com/api/openai-compatibility),
  [models with tools](https://ollama.com/search?c=tools),
  [context length](https://docs.ollama.com/context-length)
- vLLM: [tool calling](https://docs.vllm.ai/en/latest/features/tool_calling.html),
  [reasoning outputs](https://docs.vllm.ai/en/latest/features/reasoning_outputs.html),
  [Gemma 4 recipe](https://docs.vllm.ai/projects/recipes/en/stable/Google/Gemma4.html);
  Qwen: [function calling](https://qwen.readthedocs.io/en/latest/framework/function_call.html)
- Jan: [MCP](https://docs.jan.ai/desktop/mcp)

## Safety

- The tools run on your machine and contact no network service. Molecular
  arguments are SMILES or SMARTS text, never paths, URLs or commands.
- Limits are fixed: 16,384 characters per text argument and 1,024 atoms per
  molecule (256 for depiction and 128 for standardisation). A similarity or
  screening call takes at most 4,096 atoms in total, and a screen at most 64
  targets and 65,536 characters. `search_smarts` returns at most 32 matches,
  and `screen_molecules` at most 32 hits per page. A search stops at
  `timeout_ms` (default 1,000, at most 30,000) or `node_limit` (default
  100,000, at most 1,000,000). These budgets cover the search, not the whole
  call. Standardisation has no time budget. Use process or container limits for
  a hard deadline.
  The [container recipe](../packaging/mcp/README.md) runs the server with no
  network access.
- `smsd-mcp` and `call_tool` check every argument themselves, so a client that
  drops schema keywords cannot widen these limits.
- Ask the assistant to report `status`, `optimal` and `complete` with each
  answer. A search that stopped early may not be maximum, and `found: false`
  from an unfinished search does not prove absence.
- A hosted assistant sends your prompts, tool arguments and results, including
  any structures, to its provider. Use a local runtime when structures must
  stay on your machine.
- API keys come from environment variables such as `ANTHROPIC_API_KEY`,
  `OPENAI_API_KEY`, `GEMINI_API_KEY` and `DEEPSEEK_API_KEY`. Do not write keys
  into code or committed configuration files. `smsd-mcp` needs no key, so its
  client entries need no environment settings.
- Keep your assistant's tool approval prompts on unless you trust every tool it
  can call.
