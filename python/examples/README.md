# Python examples

`feature_search.ipynb` demonstrates MCS-to-SMARTS library search, result
status, SVG highlighting, optional NumPy features and a scikit-learn
pipeline, followed by a temporary mapped fingerprint corpus, top-k
neighbours, Butina clustering and MaxMin selection.

Run its setup cell first, in Jupyter or Colab. It installs SMSD 8.2.0 or
later into the current kernel when it is missing or older, then attempts to
install NumPy for feature arrays and mapped storage and scikit-learn for the
pipeline; installation needs internet access, and the cells that use an
optional package report when it is missing. To run it locally, install
Jupyter in the same environment. The saved notebook has empty outputs, and
running it leaves no corpus file behind.

Five further notebooks follow the chapters of the
[user guide](../../docs/guide/README.md): `substructure_search.ipynb` searches
drugs and benzodiazepines with molecule queries and SMARTS patterns,
`mcs_search.ipynb` compares the MCS objectives, partial results and MCS
similarity, `fingerprints_similarity.ipynb` computes every fingerprint family
and coefficient, `library_screening.ipynb` screens one library five ways and
saves the hits, and `cpu_gpu_screening.ipynb` shows which screening steps can
run on a GPU and runs every cell without one. They follow the same
conventions: a setup cell that installs SMSD 8.2.0 or later and the optional
packages each notebook uses, empty saved outputs, no files left behind, and a
notice when an optional package such as NumPy, pandas or IPython is missing.

`ai_tools.ipynb` calls the molecular tools from Python with no model, prints
their tool definitions and runs an assistant tool loop against a model served
on your own machine, so that no molecule leaves it. Model requests are off by
default; set `RUN_LOCAL_MODEL` to `True` in the setup cell to run the loop.
When the switch is off or no local model server answers, the loop cells say so
and the rest of the notebook still runs. The
[user guide](../../docs/guide/README.md) also has chapters on
[clustering and diversity](../../docs/guide/CLUSTERING.md),
[machine learning and data frames](../../docs/guide/MACHINE_LEARNING.md) and
[chemistry tools for AI assistants](../../docs/guide/AI_TOOLS.md).

`feature_adapters.py` computes binary and count feature matrices with NumPy.
`mcs_highlight.py` sends an SVG to standard output and its native proof status
to standard error. `packed_corpus.py` creates a temporary packed library and
streams similarity and diversity results. `gpu_screening.py` screens with
NVIDIA CUDA when `smsd-cuda` and a device are present and checks every result
against the CPU. `agent_workflow.py` chains the molecular tools through
`call_tool` without a model, and `mcp_stdio_client.py`, which needs
`smsd[mcp]`, checks the MCP server from Python. The
[adapter guide](../../docs/PYTHON_ADAPTERS.md) describes the APIs, optional MCP
server and macOS PyTorch runtime requirements.
