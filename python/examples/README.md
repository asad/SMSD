# Python examples

`feature_search.ipynb` demonstrates MCS-to-SMARTS library search, result
status, SVG highlighting, optional NumPy features and a scikit-learn
pipeline, followed by a temporary mapped fingerprint corpus, top-k
neighbours, Butina clustering and MaxMin selection.

Use a kernel containing the matching SMSD 8.1.0 package, installed from PyPI,
a GitHub wheel or a source checkout. Install NumPy
for feature arrays and mapped storage, and scikit-learn for the pipeline;
those cells report when an optional dependency is missing. Install Jupyter
in the same environment to open the notebook. The saved notebook has empty
outputs; running it performs no downloads and leaves no corpus file behind.

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
