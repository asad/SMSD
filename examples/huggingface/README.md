---
title: SMSD MCS
emoji: 🧪
colorFrom: green
colorTo: blue
sdk: docker
app_port: 7860
pinned: false
---

# Molecular common-graph demonstration

This is source for a small Gradio/Hugging Face Docker Space demonstration. It
is not a deployed Space. The Docker recipe explicitly uses CPython 3.12 on
Linux. Supply a reviewed matching SMSD 8.1.0 CPython 3.12 wheel for the image
architecture, from PyPI or the [GitHub release](https://github.com/asad/SMSD/releases),
and compatible Gradio 5 or 6 dependency wheels in `wheelhouse/`. SMSD is
installed from that local wheel only, and the build stops if it is not 8.1.0.

```sh
cd examples/huggingface
docker build --network=none --build-arg SMSD_WHEEL=YOUR_MATCHING_SMSD_WHEEL.whl \
  -t smsd-mcs-demo:8.1.0 .
docker run --rm --memory=1g --cpus=1 -p 127.0.0.1:7860:7860 smsd-mcs-demo:8.1.0
```

Replace the wheel placeholder with the actual `smsd-8.1.0-...whl` basename.
All installation inputs must already be in the wheelhouse. A cached base image
is needed for offline building; pin a reviewed `PYTHON_BASE` digest for a
reproducible image. `--network=none` restricts build steps, not daemon base-image
fetches. For a local Python route, create a CPython 3.12 environment, install
the matching wheel with `--no-index --find-links wheelhouse`, then install
`requirements.txt` with the same flags and run `python app.py`.

The Python script's default listener is local (`127.0.0.1`) and `share=False`;
importing it does not launch anything. The container explicitly sets
`GRADIO_SERVER_NAME=0.0.0.0` for Docker/Space routing. Review platform access
controls before deploying it.

Inputs are molecular text, not uploaded files or HTML. Calls use the bounded
`depict_mcs` tool: connected MCIS or MCES, 1,000 ms and 100,000 nodes, at most
256 atoms per molecule. The queue holds at most eight requests and processes
one at a time. These search limits start after parsing; use process limits
for a hard total deadline and memory cap.

The SVG comes only from the native renderer, is checked for active/external
content, and is embedded as a base64 image rather than raw HTML. Error messages
are escaped. Search proof/status stays in the JSON; a stopped search is not a
proved optimum. MCES selected edges are given by `bond_mapping`, while the
picture highlights mapped atoms and input bonds between them. There are no
model calls, API keys, external molecule lookups or file inputs.

See the [MCP tool contracts](../../docs/MCP_SERVER.md) for chemistry options,
limits, standardisation and result interpretation.
