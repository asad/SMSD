# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Small local MCS demonstration; deployment is an explicit operator action."""
import base64
import html
import os

from smsd.mcp_server import _safe_svg, depict_mcs


def compare(query_smiles: str, target_smiles: str, objective: str = "mcis") -> tuple:
    """Return a bounded MCS result and a non-executable SVG image element.

    The JSON retains proof/stopping flags. Molecular input never becomes raw
    HTML; only validated native SVG is encoded as an image.

    Example::

        result, image = compare("CCO", "CCN")
        assert result["found"] and image.startswith('<img ')
    """
    try:
        result = depict_mcs(query_smiles, target_smiles, objective=objective,
                            timeout_ms=1000, node_limit=100000, width=800, height=360)
        svg = result.pop("svg")
        if svg is None:
            return result, "<p>No mapping to draw. Inspect the stopping status.</p>"
        encoded = base64.b64encode(_safe_svg(svg).encode("utf-8")).decode("ascii")
        return result, ('<img alt="Mapped atoms in the two input molecules" '
                        'style="max-width:100%;height:auto" src="data:image/svg+xml;base64,' + encoded + '">')
    except (ValueError, TypeError, RuntimeError) as error:
        message = str(error)[:1000]
        return {"error": message, "complete": False}, "<p>" + html.escape(message) + "</p>"


def create_demo():
    """Create a Gradio interface without launching or publishing it.

    Example::

        demo = create_demo()  # requires gradio
    """
    import gradio as gr
    with gr.Blocks(analytics_enabled=False, title="SMSD MCS") as demo:
        gr.Markdown("# SMSD MCS\nCompare two molecular graphs. An unfinished result is not an optimum proof.")
        with gr.Row():
            query = gr.Textbox(label="Query SMILES", value="CCO", max_lines=3)
            target = gr.Textbox(label="Target SMILES", value="CCN", max_lines=3)
        objective = gr.Radio(["mcis", "mces"], value="mcis", label="Objective")
        gr.Markdown("Connected search; 1 s / 100,000 nodes after parsing; at most 256 atoms per input. "
                    "The picture marks mapped atoms and input bonds between them; MCES selected edges are in JSON.")
        run = gr.Button("Find common graph")
        picture = gr.HTML()
        result = gr.JSON(label="Result, proof flags and limits")
        run.click(compare, [query, target, objective], [result, picture], concurrency_limit=1)
    return demo.queue(max_size=8, default_concurrency_limit=1)


if __name__ == "__main__":
    create_demo().launch(server_name=os.environ.get("GRADIO_SERVER_NAME", "127.0.0.1"),
                         share=False, show_error=False, max_threads=1)
