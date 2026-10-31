# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Chain the molecular tools without a model: screen a library, compare hits, draw the best."""

from pathlib import Path
from tempfile import TemporaryDirectory

from smsd.tool_definitions import call_tool, get_tool_definitions


def main():
    print("Tools:", ", ".join(definition["name"] for definition in get_tool_definitions()))
    query = "COc1ccccc1"
    library = ["c1ccccc1O", "c1ccccc1N", "Cc1ccccc1", "CCO", "CCOc1ccccc1", "C1CCCCC1",
               "COc1ccc(O)cc1", "CC(=O)O"]
    page = call_tool("screen_molecules", {"query_smiles": query, "targets_smiles": library,
                                          "threshold": 0.3, "max_hits": 3})
    print("Screen hits:", page["total_hits"], "of", page["total_targets"])
    for hit in page["hits"]:
        target = library[hit["target_index"]]
        result = call_tool("find_mcs", {"query_smiles": query, "target_smiles": target,
                                        "timeout_ms": 1000, "node_limit": 100000})
        atoms = result["matches"][0]["atom_count"] if result["found"] else 0
        print(f"{target}: Tanimoto={hit['score']:.3f}, status={result['status']}, "
              f"optimal={result['optimal']}, complete={result['complete']}, atoms={atoms}")

    best = library[page["hits"][0]["target_index"]]
    drawing = call_tool("depict_mcs", {"query_smiles": query, "target_smiles": best,
                                       "width": 600, "height": 300})
    with TemporaryDirectory() as directory:
        path = Path(directory) / "best_hit.svg"
        path.write_text(drawing["svg"], encoding="utf-8")
        print("Saved", path.name, "for", best, f"({path.stat().st_size} bytes)")

    # A calling program returns this message to the model as the tool error.
    try:
        call_tool("find_mcs", {"query_smiles": query})
    except ValueError as error:
        print("Tool error:", error)
    # To let a model choose these calls, see docs/AI_ASSISTANTS.md.


if __name__ == "__main__":
    main()
