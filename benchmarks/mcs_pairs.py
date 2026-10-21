# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Load the curated MCS pair corpus without importing a chemistry toolkit.

benchmark_leaderboard.py measures the Java CLI against RDKit, so it must not
depend on the smsd Python package merely to read the pair list.
"""
from pathlib import Path


MCS_PAIRS_FILE = Path(__file__).resolve().parent / "data" / "mcs_pairs.tsv"


def load_mcs_pairs(path=MCS_PAIRS_FILE):
    """Return the curated (query, target, name, category) MCS pairs in file order."""
    pairs = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line and not line.startswith("#"):
            fields = tuple(line.split("\t"))
            if len(fields) != 4:
                raise ValueError("Expected four tab-separated fields: " + line)
            pairs.append(fields)
    return pairs
