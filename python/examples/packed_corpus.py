# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Create a temporary packed library and stream similarity and diversity results."""

from pathlib import Path
from tempfile import TemporaryDirectory

import smsd
from smsd.clustering import butina, maxmin
from smsd.corpus import FingerprintCorpus, stream_screen_smiles


def main():
    molecules = ["CCO", "CCN", "c1ccccc1", "CC(=O)O", "CCO"]
    options = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
    with TemporaryDirectory() as directory:
        path = Path(directory) / "library.fp"
        with FingerprintCorpus.from_smiles(path, iter(molecules), options=options,
                                           chunk_size=2, num_threads=1) as corpus:
            query = smsd.fingerprint("CCO", options=options)
            print("Query hits:", list(corpus.similarities(query, threshold=0.3)))
            print("Threshold pairs:", list(corpus.threshold_pairs(0.5)))
            print("Top neighbours:", list(corpus.top_k_pairs(1)))
            print("Butina clusters:", butina(corpus, threshold=0.5))
            print("MaxMin selection:", maxmin(corpus, count=3))
    print("Native streamed hits:", list(stream_screen_smiles(
        "CCO", iter(molecules), 0.5, chunk_size=2, num_threads=1)))


if __name__ == "__main__":
    main()
