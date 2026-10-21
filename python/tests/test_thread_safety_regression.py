# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Shared native objects must remain correct while the GIL is released."""
import subprocess
import sys
import textwrap

import pytest


def _run_concurrent_check(source):
    result = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(source)],
        text=True, capture_output=True, timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_compiled_recursive_query_shared_across_targets():
    _run_concurrent_check("""
        from concurrent.futures import ThreadPoolExecutor
        from threading import Barrier
        import smsd
        query = smsd.Substructure.smarts('[C;$(CO)]')
        targets = [smsd.parse_smiles('C' * 60 + 'O'),
                   smsd.parse_smiles('OCCO'), smsd.parse_smiles('CCCC')]
        for target in targets:
            target.prewarm()
        barrier = Barrier(6)
        def match(worker):
            index = worker % 3
            barrier.wait()
            for _ in range(150):
                assert query.match(targets[index]).found == (index != 2)
        with ThreadPoolExecutor(max_workers=6) as executor:
            list(executor.map(match, range(6)))
    """)


@pytest.mark.parametrize("operation", ["substructure", "mcs", "prewarm"])
def test_shared_cold_graphs(operation):
    _run_concurrent_check(f"""
        from concurrent.futures import ThreadPoolExecutor
        from threading import Barrier
        import smsd
        for trial in range(3):
            query = smsd.parse_smiles('c1ccc(CCO)cc1')
            target = smsd.parse_smiles('c1ccc(CCO)cc1O')
            barrier = Barrier(6)
            def match(worker):
                barrier.wait()
                for _ in range(8):
                    if '{operation}' == 'prewarm':
                        query.prewarm()
                        target.prewarm()
                    elif '{operation}' == 'mcs':
                        assert smsd.MCS.find(query).match(target).best.atom_count == query.n
                    else:
                        assert smsd.Substructure.find(query).match(target).found
                    for values in [query.morgan_rank, query.canonical_label, query.orbit]:
                        assert len(values) in (0, query.n)
            with ThreadPoolExecutor(max_workers=6) as executor:
                list(executor.map(match, range(6)))
    """)


def test_parallel_batches_share_a_cold_target_library():
    _run_concurrent_check("""
        from concurrent.futures import ThreadPoolExecutor
        from threading import Barrier
        import smsd
        targets = [smsd.parse_smiles('c1ccc(CCO)cc1O') for _ in range(12)]
        queries = [smsd.parse_smiles(text) for text in ['CCO', 'CO', 'c1ccccc1', 'O']]
        barrier = Barrier(4)
        def match(worker):
            barrier.wait()
            for _ in range(4):
                assert [result.found for result in smsd.Substructure.find(queries[worker]).match_many(targets, num_threads=1)] == [True] * 12
        with ThreadPoolExecutor(max_workers=4) as executor:
            list(executor.map(match, range(4)))
    """)


def test_corpus_additions_can_overlap_queries():
    _run_concurrent_check("""
        from concurrent.futures import ThreadPoolExecutor
        from threading import Barrier
        import smsd
        corpus = smsd.TargetCorpus()
        corpus.add_target(smsd.parse_smiles('CCO'))
        corpus.prewarm(num_threads=1)
        query = smsd.parse_smiles('CO')
        barrier = Barrier(4)
        def operate(worker):
            barrier.wait()
            if worker == 0:
                for _ in range(16):
                    corpus.add_target(smsd.parse_smiles('CCO'))
            else:
                for _ in range(12):
                    result = [match.found for match in smsd.Substructure.find(query).match_many(corpus, num_threads=1)]
                    assert 1 <= len(result) <= 17 and all(result)
        with ThreadPoolExecutor(max_workers=4) as executor:
            list(executor.map(operate, range(4)))
        assert len(corpus) == 17
    """)


def test_layout_first_use_from_multiple_threads():
    _run_concurrent_check("""
        from concurrent.futures import ThreadPoolExecutor
        from threading import Barrier
        import math
        import smsd
        graphs = [smsd.parse_smiles('c1ccc2ccccc2c1') for _ in range(8)]
        barrier = Barrier(8)
        def layout(worker):
            barrier.wait()
            points = smsd.generate_coords_2d(graphs[worker])
            assert len(points) == graphs[worker].n
            assert all(math.isfinite(value) for point in points for value in point)
        with ThreadPoolExecutor(max_workers=8) as executor:
            list(executor.map(layout, range(8)))
    """)


def test_kekulising_a_prewarmed_graph_invalidates_matching_properties():
    import smsd

    query = smsd.parse_smiles("c1ccccc1")
    target = smsd.parse_smiles("Cc1ccccc1")
    query.prewarm()
    query.kekulize()
    target.kekulize()
    options = smsd.Substructure.Options()
    options.chemistry.aromaticity_mode = smsd.AromaticityMode.STRICT
    assert smsd.Substructure.find(query).with_options(options).match_all(target).found
