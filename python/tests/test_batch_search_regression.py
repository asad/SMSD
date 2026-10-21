# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Batch searches equal one search per target for any thread count, and leave no threads behind."""
import multiprocessing
import os
import subprocess
import sys
import textwrap

import pytest
import smsd

_TARGETS = [
    "CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1", "c1ccc2ccccc2c1", "OCCO", "CCN(CC)CC",
    "COc1ccc(NS(c2csc(C(=O)Nc3ccccc3)c2)(=O)=O)cc1", "Clc1cc2nc(N3CCNCC3)n(Cc3ccccc3)c2cc1",
    "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "O=C1NC(=O)c2ccccc21", "C1CCCCC1",
]


def _run(source, threads_variable):
    environment = dict(os.environ)
    environment.pop("OMP_NUM_THREADS", None)
    if threads_variable is not None:
        environment["OMP_NUM_THREADS"] = threads_variable
    result = subprocess.run([sys.executable, "-c", textwrap.dedent(source)], env=environment,
                            text=True, capture_output=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("threads_variable", [None, "1", "3"])
def test_batches_equal_single_searches_for_every_thread_count(threads_variable):
    _run(f"""
        import smsd
        smiles = {_TARGETS!r} * 3
        targets = [text if index % 2 else smsd.parse_smiles(text) for index, text in enumerate(smiles)]
        corpus = smsd.TargetCorpus.from_smiles(smiles)
        limited = smsd.MCS.Options()
        limited.node_limit = 25
        # Generous time limits leave node counts, not the clock, to end every search.
        queries = [query.with_timeout(60000) for query in [
            smsd.Substructure.find("c1ccccc1"),
            smsd.Substructure.find("CC(=O)N").with_induced(),
            smsd.Substructure.smarts("[OX2H][#6]"),
            smsd.MCS.find("c1ccc(O)cc1"),
            smsd.MCS.find("CC(=O)Nc1ccccc1").with_objective(smsd.MCS.Objective.MCES).with_fragments(),
            smsd.MCS.find("CC(=O)Nc1ccc(O)cc1C(=O)OC").with_options(limited),
        ]]
        def fields(result):
            return (result.matches, result.status, result.optimal, result.complete, result.score,
                    result.nodes_visited)
        for query in queries:
            single = [fields(query.match(target)) for target in targets]
            for threads in (0, 1, 4):
                assert [fields(result) for result in query.match_many(targets, num_threads=threads)] == single
                assert [fields(result) for result in query.match_many(corpus, num_threads=threads)] == single
    """, threads_variable)


def test_batch_failures_propagate_and_later_batches_run():
    options = smsd.MCS.Options()
    options.excluded_target_atoms = [4]
    prepared = smsd.MCS.find("CCO").with_options(options)
    for threads in (1, 4):
        with pytest.raises(IndexError):
            prepared.match_many(["CCCCCO"] * 20 + ["CC"] + ["CCCCCCO"] * 20, num_threads=threads)
        assert [result.found for result in prepared.match_many(["CCCCCO"] * 20, num_threads=threads)] == [True] * 20


@pytest.mark.skipif("fork" not in multiprocessing.get_all_start_methods(), reason="needs the fork start method")
def test_a_forked_child_searches_after_its_parent():
    _run(f"""
        import multiprocessing
        import warnings
        import smsd
        targets = [smsd.parse_smiles(text) for text in {_TARGETS!r} * 8]
        query = smsd.MCS.find("CC(=O)Nc1ccc(O)cc1")
        def fields(results):
            return [(result.score, result.nodes_visited, [sorted(match.atom_mapping.items()) for match in result.matches])
                    for result in results]
        expected = fields(query.match_many(targets, num_threads=4))
        def child(queue):
            queue.put(fields(query.match_many(targets, num_threads=4)))
        context = multiprocessing.get_context("fork")
        queue = context.Queue()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            process = context.Process(target=child, args=(queue,))
            process.start()
        found = queue.get(timeout=60)
        process.join(60)
        assert process.exitcode == 0 and found == expected
        # No batch thread outlives its call, so the parent forked a single-threaded process.
        assert not [item for item in caught if "multi-threaded" in str(item.message)]
        assert fields(query.match_many(targets, num_threads=4)) == expected
    """, None)
