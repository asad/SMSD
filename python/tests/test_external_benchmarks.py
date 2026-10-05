# SPDX-License-Identifier: Apache-2.0
"""Opt-in feature diagnostics over the checked-in external corpora.

These inputs include MoleculeNet-derived Dalke-style pairs, a Tautobase subset,
12 curated stress pairs and 1,400 collected SMARTS patterns. They do not
reproduce the source papers' complete experiments. For independently validated
RDKit comparisons use benchmarks/run_external_benchmarks.py.

Run locally with SMSD_BENCHMARK=1. SMSD_BENCHMARK_TIMEOUT_MS defaults to 10000;
SMSD_BENCHMARK_LIMIT=0 (the default) runs complete corpora. Returned mapping
calls expose no cancellation status; wall-clock budget crossings are reported
as diagnostics, not algorithmic timeout counts or performance assertions.
"""
import csv
import os
from pathlib import Path
import statistics
import time

import pytest

smsd = pytest.importorskip("smsd")
DATA_DIR = Path(__file__).resolve().parents[2] / "benchmarks/data"
benchmark = pytest.mark.skipif(not os.environ.get("SMSD_BENCHMARK"),
                               reason="Set SMSD_BENCHMARK=1 to run external benchmarks")
TIMEOUT_MS = int(os.environ.get("SMSD_BENCHMARK_TIMEOUT_MS", "10000"))
LIMIT = int(os.environ.get("SMSD_BENCHMARK_LIMIT", "0"))


def load_tsv_pairs(filename):
    path = DATA_DIR / filename
    if not path.exists():
        pytest.skip(f"{path} not found")
    rows = [line.split("\t") for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")]
    return rows[:LIMIT] if LIMIT else rows


def load_tautobase():
    with (DATA_DIR / "chodera_tautobase_subset.txt").open() as handle:
        rows = [row[:3] for row in csv.reader(handle, skipinitialspace=True)
                if len(row) >= 3 and row[0] != "name"]
    return rows[:LIMIT] if LIMIT else rows


def timed_pairs(pairs, *, tautomer_aware=False):
    times, sizes, failures = [], [], []
    for name, smi1, smi2 in pairs:
        try:
            q, t = smsd.parse_smiles(smi1), smsd.parse_smiles(smi2)
            start = time.perf_counter_ns()
            result = smsd.find_mcs(q, t, timeout_ms=TIMEOUT_MS, tautomer_aware=tautomer_aware)
            times.append((time.perf_counter_ns()-start)/1000)
            sizes.append(len(result))
        except Exception as exc:
            failures.append((name, str(exc)))
    print(f"calls={len(times)} errors={len(failures)} "
          f"median_us={statistics.median(times) if times else None} "
          f"wall_budget_crossings={sum(v >= TIMEOUT_MS*1000 for v in times)} "
          "cancellation=not_exposed")
    if failures:
        print(f"first_errors={failures[:5]}")
    assert len(times)+len(failures) == len(pairs)
    assert times, "no successful calls"
    return sizes


class TestTautobase:
    @benchmark
    def test_tautomer_recovery(self):
        pairs = load_tautobase()
        gains, failures = [], []
        for name, smi1, smi2 in pairs:
            try:
                q, t = smsd.parse_smiles(smi1), smsd.parse_smiles(smi2)
                strict = smsd.find_mcs(q, t, timeout_ms=TIMEOUT_MS)
                tautomer = smsd.find_mcs(q, t, timeout_ms=TIMEOUT_MS, tautomer_aware=True)
                gains.append(len(tautomer)-len(strict))
            except Exception as exc:
                failures.append((name, str(exc)))
        print(f"Tautobase calls={len(gains)} errors={len(failures)} "
              f"positive_atom_deltas={sum(v > 0 for v in gains)} "
              f"negative_atom_deltas={sum(v < 0 for v in gains)} sum_delta={sum(gains)}")
        if failures:
            print(f"smsd_module={smsd.__file__} first_errors={failures[:5]}")
        assert len(gains)+len(failures) == len(pairs)
        assert gains, "no successful tautomer calls"

    @benchmark
    def test_budget_diagnostics(self):
        timed_pairs(load_tautobase())


class TestDalkeBenchmark:
    @benchmark
    def test_random_pairs(self):
        rows = load_tsv_pairs("dalke_random_pairs.tsv")
        timed_pairs([(f"random_{i}", row[0], row[1]) for i, row in enumerate(rows)])

    @benchmark
    def test_nn_pairs(self):
        rows = load_tsv_pairs("dalke_nn_pairs.tsv")
        timed_pairs([(f"nn_{i}", row[0], row[1]) for i, row in enumerate(rows)])


class TestStressPairs:
    @benchmark
    def test_budget_diagnostics(self):
        rows = load_tsv_pairs("stress_pairs.tsv")
        timed_pairs([(row[2], row[0], row[1]) for row in rows])

    @benchmark
    def test_self_match_exact(self):
        for row in load_tsv_pairs("stress_pairs.tsv"):
            if len(row) < 3 or "self" not in row[2]:
                continue
            q, t = smsd.parse_smiles(row[0]), smsd.parse_smiles(row[1])
            result = smsd.find_mcs(q, t, timeout_ms=TIMEOUT_MS)
            assert len(result) == q.n, f"{row[2]}: MCS={len(result)}, expected {q.n}"


class TestEhrlichRarey:
    @benchmark
    def test_smarts_performance(self):
        patterns = load_tsv_pairs("ehrlich_rarey_smarts.txt")
        target = smsd.parse_smiles("CC(C)Cc1ccc(CC(C)C(O)=O)cc1")
        times, hits, failures = [], 0, []
        for row in patterns:
            try:
                start = time.perf_counter_ns()
                hit = smsd.smarts_match(row[0], target)
                times.append((time.perf_counter_ns()-start)/1000)
                hits += bool(hit)
            except Exception as exc:
                failures.append((row[0], str(exc)))
        print(f"SMARTS calls={len(times)} errors={len(failures)} hits={hits} "
              f"median_us={statistics.median(times) if times else None}; includes compilation")
        assert len(times)+len(failures) == len(patterns)
        assert times, "no successful SMARTS calls"
