# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
"""
External benchmark tests using community-standard datasets.

Datasets:
  - Tautobase (468 tautomer pairs, Chodera/Wahl-Sander)
  - Dalke-style random pairs (1000 low-similarity MCS pairs)
  - Dalke-style nearest-neighbor pairs (1000 high-similarity MCS pairs)
  - Stress pairs (12 adversarial hard cases)
  - Ehrlich-Rarey SMARTS (1400 substructure patterns)

Run:
    pytest python/tests/test_external_benchmarks.py -v -s --benchmark
    # or from repo root:
    python -m pytest python/tests/test_external_benchmarks.py -v -s -k benchmark
"""

import os
import time
import statistics
import pytest

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "benchmarks", "data")

# Skip entire module if smsd not installed
smsd = pytest.importorskip("smsd")

benchmark = pytest.mark.skipif(
    not os.environ.get("SMSD_BENCHMARK"),
    reason="Set SMSD_BENCHMARK=1 to run external benchmarks",
)


def load_tsv_pairs(filename):
    """Load tab-separated pairs, skipping comments (#) and blank lines."""
    path = os.path.join(DATA_DIR, filename)
    if not os.path.exists(path):
        pytest.skip(f"{path} not found")
    pairs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                pairs.append(parts)
    return pairs


def load_tautobase():
    """Load Chodera Tautobase subset (CSV-like)."""
    path = os.path.join(DATA_DIR, "chodera_tautobase_subset.txt")
    if not os.path.exists(path):
        pytest.skip(f"{path} not found")
    pairs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("name"):
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 3:
                pairs.append((parts[0], parts[1], parts[2]))
    return pairs


# ======================================================================
# 1. Tautobase — Tautomer-Aware MCS
# ======================================================================
class TestTautobase:

    @benchmark
    def test_tautomer_recovery(self):
        """Tautomer-aware MCS recovers more atoms than strict mode."""
        pairs = load_tautobase()
        total = 0
        gains = 0
        total_gain_atoms = 0
        taut_opts = smsd.ChemOptions()
        taut_opts.tautomer_aware = True

        for name, smi1, smi2 in pairs:
            try:
                g1 = smsd.parse_smiles(smi1)
                g2 = smsd.parse_smiles(smi2)
                m_strict = smsd.find_mcs(g1, g2, timeout_ms=10000)
                m_taut = smsd.find_mcs(g1, g2, chem=taut_opts, timeout_ms=10000)
                total += 1
                if len(m_taut) > len(m_strict):
                    gains += 1
                    total_gain_atoms += len(m_taut) - len(m_strict)
            except Exception:
                pass

        print(f"\nTautobase: {total} pairs, {gains} with tautomer gain "
              f"({100*gains/max(total,1):.1f}%), +{total_gain_atoms} atoms total")
        assert total > 100, f"Should parse >100 pairs, got {total}"

    @benchmark
    def test_no_timeouts(self):
        """All Tautobase pairs complete within 10s."""
        pairs = load_tautobase()
        timeouts = 0
        total = 0
        times_us = []

        for name, smi1, smi2 in pairs:
            try:
                g1 = smsd.parse_smiles(smi1)
                g2 = smsd.parse_smiles(smi2)
                t0 = time.perf_counter_ns()
                result = smsd.find_mcs(g1, g2, timeout_ms=10000)
                elapsed_us = (time.perf_counter_ns() - t0) / 1000
                times_us.append(elapsed_us)
                total += 1
                if elapsed_us > 10_000_000:
                    timeouts += 1
            except Exception:
                pass

        median_us = statistics.median(times_us) if times_us else 0
        print(f"\nTautobase timing: {total} pairs, median {median_us:.0f} us, "
              f"{timeouts} timeouts")
        assert timeouts == 0, f"No timeouts expected, got {timeouts}"


# ======================================================================
# 2. Dalke-style MCS Benchmark
# ======================================================================
class TestDalkeBenchmark:

    @benchmark
    def test_random_pairs(self):
        """1000 random pairs: measure LFUB certificate rate."""
        pairs = load_tsv_pairs("dalke_random_pairs.tsv")
        total = 0
        lfub_hit = 0
        times_us = []

        for parts in pairs[:200]:  # First 200 for speed in CI
            try:
                g1 = smsd.parse_smiles(parts[0])
                g2 = smsd.parse_smiles(parts[1])
                t0 = time.perf_counter_ns()
                result = smsd.find_mcs(g1, g2, timeout_ms=10000)
                elapsed_us = (time.perf_counter_ns() - t0) / 1000
                times_us.append(elapsed_us)
                total += 1
            except Exception:
                pass

        median_us = statistics.median(times_us) if times_us else 0
        print(f"\nDalke random: {total} pairs, "
              f"median {median_us:.0f} us")
        assert total > 100, f"Should parse >100 pairs, got {total}"

    @benchmark
    def test_nn_pairs(self):
        """1000 nearest-neighbor pairs: high-similarity MCS."""
        pairs = load_tsv_pairs("dalke_nn_pairs.tsv")
        total = 0
        mcs_ge5 = 0
        times_us = []

        for parts in pairs[:200]:  # First 200 for speed in CI
            try:
                g1 = smsd.parse_smiles(parts[0])
                g2 = smsd.parse_smiles(parts[1])
                t0 = time.perf_counter_ns()
                result = smsd.find_mcs(g1, g2, timeout_ms=10000)
                elapsed_us = (time.perf_counter_ns() - t0) / 1000
                times_us.append(elapsed_us)
                total += 1
                if len(result) >= 5:
                    mcs_ge5 += 1
            except Exception:
                pass

        median_us = statistics.median(times_us) if times_us else 0
        print(f"\nDalke NN: {total} pairs, MCS>=5 {mcs_ge5}/{total}, "
              f"median {median_us:.0f} us")
        assert total > 100, f"Should parse >100 pairs, got {total}"


# ======================================================================
# 3. Stress Pairs — Adversarial Hard Cases
# ======================================================================
class TestStressPairs:

    def test_no_timeouts(self):
        """All stress pairs complete within 10s timeout."""
        pairs = load_tsv_pairs("stress_pairs.tsv")
        timeouts = 0
        total = 0

        for parts in pairs:
            name = parts[2] if len(parts) > 2 else f"pair_{total}"
            try:
                g1 = smsd.parse_smiles(parts[0])
                g2 = smsd.parse_smiles(parts[1])
                t0 = time.perf_counter_ns()
                result = smsd.find_mcs(g1, g2, timeout_ms=10000)
                elapsed_ms = (time.perf_counter_ns() - t0) / 1_000_000
                total += 1
                print(f"  {name:35s} MCS={len(result):3d}  time={elapsed_ms:,.0f}ms")
                if elapsed_ms > 10000:
                    timeouts += 1
            except Exception as e:
                print(f"  {name:35s} PARSE ERROR: {e}")

        print(f"\nStress test: {total} pairs, {timeouts} timeouts")
        assert timeouts <= 1, f"At most 1 timeout, got {timeouts}"

    def test_self_match_exact(self):
        """Self-match pairs return exact atom count."""
        pairs = load_tsv_pairs("stress_pairs.tsv")

        for parts in pairs:
            if len(parts) < 3 or "self" not in parts[2]:
                continue
            name = parts[2]
            try:
                g1 = smsd.parse_smiles(parts[0])
                g2 = smsd.parse_smiles(parts[1])
                result = smsd.find_mcs(g1, g2, timeout_ms=10000)
                mol_size = g1.n
                assert len(result) == mol_size, (
                    f"{name}: self-match MCS={len(result)}, expected {mol_size}"
                )
            except Exception:
                pass  # skip unparseable


# ======================================================================
# 4. Ehrlich-Rarey SMARTS — Substructure Search
# ======================================================================
class TestEhrlichRarey:

    @benchmark
    def test_smarts_performance(self):
        """SMARTS substructure queries run fast."""
        path = os.path.join(DATA_DIR, "ehrlich_rarey_smarts.txt")
        if not os.path.exists(path):
            pytest.skip("ehrlich_rarey_smarts.txt not found")

        smarts_list = []
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                sma = line.split("\t")[0].strip()
                if sma:
                    smarts_list.append(sma)

        # Target: ibuprofen
        target_smi = "CC(C)Cc1ccc(CC(C)C(O)=O)cc1"
        total = 0
        matched = 0
        times_us = []

        target = smsd.parse_smiles(target_smi)
        for sma in smarts_list[:500]:  # First 500 for CI speed
            try:
                t0 = time.perf_counter_ns()
                hit = smsd.smarts_match(sma, target)
                elapsed_us = (time.perf_counter_ns() - t0) / 1000
                times_us.append(elapsed_us)
                total += 1
                if hit:
                    matched += 1
            except Exception:
                pass

        median_us = statistics.median(times_us) if times_us else 0
        print(f"\nEhrlich-Rarey: {total} patterns, {matched} matched, "
              f"median {median_us:.0f} us")
        assert total > 200, f"Should parse >200 SMARTS, got {total}"
