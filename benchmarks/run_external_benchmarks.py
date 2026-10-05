#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
"""
Run all external community-standard benchmarks — SMSD Pro vs RDKit baseline.

Both engines are called in the same Python process with identical SMILES inputs
and the same timeout. Results are written to named files per dataset.

Usage:
    python benchmarks/run_external_benchmarks.py

Output files (one per dataset):
    results_stress_12pairs.txt          — 12 adversarial hard cases
    results_tautobase_468pairs.txt      — 468 tautomer pairs (Wahl & Sander 2020)
    results_dalke_random_1000pairs.txt  — 1,000 low-similarity MCS pairs (Dalke 2013)
    results_dalke_nn_1000pairs.txt      — 1,000 nearest-neighbor MCS pairs (Dalke 2013)
    results_ehrlich_rarey_1400smarts.txt — 1,400 SMARTS patterns (Ehrlich & Rarey 2011)
    results_external_all_combined.txt   — combined summary of all five above
"""
import os, sys, time, statistics

try:
    import smsd
except ImportError:
    print("ERROR: smsd not installed. pip install smsd")
    sys.exit(1)

try:
    from rdkit import Chem
    from rdkit.Chem import rdFMCS
    RDKIT_AVAILABLE = True
    RDKIT_VERSION = Chem.rdBase.rdkitVersion
except ImportError:
    RDKIT_AVAILABLE = False
    RDKIT_VERSION = "not installed"
    print("WARNING: rdkit not installed — RDKit baseline will be skipped.")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(SCRIPT_DIR, "data")
TIMEOUT_MS  = 10_000
TIMEOUT_SEC = TIMEOUT_MS / 1000.0

combined_lines = []


def log(msg="", sink=None):
    print(msg)
    combined_lines.append(msg)
    if sink is not None:
        sink.append(msg)


def write_result_file(filename, lines):
    path = os.path.join(SCRIPT_DIR, filename)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"  → Written: {filename}")


def rdkit_mcs(smi1, smi2):
    """Run RDKit FindMCS. Returns (mcs_atoms, elapsed_ms, timed_out)."""
    if not RDKIT_AVAILABLE:
        return -1, -1, False
    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)
    if mol1 is None or mol2 is None:
        return -1, -1, False
    t0 = time.perf_counter_ns()
    res = rdFMCS.FindMCS([mol1, mol2], timeout=int(TIMEOUT_SEC))
    elapsed_ms = (time.perf_counter_ns() - t0) / 1e6
    return res.numAtoms, elapsed_ms, res.canceled


def rdkit_smarts(sma, smi_target):
    """Run RDKit SMARTS match. Returns (hit, elapsed_us) or (None, -1) on error."""
    if not RDKIT_AVAILABLE:
        return None, -1
    try:
        pat = Chem.MolFromSmarts(sma)
        mol = Chem.MolFromSmiles(smi_target)
        if pat is None or mol is None:
            return None, -1
        t0 = time.perf_counter_ns()
        hit = mol.HasSubstructMatch(pat)
        elapsed_us = (time.perf_counter_ns() - t0) / 1000
        return hit, elapsed_us
    except Exception:
        return None, -1


# ---------------------------------------------------------------------------
# Benchmark 1 — Stress Pairs (12 adversarial hard cases)
# ---------------------------------------------------------------------------
def run_stress_pairs():
    buf = []
    log("=" * 100, buf)
    log("BENCHMARK 1: Stress Pairs — 12 Adversarial Hard Cases  |  SMSD Pro vs RDKit FindMCS", buf)
    log("Source: stress_pairs.tsv (BioInception curated)", buf)
    log("Cases: cubane/cuneane, coronene, vancomycin, macrocycles, cyclopeptides", buf)
    log("=" * 100, buf)
    log(f"  {'Pair':<35s} {'SMSD MCS':>8s} {'SMSD ms':>10s} {'RDKit MCS':>9s} {'RDKit ms':>10s} {'SMSD wins?':>10s}  Challenge", buf)
    log("-" * 100, buf)

    path = os.path.join(DATA, "stress_pairs.tsv")
    total = smsd_timeouts = rdkit_timeouts = smsd_better = rdkit_better = equal = 0
    smsd_times_ms = []
    rdkit_times_ms = []

    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            smi1, smi2 = parts[0], parts[1]
            name      = parts[2] if len(parts) > 2 else f"pair_{total}"
            challenge = parts[3] if len(parts) > 3 else ""

            # SMSD
            smsd_mcs = smsd_ms = -1
            smsd_to = False
            try:
                g1 = smsd.parse_smiles(smi1)
                g2 = smsd.parse_smiles(smi2)
                t0 = time.perf_counter_ns()
                m = smsd.find_mcs(g1, g2, timeout_ms=TIMEOUT_MS)
                smsd_ms = (time.perf_counter_ns() - t0) / 1e6
                smsd_mcs = len(m)
                smsd_to = smsd_ms > TIMEOUT_MS
                smsd_times_ms.append(smsd_ms)
                total += 1
            except Exception as e:
                log(f"  {name:<35s} SMSD PARSE ERROR: {e}", buf)
                continue

            # RDKit baseline
            r_mcs, r_ms, r_to = rdkit_mcs(smi1, smi2)
            if r_ms > 0:
                rdkit_times_ms.append(r_ms)

            if smsd_to:
                smsd_timeouts += 1
            if r_to:
                rdkit_timeouts += 1

            smsd_ms_str  = "TIMEOUT" if smsd_to else f"{smsd_ms:,.1f}"
            rdkit_ms_str = "TIMEOUT" if r_to   else (f"{r_ms:,.1f}" if r_ms >= 0 else "N/A")
            rdkit_mcs_str = "TO" if r_to else (str(r_mcs) if r_mcs >= 0 else "N/A")

            if smsd_mcs >= 0 and r_mcs >= 0:
                if smsd_mcs > r_mcs:   smsd_better += 1; win = "SMSD"
                elif r_mcs > smsd_mcs: rdkit_better += 1; win = "RDKit"
                else:                  equal += 1;        win = "equal"
            else:
                win = "-"

            log(f"  {name:<35s} {smsd_mcs:>8d} {smsd_ms_str:>10s} {rdkit_mcs_str:>9s} {rdkit_ms_str:>10s} {win:>10s}  {challenge}", buf)

    log("-" * 100, buf)
    log(f"\n  SMSD  : {total} pairs, {smsd_timeouts} timeouts"
        + (f", median {statistics.median(smsd_times_ms):.1f} ms" if smsd_times_ms else ""), buf)
    if RDKIT_AVAILABLE:
        log(f"  RDKit : {total} pairs, {rdkit_timeouts} timeouts"
            + (f", median {statistics.median(rdkit_times_ms):.1f} ms" if rdkit_times_ms else ""), buf)
        log(f"  MCS quality: SMSD better={smsd_better}, RDKit better={rdkit_better}, equal={equal}", buf)
    log("", buf)
    write_result_file("results_stress_12pairs.txt", buf)


# ---------------------------------------------------------------------------
# Benchmark 2 — Tautobase (468 tautomer pairs)
# ---------------------------------------------------------------------------
def run_tautobase():
    buf = []
    log("=" * 100, buf)
    log("BENCHMARK 2: Tautobase — Tautomer-Aware MCS  |  SMSD Pro vs RDKit FindMCS", buf)
    log("Source: chodera_tautobase_subset.txt (468 pairs)", buf)
    log("Reference: Wahl & Sander, J. Chem. Inf. Model. 2020, 60(3):1085-1089", buf)
    log("Metric: % full heavy-atom match; SMSD tautomer-aware vs SMSD strict vs RDKit", buf)
    log("=" * 100, buf)

    path = os.path.join(DATA, "chodera_tautobase_subset.txt")
    total = 0
    smsd_full = smsd_partial = smsd_nogain = 0
    rdkit_full = 0
    smsd_gain_over_rdkit = 0
    total_gain_atoms = 0
    smsd_times_us = []
    rdkit_times_us = []

    taut_opts = smsd.ChemOptions()
    taut_opts.tautomer_aware = True

    with open(path) as f:
        for line in f:
            if line.startswith("name") or not line.strip():
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 3:
                continue
            name, smi1, smi2 = parts[0], parts[1], parts[2]
            try:
                g1 = smsd.parse_smiles(smi1)
                g2 = smsd.parse_smiles(smi2)
                max_atoms = min(g1.n, g2.n)

                # SMSD tautomer-aware
                t0 = time.perf_counter_ns()
                m_taut = smsd.find_mcs(g1, g2, chem=taut_opts, timeout_ms=TIMEOUT_MS)
                smsd_times_us.append((time.perf_counter_ns() - t0) / 1000)
                taut_mcs = len(m_taut)

                # SMSD strict
                m_strict = smsd.find_mcs(g1, g2, timeout_ms=TIMEOUT_MS)
                strict_mcs = len(m_strict)

                total += 1
                gain = taut_mcs - strict_mcs
                if gain > 0:
                    total_gain_atoms += gain
                if taut_mcs >= max_atoms:
                    smsd_full += 1
                elif taut_mcs > strict_mcs:
                    smsd_partial += 1
                else:
                    smsd_nogain += 1

                # RDKit baseline
                r_mcs, r_ms, r_to = rdkit_mcs(smi1, smi2)
                if r_ms > 0:
                    rdkit_times_us.append(r_ms * 1000)
                if not r_to and r_mcs >= 0:
                    if r_mcs >= max_atoms:
                        rdkit_full += 1
                    if taut_mcs > r_mcs:
                        smsd_gain_over_rdkit += 1

            except Exception:
                pass

    smsd_med = statistics.median(smsd_times_us) if smsd_times_us else 0
    rdkit_med = statistics.median(rdkit_times_us) if rdkit_times_us else 0

    log(f"  Pairs tested:                 {total}", buf)
    log("", buf)
    log(f"  SMSD tautomer-aware MCS:", buf)
    log(f"    Full heavy-atom match:      {smsd_full} ({100*smsd_full/max(total,1):.1f}%)", buf)
    log(f"    Partial gain over strict:   {smsd_partial}", buf)
    log(f"    No gain:                    {smsd_nogain}", buf)
    log(f"    Total extra atoms gained:   +{total_gain_atoms}", buf)
    log(f"    Median time:                {smsd_med:.0f} us", buf)
    if RDKIT_AVAILABLE:
        log("", buf)
        log(f"  RDKit FindMCS (strict, baseline):", buf)
        log(f"    Full heavy-atom match:      {rdkit_full} ({100*rdkit_full/max(total,1):.1f}%)", buf)
        log(f"    Median time:                {rdkit_med:.0f} us", buf)
        log("", buf)
        log(f"  SMSD tautomer-aware > RDKit strict: {smsd_gain_over_rdkit} pairs ({100*smsd_gain_over_rdkit/max(total,1):.1f}%)", buf)
    log("", buf)
    write_result_file("results_tautobase_468pairs.txt", buf)


# ---------------------------------------------------------------------------
# Benchmark 3 — Dalke Random (1,000 low-similarity pairs)
# ---------------------------------------------------------------------------
def run_dalke_random():
    buf = []
    log("=" * 100, buf)
    log("BENCHMARK 3: Dalke Random Pairs — Low-Similarity MCS  |  SMSD Pro vs RDKit FindMCS", buf)
    log("Source: dalke_random_pairs.tsv (1,000 pairs from ChEMBL drug-like set)", buf)
    log("Reference: Dalke & Hastings, J. Cheminform. 2013, 5(Suppl 1):O6", buf)
    log("Metric: median MCS time (us), timeout count, mean MCS size", buf)
    log("=" * 100, buf)

    _run_dalke_file(
        "dalke_random_pairs.tsv",
        "results_dalke_random_1000pairs.txt",
        buf,
    )


# ---------------------------------------------------------------------------
# Benchmark 4 — Dalke Nearest-Neighbor (1,000 high-similarity pairs)
# ---------------------------------------------------------------------------
def run_dalke_nn():
    buf = []
    log("=" * 100, buf)
    log("BENCHMARK 4: Dalke Nearest-Neighbor Pairs — High-Similarity MCS  |  SMSD Pro vs RDKit FindMCS", buf)
    log("Source: dalke_nn_pairs.tsv (1,000 pairs, Tanimoto >= 0.7)", buf)
    log("Reference: Dalke & Hastings, J. Cheminform. 2013, 5(Suppl 1):O6", buf)
    log("Metric: median MCS time (us), timeout count, MCS >= 5 atoms rate", buf)
    log("=" * 100, buf)

    _run_dalke_file(
        "dalke_nn_pairs.tsv",
        "results_dalke_nn_1000pairs.txt",
        buf,
    )


def _run_dalke_file(filename, outfile, buf):
    """Shared logic for Dalke random and NN benchmarks."""
    path = os.path.join(DATA, filename)
    total = smsd_tos = rdkit_tos = 0
    smsd_times_us = []
    rdkit_times_us = []
    smsd_sizes = []
    rdkit_sizes = []
    smsd_better = rdkit_better = equal = 0

    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            smi1, smi2 = parts[0], parts[1]

            # SMSD
            try:
                g1 = smsd.parse_smiles(smi1)
                g2 = smsd.parse_smiles(smi2)
                t0 = time.perf_counter_ns()
                m = smsd.find_mcs(g1, g2, timeout_ms=TIMEOUT_MS)
                elapsed_us = (time.perf_counter_ns() - t0) / 1000
                smsd_times_us.append(elapsed_us)
                smsd_sizes.append(len(m))
                total += 1
                if elapsed_us > TIMEOUT_MS * 1000:
                    smsd_tos += 1
            except Exception:
                continue

            # RDKit
            r_mcs, r_ms, r_to = rdkit_mcs(smi1, smi2)
            if r_ms > 0:
                rdkit_times_us.append(r_ms * 1000)
            if not r_to and r_mcs >= 0:
                rdkit_sizes.append(r_mcs)
            if r_to:
                rdkit_tos += 1

            # Quality comparison
            if smsd_sizes and r_mcs >= 0:
                s = smsd_sizes[-1]
                if s > r_mcs:        smsd_better += 1
                elif r_mcs > s:      rdkit_better += 1
                else:                equal += 1

    s_med  = statistics.median(smsd_times_us)  if smsd_times_us  else 0
    s_mean = statistics.mean(smsd_sizes)        if smsd_sizes     else 0
    r_med  = statistics.median(rdkit_times_us)  if rdkit_times_us else 0
    r_mean = statistics.mean(rdkit_sizes)       if rdkit_sizes    else 0

    log(f"  Pairs tested:           {total}", buf)
    log("", buf)
    log(f"  SMSD Pro:", buf)
    log(f"    Timeouts:             {smsd_tos}", buf)
    log(f"    Median time:          {s_med:.0f} us", buf)
    log(f"    Mean MCS size:        {s_mean:.1f} atoms", buf)
    if RDKIT_AVAILABLE:
        log("", buf)
        log(f"  RDKit FindMCS (baseline):", buf)
        log(f"    Timeouts:             {rdkit_tos}", buf)
        log(f"    Median time:          {r_med:.0f} us", buf)
        log(f"    Mean MCS size:        {r_mean:.1f} atoms", buf)
        log("", buf)
        speedup = r_med / s_med if s_med > 0 else 0
        log(f"  SMSD speedup (median): {speedup:.1f}x", buf)
        log(f"  MCS quality: SMSD better={smsd_better}, RDKit better={rdkit_better}, equal={equal}", buf)
    log("", buf)
    write_result_file(outfile, buf)


# ---------------------------------------------------------------------------
# Benchmark 5 — Ehrlich-Rarey SMARTS (1,400 substructure patterns)
# ---------------------------------------------------------------------------
def run_ehrlich_rarey():
    buf = []
    log("=" * 100, buf)
    log("BENCHMARK 5: Ehrlich-Rarey SMARTS v2.0 — Substructure Search  |  SMSD Pro vs RDKit", buf)
    log("Source: ehrlich_rarey_smarts.txt (1,400 patterns, ZBH Hamburg)", buf)
    log("Reference: Ehrlich & Rarey, J. Chem. Inf. Model. 2011, 51(6):1316-1324", buf)
    log("Target: ibuprofen (CC(C)Cc1ccc(CC(C)C(O)=O)cc1)", buf)
    log("Metric: median query time (us), hit rate, hit agreement", buf)
    log("=" * 100, buf)

    path = os.path.join(DATA, "ehrlich_rarey_smarts.txt")
    smarts_list = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            sma = line.strip().split("\t")[0].strip()
            if sma:
                smarts_list.append(sma)

    TARGET_SMILES = "CC(C)Cc1ccc(CC(C)C(O)=O)cc1"
    smsd_target  = smsd.parse_smiles(TARGET_SMILES)
    rdkit_target = Chem.MolFromSmiles(TARGET_SMILES) if RDKIT_AVAILABLE else None

    smsd_total = smsd_matched = smsd_errors = 0
    rdkit_total = rdkit_matched = rdkit_errors = 0
    agree = disagree = 0
    smsd_times_us = []
    rdkit_times_us = []

    for sma in smarts_list:
        # SMSD
        s_hit = None
        try:
            t0 = time.perf_counter_ns()
            s_hit = smsd.smarts_match(sma, smsd_target)
            smsd_times_us.append((time.perf_counter_ns() - t0) / 1000)
            smsd_total += 1
            if s_hit:
                smsd_matched += 1
        except Exception:
            smsd_errors += 1

        # RDKit baseline
        r_hit, r_us = rdkit_smarts(sma, TARGET_SMILES)
        if r_hit is not None:
            rdkit_times_us.append(r_us)
            rdkit_total += 1
            if r_hit:
                rdkit_matched += 1
        else:
            rdkit_errors += 1

        # Agreement
        if s_hit is not None and r_hit is not None:
            if s_hit == r_hit:
                agree += 1
            else:
                disagree += 1

    s_med = statistics.median(smsd_times_us)  if smsd_times_us  else 0
    r_med = statistics.median(rdkit_times_us) if rdkit_times_us else 0

    log(f"  SMSD Pro:", buf)
    log(f"    Patterns tested: {smsd_total}", buf)
    log(f"    Matched:         {smsd_matched} ({100*smsd_matched/max(smsd_total,1):.1f}%)", buf)
    log(f"    Unsupported:     {smsd_errors}", buf)
    log(f"    Median time:     {s_med:.0f} us", buf)
    if RDKIT_AVAILABLE:
        log("", buf)
        log(f"  RDKit HasSubstructMatch (baseline):", buf)
        log(f"    Patterns tested: {rdkit_total}", buf)
        log(f"    Matched:         {rdkit_matched} ({100*rdkit_matched/max(rdkit_total,1):.1f}%)", buf)
        log(f"    Unsupported:     {rdkit_errors}", buf)
        log(f"    Median time:     {r_med:.0f} us", buf)
        log("", buf)
        speedup = r_med / s_med if s_med > 0 else 0
        log(f"  SMSD speedup (median):    {speedup:.1f}x", buf)
        log(f"  Hit agreement:            {agree}/{agree+disagree} patterns", buf)
        if disagree > 0:
            log(f"  Disagreements:            {disagree} (check SMARTS semantics)", buf)
    log("", buf)
    write_result_file("results_ehrlich_rarey_1400smarts.txt", buf)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    header = [
        "=" * 100,
        "SMSD Pro External Benchmark Suite — SMSD vs RDKit Baseline",
        f"SMSD version : {smsd.__version__}",
        f"RDKit version: {RDKIT_VERSION}",
        f"Date         : {time.strftime('%Y-%m-%d %H:%M')}",
        f"Timeout      : {TIMEOUT_MS} ms per pair",
        "=" * 100,
        "",
        "Five benchmarks against community-standard public datasets:",
        "  1. Stress pairs        — 12 adversarial hard cases       (SMSD vs RDKit MCS)",
        "  2. Tautobase           — 468 tautomer pairs               (SMSD tautomer-aware vs RDKit strict)",
        "  3. Dalke random        — 1,000 low-similarity MCS pairs   (SMSD vs RDKit MCS)",
        "  4. Dalke NN            — 1,000 nearest-neighbor MCS pairs (SMSD vs RDKit MCS)",
        "  5. Ehrlich-Rarey       — 1,400 SMARTS patterns            (SMSD vs RDKit substructure)",
        "",
    ]
    for line in header:
        log(line)

    run_stress_pairs()
    run_tautobase()
    run_dalke_random()
    run_dalke_nn()
    run_ehrlich_rarey()

    write_result_file("results_external_all_combined.txt", combined_lines)
    print("\nAll benchmark result files written to benchmarks/")
