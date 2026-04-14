#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""
Benchmark the present SMSD Python native extension against the lighter native one.

This is an apples-to-apples Python comparison:

* Present SMSD: native `_smsd` extension built from this repository
* Light SMSD: native `_smsd` extension from the sibling `bioinception` repo

Both are benchmarked through the same Python 3.13 driver and the same timing
loop, so we avoid JVM / CLI / wrapper differences and compare the native Python
surface directly.

Coverage:

* MCS on the maintained 20-pair curated set
* Substructure on the curated query/target set

Outputs:

* `benchmarks/results_present_vs_light_scaling_mcs.tsv`
* `benchmarks/results_present_vs_light_scaling_sub.tsv`
* `benchmarks/results_present_vs_light_scaling_summary.txt`
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import statistics
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

MCS_RESULTS_TSV = SCRIPT_DIR / "results_present_vs_light_scaling_mcs.tsv"
SUB_RESULTS_TSV = SCRIPT_DIR / "results_present_vs_light_scaling_sub.tsv"
SUMMARY_TXT = SCRIPT_DIR / "results_present_vs_light_scaling_summary.txt"
SUBSTRUCTURE_PAIRS_FILE = SCRIPT_DIR / "substructure_pairs.tsv"

DEFAULT_TIMEOUT_MS = 10_000
DEFAULT_WARMUP = 3
DEFAULT_ITERS = 10

SIZE_BUCKETS: Sequence[Tuple[int, str]] = (
    (16, "<=16"),
    (32, "17-32"),
    (64, "33-64"),
    (128, "65-128"),
    (10**9, "129+"),
)

RUNNER = r"""
import importlib.util
import json
import statistics
import sys
import time
from pathlib import Path

payload = json.loads(Path(sys.argv[1]).read_text())
ext_path = Path(payload["extension_path"])
spec = importlib.util.spec_from_file_location("_smsd", ext_path)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Could not load extension from {ext_path}")
smsd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(smsd)

chem = smsd.ChemOptions()
mcs_opts = smsd.MCSOptions()
mcs_opts.timeout_ms = int(payload["timeout_ms"])
warmup = int(payload["warmup"])
iters = int(payload["iters"])
mode = payload["mode"]

rows = []
for case in payload["cases"]:
    row = {
        "label": case["label"],
        "category": case["category"],
        "query": case["query"],
        "target": case["target"],
        "error": "",
    }
    try:
        query = smsd.parse_smiles(case["query"])
        target = smsd.parse_smiles(case["target"])
        row["query_atoms"] = len(query)
        row["target_atoms"] = len(target)
        row["combined_atoms"] = len(query) + len(target)

        if hasattr(query, "prewarm"):
            query.prewarm()
        if hasattr(target, "prewarm"):
            target.prewarm()

        times_ms = []
        if mode == "mcs":
            for _ in range(warmup):
                smsd.find_mcs(query, target, chem, mcs_opts)
            result_value = 0
            for _ in range(iters):
                t0 = time.perf_counter_ns()
                mapping = smsd.find_mcs(query, target, chem, mcs_opts)
                times_ms.append((time.perf_counter_ns() - t0) / 1_000_000.0)
                result_value = len(mapping)
            row["result"] = int(result_value)
        else:
            for _ in range(warmup):
                smsd.is_substructure(query, target, chem)
            result_value = False
            for _ in range(iters):
                t0 = time.perf_counter_ns()
                result_value = bool(smsd.is_substructure(query, target, chem))
                times_ms.append((time.perf_counter_ns() - t0) / 1_000_000.0)
            row["result"] = bool(result_value)

        row["avg_ms"] = sum(times_ms) / len(times_ms)
        row["median_ms"] = statistics.median(times_ms)
        row["min_ms"] = min(times_ms)
    except Exception as exc:
        row["error"] = repr(exc)
    rows.append(row)

print(json.dumps({"rows": rows}, separators=(",", ":")))
"""

MCS_PAIRS: List[Tuple[str, str, str, str]] = [
    ("C", "CC", "methane-ethane", "Trivial"),
    ("c1ccccc1", "Cc1ccccc1", "benzene-toluene", "Small aromatic"),
    ("c1ccccc1", "Oc1ccccc1", "benzene-phenol", "Heteroatom"),
    ("CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1", "aspirin-acetaminophen", "Drug pair"),
    ("Cn1cnc2c1c(=O)n(C)c(=O)n2C", "Cn1cnc2c1c(=O)[nH]c(=O)n2C", "caffeine-theophylline", "N-methyl diff"),
    ("CN1CCC23C4C1CC5=C(C2C(C=C4)O3)C=C(C=C5)O",
     "CN1CCC23C4C1CC5=C(C2C(C=C4)OC3)C=C(C=C5)O", "morphine-codeine", "Alkaloid"),
    ("CC(C)Cc1ccc(CC(C)C(=O)O)cc1", "COc1ccc2cc(CC(C)C(=O)O)ccc2c1", "ibuprofen-naproxen", "NSAID"),
    ("C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)O)O)O)N",
     "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
     "ATP-ADP", "Nucleotide"),
    ("C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)([O-])OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N",
     "C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C4N=CN=C5N)O)O)O)O",
     "NAD-NADH", "Cofactor"),
    ("CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
     "CC(C)C1=NC(=NC(=C1C=CC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)N(C)S(=O)(=O)C",
     "atorvastatin-rosuvastatin", "Statin"),
    ("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
     "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O",
     "paclitaxel-docetaxel", "Taxane"),
    ("CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O",
     "CCC1C(C(C(N(CC(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)C)O)(C)O",
     "erythromycin-azithromycin", "Macrolide"),
    ("C1CN2CC3=CCOC4CC(=O)N5C6C4C3CC2C61C7=CC=CC=C75",
     "COC1=CC2=C(C=CN=C2C=C1)C(C3CC4CCN3CC4C=C)O", "strychnine-quinine", "Alkaloid scaffold"),
    ("CC1C(C(CC(O1)OC2C(C(C(OC2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)C(C(C(=O)"
     "NC(C(=O)NC5C(=O)NC7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)C(NC(=O)C("
     "C(C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)C(CC(C)C)NC)"
     "O)Cl)CO)O)O)(C)N)O",
     "CC1C(C(CC(O1)OC2C(C(C(OC2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)C(C(C(=O)"
     "NC(C(=O)NC5C(=O)NC7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)C(NC(=O)C("
     "C(C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)C(CC(C)C)NC)"
     "O)Cl)CO)O)O)(C)N)O", "vancomycin-self", "Self-match large"),
    ("C1C2CC3CC1CC(C2)C3", "C1C2CC3CC1CC(C2)C3", "adamantane-self", "Symmetric"),
    ("C12C3C4C1C5C4C3C25", "C12C3C4C1C5C4C3C25", "cubane-self", "Cage"),
    ("OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO",
     "OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO", "PEG12-PEG16", "Polymer"),
    ("c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67", "c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67", "coronene-self", "PAH"),
    ("O=c1[nH]c(N)nc2[nH]cnc12", "Oc1nc(N)nc2[nH]cnc12", "guanine-keto-enol", "Tautomer"),
    ("c1cc(c(c(c1)Cl)N2c3cc(cc(c3CNC2=O)c4ccc(cc4F)F)N5CCNCC5)Cl",
     "CCNc1cc(c2c(c1)N(C(=O)NC2)c3ccc(cc3)n4ccc-5ncnc5c4)c6ccnnc6", "rdkit-1585-pair", "Known failure"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--present-ext", type=Path, default=None,
                        help="Path to the present repo native _smsd extension.")
    parser.add_argument("--light-ext", type=Path, default=None,
                        help="Path to the lighter repo native _smsd extension.")
    parser.add_argument("--light-repo", type=Path, default=None,
                        help="Optional lighter repo root used for auto-detecting --light-ext.")
    parser.add_argument("--python", default=None,
                        help="Python interpreter compatible with both extensions (auto-detected by default).")
    parser.add_argument("--timeout-ms", type=int, default=DEFAULT_TIMEOUT_MS,
                        help="Per-call MCS timeout used for both native extensions.")
    parser.add_argument("--warmup", type=int, default=DEFAULT_WARMUP,
                        help="Warmup runs for both native extensions.")
    parser.add_argument("--iters", type=int, default=DEFAULT_ITERS,
                        help="Measured runs for both native extensions.")
    parser.add_argument("--mcs-limit", type=int, default=0,
                        help="Only run the first N MCS cases (0 = all).")
    parser.add_argument("--sub-limit", type=int, default=0,
                        help="Only run the first N substructure cases (0 = all).")
    parser.add_argument("--mcs-output", type=Path, default=MCS_RESULTS_TSV,
                        help="TSV output path for MCS results.")
    parser.add_argument("--sub-output", type=Path, default=SUB_RESULTS_TSV,
                        help="TSV output path for substructure results.")
    parser.add_argument("--summary-output", type=Path, default=SUMMARY_TXT,
                        help="Summary output path.")
    return parser.parse_args()


def read_substructure_pairs(path: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    category = ""
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            if line.startswith("# SECTION 1:"):
                category = "Trivial / chain"
            elif line.startswith("# SECTION 2:"):
                category = "Ring containment"
            elif line.startswith("# SECTION 3:"):
                category = "Aromatic queries"
            elif line.startswith("# SECTION 4:"):
                category = "Heterocyclic fragments"
            elif line.startswith("# SECTION 5:"):
                category = "Drug fragment searches"
            elif line.startswith("# SECTION 6:"):
                category = "Edge cases"
            elif line.startswith("# SECTION 7:"):
                category = "Negative controls"
            continue
        query, target, label = line.split("\t", 2)
        rows.append({
            "query": query.strip(),
            "target": target.strip(),
            "label": label.strip(),
            "category": category,
        })
    return rows


def build_mcs_cases(limit: int) -> List[Dict[str, str]]:
    rows = [
        {"query": q, "target": t, "label": label, "category": category}
        for q, t, label, category in MCS_PAIRS
    ]
    return rows[:limit] if limit > 0 else rows


def build_sub_cases(limit: int) -> List[Dict[str, str]]:
    rows = read_substructure_pairs(SUBSTRUCTURE_PAIRS_FILE)
    return rows[:limit] if limit > 0 else rows


def extension_versions(path: Path) -> List[str]:
    match = re.search(r"cpython-(\d{2,3})", path.name)
    if not match:
        return []
    digits = match.group(1)
    if len(digits) == 2:
        major = digits[0]
        minor = digits[1]
    else:
        major = digits[0]
        minor = digits[1:]
    return [f"{major}.{int(minor)}"]


def find_present_extension(explicit: Optional[Path]) -> Path:
    candidates: List[Path] = []
    if explicit is not None:
        candidates.append(explicit)
    candidates.extend([
        PROJECT_DIR / "cpp" / "build-py313",
        PROJECT_DIR / "cpp" / "build",
    ])
    for base in candidates:
        if base.is_file():
            if base.exists():
                return base.resolve()
            continue
        if base.exists():
            hits = sorted(base.glob("_smsd.cpython-*.*"))
            if hits:
                return hits[0].resolve()

    raise FileNotFoundError(
        "Could not find the present repo native extension. Build it first or pass --present-ext."
    )


def find_light_extension(explicit: Optional[Path], light_repo: Optional[Path]) -> Path:
    candidates: List[Path] = []
    if explicit is not None:
        candidates.append(explicit)
    roots: List[Path] = []
    if light_repo is not None:
        roots.append(light_repo)
    roots.append(PROJECT_DIR.parent / "bioinception")

    for root in roots:
        candidates.extend([
            root / "python" / "smsd",
            root / "build-local",
            root / "build" / "smsd_py",
        ])

    for base in candidates:
        if base.is_file():
            if base.exists():
                return base.resolve()
            continue
        if base.exists():
            hits = sorted(base.glob("_smsd.cpython-*.*"))
            if hits:
                return hits[0].resolve()

    raise FileNotFoundError(
        "Could not find the lighter repo native extension. Pass --light-ext or --light-repo."
    )


def detect_python_for_extensions(present_ext: Path, light_ext: Path) -> str:
    present_versions = set(extension_versions(present_ext))
    light_versions = set(extension_versions(light_ext))
    shared = sorted(present_versions & light_versions, reverse=True)
    probe_names: List[str] = []
    for version in shared:
        probe_names.append(f"python{version}")
    probe_names.extend(["python3.13", "python3.12", "python3.11", "python3.10", "python3.9"])

    seen: set[str] = set()
    for name in probe_names:
        if name in seen:
            continue
        seen.add(name)
        resolved = shutil.which(name)
        if resolved:
            return resolved

    raise FileNotFoundError(
        "Could not find a Python interpreter compatible with both native extensions. "
        "Pass --python /path/to/python3.x."
    )


def run_native_benchmark(
    python_exe: str,
    extension_path: Path,
    cases: List[Dict[str, str]],
    mode: str,
    timeout_ms: int,
    warmup: int,
    iters: int,
) -> Dict[str, Dict[str, object]]:
    payload = {
        "extension_path": str(extension_path),
        "mode": mode,
        "timeout_ms": timeout_ms,
        "warmup": warmup,
        "iters": iters,
        "cases": cases,
    }
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as handle:
        json.dump(payload, handle)
        payload_path = Path(handle.name)
    try:
        proc = subprocess.run(
            [python_exe, "-c", RUNNER, str(payload_path)],
            capture_output=True,
            text=True,
            timeout=max(60, len(cases) * 20),
        )
    finally:
        payload_path.unlink(missing_ok=True)

    if proc.returncode != 0:
        raise RuntimeError(
            "Native benchmark failed.\n"
            f"extension={extension_path}\n"
            f"stdout:\n{proc.stdout}\n\nstderr:\n{proc.stderr}"
        )

    data = json.loads(proc.stdout)
    return {row["label"]: row for row in data["rows"]}


def atom_bucket(total_atoms: int) -> str:
    for upper, label in SIZE_BUCKETS:
        if total_atoms <= upper:
            return label
    return SIZE_BUCKETS[-1][1]


def speedup_ratio(present_ms: float, light_ms: float) -> float:
    if present_ms <= 0 or light_ms <= 0:
        return -1.0
    return present_ms / light_ms


def fmt_ms(value: float) -> str:
    if value < 0:
        return "n/a"
    if value < 1.0:
        return f"{value:.3f}"
    if value < 100.0:
        return f"{value:.2f}"
    return f"{value:.1f}"


def fmt_speedup(value: float) -> str:
    if value < 0:
        return "n/a"
    return f"{value:.1f}x"


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Sequence[Dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def summarize_rows(
    name: str,
    rows: Sequence[Dict[str, object]],
    present_key: str,
    light_key: str,
) -> List[str]:
    lines = [name]
    comparable = [
        row for row in rows
        if float(row[present_key]) > 0 and float(row[light_key]) > 0
    ]
    agrees = sum(1 for row in rows if row["agree"] is True)
    mismatches = [str(row["pair"]) for row in rows if row["agree"] is False]

    lines.append(f"  cases={len(rows)} comparable={len(comparable)} agreement={agrees}/{len(rows)}")
    if comparable:
        present_vals = [float(row[present_key]) for row in comparable]
        light_vals = [float(row[light_key]) for row in comparable]
        speedups = [float(row["speedup"]) for row in comparable if float(row["speedup"]) > 0]
        lines.append(
            "  median present={}ms light={}ms speedup={}".format(
                fmt_ms(statistics.median(present_vals)),
                fmt_ms(statistics.median(light_vals)),
                fmt_speedup(statistics.median(speedups)) if speedups else "n/a",
            )
        )
    if mismatches:
        lines.append("  mismatches=" + ", ".join(mismatches))

    lines.append("  buckets by combined atom count:")
    for _, bucket_label in SIZE_BUCKETS:
        bucket_rows = [row for row in comparable if row["size_bucket"] == bucket_label]
        if not bucket_rows:
            continue
        present_vals = [float(row[present_key]) for row in bucket_rows]
        light_vals = [float(row[light_key]) for row in bucket_rows]
        speedups = [float(row["speedup"]) for row in bucket_rows if float(row["speedup"]) > 0]
        lines.append(
            "    {} cases={} present={}ms light={}ms speedup={}".format(
                bucket_label,
                len(bucket_rows),
                fmt_ms(statistics.median(present_vals)),
                fmt_ms(statistics.median(light_vals)),
                fmt_speedup(statistics.median(speedups)) if speedups else "n/a",
            )
        )

    slowest = sorted(comparable, key=lambda row: float(row[present_key]), reverse=True)[:5]
    if slowest:
        lines.append("  slowest present-extension cases:")
        for row in slowest:
            lines.append(
                "    {} [{} atoms] present={}ms light={}ms speedup={}".format(
                    row["pair"],
                    row["combined_atoms"],
                    fmt_ms(float(row[present_key])),
                    fmt_ms(float(row[light_key])),
                    fmt_speedup(float(row["speedup"])),
                )
            )
    return lines


def main() -> int:
    args = parse_args()

    present_ext = find_present_extension(args.present_ext.resolve() if args.present_ext else None)
    light_repo = args.light_repo.resolve() if args.light_repo else None
    light_ext = find_light_extension(args.light_ext.resolve() if args.light_ext else None, light_repo)
    python_exe = args.python or detect_python_for_extensions(present_ext, light_ext)

    mcs_cases = build_mcs_cases(args.mcs_limit)
    sub_cases = build_sub_cases(args.sub_limit)

    present_mcs = run_native_benchmark(
        python_exe=python_exe,
        extension_path=present_ext,
        cases=mcs_cases,
        mode="mcs",
        timeout_ms=args.timeout_ms,
        warmup=args.warmup,
        iters=args.iters,
    )
    light_mcs = run_native_benchmark(
        python_exe=python_exe,
        extension_path=light_ext,
        cases=mcs_cases,
        mode="mcs",
        timeout_ms=args.timeout_ms,
        warmup=args.warmup,
        iters=args.iters,
    )
    present_sub = run_native_benchmark(
        python_exe=python_exe,
        extension_path=present_ext,
        cases=sub_cases,
        mode="sub",
        timeout_ms=args.timeout_ms,
        warmup=args.warmup,
        iters=args.iters,
    )
    light_sub = run_native_benchmark(
        python_exe=python_exe,
        extension_path=light_ext,
        cases=sub_cases,
        mode="sub",
        timeout_ms=args.timeout_ms,
        warmup=args.warmup,
        iters=args.iters,
    )

    mcs_rows: List[Dict[str, object]] = []
    for case in mcs_cases:
        present = present_mcs[case["label"]]
        light = light_mcs[case["label"]]
        combined_atoms = int(present.get("combined_atoms", light.get("combined_atoms", -1)))
        present_result = int(present["result"]) if not present["error"] else -1
        light_result = int(light["result"]) if not light["error"] else -1
        speedup = speedup_ratio(float(present.get("avg_ms", -1.0)), float(light.get("avg_ms", -1.0)))
        mcs_rows.append({
            "pair": case["label"],
            "category": case["category"],
            "query_atoms": int(present.get("query_atoms", light.get("query_atoms", -1))),
            "target_atoms": int(present.get("target_atoms", light.get("target_atoms", -1))),
            "combined_atoms": combined_atoms,
            "size_bucket": atom_bucket(combined_atoms),
            "present_avg_ms": float(present.get("avg_ms", -1.0)),
            "present_median_ms": float(present.get("median_ms", -1.0)),
            "present_min_ms": float(present.get("min_ms", -1.0)),
            "present_mcs": present_result,
            "light_avg_ms": float(light.get("avg_ms", -1.0)),
            "light_median_ms": float(light.get("median_ms", -1.0)),
            "light_min_ms": float(light.get("min_ms", -1.0)),
            "light_mcs": light_result,
            "speedup": speedup,
            "agree": (present_result == light_result) if not present["error"] and not light["error"] else False,
            "present_error": str(present["error"]),
            "light_error": str(light["error"]),
        })

    sub_rows: List[Dict[str, object]] = []
    for case in sub_cases:
        present = present_sub[case["label"]]
        light = light_sub[case["label"]]
        combined_atoms = int(present.get("combined_atoms", light.get("combined_atoms", -1)))
        present_result = bool(present["result"]) if not present["error"] else False
        light_result = bool(light["result"]) if not light["error"] else False
        speedup = speedup_ratio(float(present.get("avg_ms", -1.0)), float(light.get("avg_ms", -1.0)))
        sub_rows.append({
            "pair": case["label"],
            "category": case["category"],
            "query_atoms": int(present.get("query_atoms", light.get("query_atoms", -1))),
            "target_atoms": int(present.get("target_atoms", light.get("target_atoms", -1))),
            "combined_atoms": combined_atoms,
            "size_bucket": atom_bucket(combined_atoms),
            "present_avg_ms": float(present.get("avg_ms", -1.0)),
            "present_median_ms": float(present.get("median_ms", -1.0)),
            "present_min_ms": float(present.get("min_ms", -1.0)),
            "present_hit": present_result,
            "light_avg_ms": float(light.get("avg_ms", -1.0)),
            "light_median_ms": float(light.get("median_ms", -1.0)),
            "light_min_ms": float(light.get("min_ms", -1.0)),
            "light_hit": light_result,
            "speedup": speedup,
            "agree": (present_result == light_result) if not present["error"] and not light["error"] else False,
            "present_error": str(present["error"]),
            "light_error": str(light["error"]),
        })

    write_tsv(
        args.mcs_output,
        [
            "pair",
            "category",
            "query_atoms",
            "target_atoms",
            "combined_atoms",
            "size_bucket",
            "present_avg_ms",
            "present_median_ms",
            "present_min_ms",
            "present_mcs",
            "light_avg_ms",
            "light_median_ms",
            "light_min_ms",
            "light_mcs",
            "speedup",
            "agree",
            "present_error",
            "light_error",
        ],
        mcs_rows,
    )
    write_tsv(
        args.sub_output,
        [
            "pair",
            "category",
            "query_atoms",
            "target_atoms",
            "combined_atoms",
            "size_bucket",
            "present_avg_ms",
            "present_median_ms",
            "present_min_ms",
            "present_hit",
            "light_avg_ms",
            "light_median_ms",
            "light_min_ms",
            "light_hit",
            "speedup",
            "agree",
            "present_error",
            "light_error",
        ],
        sub_rows,
    )

    summary_lines = [
        "Present SMSD vs Light SMSD Scaling Benchmark",
        "============================================",
        f"Present ext:   {present_ext}",
        f"Light ext:     {light_ext}",
        f"Python:        {python_exe}",
        f"Timeout ms:    {args.timeout_ms}",
        f"Warmup:        {args.warmup}",
        f"Iters:         {args.iters}",
        f"MCS cases:     {len(mcs_rows)}",
        f"Sub cases:     {len(sub_rows)}",
        "",
    ]
    summary_lines.extend(summarize_rows("MCS", mcs_rows, "present_avg_ms", "light_avg_ms"))
    summary_lines.append("")
    summary_lines.extend(summarize_rows("Substructure", sub_rows, "present_avg_ms", "light_avg_ms"))

    summary_text = "\n".join(summary_lines) + "\n"
    args.summary_output.write_text(summary_text)

    print(summary_text, end="")
    print(f"Wrote: {args.mcs_output}")
    print(f"Wrote: {args.sub_output}")
    print(f"Wrote: {args.summary_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
