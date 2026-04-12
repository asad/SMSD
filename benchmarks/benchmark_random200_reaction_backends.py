#!/usr/bin/env python3
"""
Benchmark random sampled reaction mapping runs across backend combinations.

This script samples reactions from the Golden benchmark RDF and compares the
pure-Python mapper with different MCS/substructure backend combinations.

Outputs:
  - JSON summary
  - TSV per-reaction rows
  - Markdown comparison tables
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import random
import statistics
import sys
import time
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


SCRIPT_DIR = Path(__file__).resolve().parent
SMSD_REPO = SCRIPT_DIR.parent
GITHUB_ROOT = SMSD_REPO.parent
DEFAULT_RDT_MAPPER = GITHUB_ROOT / "bioinception" / "python" / "rdt_port" / "rdkit_mapper.py"
DEFAULT_RDF = GITHUB_ROOT / "bioinception" / "test_resources" / "benchmark" / "golden_dataset.rdf"
DEFAULT_OUTPUT_PREFIX = SCRIPT_DIR / "results_random200_reaction_backend_compare"


def load_mapper_module(path: Path):
    spec = importlib.util.spec_from_file_location("rdt_mapper_benchmark", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load mapper module from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@dataclass
class ComboSummary:
    label: str
    mcs_backend: str
    substructure_backend: str
    sample_size: int
    failures: int
    coverage_pct: float
    elapsed_seconds: float
    reactions_per_second: float
    median_ms: float
    p95_ms: float
    exact_pct: float
    molecule_pct: float
    atom_pct: float
    bond_change_pct: float
    mean_abs_bond_delta: float


@dataclass
class SampleProfile:
    sample_size: int
    issue_free_pct: float
    error_free_pct: float
    warning_free_pct: float
    charge_balanced_pct: float
    non_h_balanced_pct: float
    median_reactant_mols: float
    median_product_mols: float


def percentile(sorted_values: list[float], pct: float) -> float:
    if not sorted_values:
        return 0.0
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = (len(sorted_values) - 1) * pct
    low = int(position)
    high = min(low + 1, len(sorted_values) - 1)
    fraction = position - low
    return sorted_values[low] * (1.0 - fraction) + sorted_values[high] * fraction


def sample_reactions(reactions: list, sample_size: int, seed: int) -> list:
    eligible = [reaction for reaction in reactions if reaction.reactants and reaction.products]
    if sample_size > len(eligible):
        raise ValueError(f"Requested sample_size={sample_size}, but only {len(eligible)} eligible reactions exist")
    return random.Random(seed).sample(eligible, sample_size)


def summarize_sample_profile(mapper_module, reactions: Iterable) -> SampleProfile:
    validations = [mapper_module.validate_reaction(reaction) for reaction in reactions]
    sample = list(reactions)
    issue_free = sum(1 for validation in validations if validation.issue_count == 0)
    error_free = sum(1 for validation in validations if validation.error_count == 0)
    warning_free = sum(1 for validation in validations if validation.warning_count == 0)
    charge_balanced = sum(1 for validation in validations if validation.reactant_total_charge == validation.product_total_charge)
    non_h_balanced = sum(
        1
        for validation in validations
        if validation.reactant_non_hydrogen_element_counts == validation.product_non_hydrogen_element_counts
    )
    reactant_counts = [len(reaction.reactants) for reaction in sample]
    product_counts = [len(reaction.products) for reaction in sample]
    total = max(1, len(sample))
    return SampleProfile(
        sample_size=len(sample),
        issue_free_pct=100.0 * issue_free / total,
        error_free_pct=100.0 * error_free / total,
        warning_free_pct=100.0 * warning_free / total,
        charge_balanced_pct=100.0 * charge_balanced / total,
        non_h_balanced_pct=100.0 * non_h_balanced / total,
        median_reactant_mols=float(statistics.median(reactant_counts)),
        median_product_mols=float(statistics.median(product_counts)),
    )


def run_combo(mapper_module, reactions: list, mcs_backend: str, substructure_backend: str, timeout_sec: float, num_threads: int | None):
    per_reaction_rows: list[dict[str, object]] = []
    elapsed_values_ms: list[float] = []
    failures = 0
    exact_count = 0
    molecule_count = 0
    atom_correct_total = 0
    atom_total_total = 0
    bond_change_exact = 0
    bond_change_deltas: list[int] = []
    algorithms = Counter()
    combo_start = time.perf_counter()

    for reaction in reactions:
        row: dict[str, object] = {
            "reaction_id": str(reaction.rxn_id),
            "mcs_backend": mcs_backend,
            "substructure_backend": substructure_backend,
        }
        start = time.perf_counter()
        try:
            result = mapper_module.map_reaction(
                reaction,
                mcs_timeout_sec=timeout_sec,
                mode="strict",
                mcs_backend=mcs_backend,
                substructure_backend=substructure_backend,
                num_threads=num_threads,
            )
            elapsed_ms = (time.perf_counter() - start) * 1000.0
            elapsed_values_ms.append(elapsed_ms)
            atom_correct, atom_total, exact, molecule_exact = mapper_module.check_accuracy(
                reaction,
                result.reactant_map_numbers,
                result.product_map_numbers,
            )
            if result.bond_change_calculator is not None:
                our_total, our_formed, our_cleaved, our_order = mapper_module._bond_change_totals_from_calculator(
                    result.bond_change_calculator
                )
            else:
                our_total, our_formed, our_cleaved, our_order = mapper_module.compute_bond_changes(
                    reaction,
                    result.reactant_map_numbers,
                    result.product_map_numbers,
                )
            try:
                gold_total, gold_formed, gold_cleaved, gold_order = reaction.gold_bond_changes()
            except Exception:
                gold_total, gold_formed, gold_cleaved, gold_order = -1, -1, -1, -1
            delta_total = our_total - gold_total if our_total >= 0 and gold_total >= 0 else None
            if exact:
                exact_count += 1
            if molecule_exact:
                molecule_count += 1
            atom_correct_total += atom_correct
            atom_total_total += atom_total
            if delta_total == 0:
                bond_change_exact += 1
            if delta_total is not None:
                bond_change_deltas.append(abs(delta_total))
            algorithms[result.algorithm] += 1
            row.update(
                {
                    "algorithm": result.algorithm,
                    "elapsed_ms": round(elapsed_ms, 3),
                    "exact": int(exact),
                    "molecule_ok": int(molecule_exact),
                    "atom_correct": atom_correct,
                    "atom_total": atom_total,
                    "gold_total": gold_total,
                    "our_total": our_total,
                    "delta_total": delta_total if delta_total is not None else "",
                    "error": "",
                }
            )
        except Exception as exc:
            failures += 1
            elapsed_ms = (time.perf_counter() - start) * 1000.0
            elapsed_values_ms.append(elapsed_ms)
            row.update(
                {
                    "algorithm": "FAILED",
                    "elapsed_ms": round(elapsed_ms, 3),
                    "exact": 0,
                    "molecule_ok": 0,
                    "atom_correct": 0,
                    "atom_total": 0,
                    "gold_total": "",
                    "our_total": "",
                    "delta_total": "",
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )
        per_reaction_rows.append(row)

    elapsed_total = time.perf_counter() - combo_start
    sorted_ms = sorted(elapsed_values_ms)
    label = f"{mcs_backend}/{substructure_backend}"
    summary = ComboSummary(
        label=label,
        mcs_backend=mcs_backend,
        substructure_backend=substructure_backend,
        sample_size=len(reactions),
        failures=failures,
        coverage_pct=100.0 * (len(reactions) - failures) / max(1, len(reactions)),
        elapsed_seconds=elapsed_total,
        reactions_per_second=len(reactions) / max(elapsed_total, 1e-9),
        median_ms=statistics.median(sorted_ms) if sorted_ms else 0.0,
        p95_ms=percentile(sorted_ms, 0.95),
        exact_pct=100.0 * exact_count / max(1, len(reactions)),
        molecule_pct=100.0 * molecule_count / max(1, len(reactions)),
        atom_pct=100.0 * atom_correct_total / max(1, atom_total_total),
        bond_change_pct=100.0 * bond_change_exact / max(1, len(reactions)),
        mean_abs_bond_delta=sum(bond_change_deltas) / max(1, len(bond_change_deltas)),
    )
    return summary, per_reaction_rows, algorithms


def run_combo_with_progress(
    mapper_module,
    reactions: list,
    mcs_backend: str,
    substructure_backend: str,
    timeout_sec: float,
    num_threads: int | None,
    progress_every: int,
):
    if progress_every <= 0:
        return run_combo(
            mapper_module,
            reactions,
            mcs_backend=mcs_backend,
            substructure_backend=substructure_backend,
            timeout_sec=timeout_sec,
            num_threads=num_threads,
        )

    per_reaction_rows: list[dict[str, object]] = []
    elapsed_values_ms: list[float] = []
    failures = 0
    exact_count = 0
    molecule_count = 0
    atom_correct_total = 0
    atom_total_total = 0
    bond_change_exact = 0
    bond_change_deltas: list[int] = []
    algorithms = Counter()
    combo_start = time.perf_counter()
    label = f"{mcs_backend}/{substructure_backend}"

    for index, reaction in enumerate(reactions, start=1):
        row: dict[str, object] = {
            "reaction_id": str(reaction.rxn_id),
            "mcs_backend": mcs_backend,
            "substructure_backend": substructure_backend,
        }
        start = time.perf_counter()
        try:
            result = mapper_module.map_reaction(
                reaction,
                mcs_timeout_sec=timeout_sec,
                mode="strict",
                mcs_backend=mcs_backend,
                substructure_backend=substructure_backend,
                num_threads=num_threads,
            )
            elapsed_ms = (time.perf_counter() - start) * 1000.0
            elapsed_values_ms.append(elapsed_ms)
            atom_correct, atom_total, exact, molecule_exact = mapper_module.check_accuracy(
                reaction,
                result.reactant_map_numbers,
                result.product_map_numbers,
            )
            if result.bond_change_calculator is not None:
                our_total, _our_formed, _our_cleaved, _our_order = mapper_module._bond_change_totals_from_calculator(
                    result.bond_change_calculator
                )
            else:
                our_total, _our_formed, _our_cleaved, _our_order = mapper_module.compute_bond_changes(
                    reaction,
                    result.reactant_map_numbers,
                    result.product_map_numbers,
                )
            try:
                gold_total, _gold_formed, _gold_cleaved, _gold_order = reaction.gold_bond_changes()
            except Exception:
                gold_total = -1
            delta_total = our_total - gold_total if our_total >= 0 and gold_total >= 0 else None
            if exact:
                exact_count += 1
            if molecule_exact:
                molecule_count += 1
            atom_correct_total += atom_correct
            atom_total_total += atom_total
            if delta_total == 0:
                bond_change_exact += 1
            if delta_total is not None:
                bond_change_deltas.append(abs(delta_total))
            algorithms[result.algorithm] += 1
            row.update(
                {
                    "algorithm": result.algorithm,
                    "elapsed_ms": round(elapsed_ms, 3),
                    "exact": int(exact),
                    "molecule_ok": int(molecule_exact),
                    "atom_correct": atom_correct,
                    "atom_total": atom_total,
                    "gold_total": gold_total,
                    "our_total": our_total,
                    "delta_total": delta_total if delta_total is not None else "",
                    "error": "",
                }
            )
        except Exception as exc:
            failures += 1
            elapsed_ms = (time.perf_counter() - start) * 1000.0
            elapsed_values_ms.append(elapsed_ms)
            row.update(
                {
                    "algorithm": "FAILED",
                    "elapsed_ms": round(elapsed_ms, 3),
                    "exact": 0,
                    "molecule_ok": 0,
                    "atom_correct": 0,
                    "atom_total": 0,
                    "gold_total": "",
                    "our_total": "",
                    "delta_total": "",
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )
        per_reaction_rows.append(row)

        if index % progress_every == 0 or index == len(reactions):
            elapsed = time.perf_counter() - combo_start
            atom_pct = 100.0 * atom_correct_total / max(1, atom_total_total)
            exact_pct = 100.0 * exact_count / index
            print(
                f"[{label}] {index}/{len(reactions)}  "
                f"{index / max(elapsed, 1e-9):.2f} rxn/s  "
                f"exact={exact_pct:.2f}%  atom={atom_pct:.2f}%  failures={failures}",
                flush=True,
            )

    elapsed_total = time.perf_counter() - combo_start
    sorted_ms = sorted(elapsed_values_ms)
    summary = ComboSummary(
        label=label,
        mcs_backend=mcs_backend,
        substructure_backend=substructure_backend,
        sample_size=len(reactions),
        failures=failures,
        coverage_pct=100.0 * (len(reactions) - failures) / max(1, len(reactions)),
        elapsed_seconds=elapsed_total,
        reactions_per_second=len(reactions) / max(elapsed_total, 1e-9),
        median_ms=statistics.median(sorted_ms) if sorted_ms else 0.0,
        p95_ms=percentile(sorted_ms, 0.95),
        exact_pct=100.0 * exact_count / max(1, len(reactions)),
        molecule_pct=100.0 * molecule_count / max(1, len(reactions)),
        atom_pct=100.0 * atom_correct_total / max(1, atom_total_total),
        bond_change_pct=100.0 * bond_change_exact / max(1, len(reactions)),
        mean_abs_bond_delta=sum(bond_change_deltas) / max(1, len(bond_change_deltas)),
    )
    return summary, per_reaction_rows, algorithms


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "reaction_id",
                "mcs_backend",
                "substructure_backend",
                "algorithm",
                "elapsed_ms",
                "exact",
                "molecule_ok",
                "atom_correct",
                "atom_total",
                "gold_total",
                "our_total",
                "delta_total",
                "error",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def markdown_summary(profile: SampleProfile, summaries: list[ComboSummary], algorithms: dict[str, Counter]) -> str:
    lines: list[str] = []
    lines.append("# Random 200 Reaction Backend Comparison")
    lines.append("")
    lines.append("## Sample Profile")
    lines.append("")
    lines.append("| Metric | Value |")
    lines.append("|---|---:|")
    lines.append(f"| Sample size | {profile.sample_size} |")
    lines.append(f"| Issue-free reactions | {profile.issue_free_pct:.2f}% |")
    lines.append(f"| Error-free reactions | {profile.error_free_pct:.2f}% |")
    lines.append(f"| Warning-free reactions | {profile.warning_free_pct:.2f}% |")
    lines.append(f"| Charge-balanced reactions | {profile.charge_balanced_pct:.2f}% |")
    lines.append(f"| Non-H balanced reactions | {profile.non_h_balanced_pct:.2f}% |")
    lines.append(f"| Median reactant molecules | {profile.median_reactant_mols:.1f} |")
    lines.append(f"| Median product molecules | {profile.median_product_mols:.1f} |")
    lines.append("")
    lines.append("## Backend Summary")
    lines.append("")
    lines.append("| Backend | Failures | Coverage | Exact | Mol->Mol | Atom->Atom | Bond-change | Mean abs bond delta | Median ms | P95 ms | Rxn/s |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for summary in summaries:
        lines.append(
            f"| {summary.label} | {summary.failures} | {summary.coverage_pct:.2f}% | "
            f"{summary.exact_pct:.2f}% | {summary.molecule_pct:.2f}% | {summary.atom_pct:.2f}% | "
            f"{summary.bond_change_pct:.2f}% | {summary.mean_abs_bond_delta:.3f} | "
            f"{summary.median_ms:.3f} | {summary.p95_ms:.3f} | {summary.reactions_per_second:.2f} |"
        )
    pure = [summary for summary in summaries if summary.mcs_backend == summary.substructure_backend]
    if pure:
        lines.append("")
        lines.append("## Pure Backend Focus")
        lines.append("")
        lines.append("| Backend | Exact | Mol->Mol | Atom->Atom | Bond-change | Median ms | Rxn/s |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|")
        for summary in pure:
            lines.append(
                f"| {summary.label} | {summary.exact_pct:.2f}% | {summary.molecule_pct:.2f}% | "
                f"{summary.atom_pct:.2f}% | {summary.bond_change_pct:.2f}% | "
                f"{summary.median_ms:.3f} | {summary.reactions_per_second:.2f} |"
            )
    lines.append("")
    lines.append("## Algorithm Distribution")
    lines.append("")
    for summary in summaries:
        lines.append(f"### {summary.label}")
        lines.append("")
        lines.append("| Algorithm | Count |")
        lines.append("|---|---:|")
        for algorithm, count in algorithms[summary.label].most_common():
            lines.append(f"| {algorithm} | {count} |")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def json_summary(seed: int, timeout_sec: float, sample_ids: list[str], profile: SampleProfile, summaries: list[ComboSummary]):
    return {
        "seed": seed,
        "timeout_sec": timeout_sec,
        "sample_ids": sample_ids,
        "sample_profile": asdict(profile),
        "summaries": [asdict(summary) for summary in summaries],
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compare reaction mapping backend combinations on a random 200-reaction sample.")
    parser.add_argument("--mapper", type=Path, default=DEFAULT_RDT_MAPPER, help="Path to rdt_port/rdkit_mapper.py")
    parser.add_argument("--rdf", type=Path, default=DEFAULT_RDF, help="Path to golden benchmark RDF")
    parser.add_argument("--sample-size", type=int, default=200, help="Number of reactions to sample")
    parser.add_argument("--seed", type=int, default=20260410, help="Sampling seed")
    parser.add_argument("--timeout", type=float, default=1.0, help="Per-reaction mapping timeout in seconds")
    parser.add_argument("--threads", type=int, default=1, help="Thread count to pass into the mapper")
    parser.add_argument("--combos", choices=("all", "pure"), default="all", help="Run all 4 combinations or only pure backend rows")
    parser.add_argument("--progress-every", type=int, default=25, help="Progress report cadence within each combo; 0 disables progress logging")
    parser.add_argument("--output-prefix", type=Path, default=DEFAULT_OUTPUT_PREFIX, help="Output prefix for TSV/JSON/Markdown")
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    mapper_module = load_mapper_module(args.mapper)
    reactions = mapper_module.parse_rdf(str(args.rdf), limit=None)
    sample = sample_reactions(reactions, args.sample_size, args.seed)
    profile = summarize_sample_profile(mapper_module, sample)
    sample_ids = [str(reaction.rxn_id) for reaction in sample]

    if args.combos == "pure":
        combos = [
            ("smsd", "smsd"),
            ("rdkit", "rdkit"),
        ]
    else:
        combos = [
            ("smsd", "smsd"),
            ("smsd", "rdkit"),
            ("rdkit", "smsd"),
            ("rdkit", "rdkit"),
        ]

    summaries: list[ComboSummary] = []
    all_rows: list[dict[str, object]] = []
    algorithm_tables: dict[str, Counter] = {}

    for mcs_backend, substructure_backend in combos:
        print(f"Running {mcs_backend}/{substructure_backend} on {len(sample)} sampled reactions...", flush=True)
        summary, rows, algorithms = run_combo_with_progress(
            mapper_module,
            sample,
            mcs_backend=mcs_backend,
            substructure_backend=substructure_backend,
            timeout_sec=args.timeout,
            num_threads=args.threads,
            progress_every=args.progress_every,
        )
        summaries.append(summary)
        all_rows.extend(rows)
        algorithm_tables[summary.label] = algorithms

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    tsv_path = args.output_prefix.with_suffix(".tsv")
    json_path = args.output_prefix.with_suffix(".json")
    md_path = args.output_prefix.with_suffix(".md")

    write_tsv(tsv_path, all_rows)
    json_path.write_text(
        json.dumps(json_summary(args.seed, args.timeout, sample_ids, profile, summaries), indent=2),
        encoding="utf-8",
    )
    md_path.write_text(markdown_summary(profile, summaries, algorithm_tables), encoding="utf-8")

    print(markdown_summary(profile, summaries, algorithm_tables))
    print(f"Wrote TSV: {tsv_path}")
    print(f"Wrote JSON: {json_path}")
    print(f"Wrote Markdown: {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
