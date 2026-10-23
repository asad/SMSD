# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Compare TargetCorpus screens on the GPU (smsd-cuda) with the same screens on the CPU.

The molecules of the SMILES files are repeated to each library size. In identity mode every
query, chemistry and threshold is screened on both and the ordered hits must be equal. The
accelerator's integer overlaps are also scored with the host expression and compared byte for
byte with the CPU scores. Thresholds include scores the library reaches exactly. In timing mode
the median of at least five screens after warm-up and every raw duration are recorded. The GPU
must stay attached throughout: a screen that falls back to the CPU fails
the run, since it would compare the CPU with itself.

    python screen_compare.py --mode identity --sizes 10000,100000,1000000 --out results a.smi b.txt
"""
import argparse
import array
import collections
import ctypes
import hashlib
import json
import os
import platform
import resource
import statistics
import struct
import sys
import time

import smsd
from smsd import _accelerator, _smsd


def read_molecules(paths):
    """The parsed molecules of the files, in order: the first field of each line not blank or a comment."""
    molecules, rejected = [], 0
    for path in paths:
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                fields = line.split()
                if not fields or fields[0].startswith("#"):
                    continue
                try:
                    molecules.append(smsd.parse_smiles(fields[0]))
                except Exception:  # noqa: BLE001 - a corpus line smsd rejects is counted, not fatal
                    rejected += 1
    return molecules, rejected


def chemistries():
    strict = smsd.ChemOptions()
    strict.aromaticity_mode = smsd.AromaticityMode.STRICT
    strict.ring_matches_ring_only = True
    return {"default": smsd.ChemOptions(), "strict": strict}


def peak_rss_mb():
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / (1 << 20) if sys.platform == "darwin" else peak / 1024


class Screens:
    """Screens one corpus on the GPU or the CPU, checking that the GPU never falls back."""

    def __init__(self, corpus):
        self.corpus = corpus

    def gpu(self, query, threshold, chemistry):
        os.environ.pop("SMSD_DISABLE_GPU", None)
        before = _smsd._screen_accelerator_failures()[0]
        hits = smsd.screen_targets(query, self.corpus, threshold, chemistry)
        failures, reason = _smsd._screen_accelerator_failures()
        if failures != before:
            raise SystemExit(f"the GPU screen fell back to the CPU: {reason}")
        return hits

    def cpu(self, query, threshold, chemistry, threads=0):
        os.environ["SMSD_DISABLE_GPU"] = "1"
        try:
            return smsd.screen_targets(query, self.corpus, threshold, chemistry, num_threads=threads)
        finally:
            os.environ.pop("SMSD_DISABLE_GPU", None)


class LabelCount(ctypes.Structure):
    """The public smsd_label_count ABI layout."""

    _fields_ = [("label", ctypes.c_uint16), ("count", ctypes.c_uint16)]


DeviceInfo = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.c_size_t,
                             ctypes.c_void_p, ctypes.c_size_t)
LabelOverlaps = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(ctypes.c_int32), ctypes.c_int32,
                                ctypes.c_int64, ctypes.POINTER(ctypes.c_uint64),
                                ctypes.POINTER(LabelCount), ctypes.POINTER(ctypes.c_int32),
                                ctypes.c_void_p, ctypes.c_size_t)


class AcceleratorTable(ctypes.Structure):
    """The public smsd_screen_accelerator_v1 ABI layout."""

    _fields_ = [("abi_version", ctypes.c_uint32), ("struct_size", ctypes.c_uint32),
                ("smsd_version", ctypes.c_char_p), ("name", ctypes.c_char_p),
                ("device_info", DeviceInfo), ("label_overlaps", LabelOverlaps)]


def histogram(molecule, chemistry):
    """The labels of screening_detail::ubLabel, or None when the CPU uses the graph."""
    if molecule.n > 65535:
        return None
    if not chemistry.match_atom_type:
        return {0: molecule.n} if molecule.n else {}
    labels = collections.Counter()
    strict = chemistry.aromaticity_mode == smsd.AromaticityMode.STRICT
    for atomic_number, aromatic, ring in zip(molecule.atomic_num, molecule.aromatic, molecule.ring):
        label = (atomic_number << 2) | (2 if strict and aromatic else 0)
        if chemistry.ring_matches_ring_only and ring:
            label |= 1
        if not 0 <= label < 1024:
            return None
        labels[label] += 1
    return labels


class ScoreIdentity:
    """Check the installed integer ABI without adding a production score-reporting API.

    Actual TargetCorpus calls independently check hit order. Here every ABI output is formed
    into a score using overlapSimilarity's expression and compared with the public CPU bound.
    Targets whose label table cannot hold them remain graph-scored on the CPU, as in production.
    """

    def __init__(self, pool, size, chemistry):
        import smsd_cuda

        self.pool, self.size, self.chemistry = pool, size, chemistry
        self.capsule = smsd_cuda.capsule()
        pointer = ctypes.pythonapi.PyCapsule_GetPointer
        pointer.restype = ctypes.c_void_p
        pointer.argtypes = [ctypes.py_object, ctypes.c_char_p]
        self.table = ctypes.cast(pointer(self.capsule, b"smsd.screen_accelerator"),
                                 ctypes.POINTER(AcceleratorTable)).contents
        if self.table.abi_version != 1 or self.table.struct_size < ctypes.sizeof(AcceleratorTable):
            raise SystemExit("score identity needs screening accelerator ABI 1")
        self.histograms = [histogram(molecule, chemistry) for molecule in pool]
        encoded = [b"".join(struct.pack("=HH", label, count) for label, count in (counts or {}).items())
                   for counts in self.histograms]
        copies, remainder = divmod(size, len(pool))
        data = bytearray(b"".join(encoded) * copies + b"".join(encoded[:remainder]))
        self.entries = (LabelCount * max(1, len(data) // ctypes.sizeof(LabelCount)))()
        if data:
            ctypes.memmove(self.entries, bytes(data), len(data))
        offsets = array.array("Q", [0])
        for target_index in range(size):
            offsets.append(offsets[-1] + len(encoded[target_index % len(pool)]) // 4)
        self.offsets = (ctypes.c_uint64 * len(offsets)).from_buffer(offsets)
        self.overlaps = (ctypes.c_int32 * size)()
        self.target_atoms = array.array("i", (pool[index % len(pool)].n for index in range(size)))
        self.held = bytearray(counts is not None for counts in self.histograms) * copies
        self.held.extend(counts is not None for counts in self.histograms[:remainder])

    def compare(self, query):
        counts = histogram(query, self.chemistry)
        if counts is None:
            raise SystemExit("the score identity query cannot use the GPU label table")
        query_counts = (ctypes.c_int32 * 1024)()
        for label, count in counts.items():
            query_counts[label] = count
        error = ctypes.create_string_buffer(256)
        if self.table.label_overlaps(query_counts, len(query_counts), self.size, self.offsets,
                                     self.entries, self.overlaps, error, len(error)) != 0:
            raise SystemExit(f"direct GPU score identity failed: {error.value.decode(errors='replace')}")
        expected_pool = [smsd.similarity_upper_bound(query, target, self.chemistry) for target in self.pool]
        copies, remainder = divmod(self.size, len(self.pool))
        expected = array.array("d", expected_pool) * copies
        expected.extend(expected_pool[:remainder])
        actual = array.array("d")
        query_atoms = query.n
        for target_index, (overlap, atoms, held) in enumerate(zip(self.overlaps, self.target_atoms, self.held)):
            if overlap < 0 or overlap > (atoms if held else 0):
                raise SystemExit(f"invalid GPU overlap {overlap} at target {target_index}")
            atom_overlap = min(overlap, query_atoms, atoms)
            actual.append((atom_overlap / (query_atoms + atoms - atom_overlap)
                           if query_atoms and atoms and atom_overlap else 0.0)
                          if held else expected_pool[target_index % len(self.pool)])
        expected_bytes, actual_bytes = expected.tobytes(), actual.tobytes()
        return {"equal": expected_bytes == actual_bytes, "targets": self.size,
                "gpu_histogram_targets": sum(self.held), "cpu_graph_targets": self.size - sum(self.held),
                "cpu_sha256": hashlib.sha256(expected_bytes).hexdigest(),
                "gpu_sha256": hashlib.sha256(actual_bytes).hexdigest()}


def timed(function, repeats):
    start = time.perf_counter_ns()
    function()
    warmup = time.perf_counter_ns() - start
    times, result = [], None
    for _ in range(repeats):
        start = time.perf_counter_ns()
        result = function()
        times.append(time.perf_counter_ns() - start)
    return statistics.median(times) / 1_000_000.0, result, times, warmup


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--mode", choices=("identity", "timing"), required=True)
    parser.add_argument("--sizes", default="10000,100000,1000000")
    parser.add_argument("--queries", type=int, default=0, help="queries per size (default 12, timing 5)")
    parser.add_argument("--repeats", type=int, default=5, help="screens per timing, at least five")
    parser.add_argument("--out", required=True, help="folder for the TSV and JSON results")
    parser.add_argument("smiles", nargs="+")
    args = parser.parse_args()
    if args.mode == "timing" and args.repeats < 5:
        parser.error("timing needs at least five repetitions after warm-up")
    os.makedirs(args.out, exist_ok=True)

    os.environ.pop("SMSD_DISABLE_GPU", None)
    accelerator = _accelerator.active()
    if accelerator is None:
        raise SystemExit("no GPU is attached: install smsd-cuda for this smsd on a machine with an NVIDIA GPU")

    pool, rejected = read_molecules(args.smiles)
    if not pool:
        raise SystemExit("no molecules were read")
    count = args.queries or (12 if args.mode == "identity" else 5)
    queries = [pool[(i * len(pool)) // count] for i in range(count)]
    options = chemistries() if args.mode == "identity" else {"default": smsd.ChemOptions()}
    summary = {"mode": args.mode, "smsd": smsd.__version__, "device": smsd.gpu_device_info(),
               "accelerator": list(accelerator), "python": platform.python_version(),
               "cpus": os.cpu_count(), "molecules": len(pool), "rejected_lines": rejected,
               "files": [os.path.basename(path) for path in args.smiles], "sizes": []}
    if args.mode == "timing":
        summary.update({"repetitions": args.repeats, "warmup_screens_per_case": 1,
                        "raw_timing_unit": "ns", "timings_include_transfers": True})
    table = open(os.path.join(args.out, f"{args.mode}.tsv"), "w", encoding="utf-8")
    raw = None
    if args.mode == "identity":
        table.write("size\tchemistry\tquery\tthreshold\tcpu_hits\tgpu_hits\tequal\tequal_scores\n")
    else:
        table.write("size\tchemistry\tquery\tthreshold\thits\tcpu_ms\tcpu_1_thread_ms\tgpu_ms\tcpu_over_gpu\n")
        raw = open(os.path.join(args.out, "timing-raw.tsv"), "w", encoding="utf-8")
        raw.write("size\tchemistry\tquery\tthreshold\tbackend\trepetition\twarmup\telapsed_ns\thits\n")
    mismatches = 0

    for size in (int(value) for value in args.sizes.split(",")):
        start = time.perf_counter()
        corpus = smsd.TargetCorpus()
        corpus.add_targets([pool[i % len(pool)] for i in range(size)])
        built = time.perf_counter() - start
        screens = Screens(corpus)
        record = {"size": size, "corpus_seconds": round(built, 3), "screens": 0, "mismatches": 0}
        for name, chemistry in options.items():
            # The first screen under a label rule records the histograms both paths read.
            start = time.perf_counter()
            screens.cpu(queries[0], 0.5, chemistry)
            record[f"histograms_{name}_seconds"] = round(time.perf_counter() - start, 3)
            score_identity = ScoreIdentity(pool, size, chemistry) if args.mode == "identity" else None
            for index, query in enumerate(queries):
                if args.mode == "identity":
                    score_result = score_identity.compare(query)
                    record.setdefault("score_identity", []).append({"chemistry": name, "query": index,
                                                                     **score_result})
                    exact = sorted({smsd.similarity_upper_bound(query, pool[k], chemistry)
                                    for k in range(0, len(pool), max(1, len(pool) // 6))})
                    thresholds = sorted({0.0, 0.3, 0.5, 0.7, 0.9, 1.0, *exact})
                else:
                    thresholds = [0.5, 0.8]
                for threshold in thresholds:
                    if args.mode == "identity":
                        cpu = screens.cpu(query, threshold, chemistry)
                        gpu = screens.gpu(query, threshold, chemistry)
                        equal = cpu == gpu
                        case_equal = equal and score_result["equal"]
                        mismatches += not case_equal
                        record["mismatches"] += not case_equal
                        table.write(f"{size}\t{name}\t{index}\t{threshold!r}\t{len(cpu)}\t{len(gpu)}\t{int(equal)}"
                                    f"\t{int(score_result['equal'])}\n")
                    else:
                        measurements = {}
                        for backend, function in (
                            ("cpu", lambda: screens.cpu(query, threshold, chemistry)),
                            ("cpu_1_thread", lambda: screens.cpu(query, threshold, chemistry, 1)),
                            ("gpu", lambda: screens.gpu(query, threshold, chemistry)),
                        ):
                            elapsed, hits, durations, warmup = timed(function, args.repeats)
                            measurements[backend] = (elapsed, hits)
                            raw.write(f"{size}\t{name}\t{index}\t{threshold}\t{backend}\t0\t1\t{warmup}\t{len(hits)}\n")
                            for repetition, duration in enumerate(durations, 1):
                                raw.write(f"{size}\t{name}\t{index}\t{threshold}\t{backend}\t{repetition}"
                                          f"\t0\t{duration}\t{len(hits)}\n")
                        cpu_ms, cpu = measurements["cpu"]
                        one_ms, _ = measurements["cpu_1_thread"]
                        gpu_ms, gpu = measurements["gpu"]
                        if cpu != gpu:
                            mismatches += 1
                            record["mismatches"] += 1
                        table.write(f"{size}\t{name}\t{index}\t{threshold}\t{len(cpu)}\t{cpu_ms:.6f}\t{one_ms:.6f}"
                                    f"\t{gpu_ms:.6f}\t{cpu_ms / gpu_ms:.6f}\n")
                    record["screens"] += 1
            del score_identity
        record["peak_rss_mb"] = round(peak_rss_mb(), 1)
        summary["sizes"].append(record)
        table.flush()
        if raw is not None:
            raw.flush()
        print(json.dumps(record), flush=True)
        del corpus, screens

    table.close()
    if raw is not None:
        raw.close()
    summary["mismatches"] = mismatches
    with open(os.path.join(args.out, f"{args.mode}.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    print(f"{args.mode}: {mismatches} mismatches", flush=True)
    return 1 if mismatches else 0


if __name__ == "__main__":
    sys.exit(main())
