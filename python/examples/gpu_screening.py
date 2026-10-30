# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Screen on an NVIDIA GPU when one is ready, and get the same results on the CPU otherwise."""

import os
from pathlib import Path
from tempfile import TemporaryDirectory

import smsd
from smsd.clustering import butina, maxmin
from smsd.corpus import FingerprintCorpus

try:
    import smsd_cuda
except ImportError:
    smsd_cuda = None

# Automatic GPU use for TargetCorpus screening needs a Tesla T4 and a prepared corpus of at
# least 100,000 targets; smaller corpora and other devices screen on the CPU with the same
# hits. On other NVIDIA GPUs, ask for the device with the packed backend="cuda" or with
# match_corpus(..., gpu_prefilter=True). SMSD_DISABLE_GPU=1 keeps all screening on the CPU.


def cuda_status():
    """Return whether backend="cuda" can be used here, and the reason when it cannot."""
    if os.environ.get("SMSD_DISABLE_GPU", "") not in ("", "0"):
        return False, "disabled by SMSD_DISABLE_GPU"
    if smsd_cuda is None:
        return False, "smsd-cuda is not installed"
    if smsd_cuda.__version__ != smsd.__version__:
        return False, f"smsd-cuda {smsd_cuda.__version__} does not match smsd {smsd.__version__}"
    reason = smsd_cuda.unavailable_reason()
    return reason is None, reason


def decisions(results):
    """Each target's decision, atom bound and matches, leaving out timings."""
    return [(item.target_index, item.atom_upper_bound, item.threshold_reached,
             item.proven_below_threshold, None if item.result is None else
             (item.result.status, item.result.optimal, item.result.complete, item.result.matches))
            for item in results]


def main():
    library = ["c1ccccc1O", "c1ccccc1N", "Cc1ccccc1", "CCO", "CCN", "CC(=O)O",
               "CCOc1ccccc1", "C1CCCCC1"]
    query = "c1ccccc1O"
    cuda_ready, reason = cuda_status()
    print("CUDA screening:", "available" if cuda_ready else f"unavailable ({reason})")
    print("Compute backend:", smsd.gpu_device_info())

    options = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
    with TemporaryDirectory() as directory:
        with FingerprintCorpus.from_smiles(Path(directory) / "library.fp", library,
                                           options=options, num_threads=1) as corpus:
            bits = smsd.fingerprint(query, options=options)
            hits = list(corpus.similarities(bits, threshold=0.3, backend="cpu"))
            # "auto" may run on either device, so check its results rather than where they ran.
            assert list(corpus.similarities(bits, threshold=0.3, backend="auto")) == hits
            clusters = butina(corpus, threshold=0.3, backend="auto")
            assert clusters == butina(corpus, threshold=0.3, backend="cpu")
            selection = maxmin(corpus, count=3, backend="auto")
            assert selection == maxmin(corpus, count=3, backend="cpu")
            print("Packed hits:", hits)
            print("Butina clusters:", clusters)
            print("MaxMin selection:", selection)
            if cuda_ready:
                before = corpus.cuda_calls
                assert list(corpus.similarities(bits, threshold=0.3, backend="cuda")) == hits
                assert corpus.cuda_calls > before
                print("Explicit CUDA hits match the CPU hits")

    targets = smsd.TargetCorpus.from_smiles(library)
    targets.prewarm(num_threads=1)
    print("TargetCorpus hits:", targets.screen(query, 0.3, num_threads=1))
    # Targets sharing all seven atoms of the query; the decisions do not depend on the device.
    search = smsd.MCS.find(query)
    on_cpu = search.match_corpus(targets, 7, gpu_prefilter=False, num_threads=1)
    gpu_requested = search.match_corpus(targets, 7, gpu_prefilter=True, num_threads=1)
    assert decisions(gpu_requested) == decisions(on_cpu)
    print("All 7 atoms shared:", [item.target_index for item in on_cpu if item.threshold_reached])
    print("Proved to share fewer:", [item.target_index for item in on_cpu
                                     if item.proven_below_threshold])


if __name__ == "__main__":
    main()
