# CPU and GPU screening

Graph matching in SMSD Pro, which covers every substructure, SMARTS and MCS
search, always runs on the CPU. An NVIDIA GPU can take over screening steps:
the ECFP4 screens of a prepared `TargetCorpus`, the pre-screen of a corpus
search and the scans of a packed fingerprint index. Hits, scores and search
decisions are identical on either device. Every Python program on this page
runs on a CPU-only machine; the GPU branches run only where a device is ready.

## What runs where

| Operation | CPU | NVIDIA GPU (`smsd-cuda`) | How it is chosen |
| --- | --- | --- | --- |
| Substructure, SMARTS and MCS matching | Always | Never | Not applicable |
| Fingerprints, `batch_similarity` and similarity coefficients | Always | Never | Not applicable |
| `screen_targets`, `TargetCorpus.screen_targets_many`, and `Substructure.match_many` on a `TargetCorpus` | Always | Never with `smsd-cuda` 8.1.0 | Not applicable |
| ECFP4 screens of a prepared `TargetCorpus`: `screen`, `screen_many`, `warm_start`, `stream_screen_smiles` and `stream_screen_many_smiles` | Below the device rule, without a usable device, or after a device failure | A Tesla T4 for 100,000 targets or more, or any supported device from `SMSD_CUDA_MIN_TARGETS` targets | Automatically, for each screen |
| Pre-screen of `Substructure.match_corpus` and `MCS.match_corpus` | `gpu_prefilter=False` (the default), or without a usable device | `gpu_prefilter=True` with an attached device | By the caller, with CPU fallback |
| Packed `FingerprintCorpus` scans, `butina` and `maxmin` | `backend="cpu"` (the default), or `"auto"` without a usable device | `backend="cuda"`, or `"auto"` with a usable device | By the caller; `"cuda"` raises without a device |

The searches that follow a pre-screen, and every `match`, `match_all` and
`match_many` call, use the CPU. Python never uses Metal; see
[Metal and C++](#metal-and-c).

The `smsd` package looks for `smsd-cuda` once per process, the first time a
prepared screen, `warm_start`, a `gpu_prefilter=True` search,
`gpu_is_available()` or `gpu_device_info()` needs it. Importing `smsd` does
not load the CUDA runtime, so the first such call carries the device start-up
cost; `smsd_cuda.preload()` moves that cost to a point the program chooses.

## Check the backend

`smsd.gpu_is_available()` returns True when an attached `smsd-cuda` package
has found a usable device. `smsd.gpu_device_info()` describes the backend: on
a CPU-only machine it starts with `CPU:` and names the OpenMP version and
thread count; with an attached device it starts with `GPU:`, followed by the
device name, its compute capability, its memory and, in brackets, the CPU
threading. Both report what is available, not which device a call used.

The `cuda_status()` helper below is the one in the
[GPU screening example](../../python/examples/gpu_screening.py). It decides
whether `backend="cuda"` can be used and gives the reason when it cannot.

```python
import os

import smsd

try:
    import smsd_cuda
except ImportError:
    smsd_cuda = None


def cuda_status():
    """Return whether backend="cuda" can be used here, and the reason when it cannot."""
    if os.environ.get("SMSD_DISABLE_GPU", "") not in ("", "0"):
        return False, "disabled by SMSD_DISABLE_GPU"
    if smsd_cuda is None:
        return False, "smsd-cuda is not installed"
    if smsd_cuda.__version__ != smsd.__version__:
        return False, (f"smsd-cuda {smsd_cuda.__version__} does not match "
                       f"smsd {smsd.__version__}")
    reason = smsd_cuda.unavailable_reason()
    return reason is None, reason


ready, reason = cuda_status()
print("CUDA screening:", "ready" if ready else f"unavailable ({reason})")
assert ready or reason

available = smsd.gpu_is_available()
print("GPU available:", available)
assert isinstance(available, bool)

# 'CPU: ...' without a device; 'GPU: <device>  [...]' once smsd-cuda is attached.
backend = smsd.gpu_device_info()
assert isinstance(backend, str) and backend

if ready:
    smsd_cuda.preload()  # initialise the device and load its kernels now
    assert "compute" in smsd_cuda.device_info()
```

On a machine without `smsd-cuda`, the program prints:

```text
CUDA screening: unavailable (smsd-cuda is not installed)
GPU available: False
```

With `smsd-cuda` and a usable device, the program prints
`CUDA screening: ready` and `GPU available: True`. Otherwise
`smsd_cuda.unavailable_reason()` explains why, for example
`no CUDA device was found`, or a reason that names a device below compute
capability 7.5.

## Install the CUDA package

The `smsd-cuda` package adds NVIDIA screening to `smsd` on Linux x86_64.
Install it with the same version of `smsd`:

```bash
python -m pip install "smsd==8.1.0" "smsd-cuda==8.1.0"
```

Packed indexes and clustering also need NumPy, which neither package
installs:

```bash
python -m pip install "smsd[numpy]==8.1.0" "smsd-cuda==8.1.0"
```

- The package needs Linux x86_64 and CPython 3.12–3.15; the
  [supported-platforms table](../HOWTO-INSTALL.md#supported-platforms) gives
  the C library version.
- The device must be an NVIDIA GPU of compute capability 7.5 or later, with a
  driver compatible with the CUDA 12 runtime. The runtime is installed as a
  dependency; the CUDA toolkit is not needed.
- The versions must match exactly. `smsd-cuda` 8.1.0 requires `smsd` 8.1.0,
  and `smsd` does not attach a package built for another version.
- No `smsd-cuda` wheel is published for Windows or macOS. On Windows, the GPU
  container runs through WSL2; see [Containers](#containers).

Prepared screens, the stream functions and `match_corpus` need no NumPy.

## Automatic use in prepared screens

Prepared ECFP4 screens choose their device automatically.
`TargetCorpus.screen` and `screen_many` after `prewarm()`, `warm_start`, and
the stream functions `stream_screen_smiles` and `stream_screen_many_smiles`
use the GPU when `smsd-cuda` is attached and the corpus is large enough.
Otherwise, and whenever device work fails, they run on the CPU with the same
hits. The code is the same on every machine.

By default the GPU is used only on a Tesla T4 (compute capability 7.5), and
only for a prepared corpus of at least 100,000 targets. Setting
`SMSD_CUDA_MIN_TARGETS` to a positive whole number N lets any supported device
take corpora of N targets or more; on a Tesla T4, N replaces 100,000. The
stream functions apply the rule to each chunk of targets (4,096 by default),
so they use the device only when a chunk is large enough.

```python
import os

import smsd
from smsd.corpus import stream_screen_smiles

benzodiazepines = {
    "diazepam": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
    "nordazepam": "O=C1CN=C(c2ccccc2)c2cc(Cl)ccc2N1",
    "temazepam": "CN1C(=O)C(O)N=C(c2ccccc2)c2cc(Cl)ccc21",
    "oxazepam": "OC1N=C(c2ccccc2)c2cc(Cl)ccc2NC1=O",
    "lorazepam": "OC1N=C(c2ccccc2Cl)c2cc(Cl)ccc2NC1=O",
    "clonazepam": "O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1",
    "nitrazepam": "O=C1CN=C(c2ccccc2)c2cc([N+](=O)[O-])ccc2N1",
    "flunitrazepam": "CN1C(=O)CN=C(c2ccccc2F)c2cc([N+](=O)[O-])ccc21",
    "bromazepam": "O=C1CN=C(c2ccccn2)c2cc(Br)ccc2N1",
    "tetrazepam": "CN1C(=O)CN=C(C2=CCCCC2)c2cc(Cl)ccc21",
    "alprazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2",
    "triazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1Cl)=NC2",
    "midazolam": "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
    "clobazam": "CN1C(=O)CC(=O)N(c2ccccc2)c2cc(Cl)ccc21",
}
names = list(benzodiazepines)
library = list(benzodiazepines.values())
queries = [benzodiazepines["diazepam"], benzodiazepines["alprazolam"]]

# The same calls on every machine; each screen chooses its device.
corpus = smsd.TargetCorpus.from_smiles(library)
corpus.prewarm()
hits = corpus.screen_many(queries, 0.4)
for query, found in zip(("diazepam", "alprazolam"), hits):
    print(f"{query}:", ", ".join(names[index] for index in found))
assert hits == [[0, 1, 2, 9, 13], [10, 11, 12]]

# warm_start prepares the corpus it is given; streams prepare one chunk at a time.
fresh = smsd.TargetCorpus.from_smiles(library)
assert smsd.warm_start(fresh, queries[0], 0.4) == hits[0]
streamed = list(stream_screen_smiles(queries[0], iter(library), 0.4, chunk_size=5))
assert streamed == hits[0]

# Let any supported GPU take corpora of one target or more; the CPU ignores it.
os.environ["SMSD_CUDA_MIN_TARGETS"] = "1"
same = corpus.screen_many(queries, 0.4) == hits
print("same hits with SMSD_CUDA_MIN_TARGETS=1:", same)
assert same
```

The program prints:

```text
diazepam: diazepam, nordazepam, temazepam, tetrazepam, clobazam
alprazolam: alprazolam, triazolam, midazolam
same hits with SMSD_CUDA_MIN_TARGETS=1: True
```

A CPU-only machine ignores `SMSD_CUDA_MIN_TARGETS`. No public record shows
which device ran a prepared screen; the hits are the same either way.

## Request the GPU for a corpus pre-screen

`Substructure.match_corpus` and `MCS.match_corpus` rule out, with a proof,
the targets that cannot contain the query or cannot share `minimum_atoms`
atoms with it, and search the remaining targets on the CPU. With
`gpu_prefilter=True` they ask an attached device to run that pre-screen,
whatever the corpus size. Without a usable device, or if device work fails,
the CPU pre-screen runs instead. Decisions, bounds and matches are identical
either way, and no `prewarm()` is needed.

SMARTS queries, and searches with `match_atom_type=False`, a node limit, an
immediate timeout or a cancellation token, search every target on the CPU
without a pre-screen. MCS searches do the same for MCES, atom weights and
excluded target atoms.

```python
import smsd

benzodiazepines = {
    "diazepam": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
    "nordazepam": "O=C1CN=C(c2ccccc2)c2cc(Cl)ccc2N1",
    "temazepam": "CN1C(=O)C(O)N=C(c2ccccc2)c2cc(Cl)ccc21",
    "oxazepam": "OC1N=C(c2ccccc2)c2cc(Cl)ccc2NC1=O",
    "lorazepam": "OC1N=C(c2ccccc2Cl)c2cc(Cl)ccc2NC1=O",
    "clonazepam": "O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1",
    "nitrazepam": "O=C1CN=C(c2ccccc2)c2cc([N+](=O)[O-])ccc2N1",
    "flunitrazepam": "CN1C(=O)CN=C(c2ccccc2F)c2cc([N+](=O)[O-])ccc21",
    "bromazepam": "O=C1CN=C(c2ccccn2)c2cc(Br)ccc2N1",
    "tetrazepam": "CN1C(=O)CN=C(C2=CCCCC2)c2cc(Cl)ccc21",
    "alprazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2",
    "triazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1Cl)=NC2",
    "midazolam": "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
    "clobazam": "CN1C(=O)CC(=O)N(c2ccccc2)c2cc(Cl)ccc21",
}
names = list(benzodiazepines)
corpus = smsd.TargetCorpus.from_smiles(list(benzodiazepines.values()))


def decisions(results):
    """Each target's decision, atom bound and matches, leaving out timings."""
    return [(item.target_index, item.atom_upper_bound, item.threshold_reached,
             item.proven_below_threshold, None if item.result is None else
             (item.result.status, item.result.optimal, item.result.complete,
              item.result.matches))
            for item in results]


# 5-phenyl-1,4-benzodiazepin-2-one core
core = smsd.Substructure.find("O=C1CN=C(c2ccccc2)c2ccccc2N1")
on_cpu = core.match_corpus(corpus)
requested = core.match_corpus(corpus, gpu_prefilter=True)
assert decisions(requested) == decisions(on_cpu)
found = [names[item.target_index] for item in on_cpu if item.threshold_reached]
ruled_out = [names[item.target_index] for item in on_cpu if item.result is None]
print("core found in", len(found), "targets")
print("ruled out before matching:", ", ".join(ruled_out))
assert len(found) == 9
assert ruled_out == ["alprazolam", "triazolam", "midazolam", "clobazam"]

# Targets that share all 22 atoms of alprazolam
mcs = smsd.MCS.find(benzodiazepines["alprazolam"])
on_cpu = mcs.match_corpus(corpus, 22)
requested = mcs.match_corpus(corpus, 22, gpu_prefilter=True)
assert decisions(requested) == decisions(on_cpu)
shared = [names[item.target_index] for item in on_cpu if item.threshold_reached]
print("all 22 atoms shared:", ", ".join(shared))
assert shared == ["alprazolam", "triazolam"]
```

The program prints:

```text
core found in 9 targets
ruled out before matching: alprazolam, triazolam, midazolam, clobazam
all 22 atoms shared: alprazolam, triazolam
```

A target that was ruled out has `result` set to `None` and
`proven_below_threshold` set to True. The other `CorpusSearchResult` fields
are described in
[Screening a library](SCREENING.md#substructure-screen-with-exact-confirmation).

## Choose a backend for packed indexes

The searches of a packed `FingerprintCorpus` (`similarities`, `screen_many`,
`threshold_pairs` and `top_k_pairs`) and the clustering functions `butina`
and `maxmin` take a `backend` argument:

| Backend | Behaviour |
| --- | --- |
| `"cpu"` | The default. The scan runs on the CPU and does not load `smsd-cuda`. |
| `"cuda"` | The scan runs on the GPU. It raises `RuntimeError` when no device is usable and `ValueError` for fingerprints wider than 65,536 bits. |
| `"auto"` | The scan runs on the GPU whenever a device is usable and the width is supported, for an index of any size, and on the CPU otherwise. |

Hits, scores, their order, clusters and selections are identical for all
three. Once device work has started, an error is raised rather than retried
on the CPU, under `"auto"` as well. A packed index needs NumPy; see
[Packed fingerprint index](PACKED_INDEX.md) for building and reopening one.

```python
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


def cuda_status():
    """Return whether backend="cuda" can be used here, and the reason when it cannot."""
    if os.environ.get("SMSD_DISABLE_GPU", "") not in ("", "0"):
        return False, "disabled by SMSD_DISABLE_GPU"
    if smsd_cuda is None:
        return False, "smsd-cuda is not installed"
    if smsd_cuda.__version__ != smsd.__version__:
        return False, (f"smsd-cuda {smsd_cuda.__version__} does not match "
                       f"smsd {smsd.__version__}")
    reason = smsd_cuda.unavailable_reason()
    return reason is None, reason


benzodiazepines = {
    "diazepam": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
    "nordazepam": "O=C1CN=C(c2ccccc2)c2cc(Cl)ccc2N1",
    "temazepam": "CN1C(=O)C(O)N=C(c2ccccc2)c2cc(Cl)ccc21",
    "oxazepam": "OC1N=C(c2ccccc2)c2cc(Cl)ccc2NC1=O",
    "lorazepam": "OC1N=C(c2ccccc2Cl)c2cc(Cl)ccc2NC1=O",
    "clonazepam": "O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1",
    "nitrazepam": "O=C1CN=C(c2ccccc2)c2cc([N+](=O)[O-])ccc2N1",
    "flunitrazepam": "CN1C(=O)CN=C(c2ccccc2F)c2cc([N+](=O)[O-])ccc21",
    "bromazepam": "O=C1CN=C(c2ccccn2)c2cc(Br)ccc2N1",
    "tetrazepam": "CN1C(=O)CN=C(C2=CCCCC2)c2cc(Cl)ccc21",
    "alprazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2",
    "triazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1Cl)=NC2",
    "midazolam": "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
    "clobazam": "CN1C(=O)CC(=O)N(c2ccccc2)c2cc(Cl)ccc21",
}
names = list(benzodiazepines)
options = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
ready, reason = cuda_status()

with TemporaryDirectory() as directory:
    path = Path(directory) / "benzodiazepines.fp"
    with FingerprintCorpus.from_smiles(path, list(benzodiazepines.values()),
                                       options=options) as corpus:
        query = smsd.fingerprint(benzodiazepines["diazepam"], options=options)
        hits = list(corpus.similarities(query, threshold=0.4, backend="cpu"))
        named = [(names[index], round(score, 3)) for index, score in hits]
        for name, score in named:
            print(f"{name:<11} {score:.3f}")
        assert named == [("diazepam", 1.0), ("nordazepam", 0.528),
                         ("temazepam", 0.582), ("tetrazepam", 0.542),
                         ("clobazam", 0.5)]
        print("CUDA calls after the CPU scan:", corpus.cuda_calls)
        assert corpus.cuda_calls == 0

        # "auto" may run on either device, so check its results, not where it ran.
        automatic = list(corpus.similarities(query, threshold=0.4, backend="auto"))
        assert automatic == hits
        clusters = butina(corpus, threshold=0.4, backend="auto")
        assert clusters == butina(corpus, threshold=0.4, backend="cpu")
        picked = maxmin(corpus, count=3, backend="auto")
        assert picked == maxmin(corpus, count=3, backend="cpu")

        if ready:
            before = corpus.cuda_calls
            on_gpu = list(corpus.similarities(query, threshold=0.4, backend="cuda"))
            assert on_gpu == hits
            assert corpus.cuda_calls > before  # the device completed the scan
            assert butina(corpus, threshold=0.4, backend="cuda") == clusters
            assert maxmin(corpus, count=3, backend="cuda") == picked
```

The program prints:

```text
diazepam    1.000
nordazepam  0.528
temazepam   0.582
tetrazepam  0.542
clobazam    0.500
CUDA calls after the CPU scan: 0
```

The branch under `if ready:` runs only where `backend="cuda"` can be used;
on any other machine the program prints the same lines.

## Keep work on the CPU

Three switches keep work on the CPU:

- `SMSD_DISABLE_GPU`, set to any value other than an empty string or `0`,
  disables every GPU route in the process. Prepared screens and pre-screens
  run on the CPU, `"auto"` uses the CPU and `backend="cuda"` raises. The
  variable is read at each call, so a program can set and remove it as it
  runs. It is the only switch for prepared screens.
- `gpu_prefilter=False`, the default, keeps a corpus pre-screen on the CPU.
- `backend="cpu"`, the default, keeps a packed scan or clustering run on the
  CPU.

```python
import os
from pathlib import Path
from tempfile import TemporaryDirectory

import smsd
from smsd.corpus import FingerprintCorpus

os.environ["SMSD_DISABLE_GPU"] = "1"
available = smsd.gpu_is_available()
print("GPU available:", available)
assert available is False
marked = smsd.gpu_device_info().endswith("[GPU disabled]")
print("backend string ends with [GPU disabled]:", marked)
assert marked

# Prepared screens and pre-screens run on the CPU with their usual results.
corpus = smsd.TargetCorpus.from_smiles(["CCO", "CCN", "c1ccccc1O"])
corpus.prewarm()
assert corpus.screen("CCO", 1.0) == [0]
search = smsd.Substructure.find("c1ccccc1O")
assert search.match_corpus(corpus, gpu_prefilter=True)[2].threshold_reached

with TemporaryDirectory() as directory:
    with FingerprintCorpus.create(Path(directory) / "bits.fp", [[1], [1, 2], [3]],
                                  fp_size=8) as packed:
        hits = list(packed.similarities([1], threshold=0.5, backend="auto"))
        print("auto:", hits, "| CUDA calls:", packed.cuda_calls)
        assert hits == [(0, 1.0), (1, 0.5)] and packed.cuda_calls == 0

# The variable is read at each call, so removing it allows the GPU again.
del os.environ["SMSD_DISABLE_GPU"]
marked = smsd.gpu_device_info().endswith("[GPU disabled]")
print("after removing the variable:", marked)
assert not marked
```

The program prints:

```text
GPU available: False
backend string ends with [GPU disabled]: True
auto: [(0, 1.0), (1, 0.5)] | CUDA calls: 0
after removing the variable: False
```

## Fallback and errors

The routes treat a missing or failing device differently:

| Route | No usable device | Error during device work |
| --- | --- | --- |
| Prepared screens, `warm_start` and the stream functions | The CPU | The CPU, silently, with the same hits |
| `match_corpus(..., gpu_prefilter=True)` | The CPU pre-screen | The CPU pre-screen, with the same decisions |
| Packed `backend="auto"` | The CPU | The error is raised |
| Packed `backend="cuda"` | `RuntimeError` | The error is raised |

The packed searches are generators. A `backend="cuda"` request that cannot be
met raises when the first item is requested, not when the method is called,
so keep the iteration inside the `try` block. The program below sets
`SMSD_DISABLE_GPU` so that its output is the same on every machine.

```python
import os
from pathlib import Path
from tempfile import TemporaryDirectory

import smsd
from smsd.corpus import FingerprintCorpus

# Make the device unusable, so that every machine prints the same lines.
os.environ["SMSD_DISABLE_GPU"] = "1"

with TemporaryDirectory() as directory:
    with FingerprintCorpus.create(Path(directory) / "bits.fp", [[1], [1, 2], [3]],
                                  fp_size=8) as corpus:
        pending = corpus.similarities([1], threshold=0.5, backend="cuda")
        print("the call returns a", type(pending).__name__)
        assert type(pending).__name__ == "generator"
        try:
            next(pending)
        except RuntimeError as error:
            message = str(error)
        else:
            raise AssertionError("backend='cuda' should fail without a device")
        print("first item:", message)
        assert message.endswith("CUDA screening is disabled by SMSD_DISABLE_GPU")

    # Fingerprints wider than 65,536 bits stay on the CPU under "auto".
    with FingerprintCorpus.create(Path(directory) / "wide.fp", [[1], [70000]],
                                  fp_size=70001) as wide:
        hits = list(wide.similarities([70000], threshold=0.5, backend="auto"))
        print("auto, 70,001 bits:", hits)
        assert hits == [(1, 1.0)]
        try:
            list(wide.similarities([70000], threshold=0.5, backend="cuda"))
        except ValueError as error:
            message = str(error)
        else:
            raise AssertionError("backend='cuda' should refuse 70,001 bits")
        print("cuda, 70,001 bits:", message)
        assert message == "CUDA packed screening supports at most 65536 bits"

# Argument errors are raised at once, whatever the machine.
corpus = smsd.TargetCorpus.from_smiles(["CCO", "CCN"])
try:
    smsd.Substructure.find("CO").match_corpus(corpus, gpu_prefilter=1)
except TypeError as error:
    message = str(error)
print("gpu_prefilter=1:", message)
assert message == "GPU prefilter flag must be a bool"
```

The program prints:

```text
the call returns a generator
first item: The requested CUDA packed-screening backend is unavailable: CUDA screening is disabled by SMSD_DISABLE_GPU
auto, 70,001 bits: [(1, 1.0)]
cuda, 70,001 bits: CUDA packed screening supports at most 65536 bits
gpu_prefilter=1: GPU prefilter flag must be a bool
```

The `RuntimeError` message gives the cause after the colon:
`No module named 'smsd_cuda'` when the package is missing,
`smsd and smsd-cuda must have the same version` for a version mismatch,
`no CUDA device was found`, or a reason that names a device below compute
capability 7.5. A `backend` value other than `"cpu"`, `"cuda"` or `"auto"`
raises `ValueError`.

If the installed `smsd-cuda` serves another version of `smsd`, the first call
that looks for it issues one `RuntimeWarning`. The warning names both versions,
asks for `smsd-cuda==8.1.0` and says that screening runs on the CPU; prepared
screens and pre-screens then stay on the CPU.

## Prove that the GPU ran

`gpu_is_available()` and `gpu_device_info()` report what is available, not
what any call used. The only public evidence that work ran on the device is
`FingerprintCorpus.cuda_calls`: the number of packed CUDA calls, including
clustering passes, that the corpus has completed. Read it before and after a
call, as the packed-index example does; an increase proves that the device
ran. It is not a count of queries, hits or kernels, and it says nothing about
speed.

Prepared screens and `gpu_prefilter=True` searches have no public counter and
do not change `cuda_calls`. Their results are identical on either device, so
compare results rather than devices.

## Time and compare CPU and GPU

Timings depend on the device, the processor, the corpus, the thread count and
the load on the machine. Measure on the hardware and data that matter, and
treat any figure as specific to that set-up. A fair comparison:

- calls `smsd_cuda.preload()` first, so that device start-up is not timed;
- prepares the corpus and runs each call once before timing, so that the
  first transfer to the device is not timed;
- times with `time.perf_counter()` and keeps the best of several runs;
- asserts that the GPU results equal the CPU results.

The harness below writes its timings to standard error and prints only
results that agree on every machine. It holds the prepared screen on the CPU
with `SMSD_DISABLE_GPU`, then allows the GPU for the 2,000-target corpus with
`SMSD_CUDA_MIN_TARGETS`.

```python
import os
import sys
import time
from pathlib import Path
from tempfile import TemporaryDirectory

import smsd
from smsd.corpus import FingerprintCorpus

try:
    import smsd_cuda
except ImportError:
    smsd_cuda = None


def cuda_status():
    """Return whether backend="cuda" can be used here, and the reason when it cannot."""
    if os.environ.get("SMSD_DISABLE_GPU", "") not in ("", "0"):
        return False, "disabled by SMSD_DISABLE_GPU"
    if smsd_cuda is None:
        return False, "smsd-cuda is not installed"
    if smsd_cuda.__version__ != smsd.__version__:
        return False, (f"smsd-cuda {smsd_cuda.__version__} does not match "
                       f"smsd {smsd.__version__}")
    reason = smsd_cuda.unavailable_reason()
    return reason is None, reason


def timed(label, work, repeats=3):
    """Run once to warm up, then write the best of several timed runs to stderr."""
    result = work()
    best = float("inf")
    for _ in range(repeats):
        started = time.perf_counter()
        assert work() == result
        best = min(best, time.perf_counter() - started)
    print(f"{label}: {best:.4f} s", file=sys.stderr)
    return result


# Substituted rings with chains of 0 to 49 carbons: 5 x 8 x 50 = 2,000 molecules.
rings = ["c1ccccc1", "c1ccncc1", "C1CCCCC1", "c1ccoc1", "c1ccsc1"]
groups = ["O", "N", "Cl", "C(=O)O", "OC", "C#N", "F", "CC"]
library = [ring + group + "C" * length
           for ring in rings for group in groups for length in range(50)]
queries = ["c1ccccc1O", "c1ccncc1Cl", "C1CCCCC1C(=O)O"]
ready, reason = cuda_status()
if ready:
    smsd_cuda.preload()  # keep device start-up out of the timings

corpus = smsd.TargetCorpus.from_smiles(library)
corpus.prewarm()
os.environ["SMSD_DISABLE_GPU"] = "1"
cpu_hits = timed("prepared screen, CPU", lambda: corpus.screen_many(queries, 0.4))
del os.environ["SMSD_DISABLE_GPU"]
if ready:
    os.environ["SMSD_CUDA_MIN_TARGETS"] = str(len(library))
    gpu_hits = timed("prepared screen, GPU allowed",
                     lambda: corpus.screen_many(queries, 0.4))
    assert gpu_hits == cpu_hits
counts = [len(found) for found in cpu_hits]
print("targets:", len(library), "| hits per query:", counts)
assert len(library) == 2000 and counts == [2, 9, 2]

with TemporaryDirectory() as directory:
    path = Path(directory) / "library.fp"
    with FingerprintCorpus.from_smiles(path, library) as packed:
        probes = [packed.fingerprint(index) for index in (0, 450, 1600)]
        cpu_rows = timed("packed scan, CPU", lambda: list(
            packed.screen_many(probes, threshold=0.4, backend="cpu")))
        if ready:
            before = packed.cuda_calls
            gpu_rows = timed("packed scan, GPU", lambda: list(
                packed.screen_many(probes, threshold=0.4, backend="cuda")))
            assert gpu_rows == cpu_rows and packed.cuda_calls > before
        print("packed hits:", len(cpu_rows))
        assert len(cpu_rows) == 12
```

The program prints:

```text
targets: 2000 | hits per query: [2, 9, 2]
packed hits: 12
```

A prepared screen gives no public record of the device it used, so the
harness labels its second timing "GPU allowed". The packed timing is backed by
the `cuda_calls` assertion.

## Containers

The GPU image `smsd:8.1.0-cuda12.8` (linux/amd64) contains `smsd`, the
matching `smsd-cuda` and the CUDA 12 runtime, with Python as its entry point.
It has no Java command-line tool, no CUDA toolkit and no NumPy. As shipped,
it supports prepared screens, `warm_start`, the stream functions and
`match_corpus(..., gpu_prefilter=True)`; packed indexes and clustering need
NumPy added in an image built from it. The image needs the NVIDIA Container
Toolkit and a compatible host driver, including under WSL2 on Windows.

After loading the image archive from the release:

```bash
docker run --rm --gpus all smsd:8.1.0-cuda12.8 \
  -c "import smsd; print(smsd.__version__); print(smsd.gpu_device_info())"
```

As elsewhere, the backend string shows availability, not use.
[Command line and containers](COMMAND_LINE.md#docker-images) describes the
CPU images, and the [installation guide](../HOWTO-INSTALL.md#docker) covers
loading and building every image.

## Metal and C++

Python never screens on Metal, and the published wheels include neither
Metal nor native CUDA support. A Python source build with Metal or CUDA
switched on changes only what `gpu_is_available()` and `gpu_device_info()`
report. In C++, a program built from source can link the Metal library and
request Metal fingerprint screening explicitly; see
[Optional Metal fingerprint screening](../CPP.md#optional-metal-fingerprint-screening).
Matching stays on the CPU there as well.

## Options

| Name | Used by | Type | Default | Meaning |
| --- | --- | --- | --- | --- |
| `gpu_prefilter` | `Substructure.match_corpus`, `MCS.match_corpus` | `bool` | `False` | `True` asks an attached device to run the pre-screen; a value that is not a `bool` raises `TypeError` |
| `backend` | `FingerprintCorpus.similarities`, `screen_many`, `threshold_pairs`, `top_k_pairs`; `butina`, `maxmin` | `str` | `"cpu"` | `"cpu"`, `"cuda"` or `"auto"`; any other value raises `ValueError` |
| `num_threads` | `TargetCorpus.prewarm`, `screen`, `screen_many`; `warm_start`; the stream functions; `match_corpus` | `int` | `0` | CPU threads; `0` uses the OpenMP default |
| `threshold` | `warm_start` | `float` | `0.0` | Inclusive ECFP4 Tanimoto threshold; the default returns every target, so pass one |

| Function or attribute | Returns | Meaning |
| --- | --- | --- |
| `smsd.gpu_is_available()` | `bool` | Whether a GPU backend is available |
| `smsd.gpu_device_info()` | `str` | A description of the backend, such as `CPU: ...` or `GPU: ...` |
| `smsd.warm_start(corpus, query, threshold=0.0, *, num_threads=0)` | `list[int]` | Prepares the corpus when needed and screens one query |
| `FingerprintCorpus.cuda_calls` | `int` | Packed CUDA calls completed by this corpus |
| `smsd_cuda.is_available()` | `bool` | Whether a device of compute capability 7.5 or later is present |
| `smsd_cuda.unavailable_reason()` | `str` or `None` | Why no device can be used, or `None` |
| `smsd_cuda.device_info()` | `str` | Device name, compute capability and memory; raises `RuntimeError` without a usable device |
| `smsd_cuda.preload()` | `None` | Initialises the device and loads the screening kernels, without uploading data or changing automatic use |
| `smsd_cuda.__version__` | `str` | Must equal `smsd.__version__` |

| Limit | Value |
| --- | --- |
| Compute capability | 7.5 or later |
| Automatic device rule | A Tesla T4 and 100,000 targets, unless `SMSD_CUDA_MIN_TARGETS` is set |
| Packed fingerprint width on the GPU | 65,536 bits; wider indexes use the CPU under `"auto"` and raise `ValueError` under `"cuda"` |
| Queries in one packed `screen_many` call | 32, on either device |
| Byte order for packed GPU scans | Little-endian host |

## Environment variables

| Variable | Default | Effect |
| --- | --- | --- |
| `SMSD_DISABLE_GPU` | Unset | Unset, empty or `0` allows the GPU. Any other value keeps every route on the CPU and makes `backend="cuda"` raise. Read at each call. |
| `SMSD_CUDA_MIN_TARGETS` | Unset | A positive whole number N, written in decimal digits only, lets automatic prepared screens use any supported device for corpora of N targets or more; on a Tesla T4 it replaces 100,000. Signs, spaces, exponents, `0` and values too large for 64 bits are ignored. Read at each screen. It does not affect `gpu_prefilter` or `backend`. |
| `OMP_NUM_THREADS` | Unset | The CPU threads used when `num_threads` is `0`; when unset, OpenMP chooses. |

`SMSD_CUDA_POOL` and `SMSD_CUDA_GRAPHS` tune the device workspace; see
[Retained workspace and graph replay](../../packaging/smsd-cuda/README.md#retained-workspace-and-graph-replay).

## Notes

- Graph matching never runs on the GPU. A GPU screens or pre-screens, and
  every reported match comes from a CPU search.
- Prepared screens need `prewarm()` after construction and after every
  addition; without it, `screen` and `screen_many` return no hits rather than
  raising. `warm_start` prepares the corpus when needed. See
  [Screening a library](SCREENING.md#similarity-screen-in-memory).
- A prepared corpus keeps its fingerprints on the device between screens.
  Adding targets drops that copy, and the next screen transfers them again.
- Results do not depend on the device or the thread count, and thresholds are
  inclusive on both devices.
- `cuda_calls` counts packed calls only; it stays at 0 for CPU and fallback
  work.
- The automatic rule is narrow by default, so a small demonstration corpus
  stays on the CPU unless `SMSD_CUDA_MIN_TARGETS` is set.

## See also

- [Screening a library](SCREENING.md)
- [Packed fingerprint index](PACKED_INDEX.md)
- [Command line and containers](COMMAND_LINE.md)
- [Optional CUDA screening](../PYTHON.md#optional-cuda-screening) in the
  Python guide
- [CPU and CUDA backends](../PACKED_CORPORA.md#cpu-and-cuda-backends) for
  packed corpora
- [The smsd-cuda package](../../packaging/smsd-cuda/README.md)
- [GPU screening example](../../python/examples/gpu_screening.py) and the
  [CPU and GPU screening notebook](../../python/examples/cpu_gpu_screening.ipynb)
- [Supported platforms](../HOWTO-INSTALL.md#supported-platforms)
