# smsd-cuda

`smsd-cuda` provides optional NVIDIA CUDA screening on Linux x86_64 and needs
the same version of `smsd`:

```bash
python -m pip install "smsd==8.1.0" "smsd-cuda==8.1.0"
```

The GPU computes integer overlaps; scores use the CPU expressions and hits
keep target order. Matching runs on the CPU. The backend keeps prepared
corpus fingerprints on the device for reuse. ECFP4 Tanimoto
screening and the atom-count MCS similarity bound are distinct operations.

```python
import smsd

# Three targets: with the default settings, far too few for automatic GPU use,
# so this runs on the CPU.
corpus = smsd.TargetCorpus.from_smiles(["CCO", "CCN", "CCC"])
print(smsd.gpu_device_info())
corpus.prewarm(num_threads=1)
assert corpus.screen_many(["CCO", "CCN"], 1.0, num_threads=1) == [[0], [1]]
assert corpus.screen_targets_many(["CCO", "CCN"], 1.0, num_threads=1) == [[0], [1]]
```

`smsd` finds the package at run time. Automatic ECFP4 screening of a prepared
`TargetCorpus` of at least 100,000 targets uses CUDA on a Tesla T4 with compute
capability 7.5. Smaller corpora and other devices use the CPU unless
`SMSD_CUDA_MIN_TARGETS` is set to a positive integer N: automatic screening then
uses any NVIDIA GPU with compute capability 7.5 or later for corpora of at
least N targets; on a Tesla T4, N replaces 100,000. The
[GPU screening example](https://github.com/asad/SMSD/blob/master/python/examples/gpu_screening.py)
shows the explicit device routes.
Automatic similarity-bound screening (`screen_targets`) and substructure
pre-screening use the CPU; the explicit `match_corpus(..., gpu_prefilter=True)`
methods ask the GPU to run their pre-screen. `gpu_is_available()` and
`gpu_device_info()` report backend availability, not whether a call used the
device. Missing devices, allocation failure and device loss use the CPU.
`SMSD_DISABLE_GPU=1` disables the GPU at run time.

Corpus substructure searches can rule out targets that cannot contain the query
before CPU matching. With SMARTS queries, untyped queries or a finite node
limit, every target is searched. Target results keep their original slots, atom
numbering and search status.

Use `screen_targets` with an explicit threshold before MCS when an atom-count
Tanimoto bound is wanted; this bound does not represent bond Tanimoto or weighted
MCS scores. ECFP similarity is not a substructure test.

## Requirements

See the
[supported-platforms table](https://github.com/asad/SMSD/blob/master/docs/HOWTO-INSTALL.md#supported-platforms)
for Python versions and operating-system requirements.
The `nvidia-cuda-runtime-cu12` dependency supplies the CUDA runtime; the host
needs a compatible NVIDIA driver. The CUDA toolkit is not needed to install
the wheel.

## Building

Build from an SMSD checkout on Linux x86_64 with the CUDA 12.8 toolkit:

```bash
python -m pip wheel packaging/smsd-cuda --no-deps -w dist
```

Its version must equal the `smsd` version of the checkout; the build stops
otherwise. Building requires `nvcc`; the `smsd` package itself never needs CUDA.

## Retained workspace and graph replay

Fingerprint screening can keep up to 256 MiB of device workspace per prepared
corpus between calls, in addition to the stored corpus fingerprints.
Concurrent calls use their own temporary memory. `SMSD_CUDA_POOL=0` turns the
retained workspace off. The workspace is released with its corpus; finish
every call on a corpus before releasing it.

`SMSD_CUDA_GRAPHS=1` opts into CUDA Graph replay for screening calls of at
most 64 queries while the retained workspace is on; top-k calls are not
replayed. Errors and CPU fallback behave as without these options, and
neither option guarantees lower latency.

## Explicit initialisation and device code

`import smsd_cuda; smsd_cuda.preload()` initialises the current device and
loads the screening kernels before any data is screened. It does not upload a
corpus or change automatic dispatch.

CUDA 12.8 builds include device code for SM75, SM80, SM86, SM89, SM90, SM100
and SM120, with SM100 and SM120 PTX for forward compilation. Earlier CUDA 12
builds retain SM90 PTX.

## Top-k library export

The optional `smsd_cuda_fingerprint_top_k` library export returns at most `k`
hits per query, ranked on the device by double-precision Tanimoto score; equal
scores keep the lower target index first. Ordinary screening results and their
target order are unchanged. It does not provide device clustering or an
all-pairs matrix API.

## Licence

SMSD code is under the Apache License 2.0; see `LICENSE` and `NOTICE`. The
incorporated CUB and CCCL code keeps its own licences, whose texts ship with the
package. The package bundles no NVIDIA runtime library.
