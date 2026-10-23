# smsd-cuda

`smsd-cuda` provides optional NVIDIA CUDA screening on Linux x86_64 and needs
the same version of `smsd`:

```bash
python -m pip install "smsd==8.0.2" "smsd-cuda==8.0.2"
```

The GPU computes integer overlaps; scores use the CPU expressions and hits
keep target order. Matching runs on the CPU. The backend keeps prepared
corpus fingerprints on the device for reuse. ECFP4 Tanimoto
screening and the atom-count MCS similarity bound are distinct operations.

```python
import smsd

corpus = smsd.TargetCorpus.from_smiles(["CCO", "CCN", "CCC"])
print(smsd.gpu_device_info())
corpus.prewarm(num_threads=1)
assert corpus.screen_many(["CCO", "CCN"], 1.0, num_threads=1) == [[0], [1]]
assert corpus.screen_targets_many(["CCO", "CCN"], 1.0, num_threads=1) == [[0], [1]]
```

`smsd` finds the package at run time. Automatic ECFP4 screening uses CUDA on a
Tesla T4 with compute capability 7.5 when the prepared corpus has at least
100,000 targets. Smaller corpora and other devices use the CPU. Label-bound
and substructure candidate screening use the CPU. `gpu_is_available()` and
`gpu_device_info()` report backend availability, not whether a call used the
device. Missing devices, allocation failure and device loss use the CPU.
`SMSD_DISABLE_GPU=1` disables the GPU at run time.

Corpus substructure searches can reject impossible molecule candidates with a
conservative element-path screen before CPU matching. SMARTS queries, untyped
queries and finite node limits bypass this screen. Target results keep their
original slots, atom numbering and search status.

Use `screen_targets` with an explicit threshold before MCS when an atom-count
Tanimoto bound is wanted; this bound does not represent bond Tanimoto or weighted
MCS scores. ECFP similarity is not a substructure test.

## Requirements

See the
[supported-platforms table](https://github.com/asad/SMSD/blob/master/docs/HOWTO-INSTALL.md#supported-platforms)
for Python versions, operating-system requirements and GPU Docker support.
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

## Licence

Apache License 2.0; see `LICENSE` and `NOTICE`. The package bundles no NVIDIA library.
