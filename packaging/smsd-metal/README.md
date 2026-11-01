# smsd-metal

`smsd-metal` provides optional Apple Metal screening on Apple silicon Macs and
needs the same version of `smsd`:

```bash
python -m pip install "smsd==8.2.0" "smsd-metal==8.2.0"
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
`TargetCorpus` of at least 100,000 targets uses Metal on an Apple M5. Smaller
corpora and other Apple GPUs use the CPU unless `SMSD_METAL_MIN_TARGETS` is set
to a positive integer N: automatic screening then uses the Apple GPU for corpora
of at least N targets; on an Apple M5, N replaces 100,000.
Automatic similarity-bound screening (`screen_targets`) and substructure
pre-screening use the CPU; the explicit `match_corpus(..., gpu_prefilter=True)`
methods ask the GPU to run their pre-screen. `gpu_is_available()` and
`gpu_device_info()` report backend availability, not whether a call used the
device. Without a usable device, or when device work fails, screening runs on
the CPU with the same hits. `SMSD_DISABLE_GPU=1` disables the GPU at run time.

Corpus substructure searches can rule out targets that cannot contain the query
before CPU matching. With SMARTS queries, `match_atom_type=False` or a node
limit, every target is searched. Results keep their target order, atom
numbering and search status.

Use `screen_targets` with an explicit threshold before MCS when an atom-count
Tanimoto bound is wanted; this bound does not represent bond Tanimoto or weighted
MCS scores. ECFP similarity is not a substructure test.

Packed fingerprint corpora do not use Metal: `backend="auto"` runs on the CPU,
and `backend="cuda"` needs `smsd-cuda`.

## Requirements

macOS 14 or later on Apple silicon, with CPython 3.12–3.15. The package
uses the Metal framework that macOS provides and bundles no other library.
Intel Macs are not supported.

## Building

Build from an SMSD checkout on an Apple silicon Mac with the Xcode command-line
tools:

```bash
MACOSX_DEPLOYMENT_TARGET=14.0 python -m pip wheel packaging/smsd-metal --no-deps -w dist
```

Its version must equal the `smsd` version of the checkout; the build stops
otherwise. The `smsd` package itself never needs Metal.

## Explicit initialisation

`import smsd_metal; smsd_metal.preload()` compiles the screening kernels before
any data is screened. It does not upload a corpus or change automatic dispatch.

## Licence

SMSD code is under the Apache License 2.0; see `LICENSE` and `NOTICE`.
