# Conda packages

Two native recipes build from one frozen full-source archive:

| Recipe | Package | Intended targets |
| --- | --- | --- |
| `packaging/conda` | `smsd` | CPython 3.12–3.15; Linux x86-64/aarch64, macOS x86-64/arm64, Windows x86-64 |
| `packaging/conda/cuda` | `smsd-cuda` | CPython 3.12–3.15; Linux x86-64; CUDA 12.8 build variant |

These are recipe targets, not claims of completed builds or platform tests.
Render with the final source inputs before building.
The CPU recipe disables OpenMP, CUDA and Metal explicitly and needs no GPU.
Its compiler and C standard-library variants follow the
[conda-forge platform configuration](https://conda-forge.org/docs/maintainer/infrastructure/).
The macOS recipe targets macOS 11 or newer; package builds must verify that target.
CUDA uses GCC 14 and the Linux 2.17 sysroot with its CUDA 12.8 toolchain.
The CUDA package requires the **exact same SMSD version** and a conda CUDA
runtime. Native Windows CUDA, macOS CUDA and Linux aarch64 CUDA packages are
not produced by this recipe. Neither package is `noarch`.

## Freeze source inputs

There is no placeholder source hash and no default mutable tag. Rendering
requires the `smsd_version`, `smsd_source_url`, `smsd_source_sha256` and
`smsd_source_commit` variants. Generate them only after final source freeze:

```sh
python packaging/conda/configure_source.py \
  --archive /path/to/the/exact-downloaded-full-source.tar.gz \
  --url "https://github.com/asad/SMSD/archive/${RELEASE_COMMIT}.tar.gz" \
  --commit "$RELEASE_COMMIT" \
  --output smsd-conda-variants.yaml
```

`RELEASE_COMMIT` is the full 40-character release commit, not a branch or tag.
The local archive must be the **exact bytes served by that immutable URL**. The
configuration tool hashes those bytes, checks core/CUDA/native versions and
writes a new YAML-compatible JSON variant file. It does not download the URL or
compare the archive with Git, so check that the archive matches the commit
before building. Do not substitute a separately compressed `git archive` for the
downloaded URL payload: its compressed hash can differ.

The full archive must include `packaging/smsd-cuda`. The Python source
distribution deliberately excludes that directory and cannot supply both
recipes. Generate these inputs for the exact version and commit being built.
The configuration tool needs `tomli` when the build host runs Python 3.10;
Python 3.11 and later use `tomllib`.

## Render and build

Use native builders for each target; these commands do not establish cross
compilation. A conda-build environment and a conda-forge toolchain must
already be available. Select one Python variant initially:

```sh
conda render packaging/conda -m smsd-conda-variants.yaml --python 3.12
conda build packaging/conda -m smsd-conda-variants.yaml --python 3.12 \
  --output-folder /path/to/local-channel -c conda-forge
conda render packaging/conda/cuda -m smsd-conda-variants.yaml --python 3.12
conda build packaging/conda/cuda -m smsd-conda-variants.yaml --python 3.12 \
  --output-folder /path/to/local-channel -c /path/to/local-channel -c conda-forge
```

Build the CPU package first so the CUDA test environment can solve its exact
core dependency from the local channel. Then repeat the build for each further
Python version and platform. Both build scripts use the activated compiler,
Ninja with one worker, installed backend dependencies and pip's no-index,
no-dependency, no-build-isolation route. They check the source and package
versions before compiling. The recipes run small package smoke checks rather
than the full test suite.

The CUDA recipe's own `conda_build_config.yaml` selects `cuda-nvcc` 12.8 and
GNU C++ 14, within that toolkit's documented host-compiler range. Host
requirements use `cuda-version=12.8`, `cuda-cudart-dev` and `cuda-cccl` for CUB
headers. Runtime requirements are `cuda-cudart` and `cuda-version>=12.8,<13`.
Component packages have their own versions; select the toolkit through
`cuda-version`. Check the solved environment, the NVCC and host compiler
pairing, the device code and the linkage of each build.

### CUDA packaging adaptation

The CUDA build script makes one explicit metadata-only change to its extracted
source: it removes `nvidia-cuda-runtime-cu12` from Python `Requires-Dist`.
The conda recipe instead owns the native `cuda-cudart` dependency. It retains
the exact Python `smsd` dependency and does not alter the loader, CMake or
kernel sources. No fake NVIDIA Python distribution or runtime shim is made.
Conda-build's ordinary ELF relocation uses `lib/`; validate that the installed
library resolves libcudart from the environment's prefix without an external
`LD_LIBRARY_PATH`. The automatic CUDA smoke test checks the versions and the
loaded runtime path and version, not GPU execution.

A system NVIDIA driver is still required to execute CUDA. A GPU-free build
host may use the documented `CONDA_OVERRIDE_CUDA=12.8` **solver override** for
packaging checks; it does not emulate a device or validate GPU results. On a
GPU host, unset that override and require a device:

```sh
python -c "import smsd_cuda; assert smsd_cuda.is_available(), smsd_cuda.unavailable_reason(); print(smsd_cuda.device_info())"
python -m pytest /path/to/frozen-source/packaging/smsd-cuda/tests/test_resident_fingerprints.py -q
```

The explicit availability assertion prevents device absence being mistaken for
passing skipped GPU tests. Keep the test log, including skipped tests. Automatic
dispatch thresholds and CPU fallback are unchanged; availability alone does not
prove a normal screen ran on the GPU.

## Install the built outputs

After indexing the local build channel, use that channel first with strict
priority. These commands install local conda artefacts, not PyPI wheels:

```sh
conda index /path/to/local-channel
conda create -n smsd-cpu --strict-channel-priority \
  -c file:///path/to/local-channel -c conda-forge "python=3.12" "smsd=8.1.0"
mamba create -n smsd-gpu --strict-channel-priority \
  -c file:///path/to/local-channel -c conda-forge \
  "python=3.12" "smsd=8.1.0" "smsd-cuda=8.1.0" "cuda-version=12.8"
```

Use the frozen version in both package specifications if it changes. For Pixi:

```sh
pixi init --channel file:///path/to/local-channel --channel conda-forge \
  --platform linux-64 smsd-demo
pixi add --manifest-path smsd-demo/pixi.toml "python=3.12" "smsd=8.1.0"
pixi run --manifest-path smsd-demo/pixi.toml python -c "import smsd; print(smsd.__version__)"
```

On a compatible Linux NVIDIA host, add `"smsd-cuda=8.1.0"` and
`"cuda-version=12.8"` to the same Pixi environment. Use the native platform
for CPU installations. No lockfile is included; Pixi creates one when it
solves the environment. A separate Pixi build recipe is unnecessary: Pixi
consumes the same conda outputs.

## Install from the channel

Packages on the `bioinception` channel follow the release. Until the channel
carries a build for your platform, build from these recipes and install from a
local channel as shown above. Once it does, install SMSD from the
`bioinception` channel and the dependencies from conda-forge:

```sh
conda install -c bioinception -c conda-forge smsd
```

Where the channel also carries `smsd-cuda`, install it from the same channels on
a Linux x86-64 host with a compatible NVIDIA driver; it requires the same `smsd`
version. For mamba or Pixi, use the `bioinception` channel in place of the local
build channel in the commands above once it carries the version they name.

References: [conda-build metadata](https://docs.conda.io/projects/conda-build/en/stable/resources/define-metadata.html),
[CUDA component/version selection](https://github.com/conda-forge/cuda-feedstock/blob/main/recipe/doc/end_user_run_guide.md),
[cuda-version scheme](https://github.com/conda-forge/cuda-version-feedstock/blob/main/recipe/README.md),
and [Pixi CLI](https://pixi.prefix.dev/latest/reference/cli/pixi/init/).
