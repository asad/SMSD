# GPU screening on Windows

SMSD Pro provides no native Windows CUDA package. The Windows `smsd` wheel
searches and screens on the CPU, and an NVIDIA GPU in a Windows machine is
reached through WSL2, the Linux environment built into Windows, in one of two
ways:

- install the Linux `smsd-cuda` wheel next to `smsd` in a WSL2 distribution;
- run the GPU Docker image of the release.

Both use the NVIDIA driver installed in Windows, and neither needs the CUDA
toolkit. Graph matching stays on the CPU, and hits, scores and search
decisions are the same as on the CPU; [CPU and GPU screening](CPU_GPU.md)
describes which steps can use the GPU.

The PowerShell commands on this page run in Windows. The shell commands run in
the WSL2 distribution, where they are the same as on any Linux x86_64 machine.

## Choose a route

| Route | Use it when | Contents |
| --- | --- | --- |
| [Linux wheels in WSL2](#install-the-linux-wheels-in-wsl2) | You work in your own Python environment, with scripts, notebooks or extras such as NumPy | `smsd`, `smsd-cuda` and whatever else you install |
| [GPU Docker image](#run-the-gpu-docker-image) | You want a fixed environment that matches the release | `smsd`, `smsd-cuda` and the CUDA runtime, without NumPy |

## Requirements

- Windows 11, or Windows 10 version 21H2 or later, on x86_64.
- An NVIDIA GPU of compute capability 7.5 or later.
- A current NVIDIA driver for Windows. It serves WSL2 as well, so do not
  install an NVIDIA driver for Linux in the distribution.
- WSL2. The first version of WSL cannot reach the GPU.

## Set up WSL2

In PowerShell opened as administrator, install WSL2 with an Ubuntu 24.04
distribution, and restart Windows if asked:

```powershell
wsl --install -d Ubuntu-24.04
```

Where WSL is already installed, update it and check that the distribution
runs as version 2:

```powershell
wsl --update
wsl --list --verbose
```

`wsl --set-version Ubuntu-24.04 2` converts a distribution listed as
version 1. Then open the Ubuntu shell and check that the GPU is visible:

```bash
nvidia-smi
```

The command lists the GPU and the driver version. If it finds no device,
update the Windows driver, run `wsl --shutdown` in PowerShell and open the
shell again.

The wheels need CPython 3.12–3.15 and glibc 2.34 or later. Ubuntu 24.04
provides Python 3.12 and glibc 2.39. Ubuntu 22.04 has a suitable glibc but
Python 3.10, so it needs a newer Python first.

## Install the Linux wheels in WSL2

In the Ubuntu shell, create a virtual environment in your home directory and
install the two packages at the same version:

```bash
sudo apt update
sudo apt install -y python3-venv
python3 -m venv ~/smsd-gpu
source ~/smsd-gpu/bin/activate
python -m pip install "smsd==8.2.1" "smsd-cuda==8.2.1"
```

For packed indexes and clustering, install `"smsd[numpy]==8.2.1"` in place of
`"smsd==8.2.1"`. Then check the device and the backend:

```bash
python -c "import smsd_cuda; print(smsd_cuda.unavailable_reason() or smsd_cuda.device_info())"
python -c "import smsd; print(smsd.__version__); print(smsd.gpu_device_info())"
```

The first command prints the device name, its compute capability and its
memory, or the reason no device can be used. The second prints `8.2.1` and a
backend string that starts with `GPU:` once `smsd` has attached the device.
Both report what is available, not which device a call used.

Automatic screens of a prepared `TargetCorpus` use the GPU by default only on a
Tesla T4, and only for corpora of at least 100,000 targets, so on most Windows
machines they stay on the CPU until `SMSD_CUDA_MIN_TARGETS` is set. For
example, to apply the same minimum of 100,000 targets on any supported device:

```bash
export SMSD_CUDA_MIN_TARGETS=100000
```

[Time and compare CPU and GPU](CPU_GPU.md#time-and-compare-cpu-and-gpu) shows
how to find the corpus size from which the GPU is faster on your machine. The
requested routes, `match_corpus(..., gpu_prefilter=True)` and the packed
backend `"cuda"`, use any supported device without the variable.

Keep libraries and notebooks in the Linux file system, such as your home
directory in the distribution. Windows drives appear under `/mnt/c` and can be
read from there, but more slowly.

## Run the GPU Docker image

The [8.2.1 release](https://github.com/asad/SMSD/releases/tag/v8.2.1) has the
GPU image `smsd-8.2.1-cuda-linux-amd64.tar.xz`, with Python, `smsd`, the
matching `smsd-cuda` and the CUDA 12 runtime. Python is its entry point and
`/work` its working directory. Docker needs GPU access in WSL2, which either
of two set-ups provides:

- Docker Desktop for Windows with its WSL2 back end, which passes the GPU to
  containers with the Windows driver alone. Turn on its WSL integration for
  the Ubuntu distribution to use `docker` from the shell.
- Docker Engine installed in the distribution, with the NVIDIA Container
  Toolkit.

In the Ubuntu shell, in the directory that holds the downloaded files, check
the archive, load it and run the image:

```bash
grep ' smsd-8.2.1-cuda-linux-amd64.tar.xz$' SHA256SUMS | shasum -a 256 -c
docker load --input smsd-8.2.1-cuda-linux-amd64.tar.xz
docker run --rm --gpus all smsd:8.2.1-cuda12.8 \
  -c "import smsd; print(smsd.__version__); print(smsd.gpu_device_info())"
```

The last command prints `8.2.1` and a backend string that starts with `GPU:`.
To run a script from the current directory, mount the directory read-only at
`/work` and pass settings with `-e`:

```bash
docker run --rm --gpus all -v "$PWD:/work:ro" \
  -e SMSD_CUDA_MIN_TARGETS=100000 smsd:8.2.1-cuda12.8 screen.py
```

The image runs as an unprivileged user, so the script should write its results
to standard output, or to a separately mounted directory that this user can
write to. Packed indexes and clustering need NumPy, which the image lacks; add
it in an image built from this one. The [Docker section](../HOWTO-INSTALL.md#docker)
of the installation guide shows how to load the archive when the Docker engine
has no xz decompressor.

## Notes

- The Windows `smsd` wheel and the WSL2 environment are separate Python
  installations. `smsd-cuda` has no Windows wheel, so install it only in
  WSL2.
- `smsd` attaches only an `smsd-cuda` of the same version. With any other
  version it warns once and screens on the CPU.
- `SMSD_DISABLE_GPU=1` keeps every route on the CPU, as on Linux.
- Graph matching never runs on the GPU, and every reported match comes from a
  CPU search.

## See also

- [CPU and GPU screening](CPU_GPU.md)
- [Supported platforms](../HOWTO-INSTALL.md#supported-platforms) and
  [Docker images](../HOWTO-INSTALL.md#docker) in the installation guide
- [The smsd-cuda package](../../packaging/smsd-cuda/README.md)
- Microsoft's guide to [installing WSL](https://learn.microsoft.com/windows/wsl/install)
  and NVIDIA's [CUDA on WSL user guide](https://docs.nvidia.com/cuda/wsl-user-guide/index.html)
