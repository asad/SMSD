#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Check the CPU image or normal CUDA fingerprint dispatch and ordered CPU hits."""
import argparse
import ctypes
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import struct
import sys

SMILES = ("CCO", "CCCO", "c1ccccc1", "CC(=O)O", "CN", "CCN")
CAPSULE_NAME = b"smsd.screen_accelerator"


class AutomaticDispatchDeferred(RuntimeError):
    pass


class Hit(ctypes.Structure):
    _fields_ = [("target_index", ctypes.c_int32), ("common_bits", ctypes.c_int32)]


FingerprintScreen = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.c_int32,
    ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_int32), ctypes.c_int32,
    ctypes.POINTER(Hit), ctypes.c_uint64, ctypes.POINTER(ctypes.c_uint64), ctypes.c_void_p, ctypes.c_size_t)
MinimumTargets = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_uint32, ctypes.POINTER(ctypes.c_uint64), ctypes.c_void_p, ctypes.c_size_t)


class AcceleratorTable(ctypes.Structure):
    _fields_ = [("abi_version", ctypes.c_uint32), ("struct_size", ctypes.c_uint32),
        ("smsd_version", ctypes.c_char_p), ("name", ctypes.c_char_p),
        *[(name, ctypes.c_void_p) for name in ("device_info", "label_overlaps", "minimum_targets",
          "create_fingerprints", "fingerprint_screen", "create_labels", "resident_label_overlaps", "destroy_resident")]]


class CallbackObserver:
    def __init__(self, library):
        self.library = library
        pointer = library.smsd_cuda_screen_accelerator()
        if not pointer:
            raise RuntimeError("The installed accelerator returned no interface")
        original = AcceleratorTable.from_address(pointer)
        if (original.abi_version != 1 or original.struct_size < ctypes.sizeof(AcceleratorTable)
                or not all(getattr(original, name) for name in ("minimum_targets", "create_fingerprints", "fingerprint_screen", "destroy_resident"))):
            raise RuntimeError("The installed accelerator lacks the resident interface")
        self.table = AcceleratorTable.from_buffer_copy(ctypes.string_at(pointer, ctypes.sizeof(AcceleratorTable)))
        self.table.struct_size = ctypes.sizeof(self.table)
        self.requests = self.successes = 0
        self.original_screen = FingerprintScreen(original.fingerprint_screen)
        def screen(*arguments):
            self.requests += 1
            status = self.original_screen(*arguments)
            self.successes += status == 0
            return status
        self.callback = FingerprintScreen(screen)
        self.table.fingerprint_screen = ctypes.cast(self.callback, ctypes.c_void_p).value
        self.minimum = ctypes.c_uint64()
        error = ctypes.create_string_buffer(256)
        if MinimumTargets(original.minimum_targets)(1, ctypes.byref(self.minimum), error, len(error)):
            raise RuntimeError(error.value.decode(errors="replace"))

    def capsule(self):
        new = ctypes.pythonapi.PyCapsule_New
        new.restype = ctypes.py_object
        new.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p]
        return new(ctypes.addressof(self.table), CAPSULE_NAME, None)


def mapped_runtime():
    spec = importlib.util.find_spec("nvidia.cuda_runtime")
    expected = {str((Path(folder) / "lib/libcudart.so.12").resolve()) for folder in spec.submodule_search_locations}
    paths = {str(Path(line.split(maxsplit=5)[5]).resolve()) for line in Path("/proc/self/maps").read_text().splitlines()
             if "libcudart.so" in line and len(line.split(maxsplit=5)) == 6}
    if not paths or not paths <= expected:
        raise RuntimeError("The mapped CUDA runtime does not come from the pinned Python runtime package")
    result = []
    for name in sorted(paths):
        library = ctypes.CDLL(name)
        library.cudaRuntimeGetVersion.argtypes = [ctypes.POINTER(ctypes.c_int)]
        library.cudaRuntimeGetVersion.restype = ctypes.c_int
        version = ctypes.c_int()
        if library.cudaRuntimeGetVersion(ctypes.byref(version)) or version.value != 12080:
            raise RuntimeError("Expected an actual CUDA12.8 runtime")
        result.append({"path": name, "sha256": hashlib.sha256(Path(name).read_bytes()).hexdigest(), "version": version.value})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--cpu-reference", type=Path)
    parser.add_argument("--targets", type=int, default=100000)
    args = parser.parse_args()
    if not 100000 <= args.targets <= 1000000:
        raise ValueError("Use between 100,000 and 1,000,000 targets")
    if args.gpu:
        os.environ.pop("SMSD_DISABLE_GPU", None)
    else:
        os.environ["SMSD_DISABLE_GPU"] = "1"
    os.environ["OMP_NUM_THREADS"] = "4"
    import smsd
    from smsd import _accelerator, _smsd
    if smsd.__version__ != "8.1.0" or sys.version_info[:2] != (3, 14):
        raise RuntimeError("Expected SMSD 8.1.0 on CPython 3.14")
    artifacts = json.loads(Path("/usr/share/doc/smsd/ARTIFACTS.json").read_text())
    architecture = "amd64" if platform.machine() == "x86_64" else "arm64" if platform.machine() == "aarch64" else None
    if architecture is None:
        raise RuntimeError("Unexpected image architecture")
    mcs = smsd.MCS.find("CCO").match("CCCO")
    substructure = smsd.Substructure.find("CO").match("CCO")
    if not mcs.optimal or not mcs.complete or mcs.best.atom_count != 3 or not substructure.found or not substructure.complete:
        raise RuntimeError("MCS or substructure smoke failed")
    corpus = smsd.TargetCorpus(); observer = None
    try:
        device = "CPU"
        if args.gpu:
            import smsd_cuda
            if smsd_cuda.__version__ != "8.1.0" or not smsd_cuda.is_available():
                raise RuntimeError("The CUDA package is unavailable: " + str(smsd_cuda.unavailable_reason()))
            if not smsd.gpu_is_available():
                raise RuntimeError("SMSD did not attach the CUDA provider")
            device = smsd.gpu_device_info()
            if "T4" not in device or "compute 7.5" not in device:
                raise RuntimeError("The release smoke needs the declared T4")
            observer = CallbackObserver(smsd_cuda._library)
            if observer.minimum.value == (1 << 64)-1 or observer.minimum.value > args.targets:
                raise AutomaticDispatchDeferred("Automatic fingerprint dispatch is deferred for this target count/device")
            # Count the real callback while preserving the provider's minimum and all other callbacks.
            _smsd._attach_screen_accelerator(observer.capsule())
            if importlib.metadata.version("nvidia-cuda-runtime-cu12") != "12.8.90":
                raise RuntimeError("Expected the pinned CUDA12.8.90 package")
            if shutil.which("nvcc"):
                raise RuntimeError("The final image contains a CUDA compiler")
        graphs = [smsd.parse_smiles(text) for text in SMILES]
        corpus.add_targets([graphs[index % len(graphs)] for index in range(args.targets)])
        corpus.prewarm(num_threads=0)
        failures = _smsd._screen_accelerator_failures()[0]
        hits = corpus.screen(graphs[0], .8, num_threads=0)
        if (failures != 0 or _smsd._screen_accelerator_failures()[0] != 0
                or not hits or hits != sorted(set(hits)) or hits[0] < 0 or hits[-1] >= args.targets):
            raise RuntimeError("Screen output or fallback state is invalid")
        if args.gpu and (observer.requests < 1 or observer.requests != observer.successes):
            raise RuntimeError("No successful normal GPU fingerprint callback was observed")
        if not args.gpu and (_accelerator.active() is not None or _smsd._screen_accelerator() is not None):
            raise RuntimeError("The CPU smoke unexpectedly attached an accelerator")
        raw = b"".join(struct.pack("<i", index) for index in hits)
        result = {"status": "passed", "version": smsd.__version__, "python": sys.version,
            "source_commit": artifacts["source_commit"], "architecture": architecture,
            "core_wheel_sha256": artifacts["cpu"][architecture]["sha256"],
            "inputs": list(SMILES), "targets": args.targets, "threshold": .8,
            "fingerprint_hits": hits, "fingerprint_digest": hashlib.sha256(struct.pack("<Q", len(hits))+raw).hexdigest(),
            "mcs_atoms": mcs.best.atom_count, "substructure_found": substructure.found, "device": device,
            "gpu_requests": observer.requests if observer else 0, "fallback_failures": 0}
        if args.gpu:
            result["mapped_runtime"] = mapped_runtime()
            if args.cpu_reference is None:
                raise ValueError("The GPU smoke requires the CPU image's retained reference")
            reference = json.loads(args.cpu_reference.read_text())
            for name in ("version", "source_commit", "core_wheel_sha256", "inputs", "targets", "threshold", "fingerprint_hits", "fingerprint_digest"):
                if result[name] != reference[name]:
                    raise RuntimeError("CPU image reference differs: " + name)
            if (reference.get("status") != "passed" or reference.get("gpu_requests") != 0
                    or reference.get("fallback_failures") != 0 or reference.get("device") != "CPU"
                    or reference.get("architecture") != "amd64"):
                raise RuntimeError("Expected a completed CPU image reference")
        corpus.add_targets([])
        corpus = None
        if _smsd._screen_accelerator_failures()[0] != 0:
            raise RuntimeError("Resident cleanup failed")
        print(json.dumps(result, sort_keys=True))
    finally:
        if corpus is not None:
            corpus.add_targets([])
        corpus = None
        _accelerator._reset()


if __name__ == "__main__":
    try:
        main()
    except AutomaticDispatchDeferred as error:
        print(json.dumps({"status": "deferred", "error": str(error), "public_gpu_dispatch": False}))
        raise SystemExit(77) from error
    except Exception as error:
        print(json.dumps({"status": "failed", "error": str(error)}))
        raise SystemExit(1) from error
