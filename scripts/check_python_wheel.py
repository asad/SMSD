#!/usr/bin/env python3
"""Check that native functionality comes from the installed CPU wheel."""

import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import tomllib


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--require-rdkit", action="store_true")
    parser.add_argument("--require-openmp", action="store_true")
    args = parser.parse_args()
    source = args.source.resolve()
    expected = tomllib.loads((source / "pyproject.toml").read_text())["project"]["version"]

    import smsd

    prefix = Path(sys.prefix).resolve()
    package = Path(smsd.__file__).resolve()
    native = Path(smsd._smsd.__file__).resolve()
    if not package.is_relative_to(prefix) or not native.is_relative_to(package.parent):
        raise SystemExit("Expected package and native extension from the installed wheel")
    if package.is_relative_to(source / "python"):
        raise SystemExit("Imported the source tree instead of the installed wheel")
    if importlib.metadata.version("smsd") != expected or smsd.__version__ != expected:
        raise SystemExit("Installed wheel and source versions differ")
    backend = smsd.gpu_device_info()
    if smsd.gpu_is_available() or not backend.startswith("CPU:") or "compiled" in backend:
        raise SystemExit("Expected a CPU-only wheel")
    if args.require_openmp and ("OpenMP" not in backend or "no OpenMP" in backend):
        raise SystemExit("OpenMP was requested but is absent from this wheel")

    query = smsd.parse_smiles("c1ccccc1")
    targets = [smsd.parse_smiles(value) for value in ("Oc1ccccc1", "CCO")]
    assert len(query) == 6
    assert len(smsd.find_mcs(query, targets[0], strategy="native", timeout_ms=1000)) == 6
    assert [len(value) for value in smsd.batch_find_substructure(query, targets)] == [6, 0]
    assert smsd.smarts_match("[#8]", targets[0])
    dependencies = {"pytest": importlib.metadata.version("pytest")}
    if args.require_rdkit:
        from rdkit import Chem, rdBase
        dependencies["rdkit"] = rdBase.rdkitVersion
        rdkit_query = Chem.MolFromSmiles("c1ccccc1")
        rdkit_target = Chem.MolFromSmiles("Oc1ccccc1")
        assert len(smsd.find_mcs(rdkit_query, rdkit_target, strategy="native", timeout_ms=1000)) == 6
    print(json.dumps({"version": expected, "python": platform.python_version(),
                      "system": platform.system(), "architecture": platform.machine(),
                      "backend": backend, "dependencies": dependencies,
                      "installed_package": "smsd", "native_extension": native.name}, indent=2))


if __name__ == "__main__":
    main()
