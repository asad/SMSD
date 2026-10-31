# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Optional feature arrays preserve the native fingerprint's positions and counts."""

import json
import importlib.util
import os
import subprocess
import sys
import textwrap
from types import SimpleNamespace

import pytest
import smsd
from smsd import adapters


def test_optional_modules_do_not_import_their_dependencies(tmp_path):
    source = ("import json, sys, smsd; "
              f"smsd.__path__.insert(0, {os.path.dirname(adapters.__file__)!r}); "
              "import smsd.adapters, smsd.mcp_server, smsd.highlight; "
              "print(json.dumps(sorted(set(sys.modules) & "
              "{'numpy', 'torch', 'sklearn', 'mcp'})))")
    output = subprocess.run([sys.executable, "-c", source], cwd=tmp_path,
                            text=True, capture_output=True, timeout=30, check=True)
    assert json.loads(output.stdout) == []


def test_numpy_features_equal_native_fingerprints():
    numpy = pytest.importorskip("numpy")
    molecules = ["CCO", "c1ccccc1O", "CC(=O)O"]
    native = smsd.fingerprints(molecules, kind="ecfp", fp_size=128)
    matrix = adapters.fingerprints_numpy(iter(molecules), kind="ecfp", fp_size=128)
    assert matrix.shape == (3, 128)
    assert matrix.dtype == numpy.uint8
    assert [numpy.flatnonzero(row).tolist() for row in matrix] == native
    assert adapters.fingerprints_numpy([], fp_size=128).shape == (0, 128)


def test_numpy_count_features_equal_native_counts():
    numpy = pytest.importorskip("numpy")
    options = smsd.FingerprintOptions(family="ecfp", size=64, counts=True)
    molecules = ["CCCCCC", "CCO"]
    matrix = adapters.fingerprints_numpy(molecules, options=options, fp_size=512)
    assert matrix.shape == (2, 64)
    assert matrix.dtype == numpy.int64
    for row, molecule in zip(matrix, molecules):
        assert [(int(position), int(row[position])) for position in numpy.flatnonzero(row)] == (
            smsd.fingerprint(molecule, options=options))


def test_sparse_rows_and_empty_rows_keep_their_shape():
    numpy = pytest.importorskip("numpy")
    assert adapters.to_numpy([[], [1, 1, 3]], fp_size=4).tolist() == [[0, 0, 0, 0], [0, 1, 0, 1]]
    values = adapters.to_numpy([[], [(1, 2), (1, 3)], {3: 4}], fp_size=4, counts=True)
    assert values.tolist() == [[0, 0, 0, 0], [0, 5, 0, 0], [0, 0, 0, 4]]
    assert adapters.to_numpy([], fp_size=4, dtype="float32").dtype == numpy.float32


@pytest.mark.parametrize("fingerprints,counts,dtype", [
    ([[-1]], False, None), ([[4]], False, None), ([[1.5]], False, None),
    ([[(0, -1)]], True, None), ([[(0, 1.5)]], True, None),
    ([[(0, 256)]], True, "uint8"), ([[(0, 2)]], True, "bool"),
    ([[1]], False, "complex64"),
])
def test_invalid_sparse_input_and_overflow_raise(fingerprints, counts, dtype):
    pytest.importorskip("numpy")
    with pytest.raises(ValueError):
        adapters.to_numpy(fingerprints, fp_size=4, counts=counts, dtype=dtype)


def test_scalar_smiles_is_not_treated_as_a_batch():
    with pytest.raises(TypeError, match="single string"):
        adapters.fingerprints_numpy("CCO")


def test_numpy_boolean_positions_are_not_integer_positions():
    numpy = pytest.importorskip("numpy")
    with pytest.raises(ValueError, match="integer"):
        adapters.to_numpy([[numpy.bool_(True)]], fp_size=4)


def test_torch_features_preserve_bits_counts_dtype_and_empty_shape(tmp_path):
    if importlib.util.find_spec("torch") is None:
        pytest.skip("PyTorch is not installed")
    source = (f"import smsd, sys; smsd.__path__.insert(0, {os.path.dirname(adapters.__file__)!r})\n" +
              textwrap.dedent("""
        from smsd import adapters
        try:
            matrix = adapters.fingerprints_torch(["CCO", "CCN"], kind="ecfp", fp_size=64)
        except ImportError as error:
            if sys.platform == "darwin" and "OpenMP runtime" in str(error):
                print(str(error), file=sys.stderr)
                sys.exit(77)
            raise
        import torch
        assert matrix.shape == (2, 64)
        assert matrix.dtype == torch.uint8
        assert matrix.device.type == "cpu"
        assert [torch.nonzero(row).flatten().tolist() for row in matrix] == (
            smsd.fingerprints(["CCO", "CCN"], kind="ecfp", fp_size=64))
        counts = adapters.to_torch([[(1, 2), (1, 3)]], counts=True, fp_size=4)
        assert counts.dtype == torch.int64
        assert counts.tolist() == [[0, 5, 0, 0]]
        assert adapters.to_torch([], fp_size=4).shape == (0, 4)
        try:
            adapters.to_torch([[(0, 256)]], counts=True, fp_size=4, dtype=torch.uint8)
        except ValueError:
            pass
        else:
            raise AssertionError("Count overflow was accepted")
    """))
    output = subprocess.run([sys.executable, "-c", source], cwd=tmp_path, text=True,
                            capture_output=True, timeout=30)
    if output.returncode == 77:
        pytest.skip(output.stderr.strip())
    assert output.returncode == 0, output.stderr


def test_mac_torch_guard_rejects_distinct_runtime_before_import(monkeypatch, tmp_path):
    library = tmp_path / "torch" / "lib"
    library.mkdir(parents=True)
    (library / "libomp.dylib").touch()
    monkeypatch.setattr(sys, "platform", "darwin")
    monkeypatch.delitem(sys.modules, "torch", raising=False)
    monkeypatch.setattr(importlib.util, "find_spec", lambda name: SimpleNamespace(
        submodule_search_locations=[str(library.parent)]))
    monkeypatch.setattr(adapters, "_macos_openmp_images", lambda: {"/another/libomp.dylib"})
    with pytest.raises(ImportError, match="SMSD_BUILD_OPENMP=OFF"):
        adapters.to_torch([[1]], fp_size=4)
    assert "torch" not in sys.modules


def test_mac_torch_runtime_conflict_raises_instead_of_aborting(tmp_path):
    if sys.platform != "darwin" or importlib.util.find_spec("torch") is None:
        pytest.skip("Requires macOS and an installed PyTorch package")
    try:
        adapters._check_torch_runtime()
    except ImportError as error:
        assert "OpenMP runtime" in str(error)
    else:
        pytest.skip("This environment uses a compatible OpenMP runtime")
    source = (f"import smsd, sys; smsd.__path__.insert(0, {os.path.dirname(adapters.__file__)!r})\n" +
              textwrap.dedent("""
        from smsd import adapters
        smsd.fingerprints(["CCO", "CCN"])
        try:
            adapters.to_torch([[1]], fp_size=4)
        except ImportError as error:
            assert "SMSD_BUILD_OPENMP=OFF" in str(error)
            assert "torch" not in sys.modules
        else:
            raise AssertionError("Distinct OpenMP runtimes were accepted")
    """))
    output = subprocess.run([sys.executable, "-c", source], cwd=tmp_path, text=True,
                            capture_output=True, timeout=30)
    assert output.returncode == 0, output.stderr


def test_sklearn_transformer_clones_and_fits_in_a_pipeline():
    pytest.importorskip("sklearn")
    from sklearn.base import clone
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.pipeline import Pipeline
    from sklearn.exceptions import NotFittedError
    from smsd.adapters import FingerprintTransformer

    transformer = FingerprintTransformer(fp_size=64, counts=True)
    cloned = clone(transformer)
    assert cloned.get_params() == transformer.get_params()
    with pytest.raises(NotFittedError):
        cloned.transform(["CCO"])
    pipeline = Pipeline([("fingerprints", transformer), ("nearest", KNeighborsClassifier(n_neighbors=1))])
    molecules = ["CCO", "CCN", "c1ccccc1", "c1ccccc1O"]
    labels = [0, 0, 1, 1]
    pipeline.fit(molecules, labels)
    assert pipeline.predict(molecules).tolist() == labels
    assert transformer.transform([]).shape == (0, 64)
    assert transformer.get_feature_names_out().tolist() == [f"fp_{position}" for position in range(64)]
    matrix = FingerprintTransformer(fp_size=64).fit_transform(iter(molecules))
    assert matrix.shape == (4, 64)
