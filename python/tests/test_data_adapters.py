# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Optional data adapters retain counts, molecule identity and graph indices."""

import gc
import json
import os
import subprocess
import sys

import pytest
import smsd
from smsd import adapters


def test_new_dependencies_are_not_imported_with_adapter_module(tmp_path):
    source = ("import json, sys, smsd; "
              f"smsd.__path__.insert(0, {os.path.dirname(adapters.__file__)!r}); "
              "import smsd.adapters; print(json.dumps(sorted(set(sys.modules) & "
              "{'numpy','scipy','pandas','polars','pyarrow','torch','torch_geometric','sklearn'})))")
    result = subprocess.run([sys.executable, "-c", source], cwd=tmp_path, text=True,
                            capture_output=True, check=True, timeout=30)
    assert json.loads(result.stdout) == []


@pytest.mark.parametrize("module,extra,operation", [
    ("scipy", "scipy", lambda: adapters.to_scipy([[1]], fp_size=4)),
    ("pandas", "pandas", lambda: adapters.molecules_pandas(["CCO"])),
    ("polars", "polars", lambda: adapters.molecules_polars(["CCO"])),
    ("pyarrow", "arrow", lambda: adapters.to_arrow([[1]], fp_size=4)),
    ("sklearn", "sklearn", lambda: adapters.__getattr__("SimilarityTransformer")),
])
def test_missing_dependency_names_its_extra(monkeypatch, module, extra, operation):
    original = adapters.import_module

    def unavailable(name):
        if name == module:
            raise ModuleNotFoundError(f"No module named {module!r}", name=module)
        return original(name)

    monkeypatch.setattr(adapters, "import_module", unavailable)
    with pytest.raises(ImportError, match=rf"install smsd\[{extra}\]"):
        operation()


def test_transitive_import_error_is_not_misreported(monkeypatch):
    def broken(name):
        raise ModuleNotFoundError("Missing an internal dependency", name="internal_dependency")

    monkeypatch.setattr(adapters, "import_module", broken)
    with pytest.raises(ModuleNotFoundError, match="internal dependency"):
        adapters.to_scipy([[1]], fp_size=4)


def test_csr_is_canonical_and_preserves_empty_rows_and_duplicate_counts():
    sparse = pytest.importorskip("scipy.sparse")
    numpy = pytest.importorskip("numpy")
    matrix = adapters.to_scipy([[], [(3, 0), (2, 2), (1, 4), (2, 3)]],
                               fp_size=4, counts=True)
    assert sparse.isspmatrix_csr(matrix)
    assert matrix.dtype == numpy.int64
    assert matrix.shape == (2, 4) and matrix.has_canonical_format
    assert matrix.toarray().tolist() == [[0, 0, 0, 0], [0, 4, 5, 0]]
    assert matrix.nnz == 2
    assert adapters.to_scipy([], fp_size=4).shape == (0, 4)
    assert adapters.to_scipy([[3, 1, 1]], fp_size=4).toarray().tolist() == [[0, 1, 0, 1]]


def test_csr_retains_large_integer_counts_and_rejects_overflow():
    pytest.importorskip("scipy.sparse")
    value = (1 << 63) + 3
    matrix = adapters.to_scipy([[(0, value)]], fp_size=2, counts=True, dtype="uint64")
    assert int(matrix.data[0]) == value
    with pytest.raises(ValueError, match="range"):
        adapters.to_scipy([[(0, 255), (0, 1)]], fp_size=2, counts=True, dtype="uint8")
    with pytest.raises(ValueError, match="range"):
        adapters.to_scipy([[(0, 1 << 64)]], fp_size=2, counts=True, dtype="uint64")
    with pytest.raises(ValueError, match="integer dtype"):
        adapters.to_scipy([[(0, 2)]], fp_size=2, counts=True, dtype="float32")


def test_csr_native_fingerprints_agree_with_dense_and_count_options():
    pytest.importorskip("scipy.sparse")
    numpy = pytest.importorskip("numpy")
    molecules = ["CCO", "CCN", "c1ccccc1"]
    assert numpy.array_equal(adapters.fingerprints_scipy(molecules, fp_size=64).toarray(),
                             adapters.fingerprints_numpy(molecules, kind="ecfp", fp_size=64))
    options = smsd.FingerprintOptions(family="ecfp", size=32, counts=True)
    matrix = adapters.fingerprints_scipy(molecules, fp_size=64, options=options)
    assert matrix.shape == (3, 32)
    assert numpy.array_equal(matrix.toarray(), adapters.fingerprints_numpy(molecules, options=options))


def test_pandas_keeps_native_objects_and_exact_fingerprint_rows():
    pandas = pytest.importorskip("pandas")
    graph = smsd.parse_smiles("CCO")
    frame = adapters.molecules_pandas([graph, "CCN"], fp_size=64)
    assert isinstance(frame, pandas.DataFrame)
    assert list(frame.columns) == ["molecule", "fingerprint"]
    assert frame["molecule"].dtype == object
    assert frame["molecule"].iloc[0] is graph
    assert frame["molecule"].iloc[1].n == 3
    assert list(frame["fingerprint"]) == [tuple(row) for row in smsd.fingerprints([graph, "CCN"], kind="ecfp", fp_size=64)]
    assert adapters.molecules_pandas([]).shape == (0, 2)
    with pytest.raises(TypeError, match="MolGraph objects or SMILES"):
        adapters.molecules_pandas([object()])


def test_pandas_count_columns_do_not_coerce_pairs():
    pytest.importorskip("pandas")
    options = smsd.FingerprintOptions(family="ecfp", size=32, counts=True)
    frame = adapters.molecules_pandas(["CCCCCC"], options=options)
    expected = tuple(tuple(pair) for pair in smsd.fingerprint("CCCCCC", options=options))
    assert frame["fingerprint"].iloc[0] == expected


def test_polars_keeps_explicit_object_molecules_and_typed_fingerprints():
    polars = pytest.importorskip("polars")
    graph = smsd.parse_smiles("CCO")
    frame = adapters.molecules_polars([graph, "CCN"], fp_size=64)
    assert frame.schema == {"molecule": polars.Object, "fingerprint": polars.List(polars.UInt64)}
    assert frame["molecule"][0] is graph
    assert frame["fingerprint"].to_list() == smsd.fingerprints([graph, "CCN"], kind="ecfp", fp_size=64)
    assert adapters.molecules_polars([]).shape == (0, 2)
    with pytest.raises(TypeError, match="MolGraph objects or SMILES"):
        adapters.molecules_polars([42])


def test_polars_count_columns_retain_typed_pairs():
    polars = pytest.importorskip("polars")
    options = smsd.FingerprintOptions(family="ecfp", size=32, counts=True)
    frame = adapters.molecules_polars(["CCCCCC"], options=options)
    assert frame.schema["fingerprint"] == polars.List(polars.List(polars.UInt64))
    assert frame["fingerprint"].to_list() == [[list(pair) for pair in smsd.fingerprint("CCCCCC", options=options)]]


def test_arrow_zero_copy_records_aliasing_and_retains_buffer_owner():
    pytest.importorskip("pyarrow")
    numpy = pytest.importorskip("numpy")
    values = numpy.array([[0, 1, 0, 1]], dtype=numpy.uint8)
    table = adapters.to_arrow(values, fp_size=4, zero_copy_only=True)
    assert table.schema.metadata[b"smsd.zero_copy"] == b"true"
    assert table.column(0).chunk(0).values.buffers()[1].address == values.ctypes.data
    values[0, 0] = 1
    assert table.column(0).to_pylist() == [[1, 1, 0, 1]]
    del values
    gc.collect()
    assert table.column(0).to_pylist() == [[1, 1, 0, 1]]


def test_arrow_explicit_copy_and_non_contiguous_copy_are_honest():
    pytest.importorskip("pyarrow")
    numpy = pytest.importorskip("numpy")
    values = numpy.array([[0, 1, 0, 1]], dtype=numpy.uint8)
    copied = adapters.to_arrow(values, fp_size=4, copy=True)
    values[:] = 0
    assert copied.column(0).to_pylist() == [[0, 1, 0, 1]]
    assert copied.schema.metadata[b"smsd.copy_reason"] == b"requested"
    view = numpy.array([[0, 1, 0, 1], [1, 0, 1, 0]], dtype=numpy.uint8)[:, ::2]
    with pytest.raises(ValueError, match="requires a copy"):
        adapters.to_arrow(view, fp_size=2, zero_copy_only=True)
    table = adapters.to_arrow(view, fp_size=2)
    assert table.schema.metadata[b"smsd.zero_copy"] == b"false"
    assert table.column(0).to_pylist() == view.tolist()


def test_arrow_sparse_empty_boolean_and_endian_inputs():
    pytest.importorskip("pyarrow")
    numpy = pytest.importorskip("numpy")
    table = adapters.to_arrow([[], [(1, 2), (1, 3)]], fp_size=4, counts=True)
    assert table.column(0).to_pylist() == [[0, 0, 0, 0], [0, 5, 0, 0]]
    assert table.schema.metadata[b"smsd.copy_reason"] == b"sparse_expansion"
    assert adapters.to_arrow(numpy.zeros((0, 4), dtype=numpy.uint8), fp_size=4).num_rows == 0
    for values in (numpy.array([[True, False]]), numpy.array([[1, 0]], dtype=">i4")):
        with pytest.raises(ValueError, match="requires a copy"):
            adapters.to_arrow(values, fp_size=2, zero_copy_only=True)
        assert adapters.to_arrow(values, fp_size=2).column(0).to_pylist() == values.astype(int).tolist()


def test_arrow_counts_are_exact_and_invalid_values_rejected():
    pytest.importorskip("pyarrow")
    numpy = pytest.importorskip("numpy")
    value = (1 << 63) + 3
    array = numpy.array([[value, 0]], dtype=numpy.uint64)
    assert adapters.to_arrow(array, fp_size=2, counts=True).column(0).to_pylist() == [[value, 0]]
    assert adapters.to_arrow([[(0, value)]], fp_size=2, counts=True, dtype="uint64").column(0).to_pylist() == [[value, 0]]
    with pytest.raises(ValueError, match="range"):
        adapters.to_arrow([[(0, 1 << 64)]], fp_size=2, counts=True, dtype="uint64")
    for values, counts in ((numpy.array([[-1, 0]]), True), (numpy.array([[2, 0]]), False),
                           (numpy.array([[1.0, 0.0]]), True)):
        with pytest.raises(ValueError):
            adapters.to_arrow(values, fp_size=2, counts=counts)
    with pytest.raises(ValueError, match="Sparse fingerprints require a copy"):
        adapters.to_arrow([[1]], fp_size=2, zero_copy_only=True)


def test_similarity_transformer_matches_native_bit_set_jaccard_and_order():
    pytest.importorskip("sklearn")
    numpy = pytest.importorskip("numpy")
    from smsd.adapters import SimilarityTransformer
    references = ["CCN", "CCO", "CCN"]
    queries = ["CCO", "c1ccccc1"]
    model = SimilarityTransformer(fp_size=64).fit(iter(references))
    result = model.transform(queries)
    first = [set(row) for row in smsd.fingerprints(queries, kind="ecfp", fp_size=64)]
    second = [set(row) for row in smsd.fingerprints(references, kind="ecfp", fp_size=64)]
    expected = [[len(left & right) / len(left | right) if left | right else 1.0
                 for right in second] for left in first]
    assert numpy.array_equal(result, expected)
    assert result.dtype == numpy.float64 and result[0, 1] == 1.0
    assert model.transform([]).shape == (0, 3)
    assert model.get_feature_names_out().tolist() == ["reference_0", "reference_1", "reference_2"]


def test_similarity_transformer_fits_iterator_once_and_works_in_pipeline():
    pytest.importorskip("sklearn")
    from sklearn.base import clone
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.pipeline import make_pipeline
    from sklearn.exceptions import NotFittedError
    from smsd.adapters import SimilarityTransformer
    model = SimilarityTransformer(fp_size=64)
    assert clone(model).get_params() == model.get_params()
    with pytest.raises(NotFittedError):
        model.transform(["CCO"])
    assert model.fit_transform(iter(["CCO", "CCN"])).diagonal().tolist() == [1.0, 1.0]
    pipeline = make_pipeline(SimilarityTransformer(fp_size=64), KNeighborsClassifier(n_neighbors=1))
    pipeline.fit(["CCO", "c1ccccc1"], [0, 1])
    assert pipeline.predict(["CCO", "c1ccccc1"]).tolist() == [0, 1]
    model.set_params(fp_size=32)
    with pytest.raises(ValueError, match="refit"):
        model.transform(["CCO"])
    with pytest.raises(ValueError, match="at least one reference"):
        SimilarityTransformer().fit([])


def test_similarity_empty_graphs_have_unit_jaccard_and_inputs_are_not_coerced():
    pytest.importorskip("sklearn")
    from smsd.adapters import SimilarityTransformer
    assert SimilarityTransformer(fp_size=64).fit([smsd.MolGraph()]).transform([smsd.MolGraph()]).tolist() == [[1.0]]
    with pytest.raises(TypeError, match="MolGraph objects or SMILES"):
        SimilarityTransformer().fit([123])


def _pyg_dependencies():
    try:
        adapters._check_torch_runtime()
    except ImportError as error:
        pytest.skip(str(error))
    torch = pytest.importorskip("torch")
    pytest.importorskip("torch_geometric.data")
    return torch


def test_pyg_preserves_native_indices_atom_features_and_both_bond_directions():
    torch = _pyg_dependencies()
    molecule = smsd.parse_smiles("[13CH3][NH2+]C")
    data = adapters.to_pyg(molecule)
    columns = [molecule.atomic_num, molecule.formal_charge, molecule.mass_number,
               molecule.hydrogen_count, molecule.degree, molecule.aromatic, molecule.ring]
    assert data.x.dtype == torch.int64
    assert data.x.tolist() == [[int(column[atom]) for column in columns] for atom in range(molecule.n)]
    assert data.atom_id.tolist() == list(range(molecule.n))
    expected_edges, expected_attributes = [], []
    for first, second, order in sorted(molecule.bonds()):
        expected_edges.extend([[first, second], [second, first]])
        feature = [order, int(molecule.bond_aromatic(first, second)), int(molecule.bond_in_ring(first, second))]
        expected_attributes.extend([feature, feature])
    assert data.edge_index.t().tolist() == expected_edges
    assert data.edge_attr.tolist() == expected_attributes
    assert data.num_nodes == molecule.n and data.validate(raise_on_error=True)
    assert data.pos is None and data.y is None


@pytest.mark.parametrize("smiles,nodes", [("", 0), ("[Na+].[Cl-]", 2)])
def test_pyg_empty_and_isolated_atoms_are_not_dropped(smiles, nodes):
    _pyg_dependencies()
    data = adapters.to_pyg(smiles)
    assert data.num_nodes == nodes and tuple(data.x.shape) == (nodes, 7)
    assert tuple(data.edge_index.shape) == (2, 0) and tuple(data.edge_attr.shape) == (0, 3)
    assert data.atom_id.tolist() == list(range(nodes))


def test_pyg_preserves_recorded_aromatic_ring_features():
    _pyg_dependencies()
    molecule = smsd.parse_smiles("c1ccccc1")
    data = adapters.to_pyg(molecule)
    assert data.x[:, 5].tolist() == [int(value) for value in molecule.aromatic]
    assert data.x[:, 6].tolist() == [int(value) for value in molecule.ring]
    assert data.edge_attr[:, 1:].tolist() == [[1, 1]] * 12


def test_pyg_missing_dependency_and_torch_guard_fail_before_import(monkeypatch):
    original = adapters.import_module
    monkeypatch.setattr(adapters, "_check_torch_runtime", lambda: None)

    def unavailable(name):
        if name in ("torch", "torch_geometric"):
            raise ModuleNotFoundError("Missing optional graph dependency", name=name)
        return original(name)

    monkeypatch.setattr(adapters, "import_module", unavailable)
    with pytest.raises(ImportError, match=r"install smsd\[pyg\]"):
        adapters.to_pyg("CCO")

    def conflict():
        raise ImportError("OpenMP runtime conflict")

    monkeypatch.setattr(adapters, "_check_torch_runtime", conflict)
    monkeypatch.setattr(adapters, "import_module", lambda name: pytest.fail("Imported before guard"))
    with pytest.raises(ImportError, match="OpenMP runtime conflict"):
        adapters.to_pyg("CCO")
