# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Optional fingerprint features, molecule columns and graph-learning adapters.

Importing this module loads none of its optional data or ML dependencies. Fingerprint
rows are sparse set-bit positions, or sparse ``(position, count)`` pairs when
``counts=True``. Dense inputs are not interpreted as sparse fingerprints.
"""

from collections.abc import Mapping
from importlib import import_module
from numbers import Integral
from operator import index

from ._fingerprints import fingerprints
from ._molecule import _checked_int

__all__ = ["to_numpy", "to_torch", "fingerprints_numpy", "fingerprints_torch",
           "FingerprintTransformer", "SimilarityTransformer", "to_scipy", "fingerprints_scipy",
           "molecules_pandas", "molecules_polars", "to_arrow", "to_pyg"]


def _optional(module, extra):
    """Load an optional dependency when its adapter is called."""
    try:
        return import_module(module)
    except ModuleNotFoundError as error:
        if error.name != module:
            raise
        raise ImportError(f"This adapter requires {module}; install smsd[{extra}]") from error


def _integer(value, name):
    """An integer position or count, including integer array scalars."""
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer, not a boolean or a float")
    try:
        return index(value)
    except TypeError as error:
        raise ValueError(f"{name} must be an integer") from error


def _rows(sparse_fingerprints, fp_size, counts):
    """Validate sparse inputs and combine repeated count bins without losing counts."""
    _checked_int(fp_size, "fp_size", 1)
    if not isinstance(counts, bool):
        raise ValueError("counts must be True or False")
    rows = []
    for fingerprint in sparse_fingerprints:
        row = {}
        items = fingerprint.items() if isinstance(fingerprint, Mapping) else fingerprint
        for item in items:
            if counts:
                try:
                    position, count = item
                except (TypeError, ValueError) as error:
                    raise ValueError("Count fingerprints require (position, count) pairs") from error
                count = _integer(count, "Fingerprint count")
                if count < 0:
                    raise ValueError("Fingerprint counts must be non-negative")
            else:
                position, count = item, 1
            position = _integer(position, "Fingerprint position")
            if not 0 <= position < fp_size:
                raise ValueError(f"Fingerprint position {position} is outside [0, {fp_size})")
            row[position] = row.get(position, 0) + count if counts else 1
        rows.append(row)
    return rows


def _check_range(rows, maximum):
    """Reject count overflow instead of wrapping or silently clipping it."""
    if any(count > maximum for row in rows for count in row.values()):
        raise ValueError("Fingerprint count is outside the selected dtype's range")


def _macos_openmp_images():
    """Locate loaded OpenMP runtimes without initialising either runtime."""
    import ctypes
    import os

    process = ctypes.CDLL(None)
    count = process._dyld_image_count
    count.restype = ctypes.c_uint32
    image_name = process._dyld_get_image_name
    image_name.argtypes = [ctypes.c_uint32]
    image_name.restype = ctypes.c_char_p
    paths = set()
    for image_index in range(count()):
        name = image_name(image_index)
        if name is not None:
            path = os.fsdecode(name)
            filename = os.path.basename(path)
            if filename.startswith(("libomp", "libiomp")) and filename.endswith(".dylib"):
                paths.add(os.path.realpath(path))
    return paths


def _check_torch_runtime():
    """Reject macOS runtime combinations that can abort before Python can raise."""
    import sys
    if sys.platform != "darwin":
        return
    import os
    from importlib.util import find_spec

    specification = find_spec("torch")
    if specification is None:
        return
    try:
        loaded = _macos_openmp_images()
    except (AttributeError, OSError) as error:
        raise ImportError("Cannot check PyTorch's OpenMP runtime on macOS; use an "
                          "SMSD build with SMSD_BUILD_OPENMP=OFF") from error
    if not loaded:
        return
    if "torch" in sys.modules and len(loaded) == 1:
        return
    bundled = {os.path.realpath(os.path.join(directory, "lib", "libomp.dylib"))
               for directory in specification.submodule_search_locations or ()
               if os.path.isfile(os.path.join(directory, "lib", "libomp.dylib"))}
    if len(loaded) == 1 and loaded <= bundled:
        return
    raise ImportError("PyTorch would use a different OpenMP runtime from one already loaded "
                      "on macOS. Install an SMSD build with SMSD_BUILD_OPENMP=OFF, or use "
                      "packages linked to one shared OpenMP runtime.")


def to_numpy(sparse_fingerprints, *, fp_size=2048, counts=False, dtype=None):
    """Expand fingerprints to a NumPy array of shape ``(n_molecules, fp_size)``.

    Binary fingerprints default to ``uint8``; count fingerprints default to
    ``int64``. Empty batches have shape ``(0, fp_size)``. Repeated binary bits
    stay one, repeated count bins are summed, and invalid positions raise.
    An explicit numeric dtype may change precision, but cannot overflow.

    Example::

        from smsd.adapters import to_numpy
        assert to_numpy([[1, 3]], fp_size=4).tolist() == [[0, 1, 0, 1]]
    """
    numpy = _optional("numpy", "numpy")
    rows = _rows(sparse_fingerprints, fp_size, counts)
    dtype = numpy.dtype(dtype if dtype is not None else ("int64" if counts else "uint8"))
    if dtype.kind not in "buif":
        raise ValueError("dtype must be a boolean, integer or floating-point dtype")
    maximum = 1 if dtype.kind == "b" else (
        numpy.finfo(dtype).max if dtype.kind == "f" else numpy.iinfo(dtype).max)
    _check_range(rows, maximum)
    result = numpy.zeros((len(rows), fp_size), dtype=dtype)
    for row_index, row in enumerate(rows):
        for position, count in row.items():
            result[row_index, position] = count
    return result


def to_torch(sparse_fingerprints, *, fp_size=2048, counts=False, dtype=None, device="cpu"):
    """Expand fingerprints to a PyTorch tensor without requiring NumPy.

    Shape and sparse-input rules are those of :func:`to_numpy`. Binary rows
    default to ``torch.uint8``, counts to ``torch.int64`` and the device to CPU.
    On macOS, incompatible OpenMP runtimes raise before importing PyTorch.

    Example::

        from smsd.adapters import to_torch
        assert to_torch([[1, 3]], fp_size=4).tolist() == [[0, 1, 0, 1]]
    """
    _check_torch_runtime()
    torch = _optional("torch", "torch")
    rows = _rows(sparse_fingerprints, fp_size, counts)
    dtype = dtype if dtype is not None else (torch.int64 if counts else torch.uint8)
    if not isinstance(dtype, torch.dtype) or dtype.is_complex:
        raise ValueError("dtype must be a boolean, integer or floating-point torch dtype")
    try:
        maximum = 1 if dtype == torch.bool else (
            torch.finfo(dtype).max if dtype.is_floating_point else torch.iinfo(dtype).max)
    except TypeError as error:
        raise ValueError("dtype must be a boolean, integer or floating-point torch dtype") from error
    _check_range(rows, maximum)
    dense_rows = [[row.get(position, 0) for position in range(fp_size)] for row in rows]
    return torch.tensor(dense_rows, dtype=dtype, device=device).reshape(len(rows), fp_size)


def fingerprints_numpy(molecules, *, kind="mcs", path_length=7, fp_size=2048, radius=2,
                       options=None, num_threads=0, dtype=None):
    """Compute native fingerprints and expand them to a NumPy feature matrix.

    Fingerprint arguments have the same meaning as :func:`smsd.fingerprints`.
    When ``options`` is supplied, its ``size`` and ``counts`` define the matrix.
    A single molecule must be supplied as a one-element iterable.

    Example::

        from smsd.adapters import fingerprints_numpy
        assert fingerprints_numpy(["CCO"], fp_size=64).shape == (1, 64)
    """
    _batch(molecules)
    values = fingerprints(molecules, kind=kind, path_length=path_length, fp_size=fp_size,
                          radius=radius, options=options, num_threads=num_threads)
    return to_numpy(values, fp_size=options.size if options is not None else fp_size,
                    counts=options.counts if options is not None else False, dtype=dtype)


def fingerprints_torch(molecules, *, kind="mcs", path_length=7, fp_size=2048, radius=2,
                       options=None, num_threads=0, dtype=None, device="cpu"):
    """Compute native fingerprints and expand them to a PyTorch feature tensor.

    Example::

        from smsd.adapters import fingerprints_torch
        assert tuple(fingerprints_torch(["CCO"], fp_size=64).shape) == (1, 64)
    """
    _batch(molecules)
    _check_torch_runtime()
    values = fingerprints(molecules, kind=kind, path_length=path_length, fp_size=fp_size,
                          radius=radius, options=options, num_threads=num_threads)
    return to_torch(values, fp_size=options.size if options is not None else fp_size,
                    counts=options.counts if options is not None else False,
                    dtype=dtype, device=device)


def _batch(molecules):
    """Reject an accidental scalar SMILES string before iterating its characters."""
    if isinstance(molecules, (str, bytes)):
        raise TypeError("molecules must be an iterable of molecules, not a single string")


def __getattr__(name):
    """Create the estimator only when scikit-learn integration is requested."""
    if name == "SimilarityTransformer":
        return _similarity_transformer()
    if name != "FingerprintTransformer":
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    _optional("sklearn", "sklearn")
    from sklearn.base import BaseEstimator, TransformerMixin
    from sklearn.utils.validation import check_is_fitted

    class FingerprintTransformer(TransformerMixin, BaseEstimator):
        """A stateless molecular fingerprint transformer for scikit-learn pipelines.

        Input is a one-dimensional iterable of molecules or SMILES strings.
        ``fit`` marks the transformer fitted; ``transform`` computes native
        fingerprints. Count mode supports ECFP, FCFP, path and torsion families.

        Example::

            from smsd.adapters import FingerprintTransformer
            model = FingerprintTransformer(fp_size=64)
            assert model.fit_transform(["CCO"]).shape == (1, 64)
        """

        def __init__(self, *, kind="ecfp", fp_size=2048, radius=2, path_length=7,
                     counts=False, dtype=None, num_threads=0):
            self.kind = kind
            self.fp_size = fp_size
            self.radius = radius
            self.path_length = path_length
            self.counts = counts
            self.dtype = dtype
            self.num_threads = num_threads

        def fit(self, X, y=None):
            """Mark this stateless transformer fitted without consuming an input iterator.

            Example::

                model = FingerprintTransformer(fp_size=64).fit(["CCO"])
                assert model.is_fitted_
            """
            _batch(X)
            self.is_fitted_ = True
            return self

        def transform(self, X):
            """Return one native fingerprint vector per input molecule.

            Example::

                model = FingerprintTransformer(fp_size=64).fit(["CCO"])
                assert model.transform(["CCN"]).shape == (1, 64)
            """
            check_is_fitted(self, "is_fitted_")
            if not isinstance(self.counts, bool):
                raise ValueError("counts must be True or False")
            options = None
            if self.counts:
                from . import FingerprintOptions
                options = FingerprintOptions(family=self.kind, size=self.fp_size,
                                             radius=self.radius, path_length=self.path_length,
                                             counts=True)
            return fingerprints_numpy(X, kind=self.kind, fp_size=self.fp_size, radius=self.radius,
                                      path_length=self.path_length, options=options,
                                      num_threads=self.num_threads, dtype=self.dtype)

        def get_feature_names_out(self, input_features=None):
            """Fingerprint column names; input molecules are samples, not feature columns.

            Example::

                model = FingerprintTransformer(fp_size=64).fit(["CCO"])
                assert model.get_feature_names_out()[0] == "fp_0"
            """
            check_is_fitted(self, "is_fitted_")
            if input_features is not None:
                raise ValueError("Molecular inputs have no input feature columns")
            numpy = _optional("numpy", "numpy")
            _checked_int(self.fp_size, "fp_size", 1)
            return numpy.asarray([f"fp_{position}" for position in range(self.fp_size)], dtype=object)

    FingerprintTransformer.__qualname__ = name
    globals()[name] = FingerprintTransformer
    return FingerprintTransformer


def to_scipy(sparse_fingerprints, *, fp_size=2048, counts=False, dtype=None):
    """Return a canonical SciPy CSR matrix without expanding dense rows.

    Binary rows use uint8; counts use int64. Integer dtypes only are accepted,
    and repeated count bins are summed in Python before range checking.

    Example::

        from smsd.adapters import to_scipy
        matrix = to_scipy([[], [(1, 2), (1, 3)]], fp_size=4, counts=True)
        assert matrix.toarray().tolist() == [[0, 0, 0, 0], [0, 5, 0, 0]]
    """
    _optional("scipy", "scipy")
    sparse = import_module("scipy.sparse")
    numpy = _optional("numpy", "scipy")
    rows = _rows(sparse_fingerprints, fp_size, counts)
    dtype = numpy.dtype(dtype if dtype is not None else ("int64" if counts else "uint8"))
    if dtype.kind not in "bui":
        raise ValueError("CSR fingerprints require a boolean or integer dtype")
    _check_range(rows, 1 if dtype.kind == "b" else int(numpy.iinfo(dtype).max))
    data, indices, indptr = [], [], [0]
    for row in rows:
        for position, count in sorted(row.items()):
            if count:
                indices.append(position)
                data.append(count)
        indptr.append(len(data))
    return sparse.csr_matrix((numpy.asarray(data, dtype=dtype),
                              numpy.asarray(indices, dtype=numpy.int64),
                              numpy.asarray(indptr, dtype=numpy.int64)),
                             shape=(len(rows), fp_size), dtype=dtype)


def fingerprints_scipy(molecules, *, kind="ecfp", path_length=7, fp_size=2048,
                       radius=2, options=None, num_threads=0, dtype=None):
    """Compute native fingerprints as an integer CSR feature matrix.

    Example::

        from smsd.adapters import fingerprints_scipy
        assert fingerprints_scipy(["CCO", "CCN"], fp_size=64).shape == (2, 64)
    """
    _batch(molecules)
    values = fingerprints(molecules, kind=kind, path_length=path_length, fp_size=fp_size,
                          radius=radius, options=options, num_threads=num_threads)
    return to_scipy(values, fp_size=options.size if options is not None else fp_size,
                    counts=options.counts if options is not None else False, dtype=dtype)


def _graphs(molecules):
    """Accept native molecules and explicit SMILES, without coercing other objects."""
    _batch(molecules)
    from ._molecule import MolGraph, parse_smiles
    result = []
    for molecule in molecules:
        if isinstance(molecule, MolGraph):
            result.append(molecule)
        elif isinstance(molecule, str):
            result.append(parse_smiles(molecule))
        else:
            raise TypeError("These adapters accept MolGraph objects or SMILES strings only")
    return result


def _frame_features(molecules, kind, fp_size, radius, path_length, options, num_threads):
    graphs = _graphs(molecules)
    values = fingerprints(graphs, kind=kind, fp_size=fp_size, radius=radius,
                          path_length=path_length, options=options, num_threads=num_threads)
    counts = options.counts if options is not None else False
    rows = _rows(values, options.size if options is not None else fp_size, counts)
    values = [tuple((position, count) for position, count in sorted(row.items()) if count)
              if counts else tuple(sorted(row)) for row in rows]
    return graphs, values, counts


def molecules_pandas(molecules, *, kind="ecfp", fp_size=2048, radius=2,
                     path_length=7, options=None, num_threads=0):
    """Return a pandas DataFrame with molecule and fingerprint object columns.

    SMILES are parsed explicitly; native MolGraph references are retained, not
    stringified or copied. Fingerprints are immutable tuples of positions, or
    (position, count) pairs with count options. Other input objects are rejected.

    Example::

        from smsd.adapters import molecules_pandas
        frame = molecules_pandas(["CCO", "CCN"], fp_size=64)
        assert frame["molecule"].iloc[0].n == 3
    """
    pandas = _optional("pandas", "pandas")
    graphs, values, _ = _frame_features(molecules, kind, fp_size, radius, path_length,
                                       options, num_threads)
    return pandas.DataFrame({"molecule": pandas.Series(graphs, dtype=object),
                             "fingerprint": pandas.Series(values, dtype=object)})


def molecules_polars(molecules, *, kind="ecfp", fp_size=2048, radius=2,
                     path_length=7, options=None, num_threads=0):
    """Return a Polars DataFrame with explicit Object molecules and typed lists.

    The molecule column holds native references. Binary fingerprints are lists
    of uint64 positions; counts are lists of [position, count] uint64 pairs.
    Counts outside uint64 raise. Object molecules are not Arrow-serialisable.

    Example::

        from smsd.adapters import molecules_polars
        frame = molecules_polars(["CCO", "CCN"], fp_size=64)
        assert frame["molecule"][0].n == 3
    """
    polars = _optional("polars", "polars")
    graphs, values, counts = _frame_features(molecules, kind, fp_size, radius, path_length,
                                            options, num_threads)
    if any(value > (1 << 64) - 1 for row in values for item in row
           for value in (item if counts else (item,))):
        raise ValueError("Polars fingerprint positions and counts must fit uint64")
    dtype = polars.List(polars.List(polars.UInt64)) if counts else polars.List(polars.UInt64)
    lists = [[[position, count] for position, count in row] if counts else list(row) for row in values]
    return polars.DataFrame([
        polars.Series("molecule", graphs, dtype=polars.Object, strict=True),
        polars.Series("fingerprint", lists, dtype=dtype, strict=True),
    ])


def to_arrow(fingerprints, *, fp_size=2048, counts=False, dtype=None,
             copy=False, zero_copy_only=False):
    """Export a fixed-size-list Arrow fingerprint column with ownership metadata.

    A native-endian, C-contiguous 2-D integer NumPy array can share its buffer.
    Arrow retains that owner; callers must not mutate shared source values.
    Sparse rows, boolean buffers, other layouts and copy=True require an owned
    NumPy buffer. zero_copy_only=True rejects those cases. Counts remain exact
    non-negative integers. Schema metadata records copying and its reason.

    Example::

        import numpy as np
        from smsd.adapters import to_arrow
        bits = np.array([[0, 1, 0, 1]], dtype=np.uint8)
        table = to_arrow(bits, fp_size=4, zero_copy_only=True)
        assert table.schema.metadata[b"smsd.zero_copy"] == b"true"
    """
    arrow = _optional("pyarrow", "arrow")
    numpy = _optional("numpy", "arrow")
    _checked_int(fp_size, "fp_size", 1)
    if not all(isinstance(value, bool) for value in (counts, copy, zero_copy_only)):
        raise ValueError("counts, copy and zero_copy_only must be booleans")
    if copy and zero_copy_only:
        raise ValueError("copy=True conflicts with zero_copy_only=True")
    external = isinstance(fingerprints, numpy.ndarray)
    reason = "none"
    if external:
        if dtype is not None:
            raise ValueError("dtype is only used to expand sparse fingerprint input")
        matrix = fingerprints
        if matrix.ndim != 2 or matrix.shape[1] != fp_size:
            raise ValueError("Dense fingerprints must have shape (n_molecules, fp_size)")
    else:
        if zero_copy_only:
            raise ValueError("Sparse fingerprints require a copy")
        rows = _rows(fingerprints, fp_size, counts)
        selected_dtype = numpy.dtype(dtype if dtype is not None else ("int64" if counts else "uint8"))
        if selected_dtype.kind not in "bui":
            raise ValueError("Arrow fingerprints require a boolean or integer dtype")
        _check_range(rows, 1 if selected_dtype.kind == "b" else int(numpy.iinfo(selected_dtype).max))
        matrix = numpy.zeros((len(rows), fp_size), dtype=selected_dtype)
        for row_index, row in enumerate(rows):
            for position, count in row.items():
                matrix[row_index, position] = count
        reason = "sparse_expansion"
    if matrix.dtype.kind not in "bui":
        raise ValueError("Arrow fingerprints require a boolean or integer dtype")
    if matrix.size and (int(matrix.min()) < 0 or (not counts and int(matrix.max()) > 1)):
        raise ValueError("Fingerprint values must be non-negative counts or binary zero/one values")
    needs_layout = matrix.dtype.kind == "b" or not matrix.flags.c_contiguous or not matrix.dtype.isnative
    if zero_copy_only and needs_layout:
        raise ValueError("This dtype or layout requires a copy")
    if needs_layout or (external and copy):
        matrix = numpy.array(matrix, dtype=("uint8" if matrix.dtype.kind == "b"
                                           else matrix.dtype.newbyteorder("=")), order="C", copy=True)
        reason = "requested" if copy else "layout_or_boolean"
    shared = external and reason == "none"
    values = arrow.Array.from_buffers(arrow.from_numpy_dtype(matrix.dtype), matrix.size,
                                     [None, arrow.py_buffer(matrix.reshape(-1))])
    column = arrow.FixedSizeListArray.from_arrays(values, list_size=fp_size)
    return arrow.Table.from_arrays([column], names=["fingerprint"]).replace_schema_metadata({
        b"smsd.fp_size": str(fp_size).encode(), b"smsd.counts": str(counts).lower().encode(),
        b"smsd.zero_copy": str(shared).lower().encode(), b"smsd.copy_reason": reason.encode(),
        b"smsd.ownership": b"shared_numpy_owner_retained" if shared else b"owned_numpy_owner_retained",
    })


def to_pyg(molecule, *, device="cpu"):
    """Export native graph features as a PyTorch Geometric Data object.

    Atom rows keep native indices. Seven int64 features are atomic number,
    formal charge, mass number, recorded hydrogen count, degree, aromatic and
    ring flags. Every bond yields both directions; its three int64 features
    are recorded order, aromatic flag and ring flag. No coordinates, stereo
    labels or targets are invented. This is not a lossless molecule format.

    Example::

        from smsd.adapters import to_pyg
        graph = to_pyg("CCO")
        assert graph.num_nodes == 3 and tuple(graph.edge_index.shape) == (2, 4)
    """
    _check_torch_runtime()
    torch = _optional("torch", "pyg")
    _optional("torch_geometric", "pyg")
    data_type = import_module("torch_geometric.data").Data
    graph = _graphs([molecule])[0]
    columns = [graph.atomic_num, graph.formal_charge, graph.mass_number,
               graph.hydrogen_count, graph.degree, graph.aromatic, graph.ring]
    rows = [[int(column[atom]) for column in columns] for atom in range(graph.n)]
    edges, attributes = [], []
    for first, second, order in sorted(graph.bonds()):
        features = [int(order), int(graph.bond_aromatic(first, second)),
                    int(graph.bond_in_ring(first, second))]
        edges.extend(((first, second), (second, first)))
        attributes.extend((features, features))
    return data_type(x=torch.tensor(rows, dtype=torch.int64, device=device).reshape(graph.n, 7),
                     edge_index=torch.tensor(edges, dtype=torch.int64, device=device).reshape(-1, 2).t().contiguous(),
                     edge_attr=torch.tensor(attributes, dtype=torch.int64, device=device).reshape(-1, 3),
                     atom_id=torch.arange(graph.n, dtype=torch.int64, device=device), num_nodes=graph.n)


def _similarity_transformer():
    """Create the reference-fingerprint transformer only on explicit access."""
    _optional("sklearn", "sklearn")
    from sklearn.base import BaseEstimator, TransformerMixin
    from sklearn.utils.validation import check_is_fitted

    class SimilarityTransformer(TransformerMixin, BaseEstimator):
        """Binary native-fingerprint Tanimoto features against fitted references.

        Reference order and duplicates are retained; empty/empty similarity is
        one. fit snapshots fingerprints, not mutable input molecules. Output
        is a dense float64 matrix, not an approximate neighbour search.

        Example::

            from smsd.adapters import SimilarityTransformer
            model = SimilarityTransformer(fp_size=64).fit(["CCO", "CCN"])
            assert model.transform(["CCO"]).shape == (1, 2)
        """

        def __init__(self, *, kind="ecfp", fp_size=2048, radius=2, path_length=7, num_threads=0):
            self.kind = kind
            self.fp_size = fp_size
            self.radius = radius
            self.path_length = path_length
            self.num_threads = num_threads

        def _settings(self):
            return self.kind, self.fp_size, self.radius, self.path_length

        def _fingerprints(self, molecules):
            return tuple(frozenset(row) for row in fingerprints(
                _graphs(molecules), kind=self.kind, fp_size=self.fp_size, radius=self.radius,
                path_length=self.path_length, num_threads=self.num_threads))

        def fit(self, X, y=None):
            """Snapshot references in input order.

            Example::

                model = SimilarityTransformer(fp_size=64).fit(["CCO", "CCN"])
                assert model.n_references_ == 2
            """
            references = self._fingerprints(X)
            if not references:
                raise ValueError("SimilarityTransformer needs at least one reference molecule")
            self.reference_fingerprints_ = references
            self.n_references_ = len(references)
            self._fitted_settings_ = self._settings()
            return self

        def _matrix(self, queries):
            numpy = _optional("numpy", "sklearn")
            matrix = numpy.empty((len(queries), self.n_references_), dtype=numpy.float64)
            for row, query in enumerate(queries):
                for column, reference in enumerate(self.reference_fingerprints_):
                    union = len(query | reference)
                    matrix[row, column] = len(query & reference) / union if union else 1.0
            return matrix

        def transform(self, X):
            """Compare to fitted references.

            Example::

                model = SimilarityTransformer(fp_size=64).fit(["CCO"])
                assert model.transform(["CCO"]).tolist() == [[1.0]]
            """
            check_is_fitted(self, "reference_fingerprints_")
            if self._settings() != self._fitted_settings_:
                raise ValueError("Fingerprint parameters changed; refit SimilarityTransformer")
            return self._matrix(self._fingerprints(X))

        def fit_transform(self, X, y=None, **fit_params):
            """Consume references once, including one-shot iterators.

            Example::

                model = SimilarityTransformer(fp_size=64)
                assert model.fit_transform(iter(["CCO"])).tolist() == [[1.0]]
            """
            if fit_params:
                raise TypeError("SimilarityTransformer has no additional fit parameters")
            self.fit(X, y)
            return self._matrix(self.reference_fingerprints_)

        def get_feature_names_out(self, input_features=None):
            """Name reference columns.

            Example::

                model = SimilarityTransformer(fp_size=64).fit(["CCO"])
                assert model.get_feature_names_out().tolist() == ["reference_0"]
            """
            check_is_fitted(self, "reference_fingerprints_")
            if input_features is not None:
                raise ValueError("Molecular inputs have no input feature columns")
            numpy = _optional("numpy", "sklearn")
            return numpy.asarray([f"reference_{number}" for number in range(self.n_references_)], dtype=object)

    SimilarityTransformer.__qualname__ = "SimilarityTransformer"
    globals()["SimilarityTransformer"] = SimilarityTransformer
    return SimilarityTransformer
