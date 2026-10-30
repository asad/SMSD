# Python feature adapters and MCP

The optional modules `smsd.adapters`, `smsd.mcp_server` and `smsd.highlight`
leave the core API unchanged. They are part of `smsd` 8.1.0
(`python -m pip install "smsd==8.1.0"`). NumPy, PyTorch, scikit-learn and the
MCP SDK are loaded only by the functions that need them.

## NumPy and PyTorch

The first Python block requires both NumPy and PyTorch, supplied by the `torch`
extra. On macOS, use an OpenMP-free SMSD build or a compatible shared OpenMP
runtime before running it; the source-build command below installs both extras.
The adapters compute native molecular fingerprints and return a matrix with
one molecule per row and one fingerprint position per column:

```python
import smsd
from smsd.adapters import fingerprints_numpy, fingerprints_torch, to_numpy

molecules = ["CCO", "CCN", "c1ccccc1O"]
features = fingerprints_numpy(molecules, kind="ecfp", fp_size=2048)
assert features.shape == (3, 2048)
assert str(features.dtype) == "uint8"

tensor = fingerprints_torch(molecules, kind="ecfp", fp_size=2048)
assert tuple(tensor.shape) == (3, 2048)

options = smsd.FingerprintOptions(family="ecfp", size=512, counts=True)
counts = fingerprints_numpy(molecules, options=options)
assert counts.shape == (3, 512)
assert str(counts.dtype) == "int64"

sparse = smsd.fingerprints(molecules, kind="path", fp_size=128)
assert to_numpy(sparse, fp_size=128).shape == (3, 128)
```

`to_numpy` and `to_torch` consume sparse fingerprint rows: set-bit positions,
or `(position, count)` pairs with `counts=True`. They do not consume dense
bit vectors. Empty batches have shape `(0, fp_size)`; an empty row is a zero
vector. Repeated count bins are summed. Positions outside the declared size,
negative counts and counts that overflow the chosen dtype raise `ValueError`.
An explicitly selected floating dtype can reduce integer precision.

Binary output defaults to NumPy `uint8` or PyTorch `uint8`; count output
defaults to `int64`. PyTorch defaults to the CPU and accepts `dtype` and
`device`, such as `device="cuda"` for an NVIDIA GPU or `device="mps"` on Apple
silicon; on macOS, PyTorch may need the OpenMP-free SMSD build described below.
Its adapter does not require NumPy. Native fingerprint options,
including structural keys and count fingerprints, determine the output size
and count mode when supplied.

Packaged macOS OpenMP wheels may require an OpenMP-free SMSD build with
PyTorch. Different bundled OpenMP runtimes can abort the process; the adapters
conservatively raise `ImportError` before importing PyTorch or running their
native batch when they detect this combination. A shared runtime or an
OpenMP-free SMSD build is required. To build the latter from the release source:

```sh
python -m pip install ".[numpy,torch]" --config-settings=cmake.define.SMSD_BUILD_OPENMP=OFF
```

This is a source build requiring the documented C++ build dependencies.
Import order and a single requested thread do not resolve distinct OpenMP
runtimes. No duplicate-runtime environment override is used.

## scikit-learn

Install the `sklearn` extra. The transformer accepts an iterable of molecules
or SMILES strings and participates in ordinary pipelines and parameter cloning:

```python
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline
from smsd.adapters import FingerprintTransformer

model = make_pipeline(
    FingerprintTransformer(kind="ecfp", fp_size=512),
    KNeighborsClassifier(n_neighbors=1),
)
model.fit(["CCO", "CCN", "c1ccccc1", "c1ccccc1O"], [0, 0, 1, 1])
assert model.predict(["CCO"]).tolist() == [0]
```

The transformer is stateless: `fit` records that it is fitted, and `transform`
computes the fingerprints. Its `counts=True` mode supports ECFP, FCFP, path
and torsion fingerprints. A scalar SMILES string must be wrapped in a list;
the adapter will not treat its characters as separate molecules.

## MCP over standard input/output

The optional `smsd[mcp]` extra uses the official Python SDK over standard
input/output. Run `smsd-mcp` or `python -m smsd.mcp_server` in the environment
that has the extra installed.

The eight tools are `find_mcs`, `find_substructure`, `search_smarts`,
`molecular_fingerprint`, `molecular_similarity`, `screen_molecules`,
`standardise_molecule` and `depict_mcs`. They accept bounded molecular text,
not paths, commands or network requests. Native proof/status remains distinct
from fingerprint scores and from whether a witness was found. Default
standardisation changes representation only; chemical transformations require
explicit options and the optional backend.

See [MCP tools and function definitions](MCP_SERVER.md) for signatures, finite
budgets, original-JSON validation, result pages, timing boundaries, SVG scope,
function-definition formats, lazy framework adapters and the local-wheel
container/demo recipes. Exporting definitions makes no model request.

## MCS SVG highlight

The native SVG helper renders an existing result without running another
search or importing an optional drawing toolkit:

```python
import smsd
from smsd.highlight import mcs_svg

query, target = "c1ccccc1O", "c1ccccc1N"
result = smsd.MCS.find(query).with_timeout(1000).match(target)
svg = mcs_svg(query, target, result, width=600, height=300)
assert "<svg " in svg
```

Use the same molecules and numbering as the search. The SVG title and
`data-search-status`, `data-optimal` and `data-complete` attributes retain the
result's proof status. The pair renderer highlights mapped atoms and bonds
joining them. If an MCES result omits a bond between mapped atoms, the helper
raises rather than highlighting that extra bond; render
`smsd.depict_svg(result.best.to_graph(query))` to draw the selected common
graph in that case.

### MCP fingerprint limits

MCP inputs are limited to 16,384 text characters. The fingerprint tool accepts
at most 1,024 parsed atoms, 65,536 fingerprint positions, radius 8 and path
length 8. Path and torsion requests are also checked against a conservative
2,000,000-walk expansion bound before enumeration. Larger requests raise
`ValueError`; these service limits do not change the native fingerprint APIs.
Fingerprinting has no hard wall-clock deadline. Search budgets apply after
parsing, as described above.

## SciPy sparse features

Install `smsd[scipy]`. `to_scipy` accepts the same sparse positions or
position/count pairs as `to_numpy`, and `fingerprints_scipy` computes native
fingerprints first. Both return a canonical CSR matrix with sorted columns,
combined duplicates and no stored zero counts. They never allocate the full
dense feature matrix. Binary values default to `uint8`; counts default to
`int64`. Only boolean/integer dtypes are accepted. Counts are combined as
Python integers before checking the selected dtype; overflow raises rather
than wrapping or rounding.

```python
from smsd.adapters import to_scipy, fingerprints_scipy

counts = to_scipy([[], [(1, 2), (1, 3)]], fp_size=4, counts=True)
assert counts.toarray().tolist() == [[0, 0, 0, 0], [0, 5, 0, 0]]
assert fingerprints_scipy(["CCO", "CCN"], fp_size=64).shape == (2, 64)
```

## pandas and Polars molecule columns

Install `smsd[pandas]` or `smsd[polars]`. `molecules_pandas` and
`molecules_polars` return a DataFrame with `molecule` and `fingerprint`
columns, in input order. Inputs must be native `MolGraph` objects or explicit
SMILES strings. SMILES are parsed; native molecule references are retained.
Other objects, missing values and scalar strings supplied as a batch raise.
These functions do not stringify arbitrary objects or convert external
molecule types implicitly.

pandas uses object columns with immutable fingerprint tuples. Polars uses an
explicit Object column for molecules and typed lists of `uint64` positions.
With count options, its fingerprint rows are lists of `[position, count]`
pairs, each value checked against `uint64`. Neither frame is a copy of the
native molecules. Polars Object molecules cannot be exported as Arrow
molecular data; export the fingerprint features separately.

```python
from smsd.adapters import molecules_pandas, molecules_polars

pandas_frame = molecules_pandas(["CCO", "CCN"], fp_size=64)
polars_frame = molecules_polars(["CCO", "CCN"], fp_size=64)
assert pandas_frame["molecule"].iloc[0].n == 3
assert polars_frame["molecule"][0].n == 3
```

Both functions accept the ordinary native fingerprint arguments, including
`options` for count mode, and `num_threads`. Empty batches retain both columns
and their declared types.

## Arrow fingerprint buffers

Install `smsd[arrow]`. `to_arrow` exports one `fingerprint` column as an Arrow
fixed-size list of integer values. Sparse inputs are expanded once into an
owned dense NumPy matrix. A native-endian, C-contiguous two-dimensional
integer NumPy matrix can instead share its value buffer without copying.
Boolean NumPy storage needs conversion to byte values, since Arrow boolean
storage is bit-packed. Other layouts and byte orders also need a copy.
Floating values are rejected, including in count mode.

```python
import numpy as np
from smsd.adapters import to_arrow

features = np.array([[0, 1, 0, 1]], dtype=np.uint8)
table = to_arrow(features, fp_size=4, zero_copy_only=True)
assert table.schema.metadata[b"smsd.zero_copy"] == b"true"
assert table.column("fingerprint").to_pylist() == [[0, 1, 0, 1]]
```

Arrow retains the NumPy buffer owner after the caller drops its reference.
Shared values still alias the original array: do not mutate that array while
using the exported table. Use `copy=True` for isolation. `zero_copy_only=True`
rejects any required copy and cannot be combined with `copy=True`.

Schema metadata records `smsd.fp_size`, `smsd.counts`, `smsd.zero_copy`,
`smsd.copy_reason` and `smsd.ownership`. A false zero-copy flag means copying
from the original input was necessary, even though Arrow subsequently retains
the prepared NumPy buffer. `dtype` is an optional integer dtype for sparse
expansion only; an existing matrix retains its dtype. This is a feature-table
export, not molecular graph serialisation or a GPU transfer.

## Similarity features in scikit-learn pipelines

The `sklearn` extra also supplies `SimilarityTransformer`. `fit` snapshots
binary native fingerprints for its reference molecules, retaining order and
duplicates. `transform` returns a dense `float64` matrix of Tanimoto values,
one column per fitted reference. It computes intersection/union of native
set-bit positions; two empty fingerprints have similarity one. It does not
fit an approximate index, execute MCS, or infer molecular identity from a
fingerprint. The output matrix requires memory proportional to queries times
references.

```python
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline
from smsd.adapters import SimilarityTransformer

model = make_pipeline(
    SimilarityTransformer(kind="ecfp", fp_size=128),
    KNeighborsClassifier(n_neighbors=1),
)
model.fit(["CCO", "c1ccccc1"], [0, 1])
assert model.predict(["CCO"]).tolist() == [0]
```

`fit_transform` consumes an input iterator once. Empty query batches produce
zero rows; an empty reference fit raises. Changing fingerprint parameters
after fitting requires a refit. Native molecule objects and SMILES are the
accepted input types, as for the DataFrame adapters.

## PyTorch Geometric molecular graphs

Install `smsd[pyg]`. `to_pyg` returns a PyTorch Geometric `Data` object, using
the same macOS OpenMP guard as the existing PyTorch adapter. It defaults to
CPU; an explicitly selected tensor device does not move molecular parsing or
fingerprint generation to a GPU.

```python
from smsd.adapters import to_pyg

graph = to_pyg("CCO")
assert tuple(graph.x.shape) == (3, 7)
assert tuple(graph.edge_index.shape) == (2, 4)
assert graph.atom_id.tolist() == [0, 1, 2]
```

All feature tensors are `int64`. Atom rows preserve native atom numbering;
`atom_id` records it explicitly, and `num_nodes` retains isolated atoms.
The seven `x` columns are atomic number, formal charge, mass number,
recorded hydrogen count, degree, aromatic flag and ring flag. These are the
native graph's recorded values, not newly inferred descriptors.

Every unique bond produces both directed edges, ordered by its native endpoint
pair and then forward/reverse. The three `edge_attr` columns are recorded bond
order, aromatic flag and ring flag. Empty graphs have shapes `(0, 7)`,
`(2, 0)` and `(0, 3)`. No coordinates, stereo labels, training targets or
one-hot encodings are invented. Stereo and coordinates are not encoded by
this feature schema, so it is not a lossless molecular file format. Add
task-specific features explicitly before model training.

The external data structures follow the official
[SciPy CSR interface](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html),
[Arrow fixed-size-list interface](https://arrow.apache.org/docs/python/generated/pyarrow.FixedSizeListArray.html)
and [PyTorch Geometric Data interface](https://pytorch-geometric.readthedocs.io/en/latest/generated/torch_geometric.data.Data.html).

## Jupyter and Colab notebook

The [feature-search notebook](../python/examples/feature_search.ipynb) combines
MCS-to-SMARTS search, SVG highlighting, native fingerprint features, a small
scikit-learn pipeline and a temporary packed corpus. Its example labels are
illustrative, not a predictive-performance evaluation.

Use a local kernel with SMSD 8.1.0, NumPy and scikit-learn already installed.
The setup cell detects an actual Colab shell; ordinary local execution performs
no upload or installation. Upload the notebook file to Colab and select a fresh
Linux x86_64 runtime. There, the cell prompts for one user-uploaded SMSD 8.1.0
core wheel. The distribution name, version, CPython ABI and compatible Linux
x86_64 platform tags are checked
before installation with `sys.executable -m pip --no-index --no-deps --force-reinstall`.
The imported version is then checked. No CUDA package is needed.

Download the matching trusted wheel beforehand from PyPI or the
[GitHub releases page](https://github.com/asad/SMSD/releases); the notebook does
not download packages.
NumPy and scikit-learn must already be available; they are checked but are not
installed by this cell. Run the upload setup before importing SMSD, and restart
the Colab runtime before repeating it. A native wheel executes code when loaded.
