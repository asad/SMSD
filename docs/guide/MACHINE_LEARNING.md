# Machine learning and data frames

The `smsd.adapters` module turns SMSD Pro fingerprints into the structures that
data and machine-learning libraries expect: NumPy matrices, SciPy sparse
matrices, pandas and Polars columns, Arrow tables, PyTorch tensors and
scikit-learn transformers. This chapter shows which adapter suits which data
stack, how each one is used, and the points about data types, bit order,
copies and memory that decide whether a pipeline is correct and affordable. A
worked example at the end fingerprints a small set of molecules, trains two
scikit-learn models in pipelines and screens the same molecules by similarity.

## Choose an adapter

| Data stack | Use | Returns | Extra |
| --- | --- | --- | --- |
| NumPy arrays | `fingerprints_numpy`, or `to_numpy` on fingerprints already computed | A dense matrix, one row per molecule; `uint8` bits or `int64` counts | `numpy` |
| SciPy sparse matrices | `fingerprints_scipy`, or `to_scipy` | A CSR matrix with sorted columns; `uint8` bits or `int64` counts | `scipy` |
| pandas | `molecules_pandas` | A data frame with `molecule` and `fingerprint` object columns | `pandas` |
| Polars | `molecules_polars` | A data frame with an Object `molecule` column and a typed list `fingerprint` column | `polars` |
| Arrow, and tools that read Arrow | `to_arrow` on a NumPy matrix | A table with one fixed-size-list `fingerprint` column, shared with the matrix where possible | `arrow` |
| scikit-learn, fingerprint features | `FingerprintTransformer` | A transformer whose output is the `fingerprints_numpy` matrix | `sklearn` |
| scikit-learn, similarity features or a Tanimoto kernel | `SimilarityTransformer` | A transformer whose output is a `float64` Tanimoto matrix against fitted references | `sklearn` |
| PyTorch | `fingerprints_torch`, or `to_torch` | A tensor, on the CPU unless a device is given | `torch` |
| PyTorch Geometric | `to_pyg` | A `Data` object with atom and bond features of one molecule | `pyg` |

As a rule of thumb, take a dense NumPy matrix when the model needs one and the
library fits in memory, a CSR matrix when it does not and the estimator accepts
sparse input, and the data-frame adapters when molecules and fingerprints must
stay next to other columns. Arrow is for handing fingerprints to another
process or file format without a copy. For library-scale similarity searching,
a [packed fingerprint index](PACKED_INDEX.md) is better suited than a feature
matrix.

Install the extras that the chosen adapters need, for example
`python -m pip install "smsd[sklearn,pandas,arrow]==8.1.0"`;
[Getting started](GETTING_STARTED.md#install) lists every extra. Importing
`smsd.adapters` loads none of these packages; each function loads its own when
it is called and raises `ImportError` naming the extra when the package is
missing.

## Feature matrices with NumPy

`fingerprints_numpy(molecules, ...)` computes the fingerprints and expands them
to a matrix in one call. Column `j` is fingerprint position `j`. Bits are
`uint8` and counts are `int64` unless `dtype` says otherwise. `to_numpy` does
the same for sparse fingerprints that have already been computed, such as the
output of `smsd.fingerprints`. Both accept SMILES strings, `MolGraph` objects
and RDKit molecules.

```python
import numpy as np
import smsd
from smsd.adapters import fingerprints_numpy, to_numpy

smiles = ["CC(=O)Oc1ccccc1C(=O)O", "OC(=O)c1ccccc1O", "CC(=O)Nc1ccc(O)cc1",
          "Oc1ccccc1"]  # aspirin, salicylic acid, paracetamol, phenol
ecfp4 = smsd.FingerprintOptions.ecfp(radius=2, size=2048)

# One row per molecule, one column per bit position.
bits = fingerprints_numpy(smiles, options=ecfp4)
print("bits:", bits.shape, bits.dtype, f"{bits.nbytes:,} bytes")
print("set bits per row:", bits.sum(axis=1).tolist())
assert bits.dtype == np.uint8 and bits.sum(axis=1).tolist() == [36, 28, 24, 14]

# The same matrix from sparse positions computed earlier.
positions = smsd.fingerprints(smiles, options=ecfp4)
assert (to_numpy(positions, fp_size=2048) == bits).all()
assert bits[0].nonzero()[0].tolist() == positions[0]

# Count fingerprints: one column per bin, int64 counts.
counts_options = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
counts_options.counts = True
counts = fingerprints_numpy(smiles, options=counts_options)
print("counts:", counts.shape, counts.dtype, f"{counts.nbytes:,} bytes")
print("largest count per row:", counts.max(axis=1).tolist())
assert counts.dtype == np.int64 and counts.max(axis=1).tolist() == [3, 3, 4, 4]
assert ((counts > 0) == (bits > 0)).all()

# Without kind or options the family is "mcs", not ECFP.
default = fingerprints_numpy(smiles)
print("default kind equals ECFP4:", bool((default == bits).all()))
assert not (default == bits).all()
```

The program prints:

```text
bits: (4, 2048) uint8 8,192 bytes
set bits per row: [36, 28, 24, 14]
counts: (4, 2048) int64 65,536 bytes
largest count per row: [3, 3, 4, 4]
default kind equals ECFP4: False
```

Like `smsd.fingerprints`, `fingerprints_numpy` defaults to `kind="mcs"`; pass
`kind` or `options` every time. When `options` is given, its `size` and
`counts` decide the width and the data type of the matrix. A count matrix is
eight times the size of a bit matrix of the same width.
[NumPy arrays](FINGERPRINTS.md#numpy-arrays) in the fingerprint chapter shows a
Tanimoto matrix computed from such a matrix with one product.

## Sparse matrices with SciPy

A 2,048-bit ECFP4 row of a drug-like molecule has a few dozen set bits, so most
of a dense matrix is zeros. `fingerprints_scipy` and `to_scipy` return a SciPy
CSR matrix without ever building the dense matrix. Only integer and boolean
data types are accepted.

```python
import itertools

import smsd
from smsd.adapters import fingerprints_numpy, fingerprints_scipy

GROUPS = ["C", "CC", "OC", "N", "O", "F", "Cl", "C(=O)O"]
smiles = [f"c1c({a})cc({b})cc1{c}"
          for a, b, c in itertools.product(GROUPS, repeat=3)]
ecfp4 = smsd.FingerprintOptions.ecfp(radius=2, size=2048)

dense = fingerprints_numpy(smiles, options=ecfp4)
sparse = fingerprints_scipy(smiles, options=ecfp4)
sparse_bytes = sparse.data.nbytes + sparse.indices.nbytes + sparse.indptr.nbytes
density = sparse.nnz / (sparse.shape[0] * sparse.shape[1])
print("molecules:", sparse.shape[0], "| set bits:", sparse.nnz,
      f"| density: {density:.4f}")
print(f"dense: {dense.nbytes:,} bytes | CSR: {sparse_bytes:,} bytes")
assert (sparse.toarray() == dense).all() and sparse.dtype == dense.dtype
assert sparse.nnz == 12874
```

The program prints:

```text
molecules: 512 | set bits: 12874 | density: 0.0123
dense: 1,048,576 bytes | CSR: 66,422 bytes
```

The CSR matrix stores each set bit as a value and a 64-bit column index, so for
`uint8` bits it saves memory only while fewer than about one position in nine
is set. Most
scikit-learn linear models and nearest-neighbour estimators accept CSR input;
check the estimator before relying on it, because one that needs dense input
will convert the whole matrix.

## pandas and Polars columns

`molecules_pandas` and `molecules_polars` return a data frame with two columns,
in input order: `molecule`, holding the parsed `MolGraph`, and `fingerprint`,
holding the sorted set-bit positions. Further columns, such as names or
measured values, are added in the usual way. `to_numpy` turns the fingerprint
column into a feature matrix when a model needs one.

```python
import polars as pl
import smsd
from smsd.adapters import molecules_pandas, molecules_polars, to_numpy

names = ["aspirin", "salicylic acid", "paracetamol", "phenol"]
smiles = ["CC(=O)Oc1ccccc1C(=O)O", "OC(=O)c1ccccc1O", "CC(=O)Nc1ccc(O)cc1",
          "Oc1ccccc1"]
ecfp4 = smsd.FingerprintOptions.ecfp(radius=2, size=2048)

# pandas: object columns, one MolGraph and one tuple of positions per row.
frame = molecules_pandas(smiles, options=ecfp4)
frame.insert(0, "name", names)
frame["atoms"] = [graph.n for graph in frame["molecule"]]
frame["bits"] = frame["fingerprint"].map(len)
print(frame[["name", "atoms", "bits"]].to_string(index=False))
print("pandas columns:", frame["molecule"].dtype, frame["fingerprint"].dtype)
X = to_numpy(frame["fingerprint"], fp_size=2048)
print("pandas features:", X.shape, X.dtype)

# Polars: an Object column for molecules and a typed list of positions.
table = molecules_polars(smiles, options=ecfp4).with_columns(
    pl.Series("name", names),
    pl.col("fingerprint").list.len().alias("bits"),
)
print("polars columns:", table.schema["molecule"], table.schema["fingerprint"])
Y = to_numpy(table["fingerprint"].to_list(), fp_size=2048)
assert (X == Y).all()
assert table["bits"].to_list() == frame["bits"].tolist() == [36, 28, 24, 14]
assert isinstance(frame["molecule"].iloc[0], smsd.MolGraph)
```

The program prints:

```text
          name  atoms  bits
       aspirin     13    36
salicylic acid     10    28
   paracetamol     11    24
        phenol      7    14
pandas columns: object object
pandas features: (4, 2048) uint8
polars columns: Object List(UInt64)
```

pandas stores the fingerprints as immutable tuples; Polars stores them as lists
of `uint64` positions, or lists of `[position, count]` pairs with count
options. Both functions accept only `MolGraph` objects and SMILES strings, and
raise `TypeError` for anything else, including RDKit molecules and missing
values; convert or drop those rows first. The frames hold references to the
molecules rather than copies.

A Polars Object column cannot be exported to Arrow. Drop the `molecule` column,
or keep the SMILES in a string column, before writing a Polars frame to
Parquet or Arrow.

## Arrow tables and zero-copy export

`to_arrow` writes one `fingerprint` column as an Arrow fixed-size list of
integers. When it is given a two-dimensional, C-contiguous, native-endian
integer NumPy matrix, the Arrow column uses the matrix's memory directly, and
Arrow keeps the matrix alive for as long as the table exists. Other boolean
or integer input is copied once into a new matrix first, and floating-point
input is rejected. The schema metadata says which happened, and
`zero_copy_only=True` turns a needed copy into a `ValueError`.

```python
import numpy as np
import polars as pl
import smsd
from smsd.adapters import fingerprints_numpy, to_arrow

smiles = ["CC(=O)Oc1ccccc1C(=O)O", "OC(=O)c1ccccc1O", "CC(=O)Nc1ccc(O)cc1",
          "Oc1ccccc1"]
ecfp4 = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
bits = fingerprints_numpy(smiles, options=ecfp4)


def report(label, table):
    metadata = table.schema.metadata
    print(f"{label:18} zero copy {metadata[b'smsd.zero_copy'].decode():5}  "
          f"reason {metadata[b'smsd.copy_reason'].decode()}")
    return metadata[b"smsd.zero_copy"] == b"true"


# A C-contiguous integer matrix shares its buffer with the Arrow column.
shared = to_arrow(bits, fp_size=2048, zero_copy_only=True)
print("type:", shared.schema.field("fingerprint").type)
assert report("uint8 matrix", shared)
values = shared.column("fingerprint").chunk(0).values
assert np.shares_memory(values.to_numpy(zero_copy_only=True), bits)

# These inputs need a copy, and the metadata says why.
assert not report("boolean matrix", to_arrow(bits.astype(bool), fp_size=2048))
assert not report("column slice", to_arrow(bits[:, :1024], fp_size=1024))
positions = smsd.fingerprints(smiles, options=ecfp4)
assert not report("sparse positions", to_arrow(positions, fp_size=2048))
assert not report("copy=True", to_arrow(bits, fp_size=2048, copy=True))
try:
    to_arrow(bits.astype(bool), fp_size=2048, zero_copy_only=True)
except ValueError as error:
    print("zero_copy_only:", error)

# Polars reads the table as a fixed-width array column.
frame = pl.from_arrow(shared)
print("polars:", frame.schema["fingerprint"],
      frame["fingerprint"].to_numpy().shape)

# Shared values alias the source: a change to the matrix shows in the table.
assert bits[3, 0] == 0
bits[3, 0] = 1
seen = shared.column("fingerprint")[3].values[0].as_py()
print("table sees the change:", seen)
assert seen == 1
```

The program prints:

```text
type: fixed_size_list<item: uint8>[2048]
uint8 matrix       zero copy true   reason none
boolean matrix     zero copy false  reason layout_or_boolean
column slice       zero copy false  reason layout_or_boolean
sparse positions   zero copy false  reason sparse_expansion
copy=True          zero copy false  reason requested
zero_copy_only: This dtype or layout requires a copy
polars: Array(UInt8, shape=(2048,)) (4, 2048)
table sees the change: 1
```

The rules for a copy-free export:

- Use an integer data type. A boolean matrix is converted to `uint8`, because
  Arrow packs booleans as bits; a floating-point matrix is rejected, also in
  count mode.
- Keep the matrix C-contiguous. A column slice, a transposed matrix or a
  Fortran-ordered matrix is copied, as is a stepped slice of rows such as
  `bits[::2]`; a contiguous range of whole rows is not.
- Keep the byte order native. A big-endian matrix on a little-endian machine is
  copied.
- Do not change the matrix while the table is in use, as the last lines of the
  program show. Pass `copy=True` when the table must be independent of it.

`copy=True` and `zero_copy_only=True` cannot be combined. The `dtype` argument
applies only when sparse positions are expanded; passing it with a NumPy
matrix raises `ValueError`. The table carries fingerprint features only, not
the molecules, and the export is a memory operation, not a transfer to a GPU.

## Packed bytes and bit order

A [packed fingerprint index](PACKED_INDEX.md#file-format) stores bit `p` of a
row in byte `p // 8`, at bit `p % 8`: bit 0 is the least significant bit of the
first byte. NumPy's `packbits` and `unpackbits` use the opposite order unless
they are told otherwise, so pass `bitorder="little"` whenever packed bytes are
exchanged with an index:

```python
import tempfile
from pathlib import Path

import numpy as np
import smsd
from smsd.adapters import fingerprints_numpy

smiles = ["CC(=O)Oc1ccccc1C(=O)O", "Oc1ccccc1"]  # aspirin, phenol
bits = fingerprints_numpy(smiles, kind="ecfp", fp_size=2048, radius=2)

little = np.packbits(bits, axis=1, bitorder="little")
big = np.packbits(bits, axis=1)  # NumPy's default is bitorder="big"
print("packed:", little.shape, little.dtype)

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / "two.smsdfp"
    with smsd.FingerprintCorpus.from_smiles(path, smiles, kind="ecfp",
                                            fp_size=2048, radius=2) as index:
        (start, stored), = index.iter_packed_chunks()
        print("matches bitorder='little':", bool((stored == little).all()))
        print("matches the default order:", bool((stored == big).all()))
        assert (stored == little).all() and not (stored == big).all()
        assert (np.unpackbits(stored, axis=1, bitorder="little") == bits).all()
```

The program prints:

```text
packed: (2, 256) uint8
matches bitorder='little': True
matches the default order: False
```

A wrong bit order raises no error: the bytes are valid, but every similarity
computed from them is wrong. The matrices of `smsd.adapters` hold one column
per position and have no bit order to get wrong.

## scikit-learn transformers

`FingerprintTransformer` wraps `fingerprints_numpy` as a scikit-learn
transformer. It takes a one-dimensional sequence of molecules (a list, a
pandas Series or one column selected by a `ColumnTransformer`) and returns the
feature matrix. It learns nothing: `fit` only marks it fitted, so the same
molecules always give the same features, inside or outside cross-validation.
Its parameters are ordinary estimator parameters and work with `clone`,
`set_params` and grid searches.

`SimilarityTransformer` stores the binary fingerprints of the molecules passed
to `fit`, in order and with duplicates, and `transform` returns one `float64`
column per stored reference: the Tanimoto similarity of each input molecule to
that reference.

```python
import pandas as pd
import smsd
from sklearn.base import clone
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from smsd.adapters import FingerprintTransformer, SimilarityTransformer

smiles = ["CC(=O)Oc1ccccc1C(=O)O", "OC(=O)c1ccccc1O",
          "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "OC(=O)c1ccccc1",
          "CC(=O)Nc1ccc(O)cc1", "Oc1ccccc1", "Clc1ccccc1",
          "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"]
data = pd.DataFrame({
    "smiles": smiles,
    "heavy_atoms": [smsd.parse_smiles(text).n for text in smiles],
    "acid": [1, 1, 1, 1, 0, 0, 0, 0],  # 1: carboxylic acid
})

# A transformer on its own: fit records nothing, transform fingerprints.
fingerprint = FingerprintTransformer(kind="ecfp", fp_size=1024)
X = fingerprint.fit_transform(data["smiles"])
print("features:", X.shape, X.dtype, "| first names:",
      fingerprint.get_feature_names_out()[:3].tolist())

# A column of SMILES beside a numeric column, in one pipeline.
features = ColumnTransformer([
    ("fingerprint", FingerprintTransformer(kind="ecfp", fp_size=1024),
     "smiles"),
    ("size", "passthrough", ["heavy_atoms"]),
])
model = make_pipeline(features, LogisticRegression(max_iter=1000))
model.fit(data[["smiles", "heavy_atoms"]], data["acid"])
predictions = model.predict(data[["smiles", "heavy_atoms"]]).tolist()
print("predictions:", predictions)
assert predictions == data["acid"].tolist()

# Parameters work with clone and set_params, as for any estimator.
wider = clone(model).set_params(columntransformer__fingerprint__fp_size=2048)
wider.fit(data[["smiles", "heavy_atoms"]], data["acid"])
print("model inputs after widening:", wider[-1].n_features_in_)
assert wider[-1].n_features_in_ == 2049

# Tanimoto similarity to fitted references, one column per reference.
references = SimilarityTransformer(kind="ecfp", fp_size=2048).fit(smiles[:2])
scores = references.transform(["OC(=O)c1ccccc1"]).round(3).tolist()
print("columns:", references.get_feature_names_out().tolist())
print("benzoic acid:", scores)
assert scores == [[0.256, 0.278]]
references.set_params(fp_size=1024)
try:
    references.transform(["OC(=O)c1ccccc1"])
except ValueError as error:
    print("after set_params:", error)
```

The program prints:

```text
features: (8, 1024) uint8 | first names: ['fp_0', 'fp_1', 'fp_2']
predictions: [1, 1, 1, 1, 0, 0, 0, 0]
model inputs after widening: 2049
columns: ['reference_0', 'reference_1']
benzoic acid: [[0.256, 0.278]]
after set_params: Fingerprint parameters changed; refit SimilarityTransformer
```

The training predictions only show that the pipeline runs; they are not a
measure of model quality.

| Transformer | Learns in `fit` | Output | Output size | Molecule input |
| --- | --- | --- | --- | --- |
| `FingerprintTransformer` | Nothing | `uint8` bits, `int64` counts or the chosen `dtype` | Molecules times `fp_size` | SMILES, `MolGraph` or RDKit molecules |
| `SimilarityTransformer` | The reference fingerprints | `float64` Tanimoto values | Molecules times references | SMILES or `MolGraph` only |

Choose `FingerprintTransformer` when the model should weigh individual
structural features, as linear models, tree ensembles and neural networks do.
Choose `SimilarityTransformer` when similarity to known molecules is the
signal: as the features of a nearest-neighbour model, as a precomputed kernel
(see the worked example), or as a column of "distance to the nearest active"
beside other descriptors.

## PyTorch tensors

`fingerprints_torch` and `to_torch` return a tensor of the same shape and with
the same default types as the NumPy adapters, without needing NumPy. The
tensor is created on the CPU unless `device` names another device, and its
type follows `dtype`. Write the call so that a missing or unusable PyTorch is
handled:

```python
import smsd
from smsd.adapters import fingerprints_numpy, fingerprints_torch

smiles = ["CC(=O)Oc1ccccc1C(=O)O", "OC(=O)c1ccccc1O", "CC(=O)Nc1ccc(O)cc1",
          "Oc1ccccc1"]
ecfp4 = smsd.FingerprintOptions.ecfp(radius=2, size=2048)

try:
    # The adapter checks the runtime before it imports PyTorch.
    features = fingerprints_torch(smiles, options=ecfp4)
except ImportError as error:
    print("PyTorch features unavailable:", error)
    features = fingerprints_numpy(smiles, options=ecfp4)
    print("NumPy features instead:", features.shape, features.dtype)
else:
    import torch  # already loaded by the adapter

    if torch.cuda.is_available():
        features = features.to(device="cuda", dtype=torch.float32)
    print("tensor:", tuple(features.shape), features.dtype, features.device)
assert tuple(features.shape) == (4, 2048)
```

Without PyTorch, the program prints:

```text
PyTorch features unavailable: This adapter requires torch; install smsd[torch]
NumPy features instead: (4, 2048) uint8
```

With PyTorch installed, the program prints one `tensor:` line instead: a
`torch.uint8` tensor on the CPU, or, where PyTorch reports a usable CUDA
device, a `torch.float32` tensor on that device. Fingerprints are always
computed on CPU threads; the `device` argument and `.to()` move only the
finished tensor. On macOS, a packaged SMSD wheel built with OpenMP and PyTorch
can load two different OpenMP runtimes, which can stop the process. The
adapters detect this before importing PyTorch and raise `ImportError`, which
the program above catches; the remedy is an SMSD build without OpenMP, or
packages linked to one shared OpenMP runtime, as described in
[NumPy and PyTorch](../PYTHON_ADAPTERS.md#numpy-and-pytorch). Call the
adapter before importing PyTorch yourself, so that the check runs first.

`to_pyg` turns one molecule into a PyTorch Geometric `Data` object with seven
integer atom features and three integer bond features per directed edge;
[PyTorch Geometric molecular graphs](../PYTHON_ADAPTERS.md#pytorch-geometric-molecular-graphs)
lists them.

## A worked example: train and screen

This program uses 18 labelled molecules to train two models, one on
fingerprint features and one on Tanimoto similarity as a kernel, and scores 8
held-out molecules with both. It then screens the same held-out molecules by
their highest similarity to the known actives, and prints the three views side
by side. The label marks the benzodiazepines; it is illustrative, and a set of
this size says nothing about predictive performance.

```python
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.svm import SVC
from smsd.adapters import FingerprintTransformer, SimilarityTransformer

# Training set: label 1 for a benzodiazepine, 0 otherwise.
TRAIN = [
    ("diazepam", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21", 1),
    ("nordazepam", "O=C1CN=C(c2ccccc2)c2cc(Cl)ccc2N1", 1),
    ("temazepam", "CN1C(=O)C(O)N=C(c2ccccc2)c2cc(Cl)ccc21", 1),
    ("oxazepam", "OC1N=C(c2ccccc2)c2cc(Cl)ccc2NC1=O", 1),
    ("lorazepam", "OC1N=C(c2ccccc2Cl)c2cc(Cl)ccc2NC1=O", 1),
    ("clonazepam", "O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1", 1),
    ("nitrazepam", "O=C1CN=C(c2ccccc2)c2cc([N+](=O)[O-])ccc2N1", 1),
    ("alprazolam", "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2", 1),
    ("triazolam", "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1Cl)=NC2", 1),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O", 0),
    ("salicylic acid", "OC(=O)c1ccccc1O", 0),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", 0),
    ("naproxen", "COc1ccc2cc(ccc2c1)C(C)C(=O)O", 0),
    ("caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 0),
    ("phenol", "Oc1ccccc1", 0),
    ("benzophenone", "O=C(c1ccccc1)c1ccccc1", 0),
    ("4-chloroaniline", "Nc1ccc(Cl)cc1", 0),
    ("acetanilide", "CC(=O)Nc1ccccc1", 0),
]
# Held-out molecules, scored after training.
TEST = [
    ("flunitrazepam", "CN1C(=O)CN=C(c2ccccc2F)c2cc([N+](=O)[O-])ccc21", 1),
    ("bromazepam", "O=C1CN=C(c2ccccn2)c2cc(Br)ccc2N1", 1),
    ("midazolam", "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2", 1),
    ("tetrazepam", "CN1C(=O)CN=C(C2=CCCCC2)c2cc(Cl)ccc21", 1),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1", 0),
    ("2-aminobenzophenone", "Nc1ccccc1C(=O)c1ccccc1", 0),
    ("chlorobenzene", "Clc1ccccc1", 0),
    ("benzoic acid", "OC(=O)c1ccccc1", 0),
]
train_smiles = [smiles for _, smiles, _ in TRAIN]
y_train = np.array([label for _, _, label in TRAIN])
test_names = [name for name, _, _ in TEST]
test_smiles = [smiles for _, smiles, _ in TEST]
y_test = np.array([label for _, _, label in TEST])

# 1. Fingerprint features and a linear model.
linear = make_pipeline(
    FingerprintTransformer(kind="ecfp", fp_size=2048, radius=2),
    LogisticRegression(C=1.0, max_iter=1000),
)
linear.fit(train_smiles, y_train)
probability = linear.predict_proba(test_smiles)[:, 1]

# 2. A Tanimoto kernel: the similarity matrix to the training molecules is
#    the precomputed kernel of a support vector machine.
kernel = make_pipeline(
    SimilarityTransformer(kind="ecfp", fp_size=2048, radius=2),
    SVC(kernel="precomputed", C=1.0),
)
kernel.fit(train_smiles, y_train)
margin = kernel.decision_function(test_smiles)

# 3. A similarity screen: the nearest known active for each test molecule.
active_names = [name for name, _, label in TRAIN if label == 1]
actives = [smiles for _, smiles, label in TRAIN if label == 1]
screen = SimilarityTransformer(kind="ecfp", fp_size=2048, radius=2).fit(actives)
scores = screen.transform(test_smiles)
best = scores.max(axis=1)
nearest = scores.argmax(axis=1)

print(f"{'molecule':20} {'label':>5} {'P(active)':>9} {'margin':>7} "
      f"{'Tanimoto':>8}  nearest active")
order = np.argsort(-best, kind="stable")
for i in order:
    print(f"{test_names[i]:20} {y_test[i]:5d} {probability[i]:9.3f} "
          f"{margin[i]:7.3f} {best[i]:8.3f}  {active_names[nearest[i]]}")
hits = [test_names[i] for i in order if best[i] >= 0.4]
for label, model in (("linear", linear), ("kernel", kernel)):
    correct = int((model.predict(test_smiles) == y_test).sum())
    print(f"{label} model correct: {correct} of 8")
print("screen hits at Tanimoto >= 0.4:", hits)

assert ((probability > 0.5) == (y_test == 1)).all()
assert ((margin > 0) == (y_test == 1)).all()
assert best.round(3).tolist() == [0.465, 0.519, 0.541, 0.542, 0.031, 0.125,
                                  0.256, 0.154]
assert hits == ["tetrazepam", "midazolam", "bromazepam", "flunitrazepam"]
```

The program prints:

```text
molecule             label P(active)  margin Tanimoto  nearest active
tetrazepam               1     0.792   0.320    0.542  diazepam
midazolam                1     0.860   0.444    0.541  triazolam
bromazepam               1     0.841   0.387    0.519  nordazepam
flunitrazepam            1     0.835   0.424    0.465  clonazepam
chlorobenzene            0     0.181  -0.436    0.256  diazepam
benzoic acid             0     0.069  -0.730    0.154  oxazepam
2-aminobenzophenone      0     0.056  -0.735    0.125  nordazepam
paracetamol              0     0.031  -0.698    0.031  diazepam
linear model correct: 8 of 8
kernel model correct: 8 of 8
screen hits at Tanimoto >= 0.4: ['tetrazepam', 'midazolam', 'bromazepam', 'flunitrazepam']
```

The three columns answer different questions. The probability and the margin
come from models that weigh the training molecules of both classes; the
Tanimoto column depends only on the actives and is the easiest to explain,
since the nearest active can be shown beside each hit. The Tanimoto values are
exact for the chosen fingerprint. The probabilities and margins come from
numerical solvers and may differ in the last digit between platforms or
library versions.

The kernel pipeline is correct as written because `SimilarityTransformer`
returns the training molecules against themselves when the pipeline is
fitted, and new molecules against the training molecules when it predicts,
which are the two matrices a precomputed kernel needs. It holds one `float64`
value per pair, so it suits training sets of thousands of molecules, not
millions.

To take the hits further, confirm the shared scaffold with an exact search, as
in [Screening a library](SCREENING.md#choose-a-route). For a library too large
to hold as a feature matrix, screen it with
[`TargetCorpus` or a streamed screen](SCREENING.md#choose-a-route), or with a
[packed fingerprint index](PACKED_INDEX.md), and pass only the shortlist to
the model.

## Options

### Feature functions

| Function | Fingerprint arguments | Other arguments | Returns |
| --- | --- | --- | --- |
| `fingerprints_numpy(molecules, ...)` | `kind="mcs"`, `path_length=7`, `fp_size=2048`, `radius=2`, `options=None` | `num_threads=0`, `dtype=None` | NumPy matrix |
| `fingerprints_scipy(molecules, ...)` | `kind="ecfp"`, `path_length=7`, `fp_size=2048`, `radius=2`, `options=None` | `num_threads=0`, `dtype=None` | SciPy CSR matrix |
| `fingerprints_torch(molecules, ...)` | `kind="mcs"`, `path_length=7`, `fp_size=2048`, `radius=2`, `options=None` | `num_threads=0`, `dtype=None`, `device="cpu"` | PyTorch tensor |
| `molecules_pandas(molecules, ...)` | `kind="ecfp"`, `fp_size=2048`, `radius=2`, `path_length=7`, `options=None` | `num_threads=0` | pandas data frame |
| `molecules_polars(molecules, ...)` | `kind="ecfp"`, `fp_size=2048`, `radius=2`, `path_length=7`, `options=None` | `num_threads=0` | Polars data frame |

All arguments after `molecules` are keywords. The fingerprint arguments mean
what they mean for `smsd.fingerprints`, and `options`, when given, sets the
fingerprint completely, including its width and count mode. `num_threads=0`
uses the OpenMP default. A single SMILES string passed as `molecules` raises
`TypeError`; wrap it in a list.

### Conversion functions

| Function | Arguments | Returns |
| --- | --- | --- |
| `to_numpy(sparse_fingerprints, *, fp_size=2048, counts=False, dtype=None)` | Positions, or `(position, count)` pairs when `counts=True` | NumPy matrix; boolean, integer or floating-point `dtype` |
| `to_scipy(sparse_fingerprints, *, fp_size=2048, counts=False, dtype=None)` | As `to_numpy` | CSR matrix; boolean or integer `dtype` |
| `to_torch(sparse_fingerprints, *, fp_size=2048, counts=False, dtype=None, device="cpu")` | As `to_numpy` | Tensor; a `torch.dtype` |
| `to_arrow(fingerprints, *, fp_size=2048, counts=False, dtype=None, copy=False, zero_copy_only=False)` | A NumPy matrix of width `fp_size`, or sparse rows | Arrow table with schema metadata |
| `to_pyg(molecule, *, device="cpu")` | One SMILES string or `MolGraph` | PyTorch Geometric `Data` |

The conversion functions read sparse rows only. A dense row of zeros and ones
passed to `to_numpy`, `to_scipy` or `to_torch` raises no error: it is read as
positions, so `[[0, 1, 0, 1]]` becomes bits 0 and 1. Positions outside
`fp_size`, negative counts and counts that do not fit the chosen type raise
`ValueError`; repeated count bins are added together. An empty batch gives
shape `(0, fp_size)`.

### Arrow schema metadata

| Key | Values | Meaning |
| --- | --- | --- |
| `smsd.fp_size` | The width | Fingerprint width in bits or bins |
| `smsd.counts` | `true`, `false` | Whether the values are counts |
| `smsd.zero_copy` | `true`, `false` | Whether the column shares the caller's matrix |
| `smsd.copy_reason` | `none`, `layout_or_boolean`, `sparse_expansion`, `requested` | Why a copy was made |
| `smsd.ownership` | `shared_numpy_owner_retained`, `owned_numpy_owner_retained` | Whose matrix the table keeps alive |

### Transformers

| Transformer | Parameters, with defaults | Fitted attributes |
| --- | --- | --- |
| `FingerprintTransformer` | `kind="ecfp"`, `fp_size=2048`, `radius=2`, `path_length=7`, `counts=False`, `dtype=None`, `num_threads=0` | `is_fitted_` |
| `SimilarityTransformer` | `kind="ecfp"`, `fp_size=2048`, `radius=2`, `path_length=7`, `num_threads=0` | `reference_fingerprints_`, `n_references_` |

Both are keyword-only and both name their output columns through
`get_feature_names_out()`: `fp_0`, `fp_1` and so on, or `reference_0`,
`reference_1` and so on. `FingerprintTransformer(counts=True)` supports the
ECFP, FCFP, path and torsion families. `SimilarityTransformer` always uses
binary fingerprints, raises `ValueError` when fitted with no references, and
gives a similarity of 1 between two empty fingerprints.

## Pitfalls

- Fingerprint family. `fingerprints_numpy` and `fingerprints_torch` default to
  `kind="mcs"`, while `fingerprints_scipy`, the data-frame adapters and both
  transformers default to `kind="ecfp"`. Pass the kind or a
  `FingerprintOptions` everywhere, so that training features and screening
  features come from the same fingerprint.
- One molecule. A bare SMILES string raises `TypeError` rather than being read
  character by character; pass `["CCO"]`, not `"CCO"`.
- Input types. `molecules_pandas`, `molecules_polars`, `SimilarityTransformer`
  and `to_pyg` accept only SMILES strings and `MolGraph` objects. Missing
  values in a SMILES column raise `TypeError`; drop or fill them first.
- Data types. Bits arrive as `uint8` and counts as `int64`. An estimator that
  works in `float64` converts the matrix, which makes a copy eight times the
  size of the `uint8` matrix. Choose `dtype` once, at the adapter, rather than
  converting large matrices repeatedly. A floating-point `dtype` can lose
  precision for very large counts.
- Memory. A dense 2,048-bit matrix needs 2,048 bytes per molecule as `uint8`
  and 16,384 bytes as `int64` or `float64`: about 2 GB and 16 GB for a million
  molecules. Use `fingerprints_scipy`, smaller chunks or a
  [packed fingerprint index](PACKED_INDEX.md) beyond that. The output of
  `SimilarityTransformer` grows with queries times references.
- Copies. Arrow shares only C-contiguous, native-endian integer matrices.
  Slicing columns, transposing or converting to boolean before `to_arrow`
  forces a copy. A shared table changes when its matrix changes.
- Bit order. Matrices hold one column per position. Packed bytes in an index
  are little-endian within each byte, unlike NumPy's default `packbits` order.
- Width. Positions are hashed into `fp_size` columns, so different structures
  can share a column. A width that is too narrow for the molecules saturates
  the fingerprint; [fingerprint quality](FINGERPRINTS.md#fingerprint-quality)
  shows how to check.
- Refitting. Changing a `SimilarityTransformer` parameter after `fit` makes
  `transform` raise until it is fitted again. `FingerprintTransformer` reads
  its parameters at every `transform`.
- Feature names in a `ColumnTransformer`. A `FingerprintTransformer` inside a
  `ColumnTransformer` fits and transforms normally, but the
  `ColumnTransformer`'s `get_feature_names_out()` raises `ValueError`, because
  the transformer refuses input column names. Call `get_feature_names_out()` on
  the fingerprint transformer itself.
- Devices. A `device` argument places the finished tensor. Fingerprints and
  graph features are always computed on the CPU, and `smsd-cuda` plays no part
  in these adapters; see [What runs where](CPU_GPU.md#what-runs-where).

## See also

- [Python feature adapters and MCP](../PYTHON_ADAPTERS.md): the reference for
  every adapter, including SciPy, Arrow and PyTorch Geometric.
- [Packed corpora, clustering and optional data adapters](../API_REFERENCE.md#packed-corpora-clustering-and-optional-data-adapters)
  in the API reference.
- [Fingerprints and similarity](FINGERPRINTS.md): fingerprint families,
  configurable fingerprints and similarity coefficients.
- [Screening a library](SCREENING.md) and
  [Packed fingerprint index](PACKED_INDEX.md), for libraries too large for a
  feature matrix.
- The [feature-search notebook](../../python/examples/feature_search.ipynb),
  which builds a small scikit-learn pipeline beside an MCS search and a packed
  index.
