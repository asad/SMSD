# SMSD Python Guide

SMSD for Python exposes the native C++ engines for substructure, MCS,
fingerprinting, SMARTS matching, stereo/CIP assignment, layout, and chemistry
I/O. The default path is native SMSD with no RDKit dependency. RDKit can still
be used alongside SMSD for interop, cross-validation, or depiction workflows.

## Support

- CPython `3.9` through the latest stable release series
- Recommended default runtime: `Python 3.12`
- CPU-only execution by default
- Optional GPU acceleration:
  - macOS: Metal
  - Linux/Windows: CUDA
- CPU fallback remains available when no GPU backend is compiled

## Install

```bash
pip install smsd
```

Optional extras:

```bash
pip install "smsd[rdkit]"
pip install "smsd[gpu]"
```

## Core Search

```python
import smsd

query = smsd.parse_smiles("c1ccccc1")
target = smsd.parse_smiles("c1ccc(O)cc1")

assert smsd.is_substructure(query, target)
mapping = smsd.find_substructure(query, target)
mcs = smsd.find_mcs(query, target)
```

Convenience wrappers accept `MolGraph`, SMILES, or RDKit molecules:

```python
mcs = smsd.mcs("CC(=O)O", "CCC(=O)O")
hit = smsd.substructure_search("c1ccccc1", "c1ccc(O)cc1")
```

## SMARTS

```python
import smsd

assert smsd.smarts_match("[#6]~[#7]", "CCN")
matches = smsd.smarts_find_all("[OH]", "c1ccc(O)cc1")
largest = smsd.find_mcs_smarts("[CX3](=O)[NX3]", "CC(=O)Nc1ccc(O)cc1")
```

SMSD supports ring primitives (`R`, `r`, `x`), recursive SMARTS, logical
operators, aromatic atoms, charges, isotopes, and stereo-aware query handling
through the native matcher.

## Fingerprints

```python
import smsd

path_fp = smsd.path_fingerprint("c1ccccc1", path_length=7, fp_size=2048)
mcs_fp = smsd.mcs_fingerprint("c1ccc(O)cc1", path_length=7, fp_size=2048)
ecfp4 = smsd.circular_fingerprint("c1ccccc1", radius=2, fp_size=2048)
torsion = smsd.topological_torsion("CC(N)C(=O)O")

assert smsd.fingerprint_subset(path_fp, path_fp)
score = smsd.tanimoto(ecfp4, ecfp4)
```

## Stereo and CIP

```python
import smsd

g = smsd.parse_smiles("N[C@@H](C)C(=O)O")
rs = smsd.assign_rs(g)
ez = smsd.assign_ez(smsd.parse_smiles("C/C=C/C"))
```

## Native MOL/SDF I/O

Use SMSD’s native readers/writers by default:

```python
import smsd

g = smsd.read_molfile("molecule.mol")
mol_block = smsd.write_mol_block(g)
mol_block_v3000 = smsd.write_mol_block_v3000(g)
sdf_record = smsd.write_sdf_record(g)

smsd.write_molfile(g, "out_v2000.mol")
smsd.write_molfile(g, "out_v3000.mol", v3000=True)
smsd.write_molfile(g, "out.sdf", sdf=True)
smsd.export_sdf([g], "library.sdf")
```

Covered natively in `6.9.1`:
- MOL V2000 and V3000 graph round-trip
- names, program line, comments, and SDF properties
- charges, isotopes, atom classes, atom maps
- `R#`/`R<n>` and `M  RGP`
- practical stereo flags

Not yet a complete native replacement for every exotic MDL query feature:
- full atom lists and not-lists
- link nodes
- variable attachment semantics
- the complete Markush query feature set

## R-group Fragments

```python
import smsd

frag = smsd.parse_smiles("[R1]c1ccccc1")
assert smsd.to_smiles(frag).startswith("[R1]")

rgroups = smsd.decompose_r_groups(
    "c1ccccc1",
    ["c1ccc(O)cc1", "c1ccc(N)cc1"],
)
```

## GPU and CPU Backends

```python
import smsd

print(smsd.gpu_is_available())
print(smsd.gpu_device_info())
```

GPU support is opportunistic. If a GPU backend is unavailable, the same API
continues on CPU with no code changes.

## Third-party Interop

SMSD does not require RDKit or CDK. RDKit remains optional for:
- drawing and notebook depiction
- cross-validation
- mixed-toolkit workflows

Open Babel is also supported as an optional interop layer for users who prefer
that stack:

```python
import smsd
from openbabel import pybel

obmol = pybel.readstring("smi", "c1ccc(O)cc1")
g = smsd.from_openbabel(obmol)
back = smsd.to_openbabel(g)
```

When needed:

```python
import smsd
from rdkit import Chem

rdmol = Chem.MolFromSmiles("c1ccccc1")
g = smsd.from_rdkit(rdmol)
back = smsd.to_rdkit(g)
```
