# Clustering and diversity

Clustering and diverse selection turn a library into something a chemist can
review: a set of series, each led by a representative, or a short list of
molecules that cover the library's chemistry. SMSD Pro offers three exact
tools on a packed fingerprint index: `smsd.butina` groups rows around cluster
centres, `smsd.maxmin` picks a fixed number of diverse rows, and
`FingerprintCorpus.top_k_pairs` finds the nearest neighbours of every row.
This chapter shows which to use for which task, how the threshold, reordering,
seed and chunk size behave, what each costs in time and memory, and what runs
on a GPU. Every example uses the inline library of
[Screening a library](SCREENING.md): 14 benzodiazepines, 7 common drugs and
512 generated trisubstituted benzenes.

All three tools read a `FingerprintCorpus`, so they need NumPy, from the
`numpy` extra:

```bash
python -m pip install "smsd[numpy]==8.2.0"
```

## Choose a method

| Task | Use | Returns | Comparison work for n rows | Extra memory |
| --- | --- | --- | --- | --- |
| Group a library into series, each with a representative | `smsd.butina(corpus, threshold)` | A list of tuples, centre first | About n²/2, plus one pass over the library per cluster; `reordering=True` adds about n² | Two arrays of n entries, about 9 bytes per row |
| Pick a fixed number of diverse rows | `smsd.maxmin(corpus, count, seed=0)` | A list of rows in the order picked | n for each row picked | One array of n entries, 8 bytes per row |
| Find the k nearest neighbours of every row, near-duplicates or outliers | `corpus.top_k_pairs(k, threshold=0.0)` | Yields `(row, neighbour, Tanimoto)`, at most k per row, best first | About n² | k entries per row in progress; n × k if every result is kept |
| List every pair above a similarity | `corpus.threshold_pairs(threshold)` | Yields `(first, second, Tanimoto)`, each pair once | About n²/2 | Up to n²/2 pairs if every result is kept |

Some guidance on the choice:

- Use Butina when the question is "which series does this library contain?".
  The number of clusters follows from the threshold, not from a target count.
- Use MaxMin when the question is "which N molecules should we buy, test or
  inspect first?". The count is fixed and the threshold is implicit.
- Use `top_k_pairs` to inspect a library before clustering it: identical
  fingerprints, close analogues and molecules with no near neighbour all show
  up in its output.
- Every method compares fingerprints exactly; none samples the library or
  approximates a score. The work therefore grows with the square of the row
  count, or with the row count times the number picked for MaxMin, and the
  results are the same on every run, every chunk size and every backend.

## Prepare the input

The three tools take an open `FingerprintCorpus`, built once with
`FingerprintCorpus.from_smiles` or `FingerprintCorpus.create`, as
[Packed fingerprint index](PACKED_INDEX.md#build-an-index) describes. Four
points matter before clustering:

- The fingerprint decides the answer. A cluster threshold of 0.5 on ECFP4
  bits is a different statement from 0.5 on path bits. The index file records
  the width and row count but not the options, so keep the
  `FingerprintOptions` beside it and state them with every result.
- Rows are numbers. The index stores no names, and every returned row refers
  to the input order. Keep the identifiers in a list in the same order.
- Duplicates are counted. Identical fingerprints score 1 with each other, so a
  series entered many times has many neighbours and tends to supply the first
  and largest clusters. Remove duplicate structures first if the counts matter;
  see [Prepare the library](SCREENING.md#prepare-the-library).
- Empty fingerprints score 1 with each other. Empty graphs from malformed SD
  records therefore fall into one cluster, and they score 0 against every
  other row. Drop graphs with no atoms before building the index.

## Cluster a library with Butina

`smsd.butina(corpus, threshold=0.7)` counts, for every row, the other rows
whose Tanimoto similarity is at least the threshold. The row with the most
neighbours becomes the first centre and takes all its unassigned neighbours
as members; the next centre is the unassigned row with the most neighbours,
and so on until every row is assigned. Equal counts go to the lower row. The
result is a list of tuples, each with its centre first and the other members
in row order, and every row appears in exactly one tuple.

A member reaches the threshold with its centre, but not necessarily with the
other members of its cluster. A row also joins the first cluster whose centre
reaches it, even when a later centre would be closer. The program below
clusters the whole library at 0.5 and checks both properties:

```python
import itertools
from pathlib import Path
from tempfile import TemporaryDirectory

import smsd

NAMED = {
    "diazepam": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
    "nordazepam": "O=C1CN=C(c2ccccc2)c2cc(Cl)ccc2N1",
    "temazepam": "CN1C(=O)C(O)N=C(c2ccccc2)c2cc(Cl)ccc21",
    "oxazepam": "OC1N=C(c2ccccc2)c2cc(Cl)ccc2NC1=O",
    "lorazepam": "OC1N=C(c2ccccc2Cl)c2cc(Cl)ccc2NC1=O",
    "clonazepam": "O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1",
    "nitrazepam": "O=C1CN=C(c2ccccc2)c2cc([N+](=O)[O-])ccc2N1",
    "flunitrazepam": "CN1C(=O)CN=C(c2ccccc2F)c2cc([N+](=O)[O-])ccc21",
    "bromazepam": "O=C1CN=C(c2ccccn2)c2cc(Br)ccc2N1",
    "tetrazepam": "CN1C(=O)CN=C(C2=CCCCC2)c2cc(Cl)ccc21",
    "alprazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2",
    "triazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1Cl)=NC2",
    "midazolam": "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
    "clobazam": "CN1C(=O)CC(=O)N(c2ccccc2)c2cc(Cl)ccc21",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "salicylic acid": "OC(=O)c1ccccc1O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "phenol": "Oc1ccccc1",
}
GROUPS = ["C", "CC", "OC", "N", "O", "F", "Cl", "C(=O)O"]
library = list(NAMED.items()) + [
    (f"GEN-{number:03d}", f"c1c({a})cc({b})cc1{c}")
    for number, (a, b, c) in enumerate(itertools.product(GROUPS, repeat=3))
]
names = [name for name, _ in library]
smiles = [text for _, text in library]
assert len(library) == 533

options = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
with TemporaryDirectory() as directory:
    path = Path(directory) / "library.fp"
    with smsd.FingerprintCorpus.from_smiles(path, smiles,
                                            options=options) as corpus:
        clusters = smsd.butina(corpus, threshold=0.5)
        fingerprints = [corpus.fingerprint(row) for row in range(len(corpus))]

# Every row is in exactly one cluster ...
assert sorted(row for cluster in clusters for row in cluster) == list(range(533))
# ... and every member reaches the threshold with its centre.
for centre, *members in clusters:
    for member in members:
        assert smsd.tanimoto_coefficient(fingerprints[centre],
                                         fingerprints[member]) >= 0.5

sizes = [len(cluster) for cluster in clusters]
print("clusters:", len(clusters), "| singletons:", sizes.count(1))
print("first three:", [(names[c[0]], len(c)) for c in clusters[:3]])
print("clusters that hold the named drugs:")
for cluster in clusters:
    if cluster[0] < 21:
        print(f"  {names[cluster[0]]:<15}", [names[row] for row in cluster[1:]])

# Bromazepam reaches 0.5 with nordazepam, but nordazepam was taken first.
nordazepam, bromazepam = names.index("nordazepam"), names.index("bromazepam")
score = smsd.tanimoto_coefficient(fingerprints[nordazepam],
                                  fingerprints[bromazepam])
print(f"nordazepam and bromazepam: {score:.3f}")
assert len(clusters) == 34 and sizes.count(1) == 11
assert (bromazepam,) in clusters and round(score, 3) == 0.519
```

The program prints:

```text
clusters: 34 | singletons: 11
first three: [('GEN-511', 67), ('GEN-019', 66), ('GEN-087', 42)]
clusters that hold the named drugs:
  diazepam        ['nordazepam', 'temazepam', 'tetrazepam', 'clobazam']
  oxazepam        ['lorazepam']
  nitrazepam      ['clonazepam']
  alprazolam      ['triazolam', 'midazolam']
  bromazepam      []
  flunitrazepam   []
  aspirin         []
  salicylic acid  []
  paracetamol     []
  ibuprofen       []
  naproxen        []
  caffeine        []
  phenol          []
nordazepam and bromazepam: 0.519
```

The generated benzenes supply the first and largest clusters because each of
their 120 distinct structures occurs several times, so each has many identical
neighbours. Bromazepam shows the order effect: it reaches the threshold with
nordazepam, but nordazepam had already joined diazepam's cluster, so
bromazepam forms a cluster of its own.

## Thresholds and reordering

The threshold is a Tanimoto similarity, inclusive, from 0 to 1; a higher
threshold gives more, tighter clusters. With the default `reordering=False`
the neighbour counts are computed once, before any cluster is taken, so a
later centre may be chosen for neighbours that earlier clusters have already
claimed. `reordering=True` recounts after each cluster so that only unassigned
neighbours count, and later centres are chosen by the neighbours still
available. It costs up to about three times the comparisons. Neither mode
is more exact than the other; they answer slightly different questions, so
record which one a result used.

```python
import itertools
from pathlib import Path
from tempfile import TemporaryDirectory

import smsd

NAMED = {
    "diazepam": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
    "nordazepam": "O=C1CN=C(c2ccccc2)c2cc(Cl)ccc2N1",
    "temazepam": "CN1C(=O)C(O)N=C(c2ccccc2)c2cc(Cl)ccc21",
    "oxazepam": "OC1N=C(c2ccccc2)c2cc(Cl)ccc2NC1=O",
    "lorazepam": "OC1N=C(c2ccccc2Cl)c2cc(Cl)ccc2NC1=O",
    "clonazepam": "O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1",
    "nitrazepam": "O=C1CN=C(c2ccccc2)c2cc([N+](=O)[O-])ccc2N1",
    "flunitrazepam": "CN1C(=O)CN=C(c2ccccc2F)c2cc([N+](=O)[O-])ccc21",
    "bromazepam": "O=C1CN=C(c2ccccn2)c2cc(Br)ccc2N1",
    "tetrazepam": "CN1C(=O)CN=C(C2=CCCCC2)c2cc(Cl)ccc21",
    "alprazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2",
    "triazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1Cl)=NC2",
    "midazolam": "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
    "clobazam": "CN1C(=O)CC(=O)N(c2ccccc2)c2cc(Cl)ccc21",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "salicylic acid": "OC(=O)c1ccccc1O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "phenol": "Oc1ccccc1",
}
GROUPS = ["C", "CC", "OC", "N", "O", "F", "Cl", "C(=O)O"]
library = list(NAMED.items()) + [
    (f"GEN-{number:03d}", f"c1c({a})cc({b})cc1{c}")
    for number, (a, b, c) in enumerate(itertools.product(GROUPS, repeat=3))
]
names = [name for name, _ in library]
smiles = [text for _, text in library]
assert len(library) == 533

options = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
rows = []
with TemporaryDirectory() as directory:
    path = Path(directory) / "library.fp"
    with smsd.FingerprintCorpus.from_smiles(path, smiles,
                                            options=options) as corpus:
        for threshold in (0.3, 0.4, 0.5, 0.6, 0.7):
            initial = smsd.butina(corpus, threshold=threshold)
            recount = smsd.butina(corpus, threshold=threshold, reordering=True)
            rows.append((threshold, initial, recount))

print("threshold  clusters  singletons  | reordering: clusters  singletons")
for threshold, initial, recount in rows:
    print(f"{threshold:>9}  {len(initial):>8}  "
          f"{sum(len(c) == 1 for c in initial):>10}  |"
          f"{len(recount):>21}  {sum(len(c) == 1 for c in recount):>10}")

assert [len(initial) for _, initial, _ in rows] == [19, 28, 34, 53, 79]
assert [len(recount) for _, _, recount in rows] == [18, 29, 33, 50, 77]
```

The program prints:

```text
threshold  clusters  singletons  | reordering: clusters  singletons
      0.3        19           7  |                   18           8
      0.4        28          11  |                   29          12
      0.5        34          11  |                   33          11
      0.6        53          22  |                   50          20
      0.7        79          21  |                   77          21
```

There is no universally right threshold: the right value for a fingerprint
and a library is the one whose clusters a chemist accepts on review. Look at
a sample of clusters at two or three thresholds before fixing one, and keep it
unchanged when results from different runs are to be compared.

## Pick a diverse subset with MaxMin

`smsd.maxmin(corpus, count, seed=0)` starts at row `seed`, then repeatedly
takes the row whose closest already-picked row is the least similar, that is,
the row with the largest minimum Tanimoto distance (1 − Tanimoto) to the
selection. Equal distances go to the lower row, and no random numbers are
used, so the same index, count and seed always give the same list. The list is
in the order picked, and a shorter run is a prefix of a longer one: the first
three rows of a pick of eight are exactly the pick of three.

The seed decides where the selection starts. Row 0 is the default, which makes
the result depend on the input order. Choose a seed that means something, such
as a reference compound or the centre of the largest Butina cluster, and record
it.

```python
import itertools
from pathlib import Path
from tempfile import TemporaryDirectory

import smsd

NAMED = {
    "diazepam": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
    "nordazepam": "O=C1CN=C(c2ccccc2)c2cc(Cl)ccc2N1",
    "temazepam": "CN1C(=O)C(O)N=C(c2ccccc2)c2cc(Cl)ccc21",
    "oxazepam": "OC1N=C(c2ccccc2)c2cc(Cl)ccc2NC1=O",
    "lorazepam": "OC1N=C(c2ccccc2Cl)c2cc(Cl)ccc2NC1=O",
    "clonazepam": "O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1",
    "nitrazepam": "O=C1CN=C(c2ccccc2)c2cc([N+](=O)[O-])ccc2N1",
    "flunitrazepam": "CN1C(=O)CN=C(c2ccccc2F)c2cc([N+](=O)[O-])ccc21",
    "bromazepam": "O=C1CN=C(c2ccccn2)c2cc(Br)ccc2N1",
    "tetrazepam": "CN1C(=O)CN=C(C2=CCCCC2)c2cc(Cl)ccc21",
    "alprazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2",
    "triazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1Cl)=NC2",
    "midazolam": "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
    "clobazam": "CN1C(=O)CC(=O)N(c2ccccc2)c2cc(Cl)ccc21",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "salicylic acid": "OC(=O)c1ccccc1O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "phenol": "Oc1ccccc1",
}
GROUPS = ["C", "CC", "OC", "N", "O", "F", "Cl", "C(=O)O"]
library = list(NAMED.items()) + [
    (f"GEN-{number:03d}", f"c1c({a})cc({b})cc1{c}")
    for number, (a, b, c) in enumerate(itertools.product(GROUPS, repeat=3))
]
names = [name for name, _ in library]
smiles = [text for _, text in library]
assert len(library) == 533


def closest_pair(fingerprints, rows):
    """The highest Tanimoto similarity between any two of the given rows."""
    return max(smsd.tanimoto_coefficient(fingerprints[a], fingerprints[b])
               for a, b in itertools.combinations(rows, 2))


options = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
with TemporaryDirectory() as directory:
    path = Path(directory) / "library.fp"
    with smsd.FingerprintCorpus.from_smiles(path, smiles,
                                            options=options) as corpus:
        first = smsd.maxmin(corpus, 8)  # starts at row 0, diazepam
        phenol = smsd.maxmin(corpus, 8, seed=names.index("phenol"))
        assert smsd.maxmin(corpus, 3) == first[:3]  # a prefix of the longer pick
        assert smsd.maxmin(corpus, 0) == []
        fingerprints = [corpus.fingerprint(row) for row in range(len(corpus))]

print("from diazepam:", [names[row] for row in first])
print("from phenol:  ", [names[row] for row in phenol])
print(f"closest pair in the first pick: {closest_pair(fingerprints, first):.3f}")
print(f"closest pair in rows 0 to 7:    "
      f"{closest_pair(fingerprints, range(8)):.3f}")
assert [names[row] for row in first[:3]] == ["diazepam", "GEN-148", "caffeine"]
assert phenol[0] == names.index("phenol") and len(set(phenol)) == 8
```

The program prints:

```text
from diazepam: ['diazepam', 'GEN-148', 'caffeine', 'GEN-175', 'paracetamol', 'GEN-219', 'naproxen', 'GEN-294']
from phenol:   ['phenol', 'flunitrazepam', 'caffeine', 'GEN-222', 'GEN-164', 'ibuprofen', 'GEN-077', 'GEN-167']
closest pair in the first pick: 0.113
closest pair in rows 0 to 7:    0.755
```

Eight rows picked by MaxMin share at most 0.113 with each other, whereas the
first eight rows of the file, all benzodiazepines, include a pair at 0.755.
Two seeds give different selections of similar quality, which is expected:
MaxMin is a greedy pick, not a search for the most diverse subset of a given
size. `count` may be any number from 0 to the row count; a seed past the last
row raises `IndexError`, and a negative one `ValueError`.

MaxMin favours the edges of the library: its first picks are the molecules
least like everything already chosen, which often includes outliers and
oddities. To cover the populated regions instead, cluster first and take the
Butina centres, or run MaxMin on an index built from the centres.

## Nearest neighbours and pairs

`corpus.top_k_pairs(k, threshold=0.0)` yields the k most similar other rows of
every row as `(row, neighbour, Tanimoto)`. The output is grouped by row in row
order, best first within each row, and equal scores prefer the lower
neighbour. `threshold` drops neighbours below it, so a row with no neighbour
at or above the threshold yields nothing: that is a direct way to find
outliers. `include_self=False` leaves out each row's pairing with itself, and
k is reduced to the number of other rows when it is larger. The method is a
generator, so consume it as it streams when the library is large.

```python
import itertools
from pathlib import Path
from tempfile import TemporaryDirectory

import smsd

NAMED = {
    "diazepam": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
    "nordazepam": "O=C1CN=C(c2ccccc2)c2cc(Cl)ccc2N1",
    "temazepam": "CN1C(=O)C(O)N=C(c2ccccc2)c2cc(Cl)ccc21",
    "oxazepam": "OC1N=C(c2ccccc2)c2cc(Cl)ccc2NC1=O",
    "lorazepam": "OC1N=C(c2ccccc2Cl)c2cc(Cl)ccc2NC1=O",
    "clonazepam": "O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1",
    "nitrazepam": "O=C1CN=C(c2ccccc2)c2cc([N+](=O)[O-])ccc2N1",
    "flunitrazepam": "CN1C(=O)CN=C(c2ccccc2F)c2cc([N+](=O)[O-])ccc21",
    "bromazepam": "O=C1CN=C(c2ccccn2)c2cc(Br)ccc2N1",
    "tetrazepam": "CN1C(=O)CN=C(C2=CCCCC2)c2cc(Cl)ccc21",
    "alprazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2",
    "triazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1Cl)=NC2",
    "midazolam": "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
    "clobazam": "CN1C(=O)CC(=O)N(c2ccccc2)c2cc(Cl)ccc21",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "salicylic acid": "OC(=O)c1ccccc1O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "phenol": "Oc1ccccc1",
}
GROUPS = ["C", "CC", "OC", "N", "O", "F", "Cl", "C(=O)O"]
library = list(NAMED.items()) + [
    (f"GEN-{number:03d}", f"c1c({a})cc({b})cc1{c}")
    for number, (a, b, c) in enumerate(itertools.product(GROUPS, repeat=3))
]
names = [name for name, _ in library]
smiles = [text for _, text in library]
assert len(library) == 533

options = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
with TemporaryDirectory() as directory:
    path = Path(directory) / "library.fp"
    with smsd.FingerprintCorpus.from_smiles(path, smiles,
                                            options=options) as corpus:
        # The three nearest neighbours of the first two rows.
        nearest = [pair for pair in corpus.top_k_pairs(3) if pair[0] < 2]
        # Rows whose nearest neighbour has an identical fingerprint.
        twins = sum(score == 1.0 for _, _, score in corpus.top_k_pairs(1))
        # Rows with no neighbour at 0.4 or above.
        close = {row for row, _, _ in corpus.top_k_pairs(1, threshold=0.4)}

for row, neighbour, score in nearest:
    print(f"{names[row]:<11} -> {names[neighbour]:<11} {score:.3f}")
print("rows with an identical twin:", twins)
print("rows with no neighbour at 0.4:",
      [names[row] for row in range(len(names)) if row not in close])
assert len(nearest) == 6 and twins == 504 and len(close) == 526
```

The program prints:

```text
diazepam    -> temazepam   0.582
diazepam    -> tetrazepam  0.542
diazepam    -> nordazepam  0.528
nordazepam  -> oxazepam    0.558
nordazepam  -> diazepam    0.528
nordazepam  -> bromazepam  0.519
rows with an identical twin: 504
rows with no neighbour at 0.4: ['aspirin', 'salicylic acid', 'paracetamol', 'ibuprofen', 'naproxen', 'caffeine', 'phenol']
```

The 504 twins are the repeated generated benzenes; the seven rows with no
close neighbour are the common drugs, which have no analogue in this library.
Both findings are worth having before clustering: the first says the counts
behind Butina centres are inflated by duplicates, and the second says which
rows will end up as singletons at 0.4 or above.

`top_k_pairs` reports each relationship from both sides, as two directed
pairs. When every pair above a similarity is wanted once, use
`corpus.threshold_pairs(threshold)`, which yields `(first, second, Tanimoto)`
with `first < second`; [Neighbours and pairs within the
index](PACKED_INDEX.md#neighbours-and-pairs-within-the-index) shows it.

## Chunk size, determinism and cost

`chunk_size` sets how many rows are compared at a time. It changes the working
memory, never the clusters, selections or neighbours. The default of 4096 rows
is reduced automatically to keep each chunk's working arrays within about
64 MiB, so it rarely needs changing; lower it only to cut peak memory on a
small machine. Results do not depend on the chunk size, the backend or the
run, because every score is computed exactly and every tie goes to the lower
row.

The comparison work is what limits library size:

- Butina compares every pair once to count neighbours, about n²/2 comparisons,
  then makes one pass over the library for each cluster it forms. With
  `reordering=True` it also makes a pass for every row it assigns, which adds
  about n² comparisons. Doubling the library roughly quadruples the time.
- MaxMin makes one pass over the library for each row it picks, n × count
  comparisons, so it grows linearly with the library for a fixed count.
- `top_k_pairs` compares every row with every other row, about n²
  comparisons, because each row is ranked against the whole library.

Memory stays modest. Butina keeps two arrays of n entries, about 9 bytes per
row, and MaxMin one array of 8 bytes per row; neither builds a similarity
matrix or a neighbour list. The index file is memory-mapped, at
`ceil(fp_size / 8)` bytes per row: 256 bytes for 2048 bits. The returned
clusters and any pairs you keep are ordinary Python objects and come on top;
collecting all of `top_k_pairs` needs n × k results, and all of
`threshold_pairs` can need memory that grows with n². Time a representative
sample on your own hardware before committing to a library size, as the
[notes of Screening a library](SCREENING.md#notes) advise.

## CPU and GPU backends

`smsd.butina`, `smsd.maxmin`, `top_k_pairs` and `threshold_pairs` take the
same `backend` argument as the other searches of a packed index:

| Backend | Behaviour |
| --- | --- |
| `"cpu"` | The default. NumPy compares the packed rows on the CPU; the CUDA package is not loaded. |
| `"cuda"` | The comparisons run on an NVIDIA GPU. It needs the `smsd-cuda` package of the same version, a usable device and a little-endian host, and raises `RuntimeError` otherwise; fingerprints wider than 65,536 bits raise `ValueError`. |
| `"auto"` | Uses the GPU under the same conditions, for an index of any size, and the CPU otherwise. |

The GPU takes only the fingerprint comparisons. Every score, threshold test and
tie break is settled identically on either device, so clusters, selections,
neighbours and scores are the same on the GPU as on the CPU; the GPU does not
make any of these methods approximate. The index is not kept on the GPU between
calls, so a Butina run sends the library to the device once for every pass it
makes.

`corpus.cuda_calls` is the cumulative number of packed CUDA calls that the
corpus has completed, clustering passes included. Read it before and after a
call: an increase proves that the device did work for that call. It is not a
count of queries, hits or kernels, and it says nothing about speed. Under
`"auto"`, the CPU is used only when the CUDA package or a device is
unavailable, or the width is unsupported; once device work has started, an
error is raised rather than the call being retried on the CPU. Setting
`SMSD_DISABLE_GPU` to any value other than `0` keeps these calls on the CPU,
and on a Mac they always run on the CPU.

The program below runs on any machine. It checks that the chunk size and the
`"auto"` backend leave every result unchanged, and, where `backend="cuda"` can
run, that the GPU gives the same results:

```python
import itertools
import os
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import smsd

try:
    import smsd_cuda
except ImportError:
    smsd_cuda = None

NAMED = {
    "diazepam": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
    "nordazepam": "O=C1CN=C(c2ccccc2)c2cc(Cl)ccc2N1",
    "temazepam": "CN1C(=O)C(O)N=C(c2ccccc2)c2cc(Cl)ccc21",
    "oxazepam": "OC1N=C(c2ccccc2)c2cc(Cl)ccc2NC1=O",
    "lorazepam": "OC1N=C(c2ccccc2Cl)c2cc(Cl)ccc2NC1=O",
    "clonazepam": "O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1",
    "nitrazepam": "O=C1CN=C(c2ccccc2)c2cc([N+](=O)[O-])ccc2N1",
    "flunitrazepam": "CN1C(=O)CN=C(c2ccccc2F)c2cc([N+](=O)[O-])ccc21",
    "bromazepam": "O=C1CN=C(c2ccccn2)c2cc(Br)ccc2N1",
    "tetrazepam": "CN1C(=O)CN=C(C2=CCCCC2)c2cc(Cl)ccc21",
    "alprazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2",
    "triazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1Cl)=NC2",
    "midazolam": "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
    "clobazam": "CN1C(=O)CC(=O)N(c2ccccc2)c2cc(Cl)ccc21",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "salicylic acid": "OC(=O)c1ccccc1O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "phenol": "Oc1ccccc1",
}
GROUPS = ["C", "CC", "OC", "N", "O", "F", "Cl", "C(=O)O"]
library = list(NAMED.items()) + [
    (f"GEN-{number:03d}", f"c1c({a})cc({b})cc1{c}")
    for number, (a, b, c) in enumerate(itertools.product(GROUPS, repeat=3))
]
names = [name for name, _ in library]
smiles = [text for _, text in library]
assert len(library) == 533


def cuda_usable():
    """Whether backend="cuda" can run here."""
    return (smsd_cuda is not None and sys.byteorder == "little"
            and os.environ.get("SMSD_DISABLE_GPU", "") in ("", "0")
            and smsd_cuda.__version__ == smsd.__version__
            and smsd_cuda.is_available())


def results(corpus, **controls):
    """Butina clusters, a MaxMin pick and the two nearest neighbours of each row."""
    return (smsd.butina(corpus, threshold=0.5, **controls),
            smsd.maxmin(corpus, 8, **controls),
            list(corpus.top_k_pairs(2, **controls)))


options = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
with TemporaryDirectory() as directory:
    path = Path(directory) / "library.fp"
    with smsd.FingerprintCorpus.from_smiles(path, smiles,
                                            options=options) as corpus:
        reference = results(corpus)  # backend="cpu", chunk_size=4096
        print("CUDA calls after the CPU run:", corpus.cuda_calls)
        assert corpus.cuda_calls == 0
        for chunk_size in (7, 100):
            assert results(corpus, chunk_size=chunk_size) == reference
        print("chunk sizes 7, 100 and 4096: identical results")

        # "auto" may run on either device, so check its results, not where.
        assert results(corpus, backend="auto") == reference
        print("backend='auto': identical results")

        if cuda_usable():
            before = corpus.cuda_calls
            assert results(corpus, backend="cuda") == reference
            assert corpus.cuda_calls > before  # the device did the comparisons
```

The program prints:

```text
CUDA calls after the CPU run: 0
chunk sizes 7, 100 and 4096: identical results
backend='auto': identical results
```

On a machine without a usable GPU the `"auto"` run takes the CPU path and the
branch under `if cuda_usable():` is skipped. On a machine with `smsd-cuda` and
a supported NVIDIA GPU, the `"auto"` and `"cuda"` runs do their comparisons
on the device, `corpus.cuda_calls` rises, and the program prints the same
lines. [CPU and GPU screening](CPU_GPU.md#choose-a-backend-for-packed-indexes)
covers installation of the CUDA package and how to time the two devices
fairly.

## Options

### smsd.butina

`smsd.butina(corpus, threshold=0.7, *, reordering=False, chunk_size=4096,
backend="cpu")` returns a list of tuples of row indices, each with its centre
first and its other members in row order.

| Name | Type | Default | Meaning |
| --- | --- | --- | --- |
| `corpus` | `FingerprintCorpus` | required | An open packed index; anything else raises `TypeError`, a closed one `ValueError` |
| `threshold` | float | `0.7` | Inclusive Tanimoto similarity a member must reach with its centre; finite and from 0 to 1 |
| `reordering` | bool | `False` | `True` recounts unassigned neighbours after each cluster; must be a `bool` |
| `chunk_size` | int | `4096` | Rows compared at a time; at least 1, reduced to keep each chunk within about 64 MiB |
| `backend` | str | `"cpu"` | `"cpu"`, `"cuda"` or `"auto"`; any other value raises `ValueError` |

### smsd.maxmin

`smsd.maxmin(corpus, count, *, seed=0, chunk_size=4096, backend="cpu")`
returns a list of `count` distinct row indices in the order picked.

| Name | Type | Default | Meaning |
| --- | --- | --- | --- |
| `corpus` | `FingerprintCorpus` | required | An open packed index |
| `count` | int | required | Rows to pick, from 0 to the row count; outside that range raises `ValueError` |
| `seed` | int | `0` | The first row picked; a row past the last raises `IndexError`, a negative one `ValueError` |
| `chunk_size` | int | `4096` | As for `butina` |
| `backend` | str | `"cpu"` | As for `butina` |

### FingerprintCorpus.top_k_pairs

`corpus.top_k_pairs(k, *, threshold=0.0, include_self=False, chunk_size=4096,
backend="cpu")` yields `(row, neighbour, Tanimoto)` tuples.

| Name | Type | Default | Meaning |
| --- | --- | --- | --- |
| `k` | int | required | Neighbours per row, at least 0; reduced to the number of candidate rows; 0 yields nothing |
| `threshold` | float | `0.0` | Inclusive lower bound on the reported Tanimoto similarity |
| `include_self` | bool | `False` | `True` lets a row be its own neighbour, with a score of 1 |
| `chunk_size` | int | `4096` | As for `butina` |
| `backend` | str | `"cpu"` | As for `butina`; with `"cuda"` or `"auto"`, a loaded CUDA package must also provide its top-k operation, or `RuntimeError` is raised |

Integer arguments reject booleans and floats with `TypeError`. A threshold
outside 0 to 1, or one that is not finite, raises `ValueError`.

## Notes

- A Butina cluster is defined by its centre. Two members of one cluster can
  be less similar to each other than the threshold, and a row joins the first
  cluster whose centre reaches it, not necessarily the closest centre.
- Butina centres depend on the neighbour counts, so duplicates and the choice
  of `reordering` change them. Apart from the row that `seed=0` names,
  changing the row order changes only how equal counts and equal distances
  are resolved.
- MaxMin is greedy. It gives a well-spread selection quickly, but a different
  seed gives a different, equally valid selection.
- The scores are Tanimoto similarities on binary fingerprints. Count
  fingerprints cannot be stored in a packed index.
- These tools rank fingerprints, not structures. Confirm a series or a
  near-duplicate that matters with an exact search; see
  [Maximum common substructure](MCS.md) and
  [Substructure search](SUBSTRUCTURE.md).

## See also

- [Packed fingerprint index](PACKED_INDEX.md): building, reopening and
  searching the index these tools read.
- [Screening a library](SCREENING.md): preparing a library and removing
  duplicates.
- [CPU and GPU screening](CPU_GPU.md): installing the CUDA package and
  checking where work ran.
- [Fingerprints and similarity](FINGERPRINTS.md): fingerprint families and
  `FingerprintOptions`.
- [Packed fingerprint corpora](../PACKED_CORPORA.md): the reference contract,
  including [CPU and CUDA
  backends](../PACKED_CORPORA.md#cpu-and-cuda-backends) and [ordering, work
  and memory](../PACKED_CORPORA.md#ordering-work-and-memory).
- The [MCS to a small library
  search](../../python/examples/feature_search.ipynb) notebook runs Butina
  clustering and MaxMin selection.
