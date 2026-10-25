# SMSD search algorithms

The C++, Java and Python libraries share one search model: a prepared query
owns a copy of its graph, chemistry and resource options, and each `match`,
`matchAll` or `matchMany` call runs an independent search against one target.
The Python package calls the native search. The [Python](PYTHON.md),
[C++](CPP.md) and [Java](JAVA.md) guides document the API; this page states
what each search guarantees, what a result means and how the languages differ.
The [supported-platforms table](HOWTO-INSTALL.md#supported-platforms) lists the
release packages. Matching runs on the CPU; optional CUDA screening does not
perform the MCS search.

## What each search guarantees

- **Substructure** finds embeddings of the whole query in the target.
  Non-induced matching requires every query bond in the target; induced
  matching also forbids extra target bonds between mapped atoms. Enumeration
  returns every embedding up to the requested limit, symmetric self-mappings
  included.
- **MCIS** (maximum common induced subgraph) maximises mapped atoms, then
  bonds, or the supplied query atom weights, and preserves both edges and
  non-edges between mapped atoms.
- **MCES** (maximum common edge subgraph) maximises common bonds, then the
  atoms they touch.
- Every reported mapping satisfies the chemistry options, and connectivity and
  fragment limits apply to the reported components.
- Canonical mapping selection returns one representative mapping under the
  symmetries of both molecules, or raises an error when it cannot.

## Results and budgets

- Every search reports `status`, `optimal` and `complete` separately from its
  matches. `Completed` means the search ran to its end. `TimedOut` and
  `Cancelled` mean the time limit or a cancellation stopped it. `LimitReached`
  means the node limit stopped it, an enumeration reached its match limit, or
  an MCES pair beyond the size limits below was declined before any search;
  a declined pair returns no match and proves nothing.
- Complete-ring matching and strict ring fusion depend on perceived cycles.
  If ring perception was truncated in either input, MCS retains its selected
  mappings but reports `LimitReached`, with `optimal` and `complete` both false.
  `TimedOut` and `Cancelled` take precedence. Consensus applies the same status
  rule when any input has truncated ring perception, keeping its fragment and
  supporting mappings. A valid mapping's maximality cannot be proved in this
  case, so the check raises its existing proof-unavailable exception.
- A stopped search may return a valid mapping without proving the optimum, and
  a capped enumeration may prove the objective while leaving the enumeration
  incomplete. An empty incomplete result does not prove absence.
- A negative timeout selects the default budget: 10 s for a substructure search
  and, for an MCS, 500 ms plus 2 ms per atom pair (query atoms × target
  atoms), at most 30 s. That is 1.3 s for two 20-atom molecules, 5.5 s for two
  50-atom molecules and 30 s from 122 atoms each. Zero stops immediately, and a
  node limit of zero leaves the node count unlimited.
- Batch results do not depend on the thread count unless a search reaches its
  time limit.

## Size limits

An MCES search declines a pair, with `LimitReached` and no match, when query
atoms × target atoms exceed 2³⁰, as for two molecules of more than 32,768
atoms each; when the bonds of the molecule searched from × (the bonds of the
other / 64, rounded up) exceed 2²², as for two molecules of more than 16,384
bonds each; or when 2 × the kinds of bond both molecules share × (largest atom
degree + 1) exceeds 2²⁴, which needs an atom with at least 512 bonds. Without
excluded target atoms, the search starts from the target when it has fewer
atoms than the query, or fewer bonds when more than one fragment is allowed,
and otherwise from the query. MCIS searches are never declined. The three
libraries share these limits.

## Stereochemistry

With chirality or double-bond stereochemistry enabled, a mapped centre or
double bond is compared by the parity of its mapped neighbours, not by CIP
labels, whose priorities depend on atoms outside the mapping. The comparison is
a hard constraint. MCIS returns the largest common subgraph whose mapped
stereochemistry agrees, which can be smaller than the MCIS without it: for the
alanine enantiomers `N[C@@H](C)C(=O)O` and `N[C@H](C)C(=O)O` it keeps five of
the six atoms. MCES instead leaves out the fewest common bonds that resolve
each disagreement: for `C[C@@H]1CCCO1` and `C[C@H]1CCCO1` it keeps all six
atoms and five of the six bonds, where MCIS keeps five atoms.

## Cross-language differences

These differences are deliberate; each follows the conventions of its language.

| Topic | C++ | Java | Python |
| --- | --- | --- | --- |
| Enum constants | `SearchStatus::Completed`, `MCS::Connectivity::Connected`; `ChemOptions` enums such as `BondOrderMode::STRICT` | `SearchResult.Status.Completed`, `MCS.Connectivity.Connected`; `ChemOptions` enums such as `BondOrderMode.STRICT` | `SearchStatus.COMPLETED`, `MCS.Connectivity.CONNECTED`, `BondOrderMode.STRICT` |
| Time budget and elapsed time | `timeoutMs` option, `elapsedMillis` result | `timeoutMillis` option, `elapsedMillis()` result | `timeout_ms` option, `elapsed_ms` result |
| Per-match score | `SearchResult::score` only | `Match.score()` as well as `SearchResult.score()` | `SearchResult.score` only |
| Witness accessors | Inspect `matches` directly | `found()` | `found` and `best` properties |
| Fluent option helpers | `withOptions`, `withTimeout` | `withOptions`, `withTimeout` | `with_options`, `with_timeout`, `with_objective`, `with_fragments`, `with_atom_weights`, `with_induced` |
| SMILES convenience | `MCS::find(std::string)` and `Substructure::find(std::string)` parse SMILES | Graphs and CDK containers only | `find(str)` parses SMILES; RDKit molecules are accepted |
| Batch searches | `matchMany(targets, numThreads = 0)`; zero uses every processor | `matchMany(List)` runs sequentially; `matchMany(List, numThreads)` uses up to `numThreads` threads, zero for one per processor | `match_many(targets, *, num_threads=0)` |
| Tuning switches | None | Deprecated `ChemOptions` tuning switches that never change a result | None |
| Cancellation | Shared `CancellationToken` | `BooleanSupplier`; thread interruption also cancels | Shared `CancellationToken` |
