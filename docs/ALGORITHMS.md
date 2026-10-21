# SMSD search algorithms

The C++, Java and Python libraries share one search model: a prepared query
owns a copy of its graph, chemistry and resource options, and each `match`,
`matchAll` or `matchMany` call runs an independent search against one target.
The Python package calls the native search. The [Python](PYTHON.md),
[C++](CPP.md) and [Java](JAVA.md) guides document the API; this page states
what each search guarantees, what a result means and how the languages differ.

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
  matches. `Completed` means the search ran to its end; `TimedOut`, `Cancelled`
  and `LimitReached` say why it stopped early.
- A stopped search may return a valid mapping without proving the optimum, and
  a capped enumeration may prove the objective while leaving the enumeration
  incomplete. An empty incomplete result does not prove absence.
- A negative timeout selects a default budget from the graph sizes, zero stops
  immediately, and a node limit of zero leaves the node count unlimited.
- Batch results do not depend on the thread count unless a search reaches its
  time limit.

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
