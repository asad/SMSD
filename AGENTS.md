# Repository conventions

Use `MCS` for the maximum common substructure acronym in prose, class names,
and acronyms embedded in camelCase identifiers. Examples include
`MCSAlgorithmRegressionTest`, `findMCSSmiles`, and `nearMCSDelta`.
Conventional lowercase snake_case names such as `mcs.hpp`, `mcs_engine.py`,
and `find_mcs` retain their language-specific spelling.

When renaming an exposed identifier, update its callers, bindings, examples,
and documentation together. State API name changes in the changelog.
