# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
"""Molecule and SMARTS searches with explicit completion results."""
from . import _smsd
from dataclasses import dataclass
from types import MappingProxyType


@dataclass(frozen=True)
class Match:
    """Atom and selected-bond correspondences, with common-edge fragments."""
    atom_mapping: object
    bond_mapping: tuple
    fragments: tuple

    @property
    def atom_count(self):
        return len(self.atom_mapping)

    @property
    def bond_count(self):
        return len(self.bond_mapping)

    def to_smiles(self, query):
        """Write the selected query atoms and bonds as a capped molecule."""
        graph, indices = _convert(query)
        return _native_correspondence(self, indices).to_smiles(graph)

    def to_graph(self, query):
        """Extract selected edges while retaining original matching context."""
        graph, indices = _convert(query)
        return _native_correspondence(self, indices).to_graph(graph)


@dataclass(frozen=True)
class SearchResult:
    """Matches and the proof, enumeration and stopping status of a search."""
    matches: tuple
    status: object
    optimal: bool
    complete: bool
    score: float
    nodes_visited: int
    elapsed_ms: int

    @property
    def found(self):
        return bool(self.matches)

    @property
    def best(self):
        return self.matches[0] if self.matches else None


@dataclass(frozen=True)
class ConsensusResult:
    """A progressive common graph and its verified supporting embeddings."""
    fragment: object
    source_index: int
    graph: object
    mappings: object
    support_count: int
    found: bool
    status: object
    optimal: bool
    complete: bool
    nodes_visited: int
    elapsed_ms: int
    _smiles: str

    def to_smiles(self):
        return self._smiles


class _Query:
    def _replace(self, native, options):
        result = object.__new__(type(self))
        result.__dict__ = dict(self.__dict__)
        result._native = native
        result._options = options
        return result

    def with_options(self, options):
        copied = self.Options(options)
        return self._replace(self._native.with_options(self._native_options(copied)), copied)

    def _native_options(self, options):
        return options

    def _target_native(self, target_indices):
        return self._native

    def with_timeout(self, milliseconds):
        options = self.Options(self._options)
        options.timeout_ms = milliseconds
        return self.with_options(options)

    def _result(self, result, target_indices):
        query_indices = self._query_indices
        matches = tuple(_translate_match(match, query_indices, target_indices) for match in result.matches)
        return SearchResult(matches, result.status, result.optimal, result.complete,
                            result.score, result.nodes_visited, result.elapsed_ms)

    def match(self, target):
        """Return one best match, with its proof and stopping status."""
        graph, indices = _convert(target)
        return self._result(self._target_native(indices).match(graph), indices)

    def match_all(self, target, limit=1000):
        """Return optimal matches up to the limit, with enumeration status."""
        if not isinstance(limit, int) or isinstance(limit, bool) or limit < 1:
            raise ValueError("limit must be a positive integer")
        graph, indices = _convert(target)
        return self._result(self._target_native(indices).match_all(graph, limit), indices)

    def match_many(self, targets, *, num_threads=0):
        """Search targets independently and preserve each target's atom indices."""
        if not isinstance(num_threads, int) or isinstance(num_threads, bool) or num_threads < 0:
            raise ValueError("num_threads must be a nonnegative integer")
        from . import TargetCorpus
        converted = targets._snapshot() if isinstance(targets, TargetCorpus) else [_convert(target) for target in targets]
        if isinstance(self, MCS) and self._options.excluded_target_atoms:
            groups = {}
            for index, (graph, indices) in enumerate(converted):
                native = self._target_native(indices)
                key = tuple(indices) if indices is not None else None
                groups.setdefault(key, (native, []))[1].append((index, graph))
            results = [None] * len(converted)
            for native, members in groups.values():
                batch = native.match_many([graph for _, graph in members], num_threads)
                for result, (index, _) in zip(batch, members):
                    results[index] = self._result(result, converted[index][1])
            return results
        results = self._native.match_many([graph for graph, _ in converted], num_threads)
        return [self._result(result, indices)
                for result, (_, indices) in zip(results, converted)]


def _convert(molecule):
    from . import _ensure_mol_ex, get_index_map
    graph, source = _ensure_mol_ex(molecule)
    indices = tuple(get_index_map(graph)) if source is not None else None
    return graph, indices


def _translate_match(match, query_indices=None, target_indices=None):
    def query_index(index):
        return query_indices[index] if query_indices is not None else index
    def target_index(index):
        return target_indices[index] if target_indices is not None else index
    atom_mapping = MappingProxyType({query_index(first): target_index(second)
                                     for first, second in match.atom_mapping.items()})
    bonds = tuple((query_index(first), query_index(second), target_index(target_first), target_index(target_second))
                  for first, second, target_first, target_second in match.bond_mapping)
    fragments = tuple(tuple(query_index(index) for index in fragment) for fragment in match.fragments)
    return Match(atom_mapping, bonds, fragments)


def _native_correspondence(match, query_indices=None, target_indices=None):
    query_reverse = {original: native for native, original in enumerate(query_indices)} if query_indices is not None else None
    target_reverse = {original: native for native, original in enumerate(target_indices)} if target_indices is not None else None
    def query_index(index):
        if query_reverse is not None and index not in query_reverse:
            raise IndexError("Query atom index is outside the input molecule")
        return query_reverse[index] if query_reverse is not None else index
    def target_index(index):
        if target_reverse is not None and index not in target_reverse:
            raise IndexError("Target atom index is outside the input molecule")
        return target_reverse[index] if target_reverse is not None else index
    if not isinstance(match, Match):
        return {query_index(first): target_index(second) for first, second in match.items()}
    native = _smsd.Match()
    native.atom_mapping = {query_index(first): target_index(second) for first, second in match.atom_mapping.items()}
    native.bond_mapping = [(query_index(first), query_index(second), target_index(target_first), target_index(target_second))
                           for first, second, target_first, target_second in match.bond_mapping]
    native.fragments = [[query_index(index) for index in fragment] for fragment in match.fragments]
    return native


class MCS(_Query):
    """Find common induced or edge subgraphs, including disconnected fragments."""
    Options = _smsd.MCS.Options
    Objective = _smsd.MCS.Objective
    Connectivity = _smsd.MCS.Connectivity

    @classmethod
    def find(cls, query):
        """Prepare a molecule query for repeated and concurrent searches."""
        result = object.__new__(cls)
        result._query_graph, result._query_indices = _convert(query)
        result._source_atom_count = query.GetNumAtoms() if result._query_indices is not None else result._query_graph.n
        result._options = cls.Options()
        result._native = _smsd.MCS.find(result._query_graph)
        return result

    def with_objective(self, objective):
        options = self.Options(self._options)
        options.objective = objective
        return self.with_options(options)

    def with_fragments(self, *, min_atoms=1, max_fragments=2147483647):
        options = self.Options(self._options)
        options.connectivity = self.Connectivity.FRAGMENTS
        options.min_fragment_atoms = min_atoms
        options.max_fragments = max_fragments
        return self.with_options(options)

    def with_atom_weights(self, weights):
        options = self.Options(self._options)
        options.atom_weights = list(weights)
        return self.with_options(options)

    def _native_options(self, options):
        native = self.Options(options)
        if self._query_indices is not None and native.atom_weights:
            if len(native.atom_weights) != self._source_atom_count:
                raise ValueError("Weights must cover every source query atom")
            weights = native.atom_weights
            native.atom_weights = [weights[index] for index in self._query_indices]
        return native

    def _target_native(self, target_indices):
        if target_indices is None or not self._options.excluded_target_atoms:
            return self._native
        reverse = {original: native for native, original in enumerate(target_indices)}
        excluded = self._options.excluded_target_atoms
        if any(index not in reverse for index in excluded):
            raise IndexError("Excluded target atom index is outside the input molecule")
        options = self._native_options(self._options)
        options.excluded_target_atoms = [reverse[index] for index in excluded]
        return self._native.with_options(options)

    @staticmethod
    def validate_mapping(query, target, match, options=None):
        """Return selected-graph validation errors; an empty list is valid."""
        query_graph, query_indices = _convert(query)
        target_graph, target_indices = _convert(target)
        native = _native_correspondence(match, query_indices, target_indices)
        return _smsd.MCS.validate_mapping(query_graph, target_graph, native, options or MCS.Options())

    @staticmethod
    def is_mapping_maximal(query, target, match, options=None):
        """Check whether a valid common graph admits a one-step extension."""
        query_graph, query_indices = _convert(query)
        target_graph, target_indices = _convert(target)
        native = _native_correspondence(match, query_indices, target_indices)
        return _smsd.MCS.is_mapping_maximal(query_graph, target_graph, native, options or MCS.Options())

    @staticmethod
    def canonicalize_mapping(query, target, atom_mapping):
        """Choose a representative under the two molecules' symmetries."""
        query_graph, query_indices = _convert(query)
        target_graph, target_indices = _convert(target)
        native = _native_correspondence(atom_mapping, query_indices, target_indices)
        result = _smsd.MCS.canonicalize_mapping(query_graph, target_graph, native)
        return {(query_indices[first] if query_indices is not None else first):
                (target_indices[second] if target_indices is not None else second)
                for first, second in result.items()}

    @staticmethod
    def compute_tautomer_confidence(query, target, match, chemistry):
        """Score selected chemistry relaxations using a read-time pH/solvent heuristic."""
        if not isinstance(match, Match):
            raise TypeError("match must include selected bonds")
        query_graph, query_indices = _convert(query)
        target_graph, target_indices = _convert(target)
        native = _native_correspondence(match, query_indices, target_indices)
        return _smsd.MCS.compute_tautomer_confidence(query_graph, target_graph, native, chemistry)

    @staticmethod
    def match_constrained(queries, target, options=None):
        """Match queries in order without reusing original target atoms."""
        target_graph, target_indices = _convert(target)
        configured = MCS.Options(options) if options is not None else MCS.Options()
        excluded = set(configured.excluded_target_atoms)
        reverse = {original: native for native, original in enumerate(target_indices)} if target_indices is not None else None
        results = []
        for query in queries:
            prepared = MCS.find(query)
            current = MCS.Options(configured)
            if reverse is not None:
                if any(index not in reverse for index in excluded):
                    raise IndexError("Excluded target atom index is outside the input molecule")
                current.excluded_target_atoms = [reverse[index] for index in excluded]
            else:
                current.excluded_target_atoms = sorted(excluded)
            prepared = prepared.with_options(current)
            result = prepared._result(prepared._native.match(target_graph), target_indices)
            if result.found:
                excluded.update(result.best.atom_mapping.values())
            results.append(result)
        return results

    @staticmethod
    def consensus(molecules, options=None, threshold=1.0):
        """Find a progressive common graph supported by a fraction of inputs."""
        converted = [_convert(molecule) for molecule in molecules]
        native = _smsd.MCS.consensus([graph for graph, _ in converted], options or MCS.Options(), threshold)
        source_indices = converted[native.source_index][1] if native.found else None
        fragment = _translate_match(native.fragment, source_indices, source_indices) if native.found else None
        mappings = MappingProxyType({index: _translate_match(match, None, converted[index][1])
                                     for index, match in native.mappings.items()})
        return ConsensusResult(fragment, native.source_index, native.graph, mappings, native.support_count,
                               native.found, native.status, native.optimal, native.complete,
                               native.nodes_visited, native.elapsed_ms, native.to_smiles())


class Substructure(_Query):
    """Find molecule or SMARTS embeddings using one public query interface."""
    Options = _smsd.Substructure.Options

    @classmethod
    def find(cls, query):
        """Prepare a molecule query; SMILES strings are parsed as molecules."""
        result = object.__new__(cls)
        result._query_graph, result._query_indices = _convert(query)
        result._options = cls.Options()
        result._native = _smsd.Substructure.find(result._query_graph)
        return result

    @classmethod
    def smarts(cls, query):
        """Compile a SMARTS predicate query once for repeated searches."""
        if not isinstance(query, str) or not query:
            raise ValueError("SMARTS must be a non-empty string")
        result = object.__new__(cls)
        result._query_graph = result._query_indices = None
        result._options = cls.Options()
        result._native = _smsd.Substructure.smarts(query)
        return result

    def with_induced(self, enabled=True):
        options = self.Options(self._options)
        options.induced = enabled
        return self.with_options(options)

    @staticmethod
    def decompose_r_groups(core, targets, options=None):
        """Return core and substituent graphs alongside each search result."""
        core_graph, core_indices = _convert(core)
        converted = [_convert(target) for target in targets]
        native = _smsd.Substructure.decompose_r_groups(core_graph, [graph for graph, _ in converted],
                                                      options or Substructure.Options())
        query = Substructure.find(core)
        return [dict(core=result["core"], rgroups=result["rgroups"],
                     search=query._result(result["search"], indices))
                for result, (_, indices) in zip(native, converted)]
