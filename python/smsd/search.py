# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Molecule and SMARTS searches with explicit completion results."""
from dataclasses import dataclass
from types import MappingProxyType

from . import _accelerator, _smsd
from ._interop import _index_tuple
from ._molecule import _checked_int, _ensure_mol_ex
from ._screening import TargetCorpus

__all__ = ["CorpusSearchResult", "ConsensusResult", "MCS", "Match", "SearchResult", "SmartsResult", "Substructure"]


@dataclass(frozen=True)
class Match:
    """Atom and selected-bond correspondences, with common-edge fragments.

    Example::

        >>> import smsd
        >>> match = smsd.MCS.find("CCO").match("CCN").best
        >>> assert isinstance(match, smsd.Match) and match.atom_count == 2
    """
    atom_mapping: object
    bond_mapping: tuple
    fragments: tuple

    @property
    def atom_count(self):
        """The number of matched atoms.

        Example::

            >>> import smsd
            >>> assert smsd.MCS.find("CCO").match("CCN").best.atom_count == 2
        """
        return len(self.atom_mapping)

    @property
    def bond_count(self):
        """The number of matched bonds.

        Example::

            >>> import smsd
            >>> assert smsd.MCS.find("CCO").match("CCN").best.bond_count == 1
        """
        return len(self.bond_mapping)

    def to_smiles(self, query):
        """Write the selected query atoms and bonds as a capped molecule.

        Example::

            >>> import smsd
            >>> result = smsd.MCS.find("CCO").match("CCN")
            >>> assert smsd.parse_smiles(result.best.to_smiles("CCO")).n == 2
        """
        graph, indices = _convert(query)
        return _native_correspondence(self, indices).to_smiles(graph)

    def to_graph(self, query):
        """Extract selected edges while retaining original matching context.

        Example::

            >>> import smsd
            >>> result = smsd.MCS.find("CCO").match("CCN")
            >>> assert result.best.to_graph("CCO").n == 2
        """
        graph, indices = _convert(query)
        return _native_correspondence(self, indices).to_graph(graph)

    def to_smarts(self, query, chemistry=None):
        """Export selected bonds as SMARTS; reject rules SMARTS cannot represent.

        Example::

            >>> import smsd
            >>> result = smsd.MCS.find("CCO").match("CCN")
            >>> pattern = result.best.to_smarts("CCO")
            >>> assert smsd.Substructure.smarts(pattern).match("CCN").found
        """
        graph, indices = _convert(query)
        return _native_correspondence(self, indices).to_smarts(
            graph, chemistry if chemistry is not None else _smsd.ChemOptions())


@dataclass(frozen=True)
class SearchResult:
    """Matches and the proof, enumeration and stopping status of a search.

    Example::

        >>> import smsd
        >>> result = smsd.MCS.find("CCO").match("CCN")
        >>> assert result.found and result.optimal and result.complete
    """
    matches: tuple
    status: object
    optimal: bool
    complete: bool
    score: float
    nodes_visited: int
    elapsed_ms: int

    @property
    def found(self):
        """Whether the search found at least one match.

        Example::

            >>> import smsd
            >>> assert smsd.Substructure.smarts("[#8]").match("CCO").found
        """
        return bool(self.matches)

    @property
    def best(self):
        """The first match, or None when nothing was found.

        Example::

            >>> import smsd
            >>> result = smsd.Substructure.smarts("[#8]").match("CCN")
            >>> assert result.best is None
        """
        return self.matches[0] if self.matches else None


@dataclass(frozen=True)
class CorpusSearchResult:
    """One original corpus slot, with a threshold decision and any actual CPU result.

    ``result`` is None only when a conservative filter proves the threshold
    cannot be reached. Both decision flags being false means unknown, including
    a stopped search with no sufficient witness. ``atom_upper_bound`` is None
    when no element-count bound was used; substructure uses a path subset proof.

    Example::

        >>> import smsd
        >>> corpus = smsd.TargetCorpus.from_smiles(["CC", "N"])
        >>> decisions = smsd.MCS.find("CC").match_corpus(corpus, 2, num_threads=1)
        >>> assert decisions[0].threshold_reached and decisions[1].proven_below_threshold
        >>> assert decisions[1].result is None
    """
    target_index: int
    atom_upper_bound: int | None
    proven_below_threshold: bool
    threshold_reached: bool
    result: SearchResult | None


# Prepared searches build these immutable instances natively, in one call.
_smsd._register_search_types(Match, SearchResult)
_MolGraph = _smsd.MolGraph


@dataclass(frozen=True)
class SmartsResult:
    """A search result with a portable SMARTS query for its best selected graph.

    Example::

        >>> import smsd
        >>> exported = smsd.MCS.find("CCO").match_with_smarts("CCN")
        >>> assert exported.result.found and exported.smarts
    """
    result: SearchResult
    smarts: str


@dataclass(frozen=True)
class ConsensusResult:
    """A progressive common graph and its verified supporting embeddings.

    ``optimal`` is always False, including empty and single-input calls.
    ``complete`` reports whether the progressive pass finished; it does not
    prove a global optimum or the absence of a qualifying common graph.

    Example::

        >>> import smsd
        >>> result = smsd.MCS.consensus(["CCO", "CCN"])
        >>> assert result.found and result.support_count == 2
    """
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
        """The common graph as SMILES.

        Example::

            >>> import smsd
            >>> result = smsd.MCS.consensus(["CCO", "CCN"])
            >>> assert smsd.parse_smiles(result.to_smiles()).n == 2
        """
        return self._smiles


class _Query:
    @classmethod
    def _prepared(cls, graph, indices, native):
        """A query over ``graph`` with default options; ``indices`` number its atoms, None keeps native indices."""
        query = object.__new__(cls)
        query._query_graph, query._query_indices = graph, indices
        return query._prepare(native, cls.Options())

    def _prepare(self, native, options):
        """Hold the native query, its options and its numbered form; return self."""
        self._native = native
        self._options = options
        # Searches of graphs that need no renumbering take only their targets.
        self._search = native._with_numbering(self._query_indices)
        return self

    def _replace(self, native, options):
        result = object.__new__(type(self))
        result.__dict__ = dict(self.__dict__)
        return result._prepare(native, options)

    def with_options(self, options):
        """This query with a copy of ``options``.

        Example::

            >>> import smsd
            >>> options = smsd.MCS.Options()
            >>> options.node_limit = 1000
            >>> result = smsd.MCS.find("CCO").with_options(options).match("CCN")
            >>> assert result.found
        """
        copied = self.Options(options)
        return self._replace(self._native.with_options(self._native_options(copied)), copied)

    def _native_options(self, options):
        return options

    def _target_native(self, target_indices):
        return self._native

    def with_timeout(self, milliseconds):
        """This query with a time limit per target, in milliseconds; a negative value selects the default limit and zero stops at once.

        Example::

            >>> import smsd
            >>> result = smsd.MCS.find("CCO").with_timeout(0).match("CCN")
            >>> assert result.status == smsd.SearchStatus.TIMED_OUT and not result.optimal
        """
        options = self.Options(self._options)
        options.timeout_ms = milliseconds
        return self.with_options(options)

    def match(self, target):
        """Return one best match, with its proof and stopping status.

        Example::

            >>> import smsd
            >>> result = smsd.MCS.find("CCO").match("CCN")
            >>> assert result.found and result.best.atom_count == 2
        """
        if type(target) is _MolGraph:
            return self._search._match(target)
        graph, indices = _convert(target)
        if indices is None:
            return self._search._match(graph)
        return self._target_native(indices)._match(graph, self._query_indices, indices)

    def match_all(self, target, limit=1000):
        """Return optimal matches up to the limit, with enumeration status.

        Example::

            >>> import smsd
            >>> result = smsd.Substructure.smarts("[#8]").match_all("OCCO", limit=4)
            >>> assert len(result.matches) == 2 and result.complete
        """
        _checked_int(limit, "Match limit", 1)
        if type(target) is _MolGraph:
            return self._search._match_all(target, limit)
        graph, indices = _convert(target)
        if indices is None:
            return self._search._match_all(graph, limit)
        return self._target_native(indices)._match_all(graph, limit, self._query_indices, indices)

    def match_many(self, targets, *, num_threads=0):
        """Search targets independently and preserve each target's atom indices.

        Example::

            >>> import smsd
            >>> results = smsd.Substructure.smarts("[#8]").match_many(["CCO", "CCN"], num_threads=1)
            >>> assert [result.found for result in results] == [True, False]
        """
        _checked_int(num_threads, "Thread count", 0)
        if isinstance(self, Substructure) and isinstance(targets, TargetCorpus):
            _accelerator.discover()
            return self._native._match_corpus(targets._impl, targets._snapshot(), num_threads,
                                            self._query_indices)
        if isinstance(targets, TargetCorpus):
            converted = targets._snapshot()
        else:
            converted = [(target, None) if type(target) is _MolGraph else _convert(target) for target in targets]
        if isinstance(self, MCS) and self._options.excluded_target_atoms:
            groups = {}
            for index, (graph, indices) in enumerate(converted):
                native = self._target_native(indices)
                key = tuple(indices) if indices is not None else None
                groups.setdefault(key, (native, []))[1].append(index)
            results = [None] * len(converted)
            for native, members in groups.values():
                batch = native._match_many([converted[index] for index in members], num_threads,
                                           self._query_indices)
                for result, index in zip(batch, members):
                    results[index] = result
            return results
        return self._search._match_many(converted, num_threads)


def _check_corpus_arguments(corpus, gpu_prefilter, num_threads):
    if not isinstance(corpus, TargetCorpus):
        raise TypeError("Corpus must be a TargetCorpus")
    if not isinstance(gpu_prefilter, bool):
        raise TypeError("GPU prefilter flag must be a bool")
    _checked_int(num_threads, "Thread count", 0)


def _convert(molecule):
    """A molecule's graph and the caller's numbering of its atoms, or None when native indices apply."""
    if type(molecule) is _MolGraph:
        return molecule, None
    graph, source = _ensure_mol_ex(molecule)
    return graph, _index_tuple(graph) if source is not None else None


def _native_numbering(indices, role):
    """A function from the caller's atom numbers to native indices; the identity when ``indices`` is None."""
    if indices is None:
        return lambda index: index
    reverse = {original: native for native, original in enumerate(indices)}

    def native_index(index):
        if index not in reverse:
            raise IndexError(f"{role} atom index out of range")
        return reverse[index]
    return native_index


def _native_correspondence(match, query_indices=None, target_indices=None):
    """``match`` in native indices: a dict for an atom mapping, a native Match for a :class:`Match`."""
    query_index = _native_numbering(query_indices, "Query")
    target_index = _native_numbering(target_indices, "Target")
    if not isinstance(match, Match):
        return {query_index(first): target_index(second) for first, second in match.items()}
    native = _smsd.Match()
    native.atom_mapping = {query_index(first): target_index(second) for first, second in match.atom_mapping.items()}
    native.bond_mapping = [(query_index(first), query_index(second), target_index(target_first), target_index(target_second))
                           for first, second, target_first, target_second in match.bond_mapping]
    native.fragments = [[query_index(index) for index in fragment] for fragment in match.fragments]
    return native


class MCS(_Query):
    """Find common induced or edge subgraphs, including disconnected fragments.

    Example::

        >>> import smsd
        >>> result = smsd.MCS.find("CCO").match("CCN")
        >>> assert result.best.atom_count == 2
    """
    Options = _smsd.MCS.Options
    Objective = _smsd.MCS.Objective
    Connectivity = _smsd.MCS.Connectivity

    def match_corpus(self, corpus, minimum_atoms, *, gpu_prefilter=False, num_threads=0):
        """Return atom-count threshold decisions for original corpus slots.

        The conservative filter supports unweighted MCIS with element matching.
        Equality with ``minimum_atoms`` is searched. MCES, weights, excluded
        atoms, node caps, immediate timeouts and cancellation tokens retain the
        ordinary CPU search. An unfinished below-threshold search is unknown.
        ``gpu_prefilter=False`` uses the CPU version of the same filter; True
        attempts the optional GPU provider and falls back to it on failure.
        Ordinary :meth:`match_many` remains an exact-search-only result API.
        GPU use is opt-in and does not require :meth:`TargetCorpus.prewarm`.

        Example::

            >>> import smsd
            >>> corpus = smsd.TargetCorpus.from_smiles(["CCO", "N"])
            >>> decisions = smsd.MCS.find("CC").match_corpus(corpus, 2, num_threads=1)
            >>> assert decisions[0].threshold_reached
            >>> assert decisions[1].proven_below_threshold and decisions[1].atom_upper_bound == 0
        """
        _check_corpus_arguments(corpus, gpu_prefilter, num_threads)
        _checked_int(minimum_atoms, "Minimum atom count", 1)
        if gpu_prefilter:
            _accelerator.discover()
        with corpus._lock:
            # Exclusions are expressed in the caller's target atom numbers.
            if self._options.excluded_target_atoms:
                results = self.match_many(corpus, num_threads=num_threads)
                atom_objective = self._options.objective == self.Objective.MCIS and not self._options.atom_weights
                decisions = []
                for index, result in enumerate(results):
                    reached = any(match.atom_count >= minimum_atoms for match in result.matches)
                    below = not reached and atom_objective and result.optimal and result.status == _smsd.SearchStatus.COMPLETED
                    decisions.append(CorpusSearchResult(index, None, below, reached, result))
                return decisions
            return [CorpusSearchResult(*fields) for fields in self._native._match_corpus_prefiltered(
                corpus._impl, corpus._snapshot(), minimum_atoms, gpu_prefilter, num_threads, self._query_indices)]

    def match_with_smarts(self, target):
        """Return the match and its SMARTS under this query's chemistry.

        Unsupported stereo, tautomer or ring-export constraints raise rather
        than being dropped. The contained result retains all proof/status flags.

        Example::

            >>> import smsd
            >>> exported = smsd.MCS.find("CCO").match_with_smarts("CCN")
            >>> assert exported.result.optimal
            >>> assert smsd.Substructure.smarts(exported.smarts).match("CCN").found
        """
        result = self.match(target)
        pattern = "" if result.best is None else _native_correspondence(
            result.best, self._query_indices).to_smarts(self._query_graph, self._options.chemistry)
        return SmartsResult(result, pattern)

    @classmethod
    def find(cls, query):
        """Prepare a molecule query for repeated and concurrent searches.

        Example::

            >>> import smsd
            >>> prepared = smsd.MCS.find("CCO")
            >>> assert prepared.match("CCN").best.atom_count == 2
        """
        graph, indices = _convert(query)
        prepared = cls._prepared(graph, indices, _smsd.MCS.find(graph))
        prepared._source_atom_count = query.GetNumAtoms() if indices is not None else graph.n
        return prepared

    def with_objective(self, objective):
        """This query maximising ``objective``: common atoms (MCIS) or common bonds (MCES).

        Example::

            >>> import smsd
            >>> result = smsd.MCS.find("CC").with_objective(smsd.MCS.Objective.MCES).match("CCO")
            >>> assert result.best.bond_count == 1 and result.score == 1
        """
        options = self.Options(self._options)
        options.objective = objective
        return self.with_options(options)

    def with_fragments(self, *, min_atoms=1, max_fragments=2147483647):
        """This query allowing disconnected common fragments of at least ``min_atoms`` atoms each.

        Example::

            >>> import smsd
            >>> result = smsd.MCS.find("C.O").with_fragments(max_fragments=2).match("C.O")
            >>> assert result.best.atom_count == 2 and len(result.best.fragments) == 2
        """
        options = self.Options(self._options)
        options.connectivity = self.Connectivity.FRAGMENTS
        options.min_fragment_atoms = min_atoms
        options.max_fragments = max_fragments
        return self.with_options(options)

    def with_atom_weights(self, weights):
        """This query scoring matched query atoms by ``weights``, one per input atom.

        Example::

            >>> import smsd
            >>> result = smsd.MCS.find("CCO").with_atom_weights([0.5, 1.0, 3.0]).match("CCN")
            >>> assert result.found and result.score == 1.5
        """
        options = self.Options(self._options)
        options.atom_weights = list(weights)
        return self.with_options(options)

    def _native_options(self, options):
        native = self.Options(options)
        if self._query_indices is not None and native.atom_weights:
            if len(native.atom_weights) != self._source_atom_count:
                raise ValueError("Atom weights must cover every query atom")
            weights = native.atom_weights
            native.atom_weights = [weights[index] for index in self._query_indices]
        return native

    def _excluded(self, target_indices):
        """The excluded target atoms of the options in native indices."""
        native_index = _native_numbering(target_indices, "Excluded target")
        return [native_index(index) for index in self._options.excluded_target_atoms]

    def _target_native(self, target_indices):
        if target_indices is None or not self._options.excluded_target_atoms:
            return self._native
        options = self._native_options(self._options)
        options.excluded_target_atoms = self._excluded(target_indices)
        return self._native.with_options(options)

    @staticmethod
    def validate_mapping(query, target, match, options=None):
        """Return selected-graph validation errors; an empty list is valid.

        Example::

            >>> import smsd
            >>> match = smsd.MCS.find("CCO").match("CCN").best
            >>> assert smsd.MCS.validate_mapping("CCO", "CCN", match) == []
        """
        query_graph, query_indices = _convert(query)
        target_graph, target_indices = _convert(target)
        native = _native_correspondence(match, query_indices, target_indices)
        return _smsd.MCS.validate_mapping(query_graph, target_graph, native, options or MCS.Options())

    @staticmethod
    def is_mapping_maximal(query, target, match, options=None):
        """Check whether a valid common graph admits a compatible extension.

        Example::

            >>> import smsd
            >>> match = smsd.MCS.find("CCO").match("CCN").best
            >>> assert smsd.MCS.is_mapping_maximal("CCO", "CCN", match)
        """
        query_graph, query_indices = _convert(query)
        target_graph, target_indices = _convert(target)
        native = _native_correspondence(match, query_indices, target_indices)
        return _smsd.MCS.is_mapping_maximal(query_graph, target_graph, native, options or MCS.Options())

    @staticmethod
    def validate_tautomer_consistency(query, target, match):
        """Whether a mapping conserves mobile protons across tautomeric groups.

        Example::

            >>> import smsd
            >>> match = smsd.MCS.find("CCO").match("CCO").best
            >>> assert smsd.MCS.validate_tautomer_consistency("CCO", "CCO", match)
        """
        query_graph, query_indices = _convert(query)
        target_graph, target_indices = _convert(target)
        mapping = match.atom_mapping if isinstance(match, Match) else match
        native = _native_correspondence(mapping, query_indices, target_indices)
        return _smsd.validate_tautomer_consistency(query_graph, target_graph, native)

    @staticmethod
    def canonicalize_mapping(query, target, atom_mapping):
        """Choose a representative under the two molecules' symmetries.

        Example::

            >>> import smsd
            >>> mapping = smsd.MCS.canonicalize_mapping("CC", "CC", {0: 1, 1: 0})
            >>> assert mapping == {0: 0, 1: 1}
        """
        query_graph, query_indices = _convert(query)
        target_graph, target_indices = _convert(target)
        native = _native_correspondence(atom_mapping, query_indices, target_indices)
        result = _smsd.MCS.canonicalize_mapping(query_graph, target_graph, native)
        return {(query_indices[first] if query_indices is not None else first):
                (target_indices[second] if target_indices is not None else second)
                for first, second in result.items()}

    @staticmethod
    def compute_tautomer_confidence(query, target, match, chemistry):
        """Score selected chemistry relaxations using a read-time pH/solvent heuristic.

        Example::

            >>> import smsd
            >>> match = smsd.MCS.find("CCO").match("CCO").best
            >>> confidence = smsd.MCS.compute_tautomer_confidence("CCO", "CCO", match, smsd.ChemOptions())
            >>> assert 0 <= confidence <= 1
        """
        if not isinstance(match, Match):
            raise TypeError("match must include selected bonds")
        query_graph, query_indices = _convert(query)
        target_graph, target_indices = _convert(target)
        native = _native_correspondence(match, query_indices, target_indices)
        return _smsd.MCS.compute_tautomer_confidence(query_graph, target_graph, native, chemistry)

    @staticmethod
    def match_constrained(queries, target, options=None):
        """Match queries in order without reusing original target atoms.

        Example::

            >>> import smsd
            >>> results = smsd.MCS.match_constrained(["CC", "O"], "CCO")
            >>> assert all(result.found for result in results)
            >>> assert set(results[0].best.atom_mapping.values()).isdisjoint(results[1].best.atom_mapping.values())
        """
        target_graph, target_indices = _convert(target)
        configured = MCS.Options(options) if options is not None else MCS.Options()
        excluded = set(configured.excluded_target_atoms)
        native_index = _native_numbering(target_indices, "Excluded target")
        results = []
        for query in queries:
            current = MCS.Options(configured)
            if target_indices is not None:
                current.excluded_target_atoms = [native_index(index) for index in excluded]
            else:
                current.excluded_target_atoms = sorted(excluded)
            prepared = MCS.find(query).with_options(current)
            result = prepared._native._match(target_graph, prepared._query_indices, target_indices)
            if result.found:
                excluded.update(result.best.atom_mapping.values())
            results.append(result)
        return results

    @staticmethod
    def consensus(molecules, options=None, threshold=1.0):
        """Find a progressive common graph supported by a fraction of inputs.

        This is an approximation under one shared budget. ``optimal`` is
        always False, including empty and single-input calls. ``complete``
        and ``COMPLETED`` report that the progressive pass finished, not that
        a global optimum or the absence of a qualifying common graph was proved.

        Example::

            >>> import smsd
            >>> result = smsd.MCS.consensus(["CCO", "CCN"])
            >>> assert result.found and result.support_count == 2
        """
        converted = [_convert(molecule) for molecule in molecules]
        native = _smsd.MCS.consensus([graph for graph, _ in converted], options or MCS.Options(), threshold)
        source_indices = converted[native.source_index][1] if native.found else None
        fragment = _smsd._python_match(native.fragment, source_indices, source_indices) if native.found else None
        mappings = MappingProxyType({index: _smsd._python_match(match, None, converted[index][1])
                                     for index, match in native.mappings.items()})
        return ConsensusResult(fragment, native.source_index, native.graph, mappings, native.support_count,
                               native.found, native.status, native.optimal, native.complete,
                               native.nodes_visited, native.elapsed_ms, native.to_smiles())


class Substructure(_Query):
    """Find molecule or SMARTS embeddings using one public query interface.

    Example::

        >>> import smsd
        >>> assert smsd.Substructure.smarts("[OX2H]").match("CCO").found
    """
    Options = _smsd.Substructure.Options

    def match_corpus(self, corpus, *, gpu_prefilter=False, num_threads=0):
        """Search corpus slots after a conservative element-path subset filter.

        A rejected slot has an explicit negative proof and no invented search
        result. Survivors retain ordinary CPU mappings and proof/status fields.
        SMARTS, untyped, empty, pseudo and hydrogen-only queries bypass the
        filter, as do capped, immediately stopped or cancellable searches.
        ``gpu_prefilter=False`` uses the CPU version of the same filter; True
        attempts optional GPU screening and falls back to the CPU filter.
        ECFP similarity fingerprints are never used for subset rejection.
        No :meth:`TargetCorpus.prewarm` is required; ordinary :meth:`match_many`
        retains its existing behaviour.

        Example::

            >>> import smsd
            >>> corpus = smsd.TargetCorpus.from_smiles(["CCO", "N"])
            >>> decisions = smsd.Substructure.find("CO").match_corpus(corpus, num_threads=1)
            >>> assert decisions[0].threshold_reached and decisions[0].result.found
            >>> assert decisions[1].proven_below_threshold and decisions[1].result is None
        """
        _check_corpus_arguments(corpus, gpu_prefilter, num_threads)
        if gpu_prefilter:
            _accelerator.discover()
        with corpus._lock:
            return [CorpusSearchResult(*fields) for fields in self._native._match_corpus_prefiltered(
                corpus._impl, corpus._snapshot(), gpu_prefilter, num_threads, self._query_indices)]

    @classmethod
    def find(cls, query):
        """Prepare a molecule query; SMILES strings are parsed as molecules.

        Example::

            >>> import smsd
            >>> assert smsd.Substructure.find("CO").match("CCO").found
        """
        graph, indices = _convert(query)
        return cls._prepared(graph, indices, _smsd.Substructure.find(graph))

    @classmethod
    def smarts(cls, query):
        """Compile a SMARTS predicate query once for repeated searches.

        Example::

            >>> import smsd
            >>> query = smsd.Substructure.smarts("[OX2H]")
            >>> assert query.match("CCO").found and not query.match("CCN").found
        """
        if not isinstance(query, str) or not query:
            raise ValueError("SMARTS must be a non-empty string")
        return cls._prepared(None, None, _smsd.Substructure.smarts(query))

    def with_induced(self, enabled=True):
        """This query requiring an induced embedding, so query non-bonds are kept apart too.

        Example::

            >>> import smsd
            >>> query = smsd.Substructure.find("CCC").with_induced()
            >>> assert not query.match("C1CC1").found
        """
        options = self.Options(self._options)
        options.induced = enabled
        return self.with_options(options)

    @staticmethod
    def decompose_r_groups(core, targets, options=None):
        """Return core and substituent graphs alongside each search result.

        Example::

            >>> import smsd
            >>> results = smsd.Substructure.decompose_r_groups("CO", ["CCO"])
            >>> assert len(results) == 1 and results[0]["search"].found
        """
        core_graph, core_indices = _convert(core)
        converted = [_convert(target) for target in targets]
        options = options or Substructure.Options()
        native = _smsd.Substructure.decompose_r_groups(core_graph, [graph for graph, _ in converted], options)
        return [dict(core=result["core"], rgroups=result["rgroups"],
                     search=_smsd._search_result(result["search"], core_indices, indices, options.cancellation))
                for result, (_, indices) in zip(native, converted)]
