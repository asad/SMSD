# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Similarity bounds, ranked similarity and reusable target corpora."""

import threading
import warnings

from . import _accelerator, _smsd
from ._interop import _index_tuple
from ._molecule import _checked_int, _ensure_mol, _ensure_mol_ex

__all__ = [
    "TargetCorpus", "batch_similarity", "gpu_device_info", "gpu_is_available", "screen_targets",
    "similarity", "similarity_upper_bound",
]


def gpu_is_available():
    """True when a GPU is available at run time: through an installed smsd-cuda package that
    finds a device, or through a native module built with CUDA or Metal."""
    return _accelerator.active() is not None or _smsd.gpu_is_available()


def gpu_device_info():
    """Human-readable compute backend string.

    Examples::

        'GPU: Tesla T4 (compute 7.5, 14930 MB)  [OpenMP 4.5, 8 threads]'
        'Metal GPU: Apple M2 Pro  [OpenMP 5.0, 10 threads]'
        'CPU: OpenMP 4.5, 16 threads'
    """
    native = _smsd.gpu_device_info()
    accelerator = _accelerator.active()
    if accelerator is None or not native.startswith("CPU: "):
        return native
    return f"GPU: {accelerator[1] or accelerator[0]}  [{native[len('CPU: '):]}]"


def _chemistry(function, chemistry, chem):
    """The chemistry options given under either keyword; ``chem`` is the deprecated spelling."""
    if chem is not None:
        warnings.warn(f"the chem argument of smsd.{function}() is deprecated; use chemistry",
                      DeprecationWarning, stacklevel=3)
        if chemistry is None:
            chemistry = chem
    return _smsd.ChemOptions() if chemistry is None else chemistry


def similarity_upper_bound(g1, g2, chemistry=None, *, chem=None):
    """Linear-time upper bound on the similarity of two molecules, in [0.0, 1.0].

    Accepts MolGraph objects, SMILES strings or RDKit molecules. ``chem`` is
    the deprecated spelling of ``chemistry``.
    """
    options = _chemistry("similarity_upper_bound", chemistry, chem)
    return _smsd.similarity_upper_bound(_ensure_mol(g1), _ensure_mol(g2), options)


def similarity(mol1, mol2):
    """Upper bound on the similarity of two molecules, in [0.0, 1.0].

    Accepts MolGraph objects, SMILES strings or RDKit molecules. The bound is
    the native :func:`similarity_upper_bound`, which cannot exceed 1.0.
    """
    return similarity_upper_bound(mol1, mol2)


def screen_targets(query, targets, threshold, chemistry=None, *, chem=None, num_threads=0):
    """Indices of the targets whose similarity bound to the query reaches ``threshold``.

    Accepts MolGraph objects, SMILES strings or RDKit molecules, or a
    :class:`TargetCorpus`, whose targets are prepared once and so are not
    converted again on each call. The screen runs on ``num_threads`` threads
    (0, the default, uses the OpenMP default, ``OMP_NUM_THREADS``); the
    indices, in target order, do not depend on the thread count. ``chem`` is the
    deprecated spelling of ``chemistry``.
    """
    _checked_int(num_threads, "Thread count", 0)
    options = _chemistry("screen_targets", chemistry, chem)
    query = _ensure_mol(query)
    if isinstance(targets, TargetCorpus):
        _accelerator.discover()  # attaches an installed accelerator before the first corpus screen
        return targets._impl.screen_targets(query, threshold, options, num_threads)
    return _smsd.screen_targets(query, targets, threshold, options, num_threads, _ensure_mol)


def batch_similarity(query, targets, fp_size=2048, radius=2):
    """Rank targets by ECFP Tanimoto similarity to a query.

    Computes ECFP fingerprints for the query and each target and returns
    ``(index, Tanimoto)`` pairs sorted from most to least similar; targets
    with equal scores keep their input order. Two empty fingerprints (empty
    molecules) score 1.0, as in :func:`tanimoto_coefficient`.

    Args:
        query: Query molecule (MolGraph or SMILES string).
        targets: List of target molecules (MolGraph or SMILES strings).
        fp_size (int): Fingerprint size in bits (default 2048).
        radius (int): ECFP radius (default 2).

    Returns:
        list[tuple[int, float]]: Sorted ``(target_index, Tanimoto)`` pairs,
            highest similarity first.

    Example::

        import smsd

        results = smsd.batch_similarity("c1ccccc1", ["c1ccc(O)cc1", "CCCC"])
        for idx, sim in results:
            print(f"Target {idx}: Tanimoto = {sim:.3f}")
    """
    query_fp = _smsd.circular_fingerprint(_ensure_mol(query), radius, fp_size, "ecfp")
    scores = [(index, _smsd.tanimoto(query_fp, _smsd.circular_fingerprint(_ensure_mol(target), radius, fp_size, "ecfp"),
                                     fp_size))
              for index, target in enumerate(targets)]
    scores.sort(key=lambda pair: pair[1], reverse=True)
    return scores


class TargetCorpus:
    """Targets captured once, with their input atom numbering, for repeated searches and screening."""

    def __init__(self):
        self._impl = _smsd.TargetCorpus()
        self._indices = []
        self._lock = threading.RLock()
        # The (graph, indices) pairs that searches read, built on first use and dropped on additions.
        self._targets = None

    @classmethod
    def from_smiles(cls, targets):
        """Create a corpus from molecule objects or SMILES strings."""
        corpus = cls()
        corpus.add_targets(targets)
        return corpus

    def add_target(self, target):
        """Capture a molecule and its current input atom numbering."""
        self.add_targets([target])

    def add_targets(self, targets):
        """Capture targets together; later input changes do not alter them."""
        converted = [_ensure_mol_ex(target) for target in targets]
        indices = [_index_tuple(graph) if source is not None else None for graph, source in converted]
        with self._lock:
            self._impl.add_targets([graph for graph, _ in converted])
            self._indices.extend(indices)
            self._targets = None

    def _snapshot(self):
        """The captured targets as (graph, indices) pairs, shared by every search of this corpus."""
        with self._lock:
            if self._targets is None:
                self._targets = list(zip(self._impl._search_targets(), self._indices))
            return self._targets

    def prewarm(self, *, num_threads=0):
        """Compute graph invariants and ECFP fingerprints once; later searches read the prewarmed graphs."""
        with self._lock:
            self._impl.prewarm(num_threads)
            self._targets = None

    def __len__(self):
        return len(self._impl)

    @property
    def is_prewarmed(self):
        """Whether :meth:`prewarm` has run since the last targets were added."""
        return self._impl.is_prewarmed

    def screen(self, query, threshold, *, num_threads=0):
        """Return target indices meeting the ECFP Tanimoto threshold."""
        _accelerator.discover()
        return self._impl.screen(_ensure_mol(query), threshold, num_threads)

    def screen_many(self, queries, threshold, *, num_threads=0):
        """Return ECFP Tanimoto hit indices for each query in input order."""
        _checked_int(num_threads, "Thread count", 0)
        graphs = [_ensure_mol(query) for query in queries]
        _accelerator.discover()
        return self._impl.screen_many(graphs, threshold, num_threads)

    def screen_targets_many(self, queries, threshold, chemistry=None, *, num_threads=0):
        """Return similarity-bound hit indices for each query in input order."""
        _checked_int(num_threads, "Thread count", 0)
        graphs = [_ensure_mol(query) for query in queries]
        options = _chemistry("screen_targets_many", chemistry, None)
        _accelerator.discover()
        return self._impl.screen_targets_many(graphs, threshold, options, num_threads)
