"""Lightweight MCS engine with coverage-driven funnel.

Stages escalate from polynomial heuristics to exponential solvers
only when cheap stages leave a gap below the label-frequency upper bound:

  L0.75  Greedy atom-by-atom (Morgan rank + NLF)
  L1     Substructure containment (native C++ VF2 or RDKit fallback)
  L1.5   Seed-and-extend from heteroatom anchors
  L3     Bron-Kerbosch / C++ clique solver on modular product graph
  L4     McGregor bond-grow backtracking

Works with SMSD MolGraph, RDKit Mol, or SMILES strings.

Usage::

    import smsd
    from smsd.mcs_engine import find_mcs_lightweight, MCSConfig

    # From SMILES
    result = find_mcs_lightweight("c1ccc(O)cc1", "c1ccc(N)cc1")

    # From MolGraph
    g1 = smsd.parse_smiles("c1ccc(O)cc1")
    g2 = smsd.parse_smiles("c1ccc(N)cc1")
    result = find_mcs_lightweight(g1, g2)

    # With options
    result = find_mcs_lightweight(g1, g2, timeout=2.0,
                                  ring_matches_ring=True)

    print(result.size, result.mapping, result.candidates)
"""

from __future__ import annotations

import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Optional, Union

try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

# Element symbol table (atomic number -> symbol)
_ELEMENT_SYMBOL = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O",
    9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",
    16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 25: "Mn", 26: "Fe",
    27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 33: "As", 34: "Se", 35: "Br",
    44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 50: "Sn", 51: "Sb",
    52: "Te", 53: "I", 54: "Xe", 78: "Pt", 79: "Au", 80: "Hg", 82: "Pb",
}


# ─── Adapter: wrap SMSD MolGraph to the interface the engine expects ────────

class _Atom:
    __slots__ = ("index", "element")
    def __init__(self, index, element):
        self.index = index
        self.element = element


class MolGraphWrapper:
    """Adapt an SMSD MolGraph (C++ binding) to the dict-based interface
    that the coverage-funnel engine uses internally."""

    def __init__(self, g):
        self._g = g
        n = len(g.atomic_num)
        self.n_atoms = n

        # Atoms (1-indexed to match engine convention)
        self.atoms = []
        for i in range(n):
            sym = _ELEMENT_SYMBOL.get(g.atomic_num[i], f"#{g.atomic_num[i]}")
            self.atoms.append(_Atom(i + 1, sym))

        # Bond lookup: {(min_1idx, max_1idx): order}
        self.bond_lookup = {}
        # Adjacency: {atom_1idx: [(neighbor_1idx, order), ...]}
        self.adjacency = defaultdict(list)
        for i in range(n):
            for j in range(i + 1, n):
                bo = g.bond_order(i, j)
                if bo > 0:
                    a, b = i + 1, j + 1
                    self.bond_lookup[(a, b)] = bo
                    self.bond_lookup[(b, a)] = bo
                    self.adjacency[a].append((b, bo))
                    self.adjacency[b].append((a, bo))

        # Ring / aromatic flags: {atom_1idx: bool}
        self.ring_flags = {i + 1: bool(g.ring[i]) for i in range(n)}
        self.aromatic_flags = {i + 1: bool(g.aromatic[i]) for i in range(n)}

        # Element frequency
        self.element_counts = Counter(a.element for a in self.atoms)

        # RDKit mol (lazy, for fallback paths)
        self._rdkit_mol = None

    def atom(self, idx):
        return self.atoms[idx - 1]

    @property
    def rdkit_mol(self):
        if self._rdkit_mol is None and HAS_RDKIT:
            try:
                import smsd._smsd as _smsd
                smi = _smsd.to_smiles(self._g)
                self._rdkit_mol = Chem.MolFromSmiles(smi)
            except Exception:
                pass
        return self._rdkit_mol


@dataclass
class MCSConfig:
    """Configuration for the lightweight MCS engine."""
    ring_matches_ring_only: bool = True
    bond_compare_name: str = "CompareOrder"  # "CompareOrder" or "CompareAny"


@dataclass
class LightMCSResult:
    """Result from the lightweight MCS engine."""
    size: int = 0
    mapping: list = field(default_factory=list)
    candidates: list = field(default_factory=list)
    lfub: int = 0
    elapsed_ms: float = 0.0


def _ensure_wrapper(mol) -> MolGraphWrapper:
    """Convert SMILES, MolGraph, or RDKit Mol to MolGraphWrapper."""
    if isinstance(mol, MolGraphWrapper):
        return mol
    if isinstance(mol, str):
        import smsd._smsd as _smsd
        return MolGraphWrapper(_smsd.parse_smiles(mol))
    # Check if it's an SMSD MolGraph (has atomic_num attribute)
    if hasattr(mol, 'atomic_num') and hasattr(mol, 'bond_order'):
        return MolGraphWrapper(mol)
    # RDKit Mol
    if HAS_RDKIT and hasattr(mol, 'GetNumAtoms'):
        smi = Chem.MolToSmiles(mol)
        import smsd._smsd as _smsd
        return MolGraphWrapper(_smsd.parse_smiles(smi))
    raise TypeError(f"Cannot convert {type(mol)} to molecule wrapper")


def find_mcs_lightweight(
    mol1, mol2, *,
    timeout: float = 1.0,
    ring_matches_ring: bool = True,
    bond_any: bool = False,
) -> LightMCSResult:
    """Find MCS using the lightweight coverage-driven pipeline.

    Accepts SMILES strings, SMSD MolGraph objects, or RDKit Mol objects.
    Uses the C++ clique solver backend for the heavy combinatorial work.

    Args:
        mol1: First molecule (SMILES, MolGraph, or RDKit Mol).
        mol2: Second molecule.
        timeout: Wall-clock timeout in seconds (default 1.0).
        ring_matches_ring: Ring atom only matches ring atom.
        bond_any: If True, any bond matches any bond.

    Returns:
        LightMCSResult with size, mapping, candidates, lfub, elapsed_ms.
    """
    t0 = time.monotonic()
    w1 = _ensure_wrapper(mol1)
    w2 = _ensure_wrapper(mol2)
    config = MCSConfig(
        ring_matches_ring_only=ring_matches_ring,
        bond_compare_name="CompareAny" if bond_any else "CompareOrder",
    )

    candidates = python_mcs_candidates(w1, w2, config, timeout_sec=timeout)

    best = max(candidates, key=len) if candidates else ()
    elapsed = (time.monotonic() - t0) * 1000

    return LightMCSResult(
        size=len(best),
        mapping=list(best),
        candidates=[list(c) for c in candidates],
        lfub=lfub(w1, w2),
        elapsed_ms=elapsed,
    )


# ---------------------------------------------------------------------------
# Atom / bond compatibility
# ---------------------------------------------------------------------------

def _atoms_compatible(mol_a, idx_a: int, mol_b, idx_b: int, config) -> bool:
    """Check if two atoms can be matched under the given config."""
    a = mol_a.atom(idx_a)
    b = mol_b.atom(idx_b)
    if a.element != b.element:
        return False
    if config.ring_matches_ring_only:
        a_ring = mol_a.ring_flags.get(idx_a, False)
        b_ring = mol_b.ring_flags.get(idx_b, False)
        if a_ring != b_ring:
            return False
    return True


def _bonds_compatible(mol_a, a1: int, a2: int, mol_b, b1: int, b2: int, config) -> bool:
    """Check if two bonds can be matched under the given config."""
    key_a = (min(a1, a2), max(a1, a2))
    key_b = (min(b1, b2), max(b1, b2))
    order_a = mol_a.bond_lookup.get(key_a)
    order_b = mol_b.bond_lookup.get(key_b)
    if order_a is None or order_b is None:
        return False
    if order_a == order_b:
        return True
    # Aromatic equivalence: both endpoints aromatic in both molecules
    if (mol_a.aromatic_flags.get(a1, False) and mol_a.aromatic_flags.get(a2, False)
            and mol_b.aromatic_flags.get(b1, False) and mol_b.aromatic_flags.get(b2, False)):
        return True
    # Bond-order insensitive mode (CompareAny)
    if config.bond_compare_name == "CompareAny":
        return True
    return False


# ---------------------------------------------------------------------------
# LFUB — Label-Frequency Upper Bound (Theorem 1)
# ---------------------------------------------------------------------------

def lfub(mol_a, mol_b) -> int:
    """Maximum possible MCS size based on element frequency overlap."""
    return sum((mol_a.element_counts & mol_b.element_counts).values())


# ---------------------------------------------------------------------------
# L0.75 — Greedy Probe
# ---------------------------------------------------------------------------

def _morgan_rank_order(mol, strategy: str = "rarest") -> list[int]:
    """Return atom indices sorted by priority for greedy matching.

    Strategies:
      rarest  — rarest element first, then highest degree
      degree  — highest degree first, then rarest element
      ring    — ring atoms first, then rarest
    """
    elem_counts = mol.element_counts
    atoms = []
    for a in mol.atoms:
        rarity = elem_counts.get(a.element, 0)
        degree = len(mol.adjacency.get(a.index, []))
        is_ring = 1 if mol.ring_flags.get(a.index, False) else 0
        if strategy == "degree":
            atoms.append((-degree, rarity, a.index))
        elif strategy == "ring":
            atoms.append((-is_ring, rarity, -degree, a.index))
        else:
            atoms.append((rarity, -degree, a.index))
    atoms.sort()
    return [idx for *_, idx in atoms]


def greedy_probe(mol_a, mol_b, config, strategy: str = "rarest") -> list[tuple[int, int]]:
    """Greedy atom-by-atom matching with configurable ordering.

    Scores each candidate pair by:
      - Bond preservation to mapped neighbors (weight 4.0, +2.0 if order matches)
      - Degree similarity (weight 1.0)
      - Ring flag match (weight 1.5)
      - Aromatic flag match (weight 1.5)
    """
    order_a = _morgan_rank_order(mol_a, strategy)
    used_b: set[int] = set()
    mapping: list[tuple[int, int]] = []
    mapped_a: set[int] = set()
    a_to_b: dict[int, int] = {}

    deg_a = {a.index: len(mol_a.adjacency.get(a.index, [])) for a in mol_a.atoms}
    deg_b = {b.index: len(mol_b.adjacency.get(b.index, [])) for b in mol_b.atoms}

    for idx_a in order_a:
        if idx_a in mapped_a:
            continue
        best_b, best_score = None, -0.5
        for atom_b in mol_b.atoms:
            if atom_b.index in used_b:
                continue
            if not _atoms_compatible(mol_a, idx_a, mol_b, atom_b.index, config):
                continue

            score = 0.0
            # Bond preservation to mapped neighbors
            for nb_a, _ in mol_a.adjacency.get(idx_a, []):
                mapped_nb_b = a_to_b.get(nb_a)
                if mapped_nb_b is None:
                    continue
                if _bonds_compatible(mol_a, idx_a, nb_a, mol_b, atom_b.index, mapped_nb_b, config):
                    score += 4.0
                    # Bonus for exact bond order match
                    key_a = (min(idx_a, nb_a), max(idx_a, nb_a))
                    key_b = (min(atom_b.index, mapped_nb_b), max(atom_b.index, mapped_nb_b))
                    if mol_a.bond_lookup.get(key_a) == mol_b.bond_lookup.get(key_b):
                        score += 2.0

            # Degree similarity
            if deg_a.get(idx_a, 0) == deg_b.get(atom_b.index, 0):
                score += 1.0

            # Ring match
            if mol_a.ring_flags.get(idx_a, False) == mol_b.ring_flags.get(atom_b.index, False):
                score += 1.5

            # Aromatic match
            if mol_a.aromatic_flags.get(idx_a, False) == mol_b.aromatic_flags.get(atom_b.index, False):
                score += 1.5

            if score > best_score:
                best_score = score
                best_b = atom_b.index

        if best_b is not None:
            mapping.append((idx_a, best_b))
            mapped_a.add(idx_a)
            used_b.add(best_b)
            a_to_b[idx_a] = best_b

    return _largest_connected_mapping(mol_a, mapping)


def _largest_connected_mapping(mol_a, mapping: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Keep only the largest connected component of the mapping."""
    if len(mapping) <= 1:
        return mapping
    mapped_a = {a for a, _ in mapping}
    adj = defaultdict(set)
    for a, _ in mapping:
        for nb, _ in mol_a.adjacency.get(a, []):
            if nb in mapped_a:
                adj[a].add(nb)
                adj[nb].add(a)

    # BFS to find components
    visited: set[int] = set()
    best_component: list[int] = []
    for a, _ in mapping:
        if a in visited:
            continue
        component = []
        queue = [a]
        while queue:
            cur = queue.pop(0)
            if cur in visited:
                continue
            visited.add(cur)
            component.append(cur)
            for nb in adj.get(cur, []):
                if nb not in visited:
                    queue.append(nb)
        if len(component) > len(best_component):
            best_component = component

    best_set = set(best_component)
    return [(a, b) for a, b in mapping if a in best_set]


# ---------------------------------------------------------------------------
# L1 — Substructure Containment (native VF2 or RDKit fallback)
# ---------------------------------------------------------------------------

def _native_substructure_check(
    query, target, flip: bool, max_matches: int = 8,
) -> list[list[tuple[int, int]]]:
    """VF2 substructure match using our C++ engine.

    Compat table built in C++ from element lists — no Python O(Q×T) loop.
    """
    try:
        from smsd._smsd import match_substructure_from_elements
    except ImportError:
        return []
    elements_q = [a.element for a in query.atoms]
    elements_t = [a.element for a in target.atoms]
    bonds_q = {(min(a - 1, b - 1), max(a - 1, b - 1)): o for (a, b), o in query.bond_lookup.items()}
    bonds_t = {(min(a - 1, b - 1), max(a - 1, b - 1)): o for (a, b), o in target.bond_lookup.items()}
    try:
        embeddings = match_substructure_from_elements(
            elements_q, elements_t, bonds_q, bonds_t,
            False, 500, max_matches,
        )
    except Exception:
        return []
    results = []
    for emb in embeddings:
        if len(emb) != query.n_atoms:
            continue
        if flip:
            results.append([(t + 1, q + 1) for q, t in sorted(emb)])
        else:
            results.append([(q + 1, t + 1) for q, t in sorted(emb)])
    return results


def substructure_check(mol_a, mol_b, config, max_matches: int = 8) -> list[list[tuple[int, int]]]:
    """If smaller molecule is a substructure of larger, return mappings."""
    if mol_a.n_atoms == mol_b.n_atoms:
        return []
    if mol_a.n_atoms < mol_b.n_atoms:
        query, target = mol_a, mol_b
        flip = False
    else:
        query, target = mol_b, mol_a
        flip = True

    # Try native C++ VF2 first
    native = _native_substructure_check(query, target, flip, max_matches)
    if native:
        return native

    # RDKit fallback
    if not HAS_RDKIT:
        return []
    q_mol = query.rdkit_mol
    t_mol = target.rdkit_mol
    if q_mol is None or t_mol is None:
        return []

    try:
        matches = t_mol.GetSubstructMatches(q_mol, uniquify=False, maxMatches=max_matches)
    except Exception:
        return []

    results = []
    for match in matches:
        if flip:
            results.append([(match[i] + 1, i + 1) for i in range(len(match))])
        else:
            results.append([(i + 1, match[i] + 1) for i in range(len(match))])
    return results


# ---------------------------------------------------------------------------
# L1.5 — Seed-and-Extend from Heteroatom Anchors
# ---------------------------------------------------------------------------

_HETEROATOMS = {"N", "O", "S", "P", "F", "Cl", "Br", "I", "Si", "Se"}


def seed_and_extend(mol_a, mol_b, config, deadline: float = 0.0) -> Optional[list[tuple[int, int]]]:
    """Seed MCS from compatible heteroatom pairs, extend by BFS bond growth."""
    best_mapping: Optional[list[tuple[int, int]]] = None
    best_size = 0

    # Collect heteroatom indices
    hetero_a = [(a.index, a.element) for a in mol_a.atoms if a.element in _HETEROATOMS]
    hetero_b = [(b.index, b.element) for b in mol_b.atoms if b.element in _HETEROATOMS]

    if not hetero_a or not hetero_b:
        return None

    for idx_a, elem_a in hetero_a:
        for idx_b, elem_b in hetero_b:
            if deadline and time.monotonic() > deadline:
                return best_mapping
            if elem_a != elem_b:
                continue
            if not _atoms_compatible(mol_a, idx_a, mol_b, idx_b, config):
                continue

            # Seed and BFS-extend
            mapping = _bfs_extend(mol_a, mol_b, config, [(idx_a, idx_b)])
            if len(mapping) > best_size:
                best_size = len(mapping)
                best_mapping = mapping

    return best_mapping


def _bfs_extend(mol_a, mol_b, config, seed: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """BFS-extend a seed mapping along compatible bonds.

    For each unmapped neighbor of a mapped atom, scores ALL compatible
    product atoms and picks the best (by bond preservation + ring/aromatic).
    """
    mapping = list(seed)
    mapped_a = {a for a, _ in mapping}
    mapped_b = {b for _, b in mapping}
    a_to_b = {a: b for a, b in mapping}
    queue = list(seed)

    while queue:
        cur_a, cur_b = queue.pop(0)
        for nb_a, _ in mol_a.adjacency.get(cur_a, []):
            if nb_a in mapped_a:
                continue
            # Score all compatible product neighbors
            best_nb_b, best_score = None, -1.0
            for nb_b, _ in mol_b.adjacency.get(cur_b, []):
                if nb_b in mapped_b:
                    continue
                if not _atoms_compatible(mol_a, nb_a, mol_b, nb_b, config):
                    continue
                if not _bonds_compatible(mol_a, cur_a, nb_a, mol_b, cur_b, nb_b, config):
                    continue
                score = 1.0
                # Bonus: additional bond preservation to other mapped neighbors
                for nb2_a, _ in mol_a.adjacency.get(nb_a, []):
                    nb2_b = a_to_b.get(nb2_a)
                    if nb2_b is not None and _bonds_compatible(mol_a, nb_a, nb2_a, mol_b, nb_b, nb2_b, config):
                        score += 3.0
                # Ring/aromatic match
                if mol_a.ring_flags.get(nb_a, False) == mol_b.ring_flags.get(nb_b, False):
                    score += 1.0
                if mol_a.aromatic_flags.get(nb_a, False) == mol_b.aromatic_flags.get(nb_b, False):
                    score += 1.0
                if score > best_score:
                    best_score = score
                    best_nb_b = nb_b

            if best_nb_b is not None:
                mapping.append((nb_a, best_nb_b))
                mapped_a.add(nb_a)
                mapped_b.add(best_nb_b)
                a_to_b[nb_a] = best_nb_b
                queue.append((nb_a, best_nb_b))

    return mapping


# ---------------------------------------------------------------------------
# L3 — Bron-Kerbosch on Modular Product Graph
# ---------------------------------------------------------------------------

def _try_cpp_mcs_pipeline(
    mol_a, mol_b, config,
    deadline: float = 0.0,
) -> Optional[list[list[tuple[int, int]]]]:
    """Full MCS pipeline in C++: greedy + seed-extend + BK + McGregor."""
    try:
        from smsd._smsd import find_mcs_clique
    except ImportError:
        return None

    # Build atom compatibility table (Python — lightweight, O(n*m))
    compat: list[tuple[int, int]] = []
    for a in mol_a.atoms:
        for b in mol_b.atoms:
            if _atoms_compatible(mol_a, a.index, mol_b, b.index, config):
                compat.append((a.index - 1, b.index - 1))
    if not compat:
        return None

    # Bond lookup tables (0-based)
    bonds_a = {(min(a-1, b-1), max(a-1, b-1)): o for (a, b), o in mol_a.bond_lookup.items()}
    bonds_b = {(min(a-1, b-1), max(a-1, b-1)): o for (a, b), o in mol_b.bond_lookup.items()}
    bond_any = config.bond_compare_name == "CompareAny"

    remaining_ms = max(50, int((deadline - time.monotonic()) * 1000)) if deadline else 1000

    # Ring/aromatic flags (0-based, one per atom)
    ring_a = [mol_a.ring_flags.get(a.index, False) for a in mol_a.atoms]
    ring_b = [mol_b.ring_flags.get(b.index, False) for b in mol_b.atoms]
    arom_a = [mol_a.aromatic_flags.get(a.index, False) for a in mol_a.atoms]
    arom_b = [mol_b.aromatic_flags.get(b.index, False) for b in mol_b.atoms]

    # True LFUB from Python
    true_lfub = lfub(mol_a, mol_b)

    # Single C++ call: compat → full pipeline → mappings
    result = find_mcs_clique(
        compat, bonds_a, bonds_b,
        mol_a.n_atoms, mol_b.n_atoms,
        bond_any, remaining_ms, 8,
        ring_a, ring_b, arom_a, arom_b,
        true_lfub,
    )

    if result.best_size == 0:
        return None

    # Convert to 1-indexed
    mappings = []
    for candidate in result.candidates:
        mapping = [(q + 1, t + 1) for q, t in candidate]
        mappings.append(mapping)
    return mappings if mappings else None


def bron_kerbosch_mcs(
    mol_a, mol_b, config,
    incumbent: list[tuple[int, int]],
    deadline: float = 0.0,
) -> Optional[list[list[tuple[int, int]]]]:
    """Maximum clique in modular product graph = connected MCS."""
    # Build compatibility table
    compat: list[tuple[int, int]] = []
    for a in mol_a.atoms:
        for b in mol_b.atoms:
            if _atoms_compatible(mol_a, a.index, mol_b, b.index, config):
                compat.append((a.index, b.index))

    if not compat:
        return None

    n = len(compat)
    if n > 2000:
        return None  # Too large, skip

    # Build adjacency for modular product graph
    # Edge between (a1,b1) and (a2,b2) if:
    #   a1 != a2, b1 != b2, and bonds (a1,a2) and (b1,b2) are compatible
    idx_map = {pair: i for i, pair in enumerate(compat)}
    adj: list[set[int]] = [set() for _ in range(n)]

    for i in range(n):
        a1, b1 = compat[i]
        for j in range(i + 1, n):
            if deadline and (j & 0xFF) == 0 and time.monotonic() > deadline:
                break
            a2, b2 = compat[j]
            if a1 == a2 or b1 == b2:
                continue
            if _bonds_compatible(mol_a, a1, a2, mol_b, b1, b2, config):
                adj[i].add(j)
                adj[j].add(i)

    # k-core pruning: remove vertices with degree < incumbent size
    best_size = len(incumbent)
    alive = [True] * n
    changed = True
    while changed:
        changed = False
        for i in range(n):
            if alive[i] and len(adj[i]) < best_size:
                alive[i] = False
                for j in adj[i]:
                    adj[j].discard(i)
                adj[i].clear()
                changed = True

    active = [i for i in range(n) if alive[i]]
    if len(active) < best_size:
        return None

    # Bron-Kerbosch with pivoting — collect ALL maximum cliques
    best_clique: list[int] = []
    all_max_cliques: list[list[int]] = []
    _MAX_CLIQUES = 8

    def bk(R: list[int], P: set[int], X: set[int]) -> None:
        nonlocal best_clique
        if deadline and time.monotonic() > deadline:
            return
        if len(all_max_cliques) >= _MAX_CLIQUES:
            return
        if not P and not X:
            if len(R) > len(best_clique):
                best_clique = list(R)
                all_max_cliques.clear()
                all_max_cliques.append(list(R))
            elif len(R) == len(best_clique) and len(best_clique) > 0:
                all_max_cliques.append(list(R))
            return
        pivot = max(P | X, key=lambda v: len(adj[v] & P))
        for v in list(P - adj[pivot]):
            if deadline and time.monotonic() > deadline:
                return
            if len(all_max_cliques) >= _MAX_CLIQUES:
                return
            new_P = P & adj[v]
            new_X = X & adj[v]
            bk(R + [v], new_P, new_X)
            P.discard(v)
            X.add(v)

    bk([], set(active), set())

    if len(best_clique) <= best_size:
        return None

    if not all_max_cliques:
        if not best_clique:
            return None
        all_max_cliques = [best_clique]

    # Return ALL maximum cliques as list of lists for embedding diversity
    return [
        [(compat[i][0], compat[i][1]) for i in clique]
        for clique in all_max_cliques
    ]


# ---------------------------------------------------------------------------
# L4 — McGregor Extension (backtracking)
# ---------------------------------------------------------------------------

def mcgregor_extend(
    mol_a, mol_b, config,
    seed: list[tuple[int, int]],
    deadline: float = 0.0,
) -> Optional[list[tuple[int, int]]]:
    """Extend seed mapping by backtracking, ensuring connectivity."""
    if not seed:
        return None

    mapped_a = {a for a, _ in seed}
    mapped_b = {b for _, b in seed}
    best = list(seed)

    # Candidate atoms: unmapped neighbors of mapped atoms in mol_a
    def get_frontier(current_mapped_a: set[int]) -> list[int]:
        frontier = []
        for a in current_mapped_a:
            for nb, _ in mol_a.adjacency.get(a, []):
                if nb not in current_mapped_a and nb not in frontier:
                    frontier.append(nb)
        return frontier

    def backtrack(mapping: list[tuple[int, int]], m_a: set[int], m_b: set[int]) -> None:
        nonlocal best
        if deadline and time.monotonic() > deadline:
            return
        if len(mapping) > len(best):
            best = list(mapping)

        frontier = get_frontier(m_a)
        if not frontier:
            return

        for a in frontier:
            candidates = []
            for b_atom in mol_b.atoms:
                if b_atom.index in m_b:
                    continue
                if not _atoms_compatible(mol_a, a, mol_b, b_atom.index, config):
                    continue
                # Check at least one bond to mapped neighbor is compatible
                has_bond = False
                for nb_a, _ in mol_a.adjacency.get(a, []):
                    if nb_a not in m_a:
                        continue
                    nb_b = None
                    for ma, mb in mapping:
                        if ma == nb_a:
                            nb_b = mb
                            break
                    if nb_b is not None and _bonds_compatible(mol_a, a, nb_a, mol_b, b_atom.index, nb_b, config):
                        has_bond = True
                        break
                if has_bond:
                    candidates.append(b_atom.index)

            for b in candidates:
                mapping.append((a, b))
                m_a.add(a)
                m_b.add(b)
                backtrack(mapping, m_a, m_b)
                mapping.pop()
                m_a.discard(a)
                m_b.discard(b)

    backtrack(list(seed), set(mapped_a), set(mapped_b))
    return best if len(best) > len(seed) else None


# ---------------------------------------------------------------------------
# Main dispatcher
# ---------------------------------------------------------------------------

def python_mcs_candidates(mol_a, mol_b, config, timeout_sec: float = 1.0) -> list[tuple[tuple[int, int], ...]]:
    """BK + McGregor MCS with LFUB termination.

    Primary path: Bron-Kerbosch for optimal MCS, McGregor for refinement.
    Fast paths: substructure containment and LFUB early-exit.

    Returns list of 1-indexed atom pair tuples, compatible with
    rdkit_mapper's _mcs_candidates interface.
    """
    ub = lfub(mol_a, mol_b)
    if ub == 0:
        return []

    deadline = time.monotonic() + timeout_sec
    best: list[tuple[int, int]] = []
    candidates: list[list[tuple[int, int]]] = []

    # L0.75: Greedy probe — O(Nd), handles ~45% of pairs optimally
    greedy = greedy_probe(mol_a, mol_b, config)
    if len(greedy) > len(best):
        best = greedy
    if greedy:
        candidates.append(greedy)
    # C++ MCS pipeline: greedy + seed-extend + VF2 all in C++
    cpp_ran = False
    cpp_results = None
    if time.monotonic() < deadline:
        cpp_results = _try_cpp_mcs_pipeline(mol_a, mol_b, config, deadline)
        cpp_ran = cpp_results is not None
        if cpp_results:
            for mapping in cpp_results:
                if len(mapping) > len(best):
                    best = mapping
                candidates.append(mapping)

    if len(best) >= ub:
        return _dedup_candidates(candidates, best)

    # L1: Substructure containment — O(NM), multiple embeddings
    subs = substructure_check(mol_a, mol_b, config, max_matches=8)
    for sub in subs:
        if len(sub) > len(best):
            best = sub
        candidates.append(sub)
    if len(best) >= ub:
        return _dedup_candidates(candidates, best)

    # L1.5: Seed-and-extend from heteroatom anchors
    if time.monotonic() < deadline:
        seed = seed_and_extend(mol_a, mol_b, config, deadline)
        if seed and len(seed) > len(best):
            best = seed
        if seed:
            candidates.append(seed)
        if len(best) >= ub:
            return _dedup_candidates(candidates, best)

    # L3: Bron-Kerbosch (Python fallback only when C++ unavailable)
    if not cpp_ran and time.monotonic() < deadline:
        bk_cliques = bron_kerbosch_mcs(mol_a, mol_b, config, best, deadline)
        if bk_cliques:
            for clique in bk_cliques:
                if len(clique) > len(best):
                    best = clique
                candidates.append(clique)
        if len(best) >= ub:
            return _dedup_candidates(candidates, best)

    # L4: McGregor backtracking — refine from best seed
    if time.monotonic() < deadline and best:
        mg = mcgregor_extend(mol_a, mol_b, config, best, deadline)
        if mg and len(mg) > len(best):
            best = mg
        if mg:
            candidates.append(mg)

    return _dedup_candidates(candidates, best)


def _dedup_candidates(
    candidates: list[list[tuple[int, int]]], best: list[tuple[int, int]],
) -> list[tuple[tuple[int, int], ...]]:
    """Deduplicate candidates and return as sorted tuples."""
    seen: set[tuple[tuple[int, int], ...]] = set()
    result: list[tuple[tuple[int, int], ...]] = []
    for c in candidates:
        key = tuple(sorted(c))
        if key and key not in seen:
            seen.add(key)
            result.append(key)
    if not result and best:
        result.append(tuple(sorted(best)))
    return result
