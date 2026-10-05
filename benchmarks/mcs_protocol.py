#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Shared, explicit protocol for local SMSD and RDKit FMCS measurements.

FMCS can omit an edge between selected query atoms, whereas SMSD returns a
query-vertex mapping. Such a result is retained but excluded from speed claims.
Parsing, graph conversion, validation and warmup are outside the search timer.
"""
from dataclasses import asdict, dataclass
import statistics
import time

import smsd
import smsd._smsd as native
from rdkit import Chem
from rdkit.Chem import rdFMCS


class StrictAtoms(rdFMCS.MCSAtomCompare):
    """Match the additional aromatic-atom constraint in SMSD's strict profile.

    This uses Python callbacks in RDKit: timings must disclose this overhead.
    """
    def __call__(self, params, mol1, i, mol2, j):
        a, b = mol1.GetAtomWithIdx(i), mol2.GetAtomWithIdx(j)
        return (a.GetAtomicNum() == b.GetAtomicNum()
                and a.GetIsAromatic() == b.GetIsAromatic()
                and a.GetFormalCharge() == b.GetFormalCharge()
                and a.IsInRing() == b.IsInRing())


def chem_options(policy):
    if policy == "strict":
        return smsd.ChemOptions.profile("strict")
    chem = smsd.ChemOptions()
    chem.aromaticity_mode = smsd.AromaticityMode.FLEXIBLE
    # LOOSE in historical SMSD versions accepts every unequal bond order.
    # Never present it as RDKit CompareOrder compatibility.
    chem.match_bond_order = (smsd.BondOrderMode.ANY if policy == "any"
                             else smsd.BondOrderMode.STRICT)
    chem.ring_matches_ring_only = False
    chem.complete_rings_only = False
    chem.match_formal_charge = False
    return chem


def mcs_options(timeout_sec):
    opts = smsd.MCSOptions()
    opts.timeout_ms = int(timeout_sec * 1000)
    opts.connected_only = True
    opts.disconnected_mcs = False
    opts.maximize_bonds = False
    opts.induced = False
    opts.max_stage = 5
    return opts


def rdkit_parameters(policy, timeout_sec):
    params = rdFMCS.MCSParameters()
    params.Timeout = int(timeout_sec)
    params.MaximizeBonds = False
    params.Threshold = 1.0
    params.AtomTyper = StrictAtoms() if policy == "strict" else rdFMCS.AtomCompare.CompareElements
    params.BondTyper = {"strict": rdFMCS.BondCompare.CompareOrderExact,
                       "fmcs": rdFMCS.BondCompare.CompareOrder,
                       "any": rdFMCS.BondCompare.CompareAny}[policy]
    params.AtomCompareParameters.MatchFormalCharge = policy == "strict"
    params.AtomCompareParameters.RingMatchesRingOnly = policy == "strict"
    params.BondCompareParameters.RingMatchesRingOnly = policy == "strict"
    params.AtomCompareParameters.CompleteRingsOnly = False
    params.BondCompareParameters.CompleteRingsOnly = False
    return params


def graph_from_rdkit(mol, policy=None):
    """Preserve RDKit atom indices and perception, avoiding a SMILES round trip."""
    n = mol.GetNumAtoms()
    neighbors, orders = [[] for _ in range(n)], [[] for _ in range(n)]
    ring_bonds, aromatic_bonds = [[] for _ in range(n)], [[] for _ in range(n)]
    for bond in mol.GetBonds():
        if policy != "any" and bond.GetBondType() not in (
                Chem.BondType.SINGLE, Chem.BondType.DOUBLE,
                Chem.BondType.TRIPLE, Chem.BondType.AROMATIC):
            raise ValueError(f"unsupported bond type {bond.GetBondType()}")
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        order = 4 if bond.GetIsAromatic() else int(bond.GetBondTypeAsDouble())
        if order not in (1, 2, 3, 4):
            raise ValueError(f"unsupported bond type {bond.GetBondType()}")
        neighbors[i].append(j); orders[i].append(order)
        neighbors[j].append(i); orders[j].append(order)
        ring_bonds[i].append(bond.IsInRing()); ring_bonds[j].append(bond.IsInRing())
        aromatic_bonds[i].append(bond.GetIsAromatic()); aromatic_bonds[j].append(bond.GetIsAromatic())
    builder = (smsd.MolGraphBuilder().atom_count(n)
            .atomic_numbers([a.GetAtomicNum() for a in mol.GetAtoms()])
            .formal_charges([a.GetFormalCharge() for a in mol.GetAtoms()])
            .ring_flags([int(a.IsInRing()) for a in mol.GetAtoms()])
            .aromatic_flags([int(a.GetIsAromatic()) for a in mol.GetAtoms()])
            .neighbors(neighbors).bond_orders(orders))
    if hasattr(builder, "bond_ring_flags") and hasattr(builder, "bond_aromatic_flags"):
        graph = builder.bond_ring_flags(ring_bonds).bond_aromatic_flags(aromatic_bonds).build(False)
    else:
        # Historical bindings could not import bond flags. Re-perceive and
        # verify every property before accepting this graph for comparison.
        graph = builder.build(True)
    assert_graph_parity(mol, graph, policy)
    return graph


def assert_graph_parity(mol, graph, policy=None):
    if policy == "any":
        # Bond-any matching uses element identities and graph topology. Ring,
        # aromaticity and charge perception are not matching constraints.
        if graph.n != mol.GetNumAtoms() or graph.atomic_num != [a.GetAtomicNum() for a in mol.GetAtoms()]:
            raise ValueError("SMSD/RDKit topology or element identity differs")
        edges = {tuple(sorted((i, j))) for i, j, order in graph.bonds()}
        if edges != {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in mol.GetBonds()}:
            raise ValueError("SMSD/RDKit graph topology differs")
        return
    if (graph.n != mol.GetNumAtoms()
        or graph.atomic_num != [a.GetAtomicNum() for a in mol.GetAtoms()]
        or graph.formal_charge != [a.GetFormalCharge() for a in mol.GetAtoms()]
        or [bool(a) for a in graph.aromatic] != [a.GetIsAromatic() for a in mol.GetAtoms()]
        or [bool(a) for a in graph.ring] != [a.IsInRing() for a in mol.GetAtoms()]):
        raise ValueError("SMSD/RDKit atom perception differs; excluded before timing")
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        order = 4 if bond.GetIsAromatic() else int(bond.GetBondTypeAsDouble())
        if (graph.bond_order(i, j) != order
            or graph.bond_aromatic(i, j) != bond.GetIsAromatic()
            or graph.bond_in_ring(i, j) != bond.IsInRing()):
            raise ValueError("SMSD/RDKit bond perception differs; excluded before timing")


def atom_matches(a, b, policy):
    if a.GetAtomicNum() != b.GetAtomicNum():
        return False
    return policy != "strict" or (a.GetFormalCharge() == b.GetFormalCharge()
        and a.GetIsAromatic() == b.GetIsAromatic() and a.IsInRing() == b.IsInRing())


def bond_matches(a, b, policy):
    if b is None:
        return False
    if policy == "strict":
        return (a.GetBondType() == b.GetBondType()
                and a.GetIsAromatic() == b.GetIsAromatic()
                and a.IsInRing() == b.IsInRing())
    if policy == "any":
        return True
    if a.GetBondType() == b.GetBondType():
        return True
    # Exact historical SMSD STRICT + FLEXIBLE semantics, not LOOSE.
    return ((a.GetIsAromatic() and b.GetBondTypeAsDouble() in (1, 2))
            or (b.GetIsAromatic() and a.GetBondTypeAsDouble() in (1, 2))
            or (a.GetBondTypeAsDouble() in (1, 2) and b.GetBondTypeAsDouble() in (1, 2)
                and all(atom.GetIsAromatic() for atom in
                    (a.GetBeginAtom(), a.GetEndAtom(), b.GetBeginAtom(), b.GetEndAtom()))))


def validate_mapping(mol1, mol2, mapping, policy, connected=True):
    """Independent witness check in caller direction on RDKit-perceived graphs."""
    if len(set(mapping.values())) != len(mapping):
        return False, "target atom reused"
    for i, j in mapping.items():
        if not 0 <= i < mol1.GetNumAtoms() or not 0 <= j < mol2.GetNumAtoms():
            return False, "atom index out of range"
        if not atom_matches(mol1.GetAtomWithIdx(i), mol2.GetAtomWithIdx(j), policy):
            return False, "incompatible atom"
    adjacency = {i: [] for i in mapping}
    for bond in mol1.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i in mapping and j in mapping:
            other = mol2.GetBondBetweenAtoms(mapping[i], mapping[j])
            if not bond_matches(bond, other, policy):
                return False, "query edge missing or incompatible"
            adjacency[i].append(j); adjacency[j].append(i)
    if connected and mapping:
        visited, pending = set(), [next(iter(mapping))]
        while pending:
            i = pending.pop()
            if i not in visited:
                visited.add(i); pending.extend(adjacency[i])
        if len(visited) != len(mapping):
            return False, "disconnected query mapping"
    return True, "valid"


def rdkit_witness(mol1, mol2, result, policy, max_matches=128):
    if result.numAtoms == 0:
        return True, "empty", {}
    query = Chem.MolFromSmarts(result.smartsString)
    if query is None:
        return False, "invalid result SMARTS", {}
    left = mol1.GetSubstructMatches(query, uniquify=False, maxMatches=max_matches)
    right = mol2.GetSubstructMatches(query, uniquify=False, maxMatches=max_matches)
    for a in left:
        for b in right:
            mapping = dict(zip(a, b))
            valid, reason = validate_mapping(mol1, mol2, mapping, policy)
            if valid:
                return True, reason, mapping
    # This is a bounded witness search, not proof that every FMCS embedding is invalid.
    return False, f"no query-vertex witness in first {max_matches} embeddings", {}


@dataclass
class Observation:
    engine: str
    elapsed_us: float
    atoms: int
    bonds: int
    valid: bool
    validation: str
    timeout: object
    budget_reached: bool
    error: str = ""


def measure_pair(smi1, smi2, policy="strict", timeout_sec=10, warmup=1, repeats=3,
                 index=0, on_observation=None, strategy="native"):
    mol1, mol2 = Chem.MolFromSmiles(smi1), Chem.MolFromSmiles(smi2)
    if mol1 is None or mol2 is None:
        raise ValueError("RDKit rejected input SMILES; no cross-engine timing")
    g1, g2 = graph_from_rdkit(mol1, policy), graph_from_rdkit(mol2, policy)
    chem, opts = chem_options(policy), mcs_options(timeout_sec)
    params = rdkit_parameters(policy, timeout_sec)
    calls = {"smsd": lambda: native.find_mcs(g1, g2, chem, opts),
             "rdkit": lambda: rdFMCS.FindMCS([mol1, mol2], params)}
    if strategy == "coverage":
        if policy == "strict":
            raise ValueError("raw coverage cannot express the strict aromatic-atom profile")
        calls["smsd"] = lambda: native.find_mcs_coverage(g1, g2,
            ring_match=False, bond_any=policy == "any", timeout_ms=int(timeout_sec*1000))
    elif strategy != "native":
        if policy == "strict":
            raise ValueError("strict aromatic-atom profile is not expressible through the public auto/lightweight API")
        calls["smsd"] = lambda: smsd.find_mcs(g1, g2, strategy=strategy,
            match_bond_order="any" if policy == "any" else "strict",
            connected_only=True, maximize_bonds=False,
            timeout_ms=int(timeout_sec*1000))
    for _ in range(warmup):
        for engine in ("smsd", "rdkit"):
            calls[engine]()
    rows = []
    for trial in range(repeats):
        order = ("smsd", "rdkit") if (index + trial) % 2 == 0 else ("rdkit", "smsd")
        for engine in order:
            start = time.perf_counter_ns()
            try:
                result = calls[engine]()
                elapsed = (time.perf_counter_ns() - start) / 1000
                if engine == "smsd":
                    valid, reason = validate_mapping(mol1, mol2, result, policy)
                    bonds = sum(b.GetBeginAtomIdx() in result and b.GetEndAtomIdx() in result
                                for b in mol1.GetBonds())
                    row = Observation(engine, elapsed, len(result), bonds, valid, reason,
                                      None, elapsed >= timeout_sec * 1e6)
                else:
                    valid, reason, mapping = rdkit_witness(mol1, mol2, result, policy)
                    row = Observation(engine, elapsed, result.numAtoms, result.numBonds,
                                      valid, reason, bool(result.canceled),
                                      elapsed >= timeout_sec * 1e6)
            except Exception as exc:
                row = Observation(engine, (time.perf_counter_ns()-start)/1000,
                                  -1, -1, False, "exception", None, False, str(exc))
            data = asdict(row); data["trial"] = trial
            # A mapping-only call that ends just below its deadline may have
            # exhausted a millisecond-granularity budget. Do not rank its speed
            # as an uncanceled completion merely because the outer timer is short.
            data["near_budget"] = row.elapsed_us >= timeout_sec*990000
            if not row.error:
                data["mapping"] = result if engine == "smsd" else mapping
                if engine == "rdkit":
                    data["smarts"] = result.smartsString
            rows.append(data)
            if on_observation is not None:
                on_observation(data)
    return rows


def summarize(rows, policy):
    s, r = [[x for x in rows if x["engine"] == e] for e in ("smsd", "rdkit")]
    result = {}
    for engine, values in (("smsd", s), ("rdkit", r)):
        result[engine + "_median_us"] = statistics.median(x["elapsed_us"] for x in values)
        result[engine + "_atoms"] = sorted(set(x["atoms"] for x in values))
        result[engine + "_valid"] = all(x["valid"] for x in values)
        result[engine + "_timeouts"] = (sum(x["timeout"] is True for x in values)
                                         if engine == "rdkit" else None)
        result[engine + "_budget_reached"] = sum(x["budget_reached"] for x in values)
        result[engine + "_near_budget"] = sum(x.get("near_budget", False) for x in values)
    # FMCS CompareOrder and historical SMSD flexible aromaticity differ.
    # Strict mode adds a Python atom-comparison callback to RDKit.
    result["rdkit_python_comparator"] = policy == "strict"
    result["chemical_policy_equivalent"] = policy in ("strict", "any")
    result["speed_comparable"] = (result["chemical_policy_equivalent"]
        and result["smsd_valid"] and result["rdkit_valid"]
        and result["smsd_atoms"] == result["rdkit_atoms"]
        and len(result["smsd_atoms"]) == 1 and result["rdkit_timeouts"] == 0
        and result["smsd_budget_reached"] == 0 and result["smsd_near_budget"] == 0
        and not result["rdkit_python_comparator"])
    result["rdkit_over_smsd_ratio"] = (result["rdkit_median_us"] / result["smsd_median_us"]
        if result["speed_comparable"] and result["smsd_median_us"] > 0 else None)
    return result
