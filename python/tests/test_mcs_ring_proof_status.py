# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Ring-dependent MCS proof statuses through the public Python API."""

import pytest
import smsd


# An 80-layer hexagonal prism stack: a synthetic cage beyond the ring-perception work limit.
TRUNCATED_STACK_SMILES = (
    'C12C3C4C5C6C1C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C6'
    '8C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C4'
    '3C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C6'
    '8C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C4'
    '3C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C6'
    '8C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C4'
    '3C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C6'
    '8C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C4'
    '3C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C6'
    '8C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C4'
    '3C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C6'
    '8C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C4'
    '3C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C6'
    '8C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C4'
    '3C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C6'
    '8C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C1C6C5C4C3C21'
)


@pytest.fixture(scope="module")
def ring_graphs():
    stack = smsd.parse_smiles(TRUNCATED_STACK_SMILES)
    ring = smsd.parse_smiles("C1CCC1")
    assert smsd.rings_truncated(stack)
    assert not smsd.rings_truncated(ring)
    return stack, ring


def _options(objective="MCIS", connectivity="CONNECTED", *, complete_rings=True):
    options = smsd.MCS.Options()
    options.objective = getattr(smsd.MCS.Objective, objective)
    options.connectivity = getattr(smsd.MCS.Connectivity, connectivity)
    options.chemistry.complete_rings_only = complete_rings
    options.timeout_ms = 10_000
    return options


def _assert_ring_proof_unavailable(result):
    assert result.status == smsd.SearchStatus.LIMIT_REACHED
    assert not result.optimal
    assert not result.complete


@pytest.mark.parametrize("objective", ["MCIS", "MCES"])
@pytest.mark.parametrize("connectivity", ["CONNECTED", "FRAGMENTS"])
@pytest.mark.parametrize("reverse", [False, True])
def test_complete_rings_keep_search_results_unproved(ring_graphs, objective, connectivity, reverse):
    stack, ring = ring_graphs
    query, target = (ring, stack) if reverse else (stack, ring)
    options = _options(objective, connectivity)
    result = smsd.MCS.find(query).with_options(options).match(target)
    _assert_ring_proof_unavailable(result)


@pytest.mark.parametrize("objective", ["MCIS", "MCES"])
def test_strict_ring_fusion_also_needs_complete_perception(ring_graphs, objective):
    stack, ring = ring_graphs
    options = _options(objective, "FRAGMENTS", complete_rings=False)
    options.chemistry.ring_fusion_mode = smsd.RingFusionMode.STRICT
    result = smsd.MCS.find(ring).with_options(options).match(stack)
    _assert_ring_proof_unavailable(result)


@pytest.mark.parametrize("objective", ["MCIS", "MCES"])
def test_ring_limited_enumeration_does_not_claim_an_optimum(ring_graphs, objective):
    stack, ring = ring_graphs
    result = smsd.MCS.find(ring).with_options(_options(objective)).match_all(stack, limit=2)
    _assert_ring_proof_unavailable(result)


@pytest.mark.parametrize("objective", ["MCIS", "MCES"])
@pytest.mark.parametrize("representation", ["atoms", "match"])
def test_mapping_maximality_requires_ring_proof(ring_graphs, objective, representation):
    stack, ring = ring_graphs
    options = _options(objective)
    mapping = {} if representation == "atoms" else smsd.Match({}, (), ())
    assert smsd.MCS.validate_mapping(stack, ring, mapping, options) == []
    with pytest.raises(RuntimeError):
        smsd.MCS.is_mapping_maximal(stack, ring, mapping, options)


@pytest.mark.parametrize("objective", ["MCIS", "MCES"])
@pytest.mark.parametrize("stop", ["timeout", "cancelled"])
def test_existing_stop_status_takes_precedence(ring_graphs, objective, stop):
    stack, ring = ring_graphs
    options = _options(objective)
    if stop == "timeout":
        options.timeout_ms = 0
        expected = smsd.SearchStatus.TIMED_OUT
    else:
        token = smsd.CancellationToken()
        token.cancel()
        options.cancellation = token
        expected = smsd.SearchStatus.CANCELLED
    result = smsd.MCS.find(stack).with_options(options).match(ring)
    assert result.status == expected
    assert not result.optimal and not result.complete


def test_match_many_keeps_later_targets_after_ring_proof_limit(ring_graphs):
    stack, ring = ring_graphs
    results = smsd.MCS.find(ring).with_options(_options()).match_many([stack, ring], num_threads=1)
    assert len(results) == 2
    _assert_ring_proof_unavailable(results[0])
    assert results[1].status == smsd.SearchStatus.COMPLETED
    assert results[1].optimal and results[1].complete
    assert results[1].best.atom_count == 4


def test_consensus_does_not_hide_truncated_input_perception(ring_graphs):
    # With the ring beside it the consensus fragment stays writable: the cage alone needs more
    # simultaneously open ring closures than SMILES output allows.
    stack, ring = ring_graphs
    result = smsd.MCS.consensus([stack, ring], _options())
    _assert_ring_proof_unavailable(result)


@pytest.mark.parametrize("objective", ["MCIS", "MCES"])
def test_complete_ring_search_with_full_perception_can_prove_result(ring_graphs, objective):
    _, ring = ring_graphs
    result = smsd.MCS.find(ring).with_options(_options(objective)).match(ring)
    assert result.status == smsd.SearchStatus.COMPLETED
    assert result.optimal and result.complete
    assert result.best.atom_count == 4 and result.best.bond_count == 4


@pytest.mark.parametrize("fusion_mode", ["IGNORE", "PERMISSIVE"])
def test_truncated_cycles_do_not_limit_unconstrained_search(ring_graphs, fusion_mode):
    stack, ring = ring_graphs
    options = _options(complete_rings=False)
    options.chemistry.ring_fusion_mode = getattr(smsd.RingFusionMode, fusion_mode)
    result = smsd.MCS.find(ring).with_options(options).match(stack)
    assert result.status == smsd.SearchStatus.COMPLETED
    assert result.optimal and result.complete
    assert result.best.atom_count == 4


@pytest.fixture(scope="module")
def relevant_cycle_graphs():
    """Four-member rings linked around one larger ring."""
    ring_count = 12
    atom_count = 1 + 3 * ring_count
    neighbors = [[] for _ in range(atom_count)]

    def add_bond(first, second):
        neighbors[first].append(second)
        neighbors[second].append(first)

    previous = 0
    for ring in range(ring_count):
        first_side = 3 * ring + 1
        second_side = first_side + 1
        following = first_side + 2
        add_bond(previous, first_side)
        add_bond(previous, second_side)
        add_bond(first_side, following)
        add_bond(second_side, following)
        previous = following
    add_bond(previous, 0)
    graph = (smsd.MolGraphBuilder().atom_count(atom_count)
             .atomic_numbers([6] * atom_count).neighbors(neighbors)
             .bond_orders([[1] * len(row) for row in neighbors]).build(False))
    assert not smsd.rings_truncated(graph)
    return graph, smsd.parse_smiles("C")


@pytest.mark.parametrize("objective", ["MCIS", "MCES"])
@pytest.mark.parametrize("connectivity", ["CONNECTED", "FRAGMENTS"])
@pytest.mark.parametrize("reverse", [False, True])
def test_relevant_cycle_cut_keeps_results_unproved(
        relevant_cycle_graphs, objective, connectivity, reverse):
    chain, carbon = relevant_cycle_graphs
    query, target = (carbon, chain) if reverse else (chain, carbon)
    options = _options(objective, connectivity)
    result = smsd.MCS.find(query).with_options(options).match(target)
    _assert_ring_proof_unavailable(result)


@pytest.mark.parametrize("objective", ["MCIS", "MCES"])
@pytest.mark.parametrize("representation", ["atoms", "match"])
def test_relevant_cycle_cut_blocks_mapping_maximality(
        relevant_cycle_graphs, objective, representation):
    chain, carbon = relevant_cycle_graphs
    options = _options(objective)
    mapping = {} if representation == "atoms" else smsd.Match({}, (), ())
    assert smsd.MCS.validate_mapping(chain, carbon, mapping, options) == []
    with pytest.raises(RuntimeError):
        smsd.MCS.is_mapping_maximal(chain, carbon, mapping, options)


@pytest.mark.parametrize("objective", ["MCIS", "MCES"])
def test_relevant_cycle_cut_keeps_consensus_unproved(relevant_cycle_graphs, objective):
    chain, _ = relevant_cycle_graphs
    result = smsd.MCS.consensus([chain], _options(objective))
    _assert_ring_proof_unavailable(result)
    assert result.found and result.graph.n == chain.n


@pytest.mark.parametrize("objective", ["MCIS", "MCES"])
@pytest.mark.parametrize("fusion_mode", ["IGNORE", "STRICT"])
def test_relevant_cycle_cut_does_not_limit_other_ring_rules(
        relevant_cycle_graphs, objective, fusion_mode):
    chain, carbon = relevant_cycle_graphs
    options = _options(objective, complete_rings=False)
    options.chemistry.ring_fusion_mode = getattr(smsd.RingFusionMode, fusion_mode)
    result = smsd.MCS.find(chain).with_options(options).match(carbon)
    assert result.status == smsd.SearchStatus.COMPLETED
    assert result.optimal and result.complete
