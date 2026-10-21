# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
import pytest

import smsd


@pytest.mark.parametrize(
    "connectivity",
    [smsd.MCS.Connectivity.CONNECTED, smsd.MCS.Connectivity.FRAGMENTS],
)
def test_one_atom_omission_is_proved_within_node_budget(connectivity):
    query = smsd.parse_smiles("C1CCCCCCCCCCC1")
    target = smsd.parse_smiles("C1CCCCCCCCCCCCCC1")
    options = smsd.MCS.Options()
    options.connectivity = connectivity
    options.node_limit = 600
    options.timeout_ms = 10000

    result = smsd.MCS.find(query).with_options(options).match(target)

    assert result.status == smsd.SearchStatus.COMPLETED
    assert result.optimal and result.complete
    assert result.nodes_visited <= options.node_limit
    assert result.score == 11
    assert len(result.matches) == 1
    match = result.matches[0]
    assert len(match.atom_mapping) == 11
    assert len(match.bond_mapping) == 10
    assert not smsd.MCS.validate_mapping(query, target, match, options)
