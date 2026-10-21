# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Read-time tautomer context must leave molecules and search objectives unchanged."""
from concurrent.futures import ThreadPoolExecutor

import pytest
import smsd


def test_selected_keto_enol_context_does_not_accumulate():
    query = smsd.parse_smiles("CC(=O)C")
    target = smsd.parse_smiles("CC(O)=C")
    chemistry = smsd.ChemOptions.tautomer_profile()
    options = smsd.MCS.Options()
    options.chemistry = chemistry
    result = smsd.MCS.find(query).with_options(options).match(target)
    assert result.optimal and result.best.atom_count == 4
    baseline = smsd.MCS.compute_tautomer_confidence(query, target, result.best, chemistry)
    chemistry.with_solvent(smsd.Solvent.CHLOROFORM)
    contextual = smsd.MCS.compute_tautomer_confidence(query, target, result.best, chemistry)
    assert 0 < contextual < baseline < 1
    with ThreadPoolExecutor(max_workers=4) as workers:
        values = list(workers.map(lambda _: smsd.MCS.compute_tautomer_confidence(
            query, target, result.best, chemistry), range(8)))
    assert values == [contextual] * 8
    aqueous = smsd.ChemOptions.tautomer_profile()
    assert smsd.MCS.compute_tautomer_confidence(query, target, result.best, aqueous) == baseline
    assert smsd.MCS.find(query).with_options(options).match(target).best == result.best
    identity = smsd.MCS.find(query).match(query).best
    chemistry.with_ph(0)
    assert smsd.MCS.compute_tautomer_confidence(query, query, identity, chemistry) == 1


@pytest.mark.parametrize("pH", [float("nan"), float("inf"), -float("inf")])
def test_tautomer_context_rejects_non_finite_ph(pH):
    with pytest.raises(ValueError):
        smsd.ChemOptions().with_ph(pH)
