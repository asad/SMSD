# SPDX-License-Identifier: Apache-2.0
import pytest
import smsd


@pytest.mark.parametrize("objective", ["MCIS", "MCES"])
@pytest.mark.parametrize("molecules,threshold,atoms,support", [
    ([], 1.0, 0, 0),
    (["CCO"], 1.0, 3, 1),
    (["N", "CCO", "CCN"], 2 / 3, 2, 2),
    (["N", "O"], 1.0, 0, 0),
])
def test_completed_consensus_retains_its_completion_flag(objective, molecules, threshold, atoms, support):
    options = smsd.MCS.Options()
    options.objective = getattr(smsd.MCS.Objective, objective)
    options.timeout_ms = 10_000
    result = smsd.MCS.consensus(molecules, options, threshold)
    assert result.status == smsd.SearchStatus.COMPLETED
    assert not result.optimal and result.complete
    assert result.found == (atoms > 0)
    assert result.graph.n == atoms and result.support_count == support
    if not molecules:
        assert result.source_index == -1
    elif len(molecules) == 1:
        assert result.source_index == 0


@pytest.mark.parametrize("objective", ["MCIS", "MCES"])
@pytest.mark.parametrize("cancelled", [False, True])
def test_stopped_consensus_retains_its_status(objective, cancelled):
    options = smsd.MCS.Options()
    options.objective = getattr(smsd.MCS.Objective, objective)
    options.timeout_ms = 10_000 if cancelled else 0
    if cancelled:
        token = smsd.CancellationToken()
        token.cancel()
        options.cancellation = token
    result = smsd.MCS.consensus(["CCO"], options)
    expected = smsd.SearchStatus.CANCELLED if cancelled else smsd.SearchStatus.TIMED_OUT
    assert result.status == expected
    assert not result.optimal and not result.complete
