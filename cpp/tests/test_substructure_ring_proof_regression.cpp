/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/substructure.hpp"
#include <cassert>
#include <memory>
#include <vector>

/** Four-member rings joined around a larger ring, with truncated relevant cycles. */
static smsd::MolGraph linkedRings() {
    constexpr int ringCount = 12;
    constexpr int atomCount = 1 + 3 * ringCount;
    std::vector<std::vector<int>> neighbours(atomCount), orders(atomCount);
    const auto addBond = [&](int first, int second) {
        neighbours[first].push_back(second);
        neighbours[second].push_back(first);
        orders[first].push_back(1);
        orders[second].push_back(1);
    };
    int previous = 0;
    for (int ring = 0; ring < ringCount; ++ring) {
        const int firstSide = 3 * ring + 1, secondSide = firstSide + 1;
        const int next = firstSide + 2;
        addBond(previous, firstSide);
        addBond(previous, secondSide);
        addBond(firstSide, next);
        addBond(secondSide, next);
        previous = next;
    }
    addBond(previous, 0);
    return smsd::MolGraph::Builder().atomCount(atomCount)
        .atomicNumbers(std::vector<int>(atomCount, 6)).setNeighbors(std::move(neighbours))
        .setBondOrders(std::move(orders)).build(false);
}

int main() {
    const auto target = linkedRings();
    bool truncated = false;
    target.computeRelevantCyclesReportingTruncation(truncated);
    assert(truncated);
    smsd::SubstructureOptions options;
    options.timeoutMs = 60000;
    options.chemistry.completeRingsOnly = true;
    for (const auto& query : {smsd::Substructure::find("C"), smsd::Substructure::smarts("[#6]")}) {
        const auto result = query.withOptions(options).match(target);
        assert(result.status == smsd::SearchStatus::LimitReached);
        assert(!result.optimal && !result.complete);
        auto stopped = options;
        stopped.timeoutMs = 0;
        const auto timedOut = query.withOptions(stopped).match(target);
        assert(timedOut.status == smsd::SearchStatus::TimedOut && !timedOut.optimal && !timedOut.complete);
        stopped.cancellation = std::make_shared<smsd::CancellationToken>();
        stopped.cancellation->cancel();
        const auto cancelled = query.withOptions(stopped).match(target);
        assert(cancelled.status == smsd::SearchStatus::Cancelled && !cancelled.optimal && !cancelled.complete);
    }
    const auto exactNegative = smsd::Substructure::find("N").withOptions(options).match(target);
    assert(exactNegative.matches.empty() && exactNegative.optimal && exactNegative.complete);
    const auto exactPatternNegative = smsd::Substructure::smarts("[#7]").withOptions(options).match(target);
    assert(exactPatternNegative.matches.empty() && exactPatternNegative.optimal && exactPatternNegative.complete);
    const auto ordinary = smsd::Substructure::find("C1CCC1").withOptions(options).match("C1CCC1");
    assert(!ordinary.matches.empty() && ordinary.optimal && ordinary.complete);
}
