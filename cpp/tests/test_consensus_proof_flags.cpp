/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/mcs.hpp"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void checkApproximate(const smsd::MCS::ConsensusResult& result) {
    check(!result.optimal, "progressive consensus must remain approximate");
    check(result.complete == (result.status == smsd::SearchStatus::Completed),
          "complete must retain its progressive-pass meaning");
}
}

int main() {
    try {
        for (const auto objective : {smsd::MCSObjective::MCIS, smsd::MCSObjective::MCES}) {
            smsd::MCSOptions options;
            options.objective = objective;
            options.timeoutMs = 10000;
            const auto empty = smsd::MCS::consensus({}, options);
            checkApproximate(empty);
            check(empty.status == smsd::SearchStatus::Completed && !empty.found
                  && empty.sourceIndex == -1 && empty.supportCount == 0,
                  "empty consensus must retain its existing result");
            const auto single = smsd::MCS::consensus({smsd::parseSMILES("CCO")}, options);
            checkApproximate(single);
            check(single.status == smsd::SearchStatus::Completed && single.found
                  && single.graph.n == 3 && single.supportCount == 1 && single.sourceIndex == 0,
                  "single-input consensus must retain its fragment and support");
            const auto supported = smsd::MCS::consensus(
                {smsd::parseSMILES("N"), smsd::parseSMILES("CCO"), smsd::parseSMILES("CCN")},
                options, 2.0 / 3.0);
            checkApproximate(supported);
            check(supported.status == smsd::SearchStatus::Completed && supported.found
                  && supported.graph.n == 2 && supported.supportCount == 2,
                  "threshold consensus must retain its fragment and support");
            const auto absent = smsd::MCS::consensus(
                {smsd::parseSMILES("N"), smsd::parseSMILES("O")}, options);
            checkApproximate(absent);
            check(absent.status == smsd::SearchStatus::Completed && !absent.found,
                  "a completed reduction without a fragment is not an absence proof");
            options.timeoutMs = 0;
            const auto timedOut = smsd::MCS::consensus({smsd::parseSMILES("CCO")}, options);
            checkApproximate(timedOut);
            check(timedOut.status == smsd::SearchStatus::TimedOut, "timeout status must be retained");
            options.timeoutMs = 10000;
            options.cancellation = std::make_shared<smsd::CancellationToken>();
            options.cancellation->cancel();
            const auto cancelled = smsd::MCS::consensus({smsd::parseSMILES("CCO")}, options);
            checkApproximate(cancelled);
            check(cancelled.status == smsd::SearchStatus::Cancelled, "cancellation status must be retained");
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
