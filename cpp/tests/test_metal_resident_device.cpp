/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/metal_gpu.hpp"
#include "smsd/screen_accelerator.hpp"
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

int main() {
    try {
        if (!smsd::metal::metalAvailable()) return 77;
        const auto* provider = smsd::metal::fingerprintAccelerator();
        smsd_test::require(smsd::screen_accelerator::hasResidentFingerprints(provider),
                           "Metal must provide resident fingerprint operations");
        std::uint64_t minimum = 0;
        char error[256] = {};
        smsd_test::require(provider->minimum_targets(1, &minimum, error, sizeof error) == 0
            && minimum == std::numeric_limits<std::uint64_t>::max(),
            "Automatic resident Metal screening must remain disabled");
        const auto failures = smsd::screen_accelerator::failureCount();
        const auto score = [](int common, int queryBits, int targetBits) {
            const int unionBits = queryBits + targetBits - common;
            return unionBits == 0 ? 0.0 : static_cast<double>(common) / unionBits;
        };
        using Fingerprint = std::vector<std::uint64_t>;
        for (std::size_t words : {1u, 2u, 16u}) {
            const Fingerprint empty(words, 0), full(words, ~std::uint64_t(0));
            const Fingerprint half(words, UINT64_C(0xffffffff));
            Fingerprint decimal(words, 0), other(words, 0), lastBit(words, 0);
            decimal[0] = 0x1f;
            other[0] = 0x3e1;
            lastBit.back() = UINT64_C(1) << 63;
            const std::vector<Fingerprint> queries{decimal, empty, full, decimal, lastBit};
            const std::vector<Fingerprint> targets{other, empty, decimal, other, full, lastBit, half};
            {
                smsd::screen_accelerator::ResidentCache cache;
                for (double threshold : {0.0, 0.1, std::nextafter(0.1, 0.0),
                        std::nextafter(0.1, 1.0), 0.5, 1.0}) {
                    std::vector<std::vector<int>> expected(queries.size());
                    for (std::size_t query = 0; query < queries.size(); ++query)
                        for (std::size_t target = 0; target < targets.size(); ++target) {
                            int common = 0, queryBits = 0, targetBits = 0;
                            for (std::size_t word = 0; word < words; ++word) {
                                common += smsd::popcount64(queries[query][word] & targets[target][word]);
                                queryBits += smsd::popcount64(queries[query][word]);
                                targetBits += smsd::popcount64(targets[target][word]);
                            }
                            if (score(common, queryBits, targetBits) >= threshold)
                                expected[query].push_back(static_cast<int>(target));
                        }
                    const auto actual = smsd::screen_accelerator::fingerprintScreens(
                        queries, targets, threshold, cache, score, provider);
                    smsd_test::require(actual.has_value(), "Physical resident Metal screen must not decline or fall back");
                    smsd_test::require(*actual == expected, "Physical resident Metal hits must equal CPU hits in order");
                }
                const auto repeated = smsd::screen_accelerator::fingerprintScreens(
                    queries, targets, 1.0, cache, score, provider);
                smsd_test::require(repeated.has_value()
                    && *repeated == std::vector<std::vector<int>>{{2}, {}, {4}, {2}, {5}},
                    "Repeated resident queries must retain exact tied hits");
            }
            smsd_test::require(smsd::screen_accelerator::failureCount() == failures,
                               "Resident Metal calls and destruction must not record a failure");
        }
        std::cout << "Resident Metal ordered fingerprint hits match CPU thresholds\n";
        return 0;
    } catch (const std::exception& failure) {
        std::cerr << failure.what() << '\n';
        return 1;
    }
}
