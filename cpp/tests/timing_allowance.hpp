/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#pragma once

#include <cstdint>

#if defined(__has_feature)
#if __has_feature(address_sanitizer) || __has_feature(undefined_behavior_sanitizer)
#define SMSD_TEST_SANITISED 1
#endif
#endif
#if defined(__SANITIZE_ADDRESS__) && !defined(SMSD_TEST_SANITISED)
#define SMSD_TEST_SANITISED 1
#endif

namespace smsd_test {

/**
 * A wall-clock bound in milliseconds for this build. Bounds are set for optimised builds;
 * unoptimised and sanitised builds run several times more slowly and often share a busy
 * runner, so they get eight times the room. A quadratic or stuck search still exceeds it.
 */
constexpr std::int64_t timingAllowance(std::int64_t optimisedMillis) {
#if !defined(SMSD_TEST_SANITISED) && (defined(__OPTIMIZE__) || (defined(_MSC_VER) && !defined(_DEBUG)))
    return optimisedMillis;
#else
    return 8 * optimisedMillis;
#endif
}

}  // namespace smsd_test
