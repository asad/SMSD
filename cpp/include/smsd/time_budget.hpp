/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#pragma once

#include <chrono>
#include <cstdint>

namespace smsd {
namespace detail {

// Avoid overflow in both the milliseconds-to-clock-duration conversion and
// the addition to the current time when callers request a very large timeout.
inline std::chrono::steady_clock::time_point steadyDeadline(
    int64_t milliseconds,
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now()) {
    using Clock = std::chrono::steady_clock;
    if (milliseconds <= 0) return start;
    const auto maxMilliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
        Clock::duration::max()).count();
    if (milliseconds >= maxMilliseconds) return Clock::time_point::max();
    const auto duration = std::chrono::duration_cast<Clock::duration>(
        std::chrono::milliseconds(milliseconds));
    if (start >= Clock::time_point::max() - duration) return Clock::time_point::max();
    return start + duration;
}

} // namespace detail
} // namespace smsd
