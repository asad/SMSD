/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * Portable 64-bit bit-operation wrappers.
 * Works on GCC, Clang, and MSVC without any external dependency.
 */
#pragma once
#include <cstdint>

#if defined(_MSC_VER)
#  include <intrin.h>
#  pragma intrinsic(_BitScanForward64, __popcnt64)
#endif

namespace smsd {

/// Population count (number of set bits) in a 64-bit word.
inline int popcount64(uint64_t x) noexcept {
#if defined(_MSC_VER)
    return static_cast<int>(__popcnt64(x));
#else
    return __builtin_popcountll(x);
#endif
}

/// Index of the least-significant set bit.  Returns 64 if x == 0
/// (safe fallback — __builtin_ctzll(0) is UB).
inline int ctz64(uint64_t x) noexcept {
    if (x == 0) return 64;
#if defined(_MSC_VER)
    unsigned long idx;
    _BitScanForward64(&idx, x);
    return static_cast<int>(idx);
#else
    return __builtin_ctzll(x);
#endif
}

} // namespace smsd
