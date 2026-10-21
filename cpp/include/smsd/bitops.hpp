/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Portable 64-bit bit-operation wrappers with fused multi-word operations.
 * Works on GCC, Clang, and MSVC without any external dependency.
 * ARM NEON optimised paths when available.
 */
#pragma once
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#if defined(_MSC_VER)
#  include <intrin.h>
#  pragma intrinsic(_BitScanForward64, __popcnt64)
#endif

#if defined(__ARM_NEON)
#  include <arm_neon.h>
#endif

namespace smsd {

// =========================================================================
// Scalar primitives
// =========================================================================

/// Population count (number of set bits) in a 64-bit word.
inline int popcount64(uint64_t x) noexcept {
#if defined(_MSC_VER)
    return static_cast<int>(__popcnt64(x));
#else
    return __builtin_popcountll(x);
#endif
}

/// Index of the least-significant set bit.  Returns 64 if x == 0
/// (safe fallback -- __builtin_ctzll(0) is UB).
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

// =========================================================================
// Fused multi-word bitset operations (scalar fallbacks)
// =========================================================================

/// Popcount of (a[w] & b[w]) summed over all words.
inline int popcountAnd_scalar(const uint64_t* a, const uint64_t* b, int nWords) noexcept {
    int count = 0;
    for (int w = 0; w < nWords; ++w)
        count += popcount64(a[w] & b[w]);
    return count;
}

/// Popcount of a single bitset over nWords.
inline int popcountN_scalar(const uint64_t* a, int nWords) noexcept {
    int count = 0;
    for (int w = 0; w < nWords; ++w)
        count += popcount64(a[w]);
    return count;
}

/// True if any bit is set in (a & b).  Early exit on first hit.
inline bool anyBitAnd_scalar(const uint64_t* a, const uint64_t* b, int nWords) noexcept {
    for (int w = 0; w < nWords; ++w)
        if (a[w] & b[w]) return true;
    return false;
}

/// Popcount of (a[w] & ~b[w]) summed over all words (bits in a not in b).
inline int popcountAndNot_scalar(const uint64_t* a, const uint64_t* b, int nWords) noexcept {
    int count = 0;
    for (int w = 0; w < nWords; ++w)
        count += popcount64(a[w] & ~b[w]);
    return count;
}

// =========================================================================
// ARM NEON optimised variants
// =========================================================================
#if defined(__ARM_NEON)

/// NEON byte-level popcount of a 128-bit register, returned as scalar sum.
inline int neon_popcnt128(uint8x16_t v) noexcept {
    // vcntq_u8 gives per-byte popcount; then horizontal sum via vaddlvq_u8
    uint8x16_t cnt = vcntq_u8(v);
    return static_cast<int>(vaddlvq_u8(cnt));
}

inline int popcountAnd_neon(const uint64_t* a, const uint64_t* b, int nWords) noexcept {
    int count = 0;
    int w = 0;
    // Process 2 uint64_t (128 bits) per iteration
    for (; w + 1 < nWords; w += 2) {
        uint64x2_t va = vld1q_u64(a + w);
        uint64x2_t vb = vld1q_u64(b + w);
        uint64x2_t vand = vandq_u64(va, vb);
        count += neon_popcnt128(vreinterpretq_u8_u64(vand));
    }
    // Handle trailing word
    if (w < nWords)
        count += popcount64(a[w] & b[w]);
    return count;
}

inline int popcountN_neon(const uint64_t* a, int nWords) noexcept {
    int count = 0;
    int w = 0;
    for (; w + 1 < nWords; w += 2) {
        uint64x2_t va = vld1q_u64(a + w);
        count += neon_popcnt128(vreinterpretq_u8_u64(va));
    }
    if (w < nWords)
        count += popcount64(a[w]);
    return count;
}

inline bool anyBitAnd_neon(const uint64_t* a, const uint64_t* b, int nWords) noexcept {
    int w = 0;
    for (; w + 1 < nWords; w += 2) {
        uint64x2_t va = vld1q_u64(a + w);
        uint64x2_t vb = vld1q_u64(b + w);
        uint64x2_t vand = vandq_u64(va, vb);
        // Check if any lane is nonzero
        if (vgetq_lane_u64(vand, 0) | vgetq_lane_u64(vand, 1)) return true;
    }
    if (w < nWords && (a[w] & b[w])) return true;
    return false;
}

inline int popcountAndNot_neon(const uint64_t* a, const uint64_t* b, int nWords) noexcept {
    int count = 0;
    int w = 0;
    for (; w + 1 < nWords; w += 2) {
        uint64x2_t va = vld1q_u64(a + w);
        uint64x2_t vb = vld1q_u64(b + w);
        // NEON: BIC = AND-NOT (a & ~b)
        uint64x2_t vbic = vbicq_u64(va, vb);
        count += neon_popcnt128(vreinterpretq_u8_u64(vbic));
    }
    if (w < nWords)
        count += popcount64(a[w] & ~b[w]);
    return count;
}

#endif // __ARM_NEON

// =========================================================================
// Dispatch wrappers -- pick NEON when available, else scalar
// =========================================================================

/// @deprecated Nothing in the library uses the fused operations; they go with their
/// _scalar and _neon variants.
[[deprecated("The fused bitset operations are unused.")]]
inline int popcountAnd(const uint64_t* a, const uint64_t* b, int nWords) noexcept {
#if defined(__ARM_NEON)
    return popcountAnd_neon(a, b, nWords);
#else
    return popcountAnd_scalar(a, b, nWords);
#endif
}

[[deprecated("The fused bitset operations are unused.")]]
inline int popcountN(const uint64_t* a, int nWords) noexcept {
#if defined(__ARM_NEON)
    return popcountN_neon(a, nWords);
#else
    return popcountN_scalar(a, nWords);
#endif
}

[[deprecated("The fused bitset operations are unused.")]]
inline bool anyBitAnd(const uint64_t* a, const uint64_t* b, int nWords) noexcept {
#if defined(__ARM_NEON)
    return anyBitAnd_neon(a, b, nWords);
#else
    return anyBitAnd_scalar(a, b, nWords);
#endif
}

[[deprecated("The fused bitset operations are unused.")]]
inline int popcountAndNot(const uint64_t* a, const uint64_t* b, int nWords) noexcept {
#if defined(__ARM_NEON)
    return popcountAndNot_neon(a, b, nWords);
#else
    return popcountAndNot_scalar(a, b, nWords);
#endif
}

// =========================================================================
// Bit positions
// =========================================================================

/// Calls visit(position) for each set bit below `limit`, in ascending order.
template<class Visit>
inline void forEachSetBit(const std::vector<uint64_t>& words, int64_t limit, Visit&& visit) {
    for (size_t w = 0; w < words.size(); ++w)
        for (uint64_t word = words[w]; word != 0; word &= word - 1) {
            const int64_t position = int64_t(w) * 64 + ctz64(word);
            if (position >= limit) return;
            visit(static_cast<int>(position));
        }
}

/// Ascending positions of the set bits below `limit`.
inline std::vector<int> bitPositions(const std::vector<uint64_t>& words,
                                     int64_t limit = std::numeric_limits<int>::max()) {
    std::vector<int> positions;
    forEachSetBit(words, limit, [&](int position) { positions.push_back(position); });
    return positions;
}

/// Words of a `size`-bit set holding the given positions; positions outside [0, size) are ignored.
inline std::vector<uint64_t> bitWords(const std::vector<int>& positions, int size) {
    std::vector<uint64_t> words(size > 0 ? (size_t(size) - 1) / 64 + 1 : 0, 0);
    for (int position : positions)
        if (position >= 0 && position < size) words[size_t(position) >> 6] |= uint64_t(1) << (position & 63);
    return words;
}

} // namespace smsd
