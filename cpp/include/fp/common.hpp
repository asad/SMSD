/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * Common fingerprint constants and hash utilities shared by all FP modules.
 */
#pragma once
#ifndef FP_COMMON_HPP
#define FP_COMMON_HPP

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace smsd {
namespace fp {

/// FNV-1a 64-bit offset basis.
constexpr uint64_t FNV1A_SEED  = 0xCBF29CE484222325ULL;
/// FNV-1a 64-bit prime.
constexpr uint64_t FNV1A_PRIME = 0x100000001B3ULL;

/// FNV-1a mixing step: fold one 64-bit value into hash state.
inline uint64_t fnvMix(uint64_t h, uint64_t val) {
    h ^= val;
    h *= FNV1A_PRIME;
    return h;
}

/// FNV-1a hash over an array of ints.
inline uint64_t fnv1a(const int* data, int len) {
    uint64_t h = FNV1A_SEED;
    for (int i = 0; i < len; ++i) {
        h ^= static_cast<uint64_t>(static_cast<uint32_t>(data[i]));
        h *= FNV1A_PRIME;
    }
    return h;
}

/// Allocate a zero-initialised bit-vector with fpSize bits.
inline std::vector<uint64_t> makeBitVector(int fpSize) {
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    int numWords = (fpSize + 63) / 64;
    return std::vector<uint64_t>(numWords, 0ULL);
}

/// Set a single bit in a fingerprint bit-vector.
inline void setBit(std::vector<uint64_t>& fp, int fpSize, uint64_t hash) {
    int bit = static_cast<int>(hash % static_cast<uint64_t>(fpSize));
    fp[bit / 64] |= (1ULL << (bit % 64));
}

} // namespace fp
} // namespace smsd

#endif // FP_COMMON_HPP
