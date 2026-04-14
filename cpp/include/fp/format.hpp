/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * Fingerprint format conversions: hex, binary string, counts-to-array.
 */
#pragma once
#ifndef FP_FORMAT_HPP
#define FP_FORMAT_HPP

#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace smsd {
namespace fp {

/// Convert a fingerprint to a hexadecimal string (lowercase, zero-padded).
inline std::string toHex(const std::vector<uint64_t>& fp) {
    std::string out;
    out.reserve(fp.size() * 16);
    for (uint64_t w : fp) {
        char buf[17];
        std::snprintf(buf, sizeof(buf), "%016" PRIx64, w);
        out.append(buf, 16);
    }
    return out;
}

/// Parse a hexadecimal string back into a fingerprint.
inline std::vector<uint64_t> fromHex(const std::string& hex) {
    int words = (static_cast<int>(hex.size()) + 15) / 16;
    std::vector<uint64_t> fp(words, 0ULL);
    for (int i = 0; i < words; ++i) {
        int start = i * 16;
        int len = std::min(16, static_cast<int>(hex.size()) - start);
        std::string chunk = hex.substr(start, len);
        fp[i] = std::strtoull(chunk.c_str(), nullptr, 16);
    }
    return fp;
}

/// Convert a fingerprint to a binary string of '0' and '1' characters.
inline std::string toBinaryString(const std::vector<uint64_t>& fp, int fpSize) {
    std::string out(fpSize, '0');
    for (int i = 0; i < fpSize; ++i) {
        if (fp[i / 64] & (1ULL << (i % 64))) out[i] = '1';
    }
    return out;
}

/// Convert a sparse count map to a dense integer array of length fpSize.
inline std::vector<int> countsToArray(const std::map<int,int>& counts, int fpSize) {
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    std::vector<int> arr(fpSize, 0);
    for (auto& kv : counts) {
        if (kv.first >= 0 && kv.first < fpSize)
            arr[kv.first] = kv.second;
    }
    return arr;
}

} // namespace fp
} // namespace smsd

#endif // FP_FORMAT_HPP
