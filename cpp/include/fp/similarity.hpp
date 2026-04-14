/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * Fingerprint similarity metrics.
 * Binary: Tanimoto, Dice, Cosine, Soergel.
 * Count-vector: Tanimoto, Dice, Cosine.
 * Subset check for substructure pre-screening.
 */
#pragma once
#ifndef FP_SIMILARITY_HPP
#define FP_SIMILARITY_HPP

#include "smsd/bitops.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace smsd {
namespace fp {

// ── Binary fingerprint metrics ───────────────────────────────────────────────

/// Tanimoto coefficient: |A&B| / |A|B|.
inline double fingerprintTanimoto(
    const std::vector<uint64_t>& a,
    const std::vector<uint64_t>& b)
{
    if (a.empty() || b.empty()) return 0.0;
    int minWords = static_cast<int>(std::min(a.size(), b.size()));
    int maxWords = static_cast<int>(std::max(a.size(), b.size()));
    int andBits = 0, orBits = 0;
    for (int i = 0; i < minWords; ++i) {
        andBits += popcount64(a[i] & b[i]);
        orBits  += popcount64(a[i] | b[i]);
    }
    const auto& longer = (a.size() > b.size()) ? a : b;
    for (int i = minWords; i < maxWords; ++i)
        orBits += popcount64(longer[i]);
    return (orBits == 0) ? 1.0 : static_cast<double>(andBits) / orBits;
}

/// Dice coefficient: 2*|A&B| / (|A| + |B|).
inline double fingerprintDice(
    const std::vector<uint64_t>& a,
    const std::vector<uint64_t>& b)
{
    if (a.empty() || b.empty()) return 0.0;
    int minWords = static_cast<int>(std::min(a.size(), b.size()));
    int maxWords = static_cast<int>(std::max(a.size(), b.size()));
    int andBits = 0, aBits = 0, bBits = 0;
    for (int i = 0; i < minWords; ++i) {
        andBits += popcount64(a[i] & b[i]);
        aBits   += popcount64(a[i]);
        bBits   += popcount64(b[i]);
    }
    const auto& longer = (a.size() > b.size()) ? a : b;
    for (int i = minWords; i < maxWords; ++i) {
        if (&longer == &a) aBits += popcount64(longer[i]);
        else               bBits += popcount64(longer[i]);
    }
    int sum = aBits + bBits;
    return (sum == 0) ? 0.0 : (2.0 * andBits) / sum;
}

/// Cosine similarity: |A&B| / sqrt(|A| * |B|).
inline double fingerprintCosine(
    const std::vector<uint64_t>& a,
    const std::vector<uint64_t>& b)
{
    if (a.empty() || b.empty()) return 0.0;
    int minWords = static_cast<int>(std::min(a.size(), b.size()));
    int maxWords = static_cast<int>(std::max(a.size(), b.size()));
    int andBits = 0, aBits = 0, bBits = 0;
    for (int i = 0; i < minWords; ++i) {
        andBits += popcount64(a[i] & b[i]);
        aBits   += popcount64(a[i]);
        bBits   += popcount64(b[i]);
    }
    const auto& longer = (a.size() > b.size()) ? a : b;
    for (int i = minWords; i < maxWords; ++i) {
        if (&longer == &a) aBits += popcount64(longer[i]);
        else               bBits += popcount64(longer[i]);
    }
    double denom = std::sqrt(static_cast<double>(aBits) * bBits);
    return (denom == 0.0) ? 0.0 : static_cast<double>(andBits) / denom;
}

/// Soergel distance: 1 - Tanimoto. A proper metric for clustering/k-NN.
inline double fingerprintSoergel(
    const std::vector<uint64_t>& a,
    const std::vector<uint64_t>& b)
{
    if (a.empty() || b.empty()) return 1.0;
    return 1.0 - fingerprintTanimoto(a, b);
}

// ── Count-vector similarity metrics ──────────────────────────────────────────

/// Count-vector Tanimoto: sum(min(a,b)) / sum(max(a,b)).
inline double countTanimoto(
    const std::vector<int>& a,
    const std::vector<int>& b)
{
    if (a.empty() || b.empty() || a.size() != b.size()) return 0.0;
    long long minSum = 0, maxSum = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        minSum += std::min(a[i], b[i]);
        maxSum += std::max(a[i], b[i]);
    }
    return (maxSum == 0) ? 0.0 : static_cast<double>(minSum) / maxSum;
}

/// Count-vector Dice: 2*sum(min(a,b)) / (sum(a) + sum(b)).
inline double countDice(
    const std::vector<int>& a,
    const std::vector<int>& b)
{
    if (a.empty() || b.empty() || a.size() != b.size()) return 0.0;
    long long minSum = 0, aSum = 0, bSum = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        minSum += std::min(a[i], b[i]);
        aSum   += a[i];
        bSum   += b[i];
    }
    long long denom = aSum + bSum;
    return (denom == 0) ? 0.0 : (2.0 * minSum) / denom;
}

/// Count-vector Cosine: dot(a,b) / (|a| * |b|).
inline double countCosine(
    const std::vector<int>& a,
    const std::vector<int>& b)
{
    if (a.empty() || b.empty() || a.size() != b.size()) return 0.0;
    long long dot = 0, aSqSum = 0, bSqSum = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        dot    += static_cast<long long>(a[i]) * b[i];
        aSqSum += static_cast<long long>(a[i]) * a[i];
        bSqSum += static_cast<long long>(b[i]) * b[i];
    }
    double denom = std::sqrt(static_cast<double>(aSqSum) * bSqSum);
    return (denom == 0.0) ? 0.0 : static_cast<double>(dot) / denom;
}

// ── Subset check ─────────────────────────────────────────────────────────────

/// Check if fingerprint a is a subset of b (all ON bits of a are ON in b).
/// Used for substructure pre-screening.
inline bool fingerprintSubset(
    const std::vector<uint64_t>& query,
    const std::vector<uint64_t>& target)
{
    if (query.empty()) return true;
    int qWords = static_cast<int>(query.size());
    int tWords = static_cast<int>(target.size());
    for (int i = 0; i < qWords; ++i) {
        uint64_t qw = query[i];
        uint64_t tw = (i < tWords) ? target[i] : 0ULL;
        if ((qw & tw) != qw) return false;
    }
    return true;
}

} // namespace fp
} // namespace smsd

#endif // FP_SIMILARITY_HPP
