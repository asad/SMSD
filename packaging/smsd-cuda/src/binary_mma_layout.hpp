/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 */
#ifndef SMSD_CUDA_BINARY_MMA_LAYOUT_HPP
#define SMSD_CUDA_BINARY_MMA_LAYOUT_HPP

#include <cstdint>

#ifdef __CUDACC__
#define SMSD_MMA_HD __host__ __device__
#else
#define SMSD_MMA_HD
#endif

namespace smsd_cuda_binary_mma {

/** Register coordinates from PTX's mma.m16n8k256 fragment diagrams. */
SMSD_MMA_HD constexpr int queryRow(int lane, int fragment) {
    return (lane >> 2) + (fragment & 1) * 8;
}
SMSD_MMA_HD constexpr int queryWord(int lane, int fragment) {
    return ((lane & 3) >> 1) + (fragment >> 1) * 2;
}
SMSD_MMA_HD constexpr int targetColumn(int lane) { return lane >> 2; }
SMSD_MMA_HD constexpr int targetWord(int lane, int fragment) {
    return ((lane & 3) >> 1) + fragment * 2;
}
SMSD_MMA_HD constexpr int wordShift(int lane) { return (lane & 1) * 32; }
SMSD_MMA_HD constexpr int resultRow(int lane, int fragment) {
    return (lane >> 2) + (fragment >> 1) * 8;
}
SMSD_MMA_HD constexpr int resultColumn(int lane, int fragment) {
    return (lane & 3) * 2 + (fragment & 1);
}

/** Zero padding covers incomplete query, target and 256-bit tiles. */
SMSD_MMA_HD inline uint32_t queryFragment(const uint64_t* packed, int queries, int words,
                                         int64_t queryBegin, int wordBegin, int lane, int fragment) {
    const int64_t query = queryBegin + queryRow(lane, fragment);
    const int word = wordBegin + queryWord(lane, fragment);
    return query < queries && word < words
        ? static_cast<uint32_t>(packed[query * words + word] >> wordShift(lane)) : 0;
}
SMSD_MMA_HD inline uint32_t targetFragment(const uint64_t* packed, int targets, int words,
                                          int64_t targetBegin, int wordBegin, int lane, int fragment) {
    const int64_t target = targetBegin + targetColumn(lane);
    const int word = wordBegin + targetWord(lane, fragment);
    return target < targets && word < words
        ? static_cast<uint32_t>(packed[static_cast<int64_t>(word) * targets + target] >> wordShift(lane)) : 0;
}

} // namespace smsd_cuda_binary_mma

#undef SMSD_MMA_HD
#endif
