/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 */
#include "../src/binary_mma_layout.hpp"

#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void check(bool condition) {
    if (!condition) throw std::runtime_error("binary MMA fragment layout differs from the integer oracle");
}
int countBits(uint64_t value) {
    int result = 0;
    while (value != 0) { value &= value - 1; ++result; }
    return result;
}
uint64_t wordValue(int row, int word, int pattern) {
    if (pattern == 0) return 0;
    if (pattern == 1) return ~uint64_t(0);
    uint64_t value = uint64_t(row + 1) * UINT64_C(0x9e3779b97f4a7c15)
                   ^ uint64_t(word + 1) * UINT64_C(0xbf58476d1ce4e5b9);
    value ^= value >> 30;
    value *= UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 27);
}

void checkCoordinates() {
    using namespace smsd_cuda_binary_mma;
    std::array<int, 16 * 256> queryUses{};
    std::array<int, 8 * 256> targetUses{};
    std::array<int, 16 * 8> resultUses{};
    for (int lane = 0; lane < 32; ++lane) {
        for (int fragment = 0; fragment < 4; ++fragment) {
            const int query = queryRow(lane, fragment);
            const int firstBit = queryWord(lane, fragment) * 64 + wordShift(lane);
            check(query >= 0 && query < 16 && firstBit >= 0 && firstBit + 32 <= 256);
            for (int bit = 0; bit < 32; ++bit) ++queryUses[query * 256 + firstBit + bit];
            const int row = resultRow(lane, fragment), column = resultColumn(lane, fragment);
            check(row == lane / 4 + (fragment / 2) * 8);
            check(column == (lane % 4) * 2 + fragment % 2);
            check(row >= 0 && row < 16 && column >= 0 && column < 8);
            ++resultUses[row * 8 + column];
        }
        for (int fragment = 0; fragment < 2; ++fragment) {
            const int target = targetColumn(lane);
            const int firstBit = targetWord(lane, fragment) * 64 + wordShift(lane);
            check(target >= 0 && target < 8 && firstBit >= 0 && firstBit + 32 <= 256);
            for (int bit = 0; bit < 32; ++bit) ++targetUses[target * 256 + firstBit + bit];
        }
    }
    for (int uses : queryUses) check(uses == 1);
    for (int uses : targetUses) check(uses == 1);
    for (int uses : resultUses) check(uses == 1);
}

void checkIntersections(int queries, int targets, int words, int pattern) {
    using namespace smsd_cuda_binary_mma;
    std::vector<uint64_t> queryWords(static_cast<size_t>(queries) * words);
    std::vector<uint64_t> targetWords(static_cast<size_t>(targets) * words);
    for (int query = 0; query < queries; ++query)
        for (int word = 0; word < words; ++word)
            queryWords[static_cast<size_t>(query) * words + word] = wordValue(query, word, pattern);
    for (int target = 0; target < targets; ++target)
        for (int word = 0; word < words; ++word)
            targetWords[static_cast<size_t>(word) * targets + target] = wordValue(target + 3, word, pattern);
    for (int queryBegin = 0; queryBegin < queries; queryBegin += 16) {
        for (int targetBegin = 0; targetBegin < targets; targetBegin += 8) {
            std::array<int, 16 * 8> actual{};
            for (int wordBegin = 0; wordBegin < words; wordBegin += 4) {
                std::array<unsigned char, 16 * 256> matrixA{};
                std::array<unsigned char, 8 * 256> matrixB{};
                for (int lane = 0; lane < 32; ++lane) {
                    // PTX Figure 99: registers 0/2 carry the first row, 1/3 the row eight below.
                    for (int fragment = 0; fragment < 4; ++fragment) {
                        const uint32_t packed = queryFragment(queryWords.data(), queries, words,
                                                             queryBegin, wordBegin, lane, fragment);
                        const int row = lane / 4 + (fragment % 2) * 8;
                        const int firstBit = (lane % 4) * 32 + (fragment / 2) * 128;
                        for (int bit = 0; bit < 32; ++bit)
                            matrixA[row * 256 + firstBit + bit] = (packed >> bit) & 1;
                    }
                    for (int fragment = 0; fragment < 2; ++fragment) {
                        const uint32_t packed = targetFragment(targetWords.data(), targets, words,
                                                              targetBegin, wordBegin, lane, fragment);
                        const int column = lane / 4;
                        const int firstBit = (lane % 4) * 32 + fragment * 128;
                        for (int bit = 0; bit < 32; ++bit)
                            matrixB[column * 256 + firstBit + bit] = (packed >> bit) & 1;
                    }
                }
                for (int lane = 0; lane < 32; ++lane)
                    for (int fragment = 0; fragment < 4; ++fragment) {
                        const int row = lane / 4 + (fragment / 2) * 8;
                        const int column = (lane % 4) * 2 + fragment % 2;
                        for (int bit = 0; bit < 256; ++bit)
                            actual[row * 8 + column] += matrixA[row * 256 + bit] & matrixB[column * 256 + bit];
                    }
            }
            for (int query = 0; query < 16; ++query)
                for (int target = 0; target < 8; ++target) {
                    int expected = 0;
                    if (queryBegin + query < queries && targetBegin + target < targets)
                        for (int word = 0; word < words; ++word)
                            expected += countBits(queryWords[static_cast<size_t>(queryBegin + query) * words + word]
                                & targetWords[static_cast<size_t>(word) * targets + targetBegin + target]);
                    check(actual[query * 8 + target] == expected);
                }
        }
    }
}
} // namespace

int main() {
    try {
        checkCoordinates();
        for (int words : {1, 2, 3, 4, 5, 9, 33})
            for (int pattern : {0, 1, 2}) {
                checkIntersections(1, 1, words, pattern);
                checkIntersections(15, 7, words, pattern);
                checkIntersections(16, 8, words, pattern);
                checkIntersections(17, 9, words, pattern);
            }
        std::cout << "binary MMA CPU layout oracle passed\n";
        return 0;
    } catch (const std::exception& failure) {
        std::cerr << failure.what() << '\n';
        return 1;
    }
}
