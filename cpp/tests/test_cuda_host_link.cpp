/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#include "smsd/gpu.hpp"
#include <iostream>

int main() {
    const int labels[] = {6 << 2};
    const uint64_t bits[] = {1};
    const auto query = smsd::cuda::toGpuMolecule(1, labels);
    const auto fingerprint = smsd::cuda::toGpuFingerprint(bits, 1);
    if (!smsd::cuda::batchRascalScreen(query, {}, 0.0).empty()
        || !smsd::cuda::batchFingerprintScreen(fingerprint, {}).empty()) return 1;
    std::cout << smsd::cuda::cudaDeviceInfo() << '\n';
    std::cout << "CUDA device available: " << smsd::cuda::cudaAvailable() << '\n';
}
