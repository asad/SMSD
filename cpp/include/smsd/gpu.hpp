/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * GPU runtime query and auto-dispatch helpers for SMSD batch operations.
 *
 * Usage:
 *   #include "smsd/gpu.hpp"
 *
 *   if (smsd::gpu::isAvailable()) {
 *       std::cout << smsd::gpu::deviceInfo() << "\n";
 *       // batchSimilarityScreenGpu is available
 *   }
 *
 * When SMSD is built with SMSD_BUILD_CUDA=ON/AUTO and nvcc is found,
 * SMSD_ENABLE_CUDA is defined and this header bridges the smsd_cuda library.
 *
 * When SMSD is built on macOS with SMSD_BUILD_METAL=ON/AUTO and
 * Metal.framework is found, SMSD_ENABLE_METAL is defined and this header
 * bridges the smsd_metal library (metal/metal_batch.mm).
 *
 * Without either, isAvailable() returns false and all helpers delegate to
 * the CPU/OpenMP path.
 *
 * Dispatch priority:  CUDA  >  Metal  >  CPU/OpenMP
 *
 * Link targets:
 *   CUDA:   target_link_libraries(my_target PRIVATE smsd smsd_cuda)
 *   Metal:  target_link_libraries(my_target PRIVATE smsd smsd_metal)
 */
#pragma once
#ifndef SMSD_GPU_HPP
#define SMSD_GPU_HPP

#include "smsd/batch.hpp"
#include "smsd/similarity_screen.hpp"
#include "smsd/gpu_kernels.hpp"

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// Metal backend (pure-C++ declarations; implementation in metal_batch.mm)
#ifdef SMSD_ENABLE_METAL
#include "smsd/metal_gpu.hpp"
#endif

// Forward-declare CUDA host types when the CUDA library is compiled in.
// The full definitions live in cuda/batch_screen.cu; we only need the
// lightweight host-side structs and function signatures here.
#ifdef SMSD_ENABLE_CUDA
namespace smsd { namespace cuda {

static constexpr int LABEL_BINS = 512;
static constexpr int FP_WORDS   = 16;

struct alignas(16) GpuMolecule {
    int n;
    int labelCounts[LABEL_BINS];
};

struct alignas(16) GpuFingerprint {
    uint64_t bits[FP_WORDS];
};

struct ScreenResult {
    int    index;
    double score;
};

bool        cudaAvailable();
std::string cudaDeviceInfo();

GpuMolecule   toGpuMolecule(int atomCount, const int* labels);
GpuFingerprint toGpuFingerprint(const uint64_t* fp, int numWords);

std::vector<ScreenResult> batchSimilarityScreen(
    const GpuMolecule& query,
    const std::vector<GpuMolecule>& targets,
    double threshold);

std::vector<ScreenResult> batchFingerprintScreen(
    const GpuFingerprint& query,
    const std::vector<GpuFingerprint>& targets);

}} // smsd::cuda
#endif // SMSD_ENABLE_CUDA


namespace smsd {
namespace gpu {

namespace detail {
inline std::string openMPVersion(long date) {
    switch (date) {
        case 200203: return "2.0";
        case 200505: return "2.5";
        case 200805: return "3.0";
        case 201107: return "3.1";
        case 201307: return "4.0";
        case 201511: return "4.5";
        case 201811: return "5.0";
        case 202011: return "5.1";
        case 202111: return "5.2";
        case 202411: return "6.0";
        default: return std::to_string(date);
    }
}

/** Pack the atom labels permitted by the requested matching options. */
template<class Descriptor>
inline Descriptor screeningDescriptor(const MolGraph& molecule, const ChemOptions& chemistry) {
    Descriptor descriptor{};
    descriptor.n = molecule.n;
    constexpr unsigned bins = sizeof(descriptor.labelCounts) / sizeof(descriptor.labelCounts[0]);
    for (int atom = 0; atom < molecule.n; ++atom)
        ++descriptor.labelCounts[static_cast<unsigned>(smsd::screening_detail::ubLabel(molecule, atom, chemistry)) & (bins - 1)];
    return descriptor;
}
}

// ============================================================================
// Runtime detection
// ============================================================================

/// True when SMSD was compiled with GPU support (CUDA or Metal) AND a device
/// is present at runtime.  Dispatch priority: CUDA > Metal > CPU.
inline bool isAvailable() noexcept {
    if (gpu_kern::detail::disabledByEnvironment()) return false;
#ifdef SMSD_ENABLE_CUDA
    return cuda::cudaAvailable();
#elif defined(SMSD_ENABLE_METAL)
    return metal::metalAvailable();
#else
    return false;
#endif
}

/// Human-readable backend description.
/// Examples:
///   "GPU: Tesla T4 (compute 7.5, 15109 MB)  [OpenMP 4.5, 8 threads]"
///   "Metal GPU: Apple M2 Pro  [OpenMP 5.0, 10 threads]"
///   "CPU: OpenMP 4.5, 16 threads"
///   "CPU: sequential (no OpenMP)"
inline std::string deviceInfo() {
    std::string cpu;
#ifdef _OPENMP
    cpu = "OpenMP " + detail::openMPVersion(_OPENMP)
        + ", " + std::to_string(omp_get_max_threads()) + " threads";
#else
    cpu = "sequential (no OpenMP)";
#endif
    if (gpu_kern::detail::disabledByEnvironment()) return "CPU: " + cpu + "  [GPU disabled]";

#ifdef SMSD_ENABLE_CUDA
    if (cuda::cudaAvailable())
        return "GPU: " + cuda::cudaDeviceInfo() + "  [" + cpu + "]";
    return "CPU: " + cpu + "  [CUDA compiled but no device found]";
#elif defined(SMSD_ENABLE_METAL)
    if (metal::metalAvailable())
        return metal::metalDeviceInfo() + "  [" + cpu + "]";
    return "CPU: " + cpu + "  [Metal compiled but no device found]";
#else
    return "CPU: " + cpu;
#endif
}

// ============================================================================
// Batch similarity-bound screening — auto-dispatch GPU → CPU
// ============================================================================

/// Result of a single screen hit: molecule index + similarity score.
struct ScreenHit {
    int    index;
    double score;
};

/**
 * Screen one query molecule against a target library using similarity upper bound.
 *
 * Automatically dispatches to:
 *   - GPU (CUDA)   when compiled with CUDA and a device is present
 *   - GPU (Metal)  when compiled with Metal on macOS and a device is present
 *   - CPU (OpenMP) otherwise
 *
 * @param query      Query MolGraph
 * @param targets    Vector of target MolGraphs
 * @param chemistry  Chemical matching options for the atom-label bound
 * @param threshold  Minimum similarity to report (default 0.0 = return all)
 * @param numThreads CPU thread count when falling back (0 = auto)
 * @return           Sorted vector of {index, score} above threshold
 */
inline std::vector<ScreenHit> batchSimilarityScreenAuto(
    const MolGraph& query,
    const std::vector<MolGraph>& targets,
    const ChemOptions& chemistry,
    double threshold = 0.0,
    int numThreads = 0)
{
    std::vector<ScreenHit> hits;

#ifdef SMSD_ENABLE_CUDA
    if (!gpu_kern::detail::disabledByEnvironment() && cuda::cudaAvailable()) {
        // Pack query
        cuda::GpuMolecule gq = detail::screeningDescriptor<cuda::GpuMolecule>(query, chemistry);
        // Pack targets
        std::vector<cuda::GpuMolecule> gt;
        gt.reserve(targets.size());
        for (auto& t : targets)
            gt.push_back(detail::screeningDescriptor<cuda::GpuMolecule>(t, chemistry));
        // Run GPU screen
        auto raw = cuda::batchSimilarityScreen(gq, gt, threshold);
        hits.reserve(raw.size());
        for (auto& r : raw)
            hits.push_back({r.index, r.score});
        return hits;
    }
#elif defined(SMSD_ENABLE_METAL)
    if (metal::metalAvailable()) {
        // Pack query
        metal::MetalMolecule mq = detail::screeningDescriptor<metal::MetalMolecule>(query, chemistry);
        // Pack targets
        std::vector<metal::MetalMolecule> mt;
        mt.reserve(targets.size());
        for (auto& t : targets)
            mt.push_back(detail::screeningDescriptor<metal::MetalMolecule>(t, chemistry));
        // Run Metal GPU screen
        auto raw = metal::batchSimilarityScreen(mq, mt,
                                            static_cast<float>(threshold));
        hits.reserve(raw.size());
        for (auto& r : raw)
            hits.push_back({r.index, static_cast<double>(r.score)});
        return hits;
    }
#endif
    // CPU fallback: the screen that smsd::screenTargets runs, on up to numThreads threads.
    const std::vector<double> scores = screening_detail::similarityScores(query, targets, chemistry, numThreads);
    const int N = static_cast<int>(scores.size());
    hits.reserve(N / 4);
    for (int i = 0; i < N; ++i)
        if (scores[i] >= threshold)
            hits.push_back({i, scores[i]});
    return hits;
}

inline std::vector<ScreenHit> batchSimilarityScreenAuto(
    const MolGraph& query, const std::vector<MolGraph>& targets,
    double threshold = 0.0, int numThreads = 0)
{
    return batchSimilarityScreenAuto(query, targets, ChemOptions(), threshold, numThreads);
}

} // namespace gpu
} // namespace smsd

#endif // SMSD_GPU_HPP
