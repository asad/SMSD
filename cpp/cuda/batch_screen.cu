/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * GPU-accelerated batch screening for SMSD.
 *
 * Two kernels:
 *   1. RASCAL screening  -- label-histogram Tanimoto upper bound
 *   2. Fingerprint subset -- path-fingerprint subset check for substructure
 *
 * Design:
 *   - One query vs N targets on GPU
 *   - Supports molecules up to 200 atoms (MAX_ATOMS)
 *   - Memory-efficient streaming: targets streamed in chunks when N > 1M
 *   - Full CUDA error handling with SMSD_CUDA_CHECK macro
 *   - Returns indices + scores to host
 */

#include <cuda_runtime.h>
#include <vector>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <utility>

namespace smsd {
namespace cuda {

// ============================================================================
// Constants
// ============================================================================

static constexpr int MAX_ATOMS        = 200;
static constexpr int LABEL_BINS       = 256;
static constexpr int FP_WORDS         = 16;    // 16 * 64 = 1024-bit fingerprint
static constexpr int BLOCK_SIZE       = 256;
static constexpr int CHUNK_SIZE       = 100000; // stream 100K targets per chunk

// ============================================================================
// Error handling
// ============================================================================

#define SMSD_CUDA_CHECK(call)                                                 \
    do {                                                                       \
        cudaError_t err = (call);                                              \
        if (err != cudaSuccess) {                                              \
            throw std::runtime_error(                                          \
                std::string("CUDA error at ") + __FILE__ + ":" +              \
                std::to_string(__LINE__) + ": " +                             \
                cudaGetErrorString(err));                                       \
        }                                                                      \
    } while (0)

// ============================================================================
// GPU data structures -- compact, fixed-size for GPU transfer
// ============================================================================

/// Compact molecule representation for RASCAL screening.
struct alignas(16) GpuMolecule {
    int n;                        // atom count
    int labelCounts[LABEL_BINS];  // histogram of atom labels (label & 0xFF)
};

/// Compact fingerprint for substructure pre-screening.
struct alignas(16) GpuFingerprint {
    uint64_t bits[FP_WORDS];      // 1024-bit path fingerprint
};

/// Result entry: index + score.
struct ScreenResult {
    int    index;
    double score;
};

// ============================================================================
// RASCAL Screening Kernel
// ============================================================================

/// CUDA kernel: compute RASCAL upper bound (Tanimoto-like) for each target.
/// Each thread handles one target molecule.
__global__ void rascalScreenKernel(
    const GpuMolecule* __restrict__ query,
    const GpuMolecule* __restrict__ targets,
    int numTargets,
    double threshold,
    int*    results,    // 1 = pass, 0 = fail
    double* scores      // similarity scores
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    // Cooperatively load query histogram into shared memory (L1 cache speed)
    __shared__ int sharedQueryLabels[LABEL_BINS];
    if (threadIdx.x < LABEL_BINS) {
        sharedQueryLabels[threadIdx.x] = query->labelCounts[threadIdx.x];
    }
    __syncthreads();

    if (tid >= numTargets) return;

    // Load query atom count into register
    const int qn = query->n;
    const int tn = targets[tid].n;

    // Compute label histogram overlap from shared memory
    int overlap = 0;
    for (int i = 0; i < LABEL_BINS; ++i) {
        int qc = sharedQueryLabels[i];
        int tc = targets[tid].labelCounts[i];
        overlap += min(qc, tc);
    }

    int total = qn + tn;
    double sim = (total == 0) ? 1.0 : static_cast<double>(overlap) /
                                       static_cast<double>(total - overlap);
    scores[tid]  = sim;
    results[tid] = (sim >= threshold) ? 1 : 0;
}

// ============================================================================
// Fingerprint Subset Kernel
// ============================================================================

/// CUDA kernel: check if query fingerprint is a subset of each target fingerprint.
/// For substructure: all ON bits in query must be ON in target.
/// Also computes Tanimoto similarity of fingerprints.
__global__ void fpSubsetKernel(
    const GpuFingerprint* __restrict__ query,
    const GpuFingerprint* __restrict__ targets,
    int numTargets,
    int*    results,    // 1 = subset (potential substructure), 0 = fail
    double* scores      // Tanimoto similarity of fingerprints
) {
    // Load query fingerprint into shared memory (128 bytes = 16 uint64 words).
    // All threads in the block read from fast on-chip SRAM instead of global DRAM.
    __shared__ uint64_t sharedQuery[FP_WORDS];
    if (threadIdx.x < FP_WORDS) {
        sharedQuery[threadIdx.x] = query->bits[threadIdx.x];
    }
    __syncthreads();

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= numTargets) return;

    int qBits = 0;     // popcount of query
    int tBits = 0;     // popcount of target
    int andBits = 0;   // popcount of query & target
    bool isSubset = true;

    for (int w = 0; w < FP_WORDS; ++w) {
        uint64_t qw = sharedQuery[w];
        uint64_t tw = targets[tid].bits[w];
        uint64_t aw = qw & tw;

        qBits   += __popcll(qw);
        tBits   += __popcll(tw);
        andBits += __popcll(aw);

        // Subset check: all query bits must be in target
        if (aw != qw) isSubset = false;
    }

    int orBits = qBits + tBits - andBits;
    double sim = (orBits == 0) ? 1.0 : static_cast<double>(andBits) /
                                        static_cast<double>(orBits);
    scores[tid]  = sim;
    results[tid] = isSubset ? 1 : 0;
}

// ============================================================================
// RAII wrappers for device memory
// ============================================================================

template<typename T>
struct DeviceBuffer {
    T* ptr = nullptr;
    size_t count = 0;

    DeviceBuffer() = default;

    explicit DeviceBuffer(size_t n) : count(n) {
        if (n > 0) {
            SMSD_CUDA_CHECK(cudaMalloc(&ptr, n * sizeof(T)));
        }
    }

    ~DeviceBuffer() {
        if (ptr) cudaFree(ptr);
    }

    // Non-copyable, movable
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    DeviceBuffer(DeviceBuffer&& o) noexcept : ptr(o.ptr), count(o.count) {
        o.ptr = nullptr; o.count = 0;
    }

    DeviceBuffer& operator=(DeviceBuffer&& o) noexcept {
        if (ptr) cudaFree(ptr);
        ptr = o.ptr; count = o.count;
        o.ptr = nullptr; o.count = 0;
        return *this;
    }

    void copyFromHost(const T* src, size_t n) {
        SMSD_CUDA_CHECK(cudaMemcpy(ptr, src, n * sizeof(T), cudaMemcpyHostToDevice));
    }

    void copyToHost(T* dst, size_t n) const {
        SMSD_CUDA_CHECK(cudaMemcpy(dst, ptr, n * sizeof(T), cudaMemcpyDeviceToHost));
    }
};

// ============================================================================
// Host API: Batch RASCAL Screening
// ============================================================================

/// Screen query against all targets on GPU using RASCAL upper bound.
/// Returns vector of (index, score) for targets passing the threshold.
/// For large N (> CHUNK_SIZE), streams targets in chunks to limit GPU memory.
inline std::vector<ScreenResult> batchRascalScreen(
    const GpuMolecule& query,
    const std::vector<GpuMolecule>& targets,
    double threshold
) {
    const int N = static_cast<int>(targets.size());
    if (N == 0) return {};

    std::vector<ScreenResult> allHits;

    // Allocate query on device (single molecule, reused across chunks)
    DeviceBuffer<GpuMolecule> d_query(1);
    d_query.copyFromHost(&query, 1);

    // Allocate device buffers ONCE for the maximum chunk size (avoid
    // costly cudaMalloc/cudaFree per chunk iteration)
    int maxChunk = std::min(CHUNK_SIZE, N);
    DeviceBuffer<GpuMolecule> d_targets(maxChunk);
    DeviceBuffer<int>         d_results(maxChunk);
    DeviceBuffer<double>      d_scores(maxChunk);

    // Process in chunks
    for (int offset = 0; offset < N; offset += CHUNK_SIZE) {
        int chunkN = std::min(CHUNK_SIZE, N - offset);

        // Copy chunk to device (reuse pre-allocated buffers)
        d_targets.copyFromHost(targets.data() + offset, chunkN);

        // Launch kernel
        int numBlocks = (chunkN + BLOCK_SIZE - 1) / BLOCK_SIZE;
        rascalScreenKernel<<<numBlocks, BLOCK_SIZE>>>(
            d_query.ptr, d_targets.ptr, chunkN, threshold,
            d_results.ptr, d_scores.ptr);

        // Check for kernel launch errors
        SMSD_CUDA_CHECK(cudaGetLastError());
        SMSD_CUDA_CHECK(cudaDeviceSynchronize());

        // Copy results back
        std::vector<int>    h_results(chunkN);
        std::vector<double> h_scores(chunkN);
        d_results.copyToHost(h_results.data(), chunkN);
        d_scores.copyToHost(h_scores.data(), chunkN);

        // Collect passing indices (adjusted for chunk offset)
        for (int i = 0; i < chunkN; ++i) {
            if (h_results[i]) {
                allHits.push_back({offset + i, h_scores[i]});
            }
        }
    }

    return allHits;
}

/// Convenience: returns just the indices (backward compatible).
inline std::vector<int> batchScreen(
    const GpuMolecule& query,
    const std::vector<GpuMolecule>& targets,
    double threshold
) {
    auto hits = batchRascalScreen(query, targets, threshold);
    std::vector<int> indices;
    indices.reserve(hits.size());
    for (auto& h : hits) indices.push_back(h.index);
    return indices;
}

// ============================================================================
// Host API: Batch Fingerprint Subset Screening
// ============================================================================

/// Screen query fingerprint against all target fingerprints on GPU.
/// Returns vector of (index, tanimotoSimilarity) for targets where
/// the query fingerprint is a subset of the target fingerprint.
inline std::vector<ScreenResult> batchFingerprintScreen(
    const GpuFingerprint& query,
    const std::vector<GpuFingerprint>& targets
) {
    const int N = static_cast<int>(targets.size());
    if (N == 0) return {};

    std::vector<ScreenResult> allHits;

    DeviceBuffer<GpuFingerprint> d_query(1);
    d_query.copyFromHost(&query, 1);

    int maxChunkFP = std::min(CHUNK_SIZE, N);
    DeviceBuffer<GpuFingerprint> d_fp_targets(maxChunkFP);
    DeviceBuffer<int>            d_fp_results(maxChunkFP);
    DeviceBuffer<double>         d_fp_scores(maxChunkFP);

    for (int offset = 0; offset < N; offset += CHUNK_SIZE) {
        int chunkN = std::min(CHUNK_SIZE, N - offset);

        d_fp_targets.copyFromHost(targets.data() + offset, chunkN);

        int numBlocks = (chunkN + BLOCK_SIZE - 1) / BLOCK_SIZE;
        fpSubsetKernel<<<numBlocks, BLOCK_SIZE>>>(
            d_query.ptr, d_fp_targets.ptr, chunkN,
            d_fp_results.ptr, d_fp_scores.ptr);

        SMSD_CUDA_CHECK(cudaGetLastError());
        SMSD_CUDA_CHECK(cudaDeviceSynchronize());

        std::vector<int>    h_results(chunkN);
        std::vector<double> h_scores(chunkN);
        d_fp_results.copyToHost(h_results.data(), chunkN);
        d_fp_scores.copyToHost(h_scores.data(), chunkN);

        for (int i = 0; i < chunkN; ++i) {
            if (h_results[i]) {
                allHits.push_back({offset + i, h_scores[i]});
            }
        }
    }

    return allHits;
}

// ============================================================================
// Conversion utilities: MolGraph -> GpuMolecule / GpuFingerprint
// ============================================================================

/// Convert label vector to GpuMolecule histogram.
/// Designed to work with MolGraph::label vector.
inline GpuMolecule toGpuMolecule(int atomCount, const int* labels) {
    GpuMolecule gm;
    gm.n = atomCount;
    std::memset(gm.labelCounts, 0, sizeof(gm.labelCounts));
    for (int i = 0; i < atomCount; ++i) {
        int bin = labels[i] & 0xFF;
        gm.labelCounts[bin]++;
    }
    return gm;
}

/// Convert a 1024-bit fingerprint (as vector<uint64_t>) to GpuFingerprint.
inline GpuFingerprint toGpuFingerprint(const uint64_t* fp, int numWords) {
    GpuFingerprint gf;
    std::memset(gf.bits, 0, sizeof(gf.bits));
    int copyN = std::min(numWords, FP_WORDS);
    std::memcpy(gf.bits, fp, copyN * sizeof(uint64_t));
    return gf;
}

// ============================================================================
// Query helpers: check CUDA availability
// ============================================================================

/// Returns true if at least one CUDA device is available.
inline bool cudaAvailable() {
    int deviceCount = 0;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);
    if (err != cudaSuccess) return false;
    return deviceCount > 0;
}

/// Returns the name and compute capability of the current CUDA device.
inline std::string cudaDeviceInfo() {
    int device = 0;
    cudaError_t err = cudaGetDevice(&device);
    if (err != cudaSuccess) return "No CUDA device";

    cudaDeviceProp prop;
    err = cudaGetDeviceProperties(&prop, device);
    if (err != cudaSuccess) return "Error querying device";

    return std::string(prop.name) + " (compute " +
           std::to_string(prop.major) + "." +
           std::to_string(prop.minor) + ", " +
           std::to_string(prop.totalGlobalMem / (1024 * 1024)) + " MB)";
}

} // namespace cuda
} // namespace smsd
