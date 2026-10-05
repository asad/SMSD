/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * CUDA kernels for graph matching acceleration:
 *   1. Domain initialization (Nq × Nt atom compatibility)
 *   2. NLF batch compute (per-atom neighbor label histograms)
 *   3. Batch feasibility checks (degree + adjacency)
 */

#include <cuda_runtime.h>
#include <cstdint>
#include <cstdio>
#include <algorithm>

// Error checking macro (matches batch_screen.cu pattern)
#define SMSD_CUDA_CHECK(call) do { \
    cudaError_t err = (call); \
    if (err != cudaSuccess) { return false; } \
} while (0)

namespace smsd {
namespace gpu_kern {

// ============================================================================
// GPU workspace: single buffer, sub-allocated via offset arithmetic.
// Avoids repeated cudaMalloc/cudaFree per kernel launch.
// ============================================================================

struct CudaWorkspace {
    void* ptr = nullptr;
    size_t capacity = 0;

    void ensure(size_t needed) {
        if (needed <= capacity) return;
        if (ptr) cudaFree(ptr);
        cudaMalloc(&ptr, needed);
        capacity = needed;
    }

    ~CudaWorkspace() { if (ptr) cudaFree(ptr); }

    CudaWorkspace() = default;
    CudaWorkspace(const CudaWorkspace&) = delete;
    CudaWorkspace& operator=(const CudaWorkspace&) = delete;
};

static thread_local CudaWorkspace workspace;

/// Return a device pointer at `offset` bytes into the workspace buffer.
static inline char* wsAt(size_t offset) {
    return static_cast<char*>(workspace.ptr) + offset;
}

// ============================================================================
// Runtime check
// ============================================================================

bool graphKernelsAvailable() noexcept {
    int count = 0;
    return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}

// ============================================================================
// Kernel 1: Domain Initialization
// ============================================================================

__global__ void domainInitKernel(
    int Nq, int Nt, int tWords,
    const int* __restrict__ qAtomicNum,
    const int* __restrict__ qCharge,
    const int* __restrict__ qRing,
    const int* __restrict__ qAromatic,
    const int* __restrict__ tAtomicNum,
    const int* __restrict__ tCharge,
    const int* __restrict__ tRing,
    const int* __restrict__ tAromatic,
    bool matchAtomType, bool matchCharge, bool ringOnly, bool strictArom,
    unsigned long long* domainOut)  // uint64_t = unsigned long long on CUDA
{
    // 2D grid: blockIdx.y = query atom group, blockIdx.x = target atom group
    int qi = blockIdx.y * blockDim.y + threadIdx.y;
    int tj = blockIdx.x * blockDim.x + threadIdx.x;

    if (qi >= Nq || tj >= Nt) return;

    // Common-case atom compatibility (matches atomsCompatFast common path)
    bool ok = true;
    if (matchAtomType && qAtomicNum[qi] != tAtomicNum[tj]) ok = false;
    if (ok && matchCharge && qCharge[qi] != tCharge[tj]) ok = false;
    if (ok && strictArom && qAromatic[qi] != tAromatic[tj]) ok = false;
    if (ok && ringOnly && qRing[qi] && !tRing[tj]) ok = false;

    if (ok) {
        int word = tj >> 6;
        unsigned long long bit = 1ULL << (tj & 63);
        atomicOr(&domainOut[qi * tWords + word], bit);
    }
}

bool domainInit(
    int Nq, int Nt, int tWords,
    const int* qAtomicNum, const int* qCharge, const int* qRing, const int* qAromatic,
    const int* tAtomicNum, const int* tCharge, const int* tRing, const int* tAromatic,
    bool matchAtomType, bool matchCharge, bool ringOnly, bool strictArom,
    uint64_t* domainOut)
{
    // Compute sizes and carve sub-buffers from a single workspace allocation.
    // Layout: [qAN | qCh | qRi | qAr | tAN | tCh | tRi | tAr | dom]
    size_t qBytes = Nq * sizeof(int);
    size_t tBytes = Nt * sizeof(int);
    size_t domBytes = static_cast<size_t>(Nq) * tWords * sizeof(unsigned long long);

    size_t off = 0;
    size_t o_qAN = off; off += qBytes;
    size_t o_qCh = off; off += qBytes;
    size_t o_qRi = off; off += qBytes;
    size_t o_qAr = off; off += qBytes;
    size_t o_tAN = off; off += tBytes;
    size_t o_tCh = off; off += tBytes;
    size_t o_tRi = off; off += tBytes;
    size_t o_tAr = off; off += tBytes;
    // Align domain output to 8 bytes for unsigned long long
    off = (off + 7) & ~size_t(7);
    size_t o_dom = off; off += domBytes;

    workspace.ensure(off);

    int* d_qAN = reinterpret_cast<int*>(wsAt(o_qAN));
    int* d_qCh = reinterpret_cast<int*>(wsAt(o_qCh));
    int* d_qRi = reinterpret_cast<int*>(wsAt(o_qRi));
    int* d_qAr = reinterpret_cast<int*>(wsAt(o_qAr));
    int* d_tAN = reinterpret_cast<int*>(wsAt(o_tAN));
    int* d_tCh = reinterpret_cast<int*>(wsAt(o_tCh));
    int* d_tRi = reinterpret_cast<int*>(wsAt(o_tRi));
    int* d_tAr = reinterpret_cast<int*>(wsAt(o_tAr));
    unsigned long long* d_dom = reinterpret_cast<unsigned long long*>(wsAt(o_dom));

    SMSD_CUDA_CHECK(cudaMemcpy(d_qAN, qAtomicNum, qBytes, cudaMemcpyHostToDevice));
    SMSD_CUDA_CHECK(cudaMemcpy(d_qCh, qCharge, qBytes, cudaMemcpyHostToDevice));
    SMSD_CUDA_CHECK(cudaMemcpy(d_qRi, qRing, qBytes, cudaMemcpyHostToDevice));
    SMSD_CUDA_CHECK(cudaMemcpy(d_qAr, qAromatic, qBytes, cudaMemcpyHostToDevice));
    SMSD_CUDA_CHECK(cudaMemcpy(d_tAN, tAtomicNum, tBytes, cudaMemcpyHostToDevice));
    SMSD_CUDA_CHECK(cudaMemcpy(d_tCh, tCharge, tBytes, cudaMemcpyHostToDevice));
    SMSD_CUDA_CHECK(cudaMemcpy(d_tRi, tRing, tBytes, cudaMemcpyHostToDevice));
    SMSD_CUDA_CHECK(cudaMemcpy(d_tAr, tAromatic, tBytes, cudaMemcpyHostToDevice));
    SMSD_CUDA_CHECK(cudaMemset(d_dom, 0, domBytes));

    dim3 block(32, 8);  // 256 threads
    dim3 grid((Nt + 31) / 32, (Nq + 7) / 8);

    domainInitKernel<<<grid, block>>>(
        Nq, Nt, tWords,
        d_qAN, d_qCh, d_qRi, d_qAr,
        d_tAN, d_tCh, d_tRi, d_tAr,
        matchAtomType, matchCharge, ringOnly, strictArom,
        d_dom);

    SMSD_CUDA_CHECK(cudaGetLastError());
    SMSD_CUDA_CHECK(cudaDeviceSynchronize());
    SMSD_CUDA_CHECK(cudaMemcpy(domainOut, d_dom, domBytes, cudaMemcpyDeviceToHost));

    // No per-call cudaFree — workspace persists across calls.
    return true;
}

// ============================================================================
// Kernel 2: NLF Batch Build
// ============================================================================

static constexpr int MAX_NLF_LABELS = 32;  // max distinct neighbor labels per atom

__global__ void nlfBuildKernel(
    int totalAtoms,
    const int* __restrict__ neighborOffsets,
    const int* __restrict__ neighborList,
    const int* __restrict__ atomicNum,
    const int* __restrict__ aromatic,
    int* nlfPairsOut,
    int* nlfLengthsOut,
    int maxPairs)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= totalAtoms) return;

    int start = neighborOffsets[idx];
    int end = neighborOffsets[idx + 1];
    int deg = end - start;

    // Build histogram in registers (small — typical degree 1-6)
    int labels[MAX_NLF_LABELS];
    int counts[MAX_NLF_LABELS];
    int nDistinct = 0;

    for (int e = start; e < end; ++e) {
        int nb = neighborList[e];
        int lab = (atomicNum[nb] << 1) | (aromatic[nb] ? 1 : 0);
        // Linear scan (degree is tiny)
        bool found = false;
        for (int k = 0; k < nDistinct; ++k) {
            if (labels[k] == lab) { counts[k]++; found = true; break; }
        }
        if (!found && nDistinct < MAX_NLF_LABELS) {
            labels[nDistinct] = lab;
            counts[nDistinct] = 1;
            nDistinct++;
        }
    }

    // Insertion sort by label
    for (int i = 1; i < nDistinct; ++i) {
        int kl = labels[i], kc = counts[i], j = i - 1;
        while (j >= 0 && labels[j] > kl) {
            labels[j + 1] = labels[j]; counts[j + 1] = counts[j]; --j;
        }
        labels[j + 1] = kl; counts[j + 1] = kc;
    }

    // Write output
    int outBase = idx * maxPairs * 2;
    int pairsWritten = min(nDistinct, maxPairs);
    for (int k = 0; k < pairsWritten; ++k) {
        nlfPairsOut[outBase + k * 2] = labels[k];
        nlfPairsOut[outBase + k * 2 + 1] = counts[k];
    }
    nlfLengthsOut[idx] = pairsWritten * 2;  // length in ints (pairs × 2)
}

bool nlfBatchBuild(
    int totalAtoms, int totalEdges,
    const int* neighborOffsets, const int* neighborList,
    const int* atomicNum, const int* aromatic,
    int* nlfPairsOut, int* nlfLengthsOut, int maxPairsPerAtom)
{
    // Compute sizes and carve sub-buffers from the workspace.
    // Layout: [offsets | nbList | atomNum | arom | nlfPairs | nlfLens]
    size_t offsetBytes = (totalAtoms + 1) * sizeof(int);
    size_t edgeBytes = totalEdges * sizeof(int);
    size_t atomBytes = totalAtoms * sizeof(int);
    size_t pairBytes = static_cast<size_t>(totalAtoms) * maxPairsPerAtom * 2 * sizeof(int);

    size_t off = 0;
    size_t o_offsets  = off; off += offsetBytes;
    size_t o_nbList   = off; off += edgeBytes;
    size_t o_atomNum  = off; off += atomBytes;
    size_t o_arom     = off; off += atomBytes;
    size_t o_nlfPairs = off; off += pairBytes;
    size_t o_nlfLens  = off; off += atomBytes;

    workspace.ensure(off);

    int* d_offsets  = reinterpret_cast<int*>(wsAt(o_offsets));
    int* d_nbList   = reinterpret_cast<int*>(wsAt(o_nbList));
    int* d_atomNum  = reinterpret_cast<int*>(wsAt(o_atomNum));
    int* d_arom     = reinterpret_cast<int*>(wsAt(o_arom));
    int* d_nlfPairs = reinterpret_cast<int*>(wsAt(o_nlfPairs));
    int* d_nlfLens  = reinterpret_cast<int*>(wsAt(o_nlfLens));

    SMSD_CUDA_CHECK(cudaMemcpy(d_offsets, neighborOffsets, offsetBytes, cudaMemcpyHostToDevice));
    SMSD_CUDA_CHECK(cudaMemcpy(d_nbList, neighborList, edgeBytes, cudaMemcpyHostToDevice));
    SMSD_CUDA_CHECK(cudaMemcpy(d_atomNum, atomicNum, atomBytes, cudaMemcpyHostToDevice));
    SMSD_CUDA_CHECK(cudaMemcpy(d_arom, aromatic, atomBytes, cudaMemcpyHostToDevice));
    SMSD_CUDA_CHECK(cudaMemset(d_nlfPairs, 0, pairBytes));

    int blocks = (totalAtoms + 255) / 256;
    nlfBuildKernel<<<blocks, 256>>>(
        totalAtoms, d_offsets, d_nbList, d_atomNum, d_arom,
        d_nlfPairs, d_nlfLens, maxPairsPerAtom);

    SMSD_CUDA_CHECK(cudaGetLastError());
    SMSD_CUDA_CHECK(cudaDeviceSynchronize());
    SMSD_CUDA_CHECK(cudaMemcpy(nlfPairsOut, d_nlfPairs, pairBytes, cudaMemcpyDeviceToHost));
    SMSD_CUDA_CHECK(cudaMemcpy(nlfLengthsOut, d_nlfLens, atomBytes, cudaMemcpyDeviceToHost));

    // No per-call cudaFree — workspace persists across calls.
    return true;
}

// ============================================================================
// Kernel 3: Batch Feasibility Check
// ============================================================================

__global__ void feasibilityKernel(
    int numCandidates, int Nq, int Nt, int tWords,
    const int* __restrict__ qiList,
    const int* __restrict__ tjList,
    const int* __restrict__ qDegree,
    const int* __restrict__ tDegree,
    const unsigned long long* __restrict__ tAdjLong,
    const int* __restrict__ q2t,
    int* feasibleOut)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numCandidates) return;

    int qi = qiList[idx];
    int tj = tjList[idx];

    // Degree check
    if (tDegree[tj] < qDegree[qi]) { feasibleOut[idx] = 0; return; }

    // Mapped-neighbor adjacency: for each mapped neighbor of qi,
    // verify that tj is adjacent to the mapped target atom
    const unsigned long long* tAdj = &tAdjLong[tj * tWords];
    // Note: we can't iterate query neighbors on GPU without the neighbor list.
    // This kernel handles the degree check; full adjacency stays on CPU.
    feasibleOut[idx] = 1;
}

bool feasibilityBatch(
    int numCandidates,
    const int* qiList, const int* tjList,
    int Nq, int Nt, int tWords,
    const int* qDegree, const int* tDegree,
    const uint64_t* tAdjLong,
    const int* q2t,
    int* feasibleOut)
{
    // Compute sizes and carve sub-buffers from the workspace.
    // Layout: [qi | tj | qDeg | tDeg | q2t | result | tAdj(aligned)]
    size_t candBytes = numCandidates * sizeof(int);
    size_t adjBytes = static_cast<size_t>(Nt) * tWords * sizeof(unsigned long long);

    size_t off = 0;
    size_t o_qi     = off; off += candBytes;
    size_t o_tj     = off; off += candBytes;
    size_t o_qDeg   = off; off += Nq * sizeof(int);
    size_t o_tDeg   = off; off += Nt * sizeof(int);
    size_t o_q2t    = off; off += Nq * sizeof(int);
    size_t o_result = off; off += candBytes;
    // Align adjacency matrix to 8 bytes for unsigned long long
    off = (off + 7) & ~size_t(7);
    size_t o_tAdj   = off; off += adjBytes;

    workspace.ensure(off);

    int* d_qi     = reinterpret_cast<int*>(wsAt(o_qi));
    int* d_tj     = reinterpret_cast<int*>(wsAt(o_tj));
    int* d_qDeg   = reinterpret_cast<int*>(wsAt(o_qDeg));
    int* d_tDeg   = reinterpret_cast<int*>(wsAt(o_tDeg));
    int* d_q2t    = reinterpret_cast<int*>(wsAt(o_q2t));
    int* d_result = reinterpret_cast<int*>(wsAt(o_result));
    unsigned long long* d_tAdj = reinterpret_cast<unsigned long long*>(wsAt(o_tAdj));

    SMSD_CUDA_CHECK(cudaMemcpy(d_qi, qiList, candBytes, cudaMemcpyHostToDevice));
    SMSD_CUDA_CHECK(cudaMemcpy(d_tj, tjList, candBytes, cudaMemcpyHostToDevice));
    SMSD_CUDA_CHECK(cudaMemcpy(d_qDeg, qDegree, Nq * sizeof(int), cudaMemcpyHostToDevice));
    SMSD_CUDA_CHECK(cudaMemcpy(d_tDeg, tDegree, Nt * sizeof(int), cudaMemcpyHostToDevice));
    SMSD_CUDA_CHECK(cudaMemcpy(d_tAdj, tAdjLong, adjBytes, cudaMemcpyHostToDevice));
    SMSD_CUDA_CHECK(cudaMemcpy(d_q2t, q2t, Nq * sizeof(int), cudaMemcpyHostToDevice));

    int blocks = (numCandidates + 255) / 256;
    feasibilityKernel<<<blocks, 256>>>(
        numCandidates, Nq, Nt, tWords,
        d_qi, d_tj, d_qDeg, d_tDeg, d_tAdj, d_q2t, d_result);

    SMSD_CUDA_CHECK(cudaGetLastError());
    SMSD_CUDA_CHECK(cudaDeviceSynchronize());
    SMSD_CUDA_CHECK(cudaMemcpy(feasibleOut, d_result, candBytes, cudaMemcpyDeviceToHost));

    // No per-call cudaFree — workspace persists across calls.
    return true;
}

} // namespace gpu_kern
} // namespace smsd
