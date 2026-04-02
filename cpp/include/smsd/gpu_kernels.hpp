/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * GPU kernel dispatch API for graph matching acceleration.
 * Pure C++ header — no CUDA or Metal types.  Algorithm headers
 * (vf2pp.hpp, mcs.hpp) include this freely.
 *
 * When SMSD_ENABLE_CUDA or SMSD_ENABLE_METAL is defined, the functions
 * are declared extern and resolved at link time from graph_kernels.cu
 * or metal_graph_kernels.mm.  Otherwise, inline stubs return false
 * and the caller falls back to the existing CPU path.
 */
#pragma once
#include <cstdint>

namespace smsd {
namespace gpu_kern {

// ============================================================================
// Kernel 1: Domain Initialization (Nq × Nt atom compatibility matrix)
//
// Computes the common-case atom compatibility checks (atomicNum, charge,
// aromaticity, ringMatchesRingOnly) for all (qi, tj) pairs.  Output is
// a flat bitmatrix: domain[qi * tWords + w] has bit j set if tj is
// compatible with qi.
//
// The GPU domain is a SUPERSET of the true domain — advanced checks
// (tautomer, isotope, chirality, ring fusion) are omitted.  The CPU
// feasibility check catches any false positives during backtracking.
//
// Returns true if GPU succeeded; false means use the CPU path.
// ============================================================================

#if defined(SMSD_ENABLE_CUDA) || defined(SMSD_ENABLE_METAL)

bool graphKernelsAvailable() noexcept;

bool domainInit(
    int Nq, int Nt, int tWords,
    const int* qAtomicNum, const int* qCharge, const int* qRing, const int* qAromatic,
    const int* tAtomicNum, const int* tCharge, const int* tRing, const int* tAromatic,
    bool matchAtomType, bool matchCharge, bool ringOnly, bool strictArom,
    uint64_t* domainOut);

// ============================================================================
// Kernel 2: NLF Batch Compute
//
// Pre-computes NLF-1 histograms for all atoms across multiple molecules.
// Input: concatenated CSR adjacency (offsets + neighbor list + labels).
// Output: flat NLF pairs array.
// Returns true if GPU succeeded.
// ============================================================================

bool nlfBatchBuild(
    int totalAtoms, int totalEdges,
    const int* neighborOffsets,   // [totalAtoms + 1]
    const int* neighborList,      // [totalEdges]
    const int* atomicNum,         // [totalAtoms]
    const int* aromatic,          // [totalAtoms]
    int* nlfPairsOut,             // [totalAtoms * maxPairsPerAtom * 2]
    int* nlfLengthsOut,           // [totalAtoms] — number of pairs per atom
    int maxPairsPerAtom);

// ============================================================================
// Kernel 3: Batch Feasibility Check
//
// Evaluates feasibility for multiple (qi, tj) candidate pairs in parallel.
// Checks: degree, NLF-1 subset, mapped-neighbor adjacency (bit-parallel).
// Bond compatibility stays on CPU.
// Returns true if GPU succeeded.
// ============================================================================

bool feasibilityBatch(
    int numCandidates,
    const int* qiList,            // [numCandidates]
    const int* tjList,            // [numCandidates]
    int Nq, int Nt, int tWords,
    const int* qDegree, const int* tDegree,
    const uint64_t* tAdjLong,     // [Nt * tWords] flat adjacency bitmatrix
    const int* q2t,               // [Nq] current mapping (-1 = unmapped)
    int* feasibleOut);            // [numCandidates] 0/1 results

#else
// No GPU — inline stubs that signal "use CPU path"
inline bool graphKernelsAvailable() noexcept { return false; }
inline bool domainInit(int,int,int,const int*,const int*,const int*,const int*,
    const int*,const int*,const int*,const int*,bool,bool,bool,bool,uint64_t*)
    { return false; }
inline bool nlfBatchBuild(int,int,const int*,const int*,const int*,const int*,
    int*,int*,int) { return false; }
inline bool feasibilityBatch(int,const int*,const int*,int,int,int,
    const int*,const int*,const uint64_t*,const int*,int*) { return false; }
#endif

} // namespace gpu_kern
} // namespace smsd
