/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * Metal/MPS batch screening backend — macOS 10.14+ / Apple Silicon.
 *
 * This header is PURE C++ — no Objective-C types, no #import.
 * It can be included from any compilation unit (.cpp / .mm / pybind11).
 *
 * The implementation lives in metal/metal_batch.mm (compiled as Obj-C++).
 * Link against smsd_metal to activate this backend:
 *
 *   target_link_libraries(my_target PRIVATE smsd smsd_metal)
 *
 * When smsd_metal is linked, SMSD_ENABLE_METAL is defined as a PUBLIC
 * compile definition, and gpu.hpp automatically dispatches to Metal.
 *
 * On Apple Silicon the Metal buffers use MTLResourceStorageModeShared —
 * the CPU and GPU share the same physical memory, so data is never copied.
 * On Intel Macs with a discrete GPU the runtime falls back gracefully.
 */
#pragma once
#ifndef SMSD_METAL_GPU_HPP
#define SMSD_METAL_GPU_HPP

#include <string>
#include <vector>

namespace smsd {
namespace metal {

// ---------------------------------------------------------------------------
// Molecule descriptor (mirrors cuda::GpuMolecule for code-path symmetry)
// ---------------------------------------------------------------------------

/// Number of atom-label histogram bins (element atomic number 0–255).
static constexpr int METAL_LABEL_BINS = 256;

/**
 * Compact molecule representation for GPU screening.
 *
 * n             — heavy-atom count
 * labelCounts   — per-element histogram: labelCounts[atomicNumber] = count
 *
 * The struct is 4 + 256*4 + 12 = 1040 bytes (16-byte aligned).
 * Explicit padding ensures the C++ stride matches the Metal shader stride.
 * On Apple Silicon (unified memory) these are passed as zero-copy shared buffers.
 */
struct alignas(16) MetalMolecule {
    int n = 0;
    int labelCounts[METAL_LABEL_BINS] = {};
    int _padding[3] = {};  // explicit 12-byte pad to reach 1040 (16-byte multiple)
};

/**
 * Single screening result: target index + RASCAL upper-bound score.
 */
struct MetalScreenResult {
    int   index = 0;
    float score = 0.f;
};

// ---------------------------------------------------------------------------
// Runtime API (implemented in metal/metal_batch.mm)
// ---------------------------------------------------------------------------

/// Returns true if a Metal-capable GPU device is available at runtime.
/// Always false when SMSD is built without Metal support.
bool metalAvailable() noexcept;

/**
 * Human-readable GPU description.
 * Examples:
 *   "Metal GPU: Apple M2 Pro"
 *   "Metal GPU: AMD Radeon Pro 5500M"
 *   "Metal: not available"
 */
std::string metalDeviceInfo();

/**
 * Pack a MolGraph-style atom array into a MetalMolecule descriptor.
 *
 * @param atomCount  Number of heavy atoms.
 * @param labels     Array of length atomCount; each entry is the atomic number.
 */
MetalMolecule toMetalMolecule(int atomCount, const int* labels);

/**
 * Batch RASCAL upper-bound screen on the Metal GPU.
 *
 * Dispatches a parallel Metal compute kernel; each GPU thread evaluates one
 * (query, target[i]) pair using the label-histogram bound:
 *
 *   score[i] = sum_Z  min(q.labelCounts[Z], t.labelCounts[Z])
 *              ÷  max(q.n, t.n)
 *
 * This is a conservative upper bound — any pair whose true RASCAL similarity
 * exceeds `threshold` is guaranteed to appear in the output.
 *
 * @param query      Packed query molecule.
 * @param targets    Packed target molecules.
 * @param threshold  Minimum score to include in results (default 0 = all).
 * @return           Unsorted list of hits above threshold.
 */
std::vector<MetalScreenResult> batchRascalScreen(
    const MetalMolecule&              query,
    const std::vector<MetalMolecule>& targets,
    float                             threshold = 0.f);

} // namespace metal
} // namespace smsd

#endif // SMSD_METAL_GPU_HPP
