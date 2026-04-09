/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * Metal compute backend for SMSD batch RASCAL screening.
 *
 * Compiled as Objective-C++ (.mm) with ARC enabled (-fobjc-arc).
 * Only compiled when SMSD_BUILD_METAL is AUTO (macOS, Metal.framework found)
 * or explicitly ON.  The public API is declared in smsd/metal_gpu.hpp (pure C++).
 *
 * Architecture
 * ------------
 * A lazy singleton (MetalContext) initialises once per process:
 *   1. Obtain the system default MTLDevice.
 *   2. Compile the rascal_screen kernel from an embedded source string
 *      using [MTLDevice newLibraryWithSource:options:error:].
 *      (No separate .metallib file needed; the Metal compiler ships with macOS.)
 *   3. Create an MTLComputePipelineState and MTLCommandQueue.
 *
 * Per call to batchRascalScreen():
 *   1. Pack query + targets into MTLBuffers (shared memory on Apple Silicon
 *      — zero-copy; on discrete GPUs the runtime copies to VRAM automatically).
 *   2. Encode a single compute pass: one thread per target molecule.
 *   3. Commit, wait for completion, read float results back.
 *   4. Filter results above threshold and return.
 *
 * Compatibility
 * -------------
 *   macOS 10.14+  (Metal Shading Language 2.1, dispatchThreadgroups API)
 *   Apple Silicon  (M1 / M2 / M3 / M4) — unified memory, zero-copy
 *   Intel Mac      (integrated/discrete AMD or Intel GPU) — normal copy path
 *   iOS / tvOS     — header-compatible; not tested
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include "smsd/metal_gpu.hpp"

#include <mutex>
#include <stdexcept>

// ---------------------------------------------------------------------------
// Embedded Metal shader source
// ---------------------------------------------------------------------------

// The kernel computes a per-target RASCAL label-histogram upper bound.
// Layout of MetalMolecule must match the C++ struct in metal_gpu.hpp:
//   int n;                  // offset 0,  4 bytes
//   int labelCounts[256];   // offset 4, 1024 bytes
//   int _padding[3];        // offset 1028, 12 bytes
//   total: 1040 bytes (16-byte aligned, matches C++ alignas(16) stride)

static NSString* const kShaderSource = @
"#include <metal_stdlib>\n"
"using namespace metal;\n"
"\n"
"struct MetalMolecule {\n"
"    int n;\n"
"    int labelCounts[256];\n"
"    int _padding[3];\n"
"};\n"
"\n"
"kernel void rascal_screen(\n"
"    device const MetalMolecule* targets [[buffer(0)]],\n"
"    device const MetalMolecule* query   [[buffer(1)]],\n"
"    device       float*         results [[buffer(2)]],\n"
"    constant     uint&          N       [[buffer(3)]],\n"
"    uint idx [[thread_position_in_grid]])\n"
"{\n"
"    if (idx >= N) return;\n"
"    int qn = query[0].n;\n"
"    int tn = targets[idx].n;\n"
"    if (qn == 0 || tn == 0) { results[idx] = 0.0f; return; }\n"
"\n"
"    // Label-histogram compatibility bound:\n"
"    //   score = sum_Z min(q.count[Z], t.count[Z]) / max(qn, tn)\n"
"    int compatible = 0;\n"
"    for (int l = 0; l < 256; ++l) {\n"
"        int a = query[0].labelCounts[l];\n"
"        int b = targets[idx].labelCounts[l];\n"
"        compatible += (a < b) ? a : b;\n"
"    }\n"
"    int denom = (qn > tn) ? qn : tn;\n"
"    results[idx] = (denom > 0) ? float(compatible) / float(denom) : 0.0f;\n"
"}\n";

// ---------------------------------------------------------------------------
// Singleton MetalContext
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Two-phase initialisation
//
// Phase 1 — DeviceProbe  (fast, < 1ms, safe to call from any context)
//   Just calls MTLCreateSystemDefaultDevice().  Used by metalAvailable()
//   and metalDeviceInfo() so those functions never trigger shader compilation.
//
// Phase 2 — MetalContext (slow: compiles Metal shader via newLibraryWithSource:)
//   Lazy-initialised only when batchRascalScreen() is first called.
//   newLibraryWithSource: invokes the Metal shader compiler and can take
//   100ms–several seconds on first call, or hang indefinitely in headless
//   CI environments.  Keeping it out of the availability probe means
//   import smsd; smsd.gpu_device_info() returns instantly everywhere.
// ---------------------------------------------------------------------------

namespace {

// --- Phase 1: lightweight device probe ---

struct DeviceProbe {
    id<MTLDevice> device = nil;
    bool          found  = false;

    static DeviceProbe& instance() {
        static DeviceProbe p;
        static std::once_flag flag;
        std::call_once(flag, [] {
#if !defined(__arm64__) && !defined(__aarch64__)
            // x86_64 macOS: Metal.framework present but no real GPU in CI.
            // newLibraryWithSource: blocks; skip entirely.
            (void)0;
#else
            @autoreleasepool {
                p.device = MTLCreateSystemDefaultDevice();
                p.found  = (p.device != nil);
            }
#endif
        });
        return p;
    }
    DeviceProbe() = default;
    DeviceProbe(const DeviceProbe&) = delete;
    DeviceProbe& operator=(const DeviceProbe&) = delete;
};

// --- Phase 2: full pipeline (shader compilation, lazy) ---

struct MetalContext {
    id<MTLCommandQueue>         queue    = nil;
    id<MTLComputePipelineState> pipeline = nil;
    bool valid = false;

    static MetalContext& instance() {
        static MetalContext ctx;
        static std::once_flag flag;
        MetalContext* p = &ctx;
        std::call_once(flag, [p] { p->init(); });
        return ctx;
    }

    void init() {
        auto& probe = DeviceProbe::instance();
        if (!probe.found) return;
        @autoreleasepool {
            queue = [probe.device newCommandQueue];
            if (!queue) return;

            NSError* error = nil;
            id<MTLLibrary> lib =
                [probe.device newLibraryWithSource:kShaderSource
                                           options:nil
                                             error:&error];
            if (!lib) {
                NSLog(@"SMSD Metal: shader compilation failed: %@", error);
                return;
            }

            id<MTLFunction> fn = [lib newFunctionWithName:@"rascal_screen"];
            if (!fn) {
                NSLog(@"SMSD Metal: function 'rascal_screen' not found");
                return;
            }

            pipeline = [probe.device newComputePipelineStateWithFunction:fn
                                                                   error:&error];
            if (!pipeline) {
                NSLog(@"SMSD Metal: pipeline creation failed: %@", error);
                return;
            }

            valid = true;
        }
    }

    MetalContext() = default;
    MetalContext(const MetalContext&) = delete;
    MetalContext& operator=(const MetalContext&) = delete;
};

} // anonymous namespace

// ---------------------------------------------------------------------------
// Public API implementation
// ---------------------------------------------------------------------------

namespace smsd {
namespace metal {

bool metalAvailable() noexcept {
    // Phase 1 only — no shader compilation.
    return DeviceProbe::instance().found;
}

std::string metalDeviceInfo() {
    // Phase 1 only — returns instantly, no shader compilation.
    @autoreleasepool {
        auto& probe = DeviceProbe::instance();
        if (!probe.found) return "Metal: not available";
        NSString* name = probe.device.name;
        return "Metal GPU: " + std::string(name ? [name UTF8String] : "unknown");
    }
}

MetalMolecule toMetalMolecule(int atomCount, const int* labels) {
    MetalMolecule m{};
    m.n = atomCount;
    for (int i = 0; i < atomCount; ++i) {
        int z = labels[i];
        if (z >= 0 && z < METAL_LABEL_BINS)
            ++m.labelCounts[z];
    }
    return m;
}

std::vector<MetalScreenResult> batchRascalScreen(
    const MetalMolecule&              query,
    const std::vector<MetalMolecule>& targets,
    float                             threshold)
{
    std::vector<MetalScreenResult> hits;
    if (targets.empty()) return hits;

    @autoreleasepool {
        // Phase 2: trigger lazy shader compilation on first real batch call.
        auto& ctx   = MetalContext::instance();
        auto& probe = DeviceProbe::instance();
        if (!ctx.valid || !probe.found) return hits;

        id<MTLDevice> dev = probe.device;
        const NSUInteger N = static_cast<NSUInteger>(targets.size());

        // ---- Allocate buffers -----------------------------------------------
        // MTLResourceStorageModeShared = CPU+GPU share the same physical memory
        // on Apple Silicon (zero-copy).  On Intel Macs with discrete GPUs the
        // runtime transparently copies to VRAM and back.

        id<MTLBuffer> targBuf =
            [dev newBufferWithBytes:targets.data()
                             length:N * sizeof(MetalMolecule)
                            options:MTLResourceStorageModeShared];

        id<MTLBuffer> queryBuf =
            [dev newBufferWithBytes:&query
                             length:sizeof(MetalMolecule)
                            options:MTLResourceStorageModeShared];

        id<MTLBuffer> resBuf =
            [dev newBufferWithLength:N * sizeof(float)
                             options:MTLResourceStorageModeShared];

        if (!targBuf || !queryBuf || !resBuf) return hits;

        // ---- Encode compute pass -------------------------------------------
        id<MTLCommandBuffer>        cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

        [enc setComputePipelineState:ctx.pipeline];
        [enc setBuffer:targBuf  offset:0 atIndex:0];
        [enc setBuffer:queryBuf offset:0 atIndex:1];
        [enc setBuffer:resBuf   offset:0 atIndex:2];
        uint nVal = static_cast<uint>(N);
        [enc setBytes:&nVal length:sizeof(uint) atIndex:3];

        // Thread-group size: clamp to pipeline max, prefer 256 for efficiency.
        NSUInteger tgsz = ctx.pipeline.maxTotalThreadsPerThreadgroup;
        if (tgsz > 256) tgsz = 256;
        if (tgsz < 1)   tgsz = 1;

        // Ceiling-division dispatch (works on macOS 10.11+ via
        // dispatchThreadgroups:threadsPerThreadgroup:).
        NSUInteger groups = (N + tgsz - 1) / tgsz;
        [enc dispatchThreadgroups:MTLSizeMake(groups, 1, 1)
           threadsPerThreadgroup:MTLSizeMake(tgsz, 1, 1)];
        [enc endEncoding];

        [cmd commit];
        [cmd waitUntilCompleted];

        // ---- Collect hits ---------------------------------------------------
        const float* scores = static_cast<const float*>(resBuf.contents);
        hits.reserve(N / 4 + 1);
        for (NSUInteger i = 0; i < N; ++i)
            if (scores[i] >= threshold)
                hits.push_back({static_cast<int>(i), scores[i]});
    }
    return hits;
}

} // namespace metal
} // namespace smsd
