/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * Metal compute kernels for graph matching acceleration on Apple Silicon.
 * Mirrors the CUDA kernels in cuda/graph_kernels.cu.
 */

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <cstdint>
#include <vector>

namespace smsd {
namespace gpu_kern {

// ============================================================================
// Metal shader sources
// ============================================================================

static NSString* const kDomainInitShader = @
"#include <metal_stdlib>\n"
"using namespace metal;\n"
"\n"
"kernel void domain_init(\n"
"    device const int* qAtomicNum [[buffer(0)]],\n"
"    device const int* qCharge    [[buffer(1)]],\n"
"    device const int* qRing      [[buffer(2)]],\n"
"    device const int* qAromatic  [[buffer(3)]],\n"
"    device const int* tAtomicNum [[buffer(4)]],\n"
"    device const int* tCharge    [[buffer(5)]],\n"
"    device const int* tRing      [[buffer(6)]],\n"
"    device const int* tAromatic  [[buffer(7)]],\n"
"    device atomic_uint* domainOut [[buffer(8)]],\n"
"    constant int& Nq       [[buffer(9)]],\n"
"    constant int& Nt       [[buffer(10)]],\n"
"    constant int& tWords2  [[buffer(11)]],\n"
"    constant int& flags    [[buffer(12)]],\n"
"    uint2 tid [[thread_position_in_grid]])\n"
"{\n"
"    int qi = tid.y;\n"
"    int tj = tid.x;\n"
"    if (qi >= Nq || tj >= Nt) return;\n"
"\n"
"    bool matchType = (flags & 1) != 0;\n"
"    bool matchChrg = (flags & 2) != 0;\n"
"    bool ringOnly  = (flags & 4) != 0;\n"
"    bool strictAr  = (flags & 8) != 0;\n"
"\n"
"    bool ok = true;\n"
"    if (matchType && qAtomicNum[qi] != tAtomicNum[tj]) ok = false;\n"
"    if (ok && matchChrg && qCharge[qi] != tCharge[tj]) ok = false;\n"
"    if (ok && strictAr && qAromatic[qi] != tAromatic[tj]) ok = false;\n"
"    if (ok && ringOnly && qRing[qi] != 0 && tRing[tj] == 0) ok = false;\n"
"\n"
"    if (ok) {\n"
"        // Pack into uint32 pairs (Metal atomic_uint is 32-bit)\n"
"        int word32 = tj >> 5;\n"
"        uint bit = 1u << (tj & 31);\n"
"        atomic_fetch_or_explicit(&domainOut[qi * tWords2 + word32], bit, memory_order_relaxed);\n"
"    }\n"
"}\n";

// ============================================================================
// Metal context (reuses DeviceProbe from metal_batch.mm)
// ============================================================================

struct GraphKernelContext {
    id<MTLDevice>               device   = nil;
    id<MTLCommandQueue>         queue    = nil;
    id<MTLComputePipelineState> domainPipeline = nil;
    bool valid = false;

    static GraphKernelContext& instance() {
        static GraphKernelContext ctx;
        static std::once_flag flag;
        GraphKernelContext* p = &ctx;
        std::call_once(flag, [p] { p->init(); });
        return ctx;
    }

    void init() {
        @autoreleasepool {
            device = MTLCreateSystemDefaultDevice();
            if (!device) return;
            queue = [device newCommandQueue];
            if (!queue) return;

            NSError* err = nil;
            id<MTLLibrary> lib = [device newLibraryWithSource:kDomainInitShader
                                                      options:nil error:&err];
            if (!lib) return;
            id<MTLFunction> fn = [lib newFunctionWithName:@"domain_init"];
            if (!fn) return;
            domainPipeline = [device newComputePipelineStateWithFunction:fn error:&err];
            valid = (domainPipeline != nil);
        }
    }
};

// ============================================================================
// API implementations
// ============================================================================

bool graphKernelsAvailable() noexcept {
#if defined(__arm64__) || defined(__aarch64__)
    return GraphKernelContext::instance().valid;
#else
    return false;
#endif
}

bool domainInit(
    int Nq, int Nt, int tWords,
    const int* qAtomicNum, const int* qCharge, const int* qRing, const int* qAromatic,
    const int* tAtomicNum, const int* tCharge, const int* tRing, const int* tAromatic,
    bool matchAtomType, bool matchCharge, bool ringOnly, bool strictArom,
    uint64_t* domainOut)
{
#if !defined(__arm64__) && !defined(__aarch64__)
    return false;
#else
    auto& ctx = GraphKernelContext::instance();
    if (!ctx.valid) return false;

    @autoreleasepool {
        id<MTLDevice> dev = ctx.device;

        // Use 32-bit words for Metal (atomic_uint is 32-bit)
        int tWords32 = (Nt + 31) / 32;
        size_t domBytes32 = static_cast<size_t>(Nq) * tWords32 * sizeof(uint32_t);

        id<MTLBuffer> bQAN = [dev newBufferWithBytes:qAtomicNum length:Nq*sizeof(int) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bQCh = [dev newBufferWithBytes:qCharge    length:Nq*sizeof(int) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bQRi = [dev newBufferWithBytes:qRing      length:Nq*sizeof(int) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bQAr = [dev newBufferWithBytes:qAromatic  length:Nq*sizeof(int) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bTAN = [dev newBufferWithBytes:tAtomicNum length:Nt*sizeof(int) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bTCh = [dev newBufferWithBytes:tCharge    length:Nt*sizeof(int) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bTRi = [dev newBufferWithBytes:tRing      length:Nt*sizeof(int) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bTAr = [dev newBufferWithBytes:tAromatic  length:Nt*sizeof(int) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bDom = [dev newBufferWithLength:domBytes32 options:MTLResourceStorageModeShared];
        memset([bDom contents], 0, domBytes32);

        int flags = (matchAtomType ? 1 : 0) | (matchCharge ? 2 : 0) | (ringOnly ? 4 : 0) | (strictArom ? 8 : 0);

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:ctx.domainPipeline];
        [enc setBuffer:bQAN offset:0 atIndex:0];
        [enc setBuffer:bQCh offset:0 atIndex:1];
        [enc setBuffer:bQRi offset:0 atIndex:2];
        [enc setBuffer:bQAr offset:0 atIndex:3];
        [enc setBuffer:bTAN offset:0 atIndex:4];
        [enc setBuffer:bTCh offset:0 atIndex:5];
        [enc setBuffer:bTRi offset:0 atIndex:6];
        [enc setBuffer:bTAr offset:0 atIndex:7];
        [enc setBuffer:bDom offset:0 atIndex:8];
        [enc setBytes:&Nq       length:sizeof(int) atIndex:9];
        [enc setBytes:&Nt       length:sizeof(int) atIndex:10];
        [enc setBytes:&tWords32 length:sizeof(int) atIndex:11];
        [enc setBytes:&flags    length:sizeof(int) atIndex:12];

        NSUInteger tgsz = ctx.domainPipeline.maxTotalThreadsPerThreadgroup;
        if (tgsz > 256) tgsz = 256;
        MTLSize threadgroupSize = MTLSizeMake(16, 16, 1);
        MTLSize gridSize = MTLSizeMake(((Nt + 15) / 16) * 16, ((Nq + 15) / 16) * 16, 1);
        [enc dispatchThreads:gridSize threadsPerThreadgroup:threadgroupSize];
        [enc endEncoding];
        [cmd commit];
        [cmd waitUntilCompleted];

        // Convert 32-bit domain words to 64-bit
        const uint32_t* dom32 = static_cast<const uint32_t*>([bDom contents]);
        memset(domainOut, 0, static_cast<size_t>(Nq) * tWords * sizeof(uint64_t));
        for (int qi = 0; qi < Nq; ++qi) {
            for (int w32 = 0; w32 < tWords32; ++w32) {
                uint32_t val = dom32[qi * tWords32 + w32];
                if (val == 0) continue;
                int w64 = w32 / 2;
                int shift = (w32 & 1) * 32;
                domainOut[qi * tWords + w64] |= static_cast<uint64_t>(val) << shift;
            }
        }
    }
    return true;
#endif
}

// NLF and feasibility batch: delegate to CPU for now (Metal implementation
// follows the same pattern as domainInit but with additional shader sources)
bool nlfBatchBuild(int, int, const int*, const int*, const int*, const int*,
                   int*, int*, int) { return false; }
bool feasibilityBatch(int, const int*, const int*, int, int, int,
                      const int*, const int*, const uint64_t*, const int*, int*) { return false; }

} // namespace gpu_kern
} // namespace smsd
