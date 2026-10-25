/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#define SMSD_ENABLE_METAL
#include "smsd/gpu_kernels.hpp"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <stdexcept>

static std::string failure;

/** Deterministic Metal protocol double; no command reaches a physical device. */
@interface TestMetalObject : NSObject {
    NSMutableData* storage_;
}
- (instancetype)initWithBytes:(const void*)bytes length:(NSUInteger)length;
- (id<MTLCommandQueue>)newCommandQueue;
- (id<MTLLibrary>)newLibraryWithSource:(NSString*)source options:(MTLCompileOptions*)options error:(NSError**)error;
- (id<MTLFunction>)newFunctionWithName:(NSString*)name;
- (id<MTLComputePipelineState>)newComputePipelineStateWithFunction:(id<MTLFunction>)function error:(NSError**)error;
- (id<MTLBuffer>)newBufferWithBytes:(const void*)bytes length:(NSUInteger)length options:(MTLResourceOptions)options;
- (id<MTLBuffer>)newBufferWithLength:(NSUInteger)length options:(MTLResourceOptions)options;
- (id<MTLCommandBuffer>)commandBuffer;
- (id<MTLComputeCommandEncoder>)computeCommandEncoder;
- (void)setComputePipelineState:(id<MTLComputePipelineState>)pipeline;
- (void)setBuffer:(id<MTLBuffer>)buffer offset:(NSUInteger)offset atIndex:(NSUInteger)index;
- (void)setBytes:(const void*)bytes length:(NSUInteger)length atIndex:(NSUInteger)index;
- (void)dispatchThreads:(MTLSize)grid threadsPerThreadgroup:(MTLSize)group;
- (void)dispatchThreadgroups:(MTLSize)grid threadsPerThreadgroup:(MTLSize)group;
- (void)endEncoding;
- (void)commit;
- (void)waitUntilCompleted;
@property(readonly) void* contents;
@property(readonly) NSUInteger maxTotalThreadsPerThreadgroup;
@property(readonly) MTLCommandBufferStatus status;
@property(readonly) NSString* name;
@end

@implementation TestMetalObject
- (instancetype)initWithBytes:(const void*)bytes length:(NSUInteger)length {
    self = [super init];
    if (self) { storage_ = [NSMutableData dataWithLength:length]; if (bytes) std::memcpy(storage_.mutableBytes, bytes, length); }
    return self;
}
- (id<MTLCommandQueue>)newCommandQueue { return (id<MTLCommandQueue>)self; }
- (id<MTLLibrary>)newLibraryWithSource:(NSString*)source options:(MTLCompileOptions*)options error:(NSError**)error {
    if (!options || options.fastMathEnabled) throw std::runtime_error("Metal fast math must be disabled");
    return failure == "pipeline" ? nil : (id<MTLLibrary>)self;
}
- (id<MTLFunction>)newFunctionWithName:(NSString*)name { return (id<MTLFunction>)self; }
- (id<MTLComputePipelineState>)newComputePipelineStateWithFunction:(id<MTLFunction>)function error:(NSError**)error { return (id<MTLComputePipelineState>)self; }
- (id<MTLBuffer>)newBufferWithBytes:(const void*)bytes length:(NSUInteger)length options:(MTLResourceOptions)options {
    return failure == "buffer" ? nil : (id<MTLBuffer>)[[TestMetalObject alloc] initWithBytes:bytes length:length];
}
- (id<MTLBuffer>)newBufferWithLength:(NSUInteger)length options:(MTLResourceOptions)options {
    return [self newBufferWithBytes:nullptr length:length options:options];
}
- (id<MTLCommandBuffer>)commandBuffer { return failure == "command" ? nil : (id<MTLCommandBuffer>)self; }
- (id<MTLComputeCommandEncoder>)computeCommandEncoder { return failure == "encoder" ? nil : (id<MTLComputeCommandEncoder>)self; }
- (void)setComputePipelineState:(id<MTLComputePipelineState>)pipeline {}
- (void)setBuffer:(id<MTLBuffer>)buffer offset:(NSUInteger)offset atIndex:(NSUInteger)index {}
- (void)setBytes:(const void*)bytes length:(NSUInteger)length atIndex:(NSUInteger)index {}
- (void)dispatchThreads:(MTLSize)grid threadsPerThreadgroup:(MTLSize)group {}
- (void)dispatchThreadgroups:(MTLSize)grid threadsPerThreadgroup:(MTLSize)group {}
- (void)endEncoding {}
- (void)commit {}
- (void)waitUntilCompleted {}
- (void*)contents { return storage_.mutableBytes; }
- (NSUInteger)maxTotalThreadsPerThreadgroup { return 256; }
- (MTLCommandBufferStatus)status { return MTLCommandBufferStatusError; }
- (NSString*)name { return @"test device"; }
@end

static id<MTLDevice> testDevice() { return (id<MTLDevice>)[TestMetalObject new]; }
#define MTLCreateSystemDefaultDevice testDevice
#include "../metal/metal_batch.mm"
#include "../metal/metal_graph_kernels.mm"
#undef MTLCreateSystemDefaultDevice

int main(int argc, char** argv) {
    failure = argc > 1 ? argv[1] : "status";
    if (failure == "disabled") setenv("SMSD_DISABLE_GPU", "1", 1);
    const int label = 6;
    const auto query = smsd::metal::toMetalMolecule(1, &label);
    const std::vector<smsd::metal::MetalMolecule> targets{query, {}};
    if (failure == "disabled" && smsd::metal::metalAvailable()) {
        std::cerr << "GPU disable flag was ignored\n"; return 1;
    }
    const auto hits = smsd::metal::batchSimilarityScreen(query, targets, 0.0f);
    if (hits.size() != 2 || hits[0].index != 0 || hits[0].score != 1.0f || hits[1].score != 0.0f) {
        std::cerr << "Failed Metal screen did not return the CPU result\n"; return 1;
    }
    smsd::metal::MetalMolecule decimalQuery{}, decimalTarget{};
    decimalQuery.n = 5; decimalQuery.labelCounts[0] = 1; decimalQuery.labelCounts[1] = 4;
    decimalTarget.n = 6; decimalTarget.labelCounts[0] = 1; decimalTarget.labelCounts[2] = 5;
    const auto decimalHits = smsd::metal::batchSimilarityScreen(decimalQuery, {decimalTarget}, 0.1f);
    if (decimalHits.size() != 1 || decimalHits[0].score != 1.0f / 10.0f) {
        std::cerr << "The direct float API must compare its rounded float score\n"; return 1;
    }
    smsd::metal::MetalMolecule largeQuery{}, largeTarget{};
    largeQuery.n = (1 << 24) + 1; largeQuery.labelCounts[0] = largeQuery.n;
    largeTarget.n = (1 << 24) + 3; largeTarget.labelCounts[0] = largeTarget.n;
    const float expectedFloat = static_cast<float>(largeQuery.n) / static_cast<float>(largeTarget.n);
    const auto largeHits = smsd::metal::batchSimilarityScreen(largeQuery, {largeTarget}, 0.0f);
    if (expectedFloat == static_cast<float>(static_cast<double>(largeQuery.n) / largeTarget.n)
        || largeHits.size() != 1 || largeHits[0].score != expectedFloat) {
        std::cerr << "The direct float API must preserve float operand division\n"; return 1;
    }
    const int atomicNumber = 6, zero = 0;
    uint64_t domain = 0;
    if (smsd::gpu_kern::domainInit(1, 1, 1, &atomicNumber, &zero, &zero, &zero,
                                  &atomicNumber, &zero, &zero, &zero, true, false, false, false, &domain)) {
        std::cerr << "Failed Metal command was reported as a valid domain\n"; return 1;
    }
    std::cout << "PASS simulated Metal " << failure << " failure fallback\n";
}
