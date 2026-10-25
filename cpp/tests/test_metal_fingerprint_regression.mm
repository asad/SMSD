/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "smsd/metal_gpu.hpp"
#include "smsd/batch.hpp"
#include "smsd/screen_accelerator.hpp"
#include <array>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <thread>

static std::string failure;
static std::atomic<int> allocations{0}, liveBuffers{0}, deviceProbes{0}, dispatches{0};

static void check(bool condition, const char* reason) {
    if (!condition) throw std::runtime_error(reason);
}

/** Protocol double with a CPU intersection oracle; no physical device is used. */
@interface FingerprintMetalObject : NSObject {
    NSMutableData* storage_;
    NSMutableArray* buffers_;
    std::array<uint32_t, 4> shape_;
    bool buffer_;
}
- (instancetype)initWithLength:(NSUInteger)length buffer:(bool)buffer;
- (id<MTLCommandQueue>)newCommandQueue;
- (id<MTLLibrary>)newLibraryWithSource:(NSString*)source options:(MTLCompileOptions*)options error:(NSError**)error;
- (id<MTLFunction>)newFunctionWithName:(NSString*)name;
- (id<MTLComputePipelineState>)newComputePipelineStateWithFunction:(id<MTLFunction>)function error:(NSError**)error;
- (id<MTLBuffer>)newBufferWithLength:(NSUInteger)length options:(MTLResourceOptions)options;
- (id<MTLCommandBuffer>)commandBuffer;
- (id<MTLComputeCommandEncoder>)computeCommandEncoder;
- (void)setComputePipelineState:(id<MTLComputePipelineState>)pipeline;
- (void)setBuffer:(id<MTLBuffer>)buffer offset:(NSUInteger)offset atIndex:(NSUInteger)index;
- (void)setBytes:(const void*)bytes length:(NSUInteger)length atIndex:(NSUInteger)index;
- (void)dispatchThreadgroups:(MTLSize)grid threadsPerThreadgroup:(MTLSize)group;
- (void)endEncoding;
- (void)commit;
- (void)waitUntilCompleted;
@property(readonly) void* contents;
@property(readonly) NSUInteger maxBufferLength;
@property(readonly) NSUInteger maxTotalThreadsPerThreadgroup;
@property(readonly) MTLCommandBufferStatus status;
@end

@implementation FingerprintMetalObject
- (instancetype)initWithLength:(NSUInteger)length buffer:(bool)buffer {
    self = [super init];
    if (self) {
        storage_ = [NSMutableData dataWithLength:length];
        buffers_ = [NSMutableArray arrayWithArray:@[[NSNull null], [NSNull null], [NSNull null], [NSNull null], [NSNull null]]];
        buffer_ = buffer;
        if (buffer_) ++liveBuffers;
    }
    return self;
}
- (void)dealloc { if (buffer_) --liveBuffers; }
- (id<MTLCommandQueue>)newCommandQueue { return failure == "queue" ? nil : (id<MTLCommandQueue>)self; }
- (id<MTLLibrary>)newLibraryWithSource:(NSString*)source options:(MTLCompileOptions*)options error:(NSError**)error {
    check([source containsString:@"kernel void fingerprint_screen"], "Wrong shader source");
    check(options && !options.fastMathEnabled, "Metal fast math must be disabled");
    return failure == "pipeline" ? nil : (id<MTLLibrary>)self;
}
- (id<MTLFunction>)newFunctionWithName:(NSString*)name {
    check([name isEqualToString:@"fingerprint_screen"], "Wrong kernel name");
    return failure == "function" ? nil : (id<MTLFunction>)self;
}
- (id<MTLComputePipelineState>)newComputePipelineStateWithFunction:(id<MTLFunction>)function error:(NSError**)error {
    return (id<MTLComputePipelineState>)self;
}
- (id<MTLBuffer>)newBufferWithLength:(NSUInteger)length options:(MTLResourceOptions)options {
    check((options & MTLResourceStorageModeMask) == MTLResourceStorageModeShared, "Buffers must use shared storage");
    if (failure == "exception") [NSException raise:@"MetalAllocation" format:@"Simulated allocation exception"];
    if (failure == "buffer") return nil;
    ++allocations;
    return (id<MTLBuffer>)[[FingerprintMetalObject alloc] initWithLength:length buffer:true];
}
- (id<MTLCommandBuffer>)commandBuffer {
    return failure == "command" ? nil : (id<MTLCommandBuffer>)[[FingerprintMetalObject alloc] initWithLength:0 buffer:false];
}
- (id<MTLComputeCommandEncoder>)computeCommandEncoder { return failure == "encoder" ? nil : (id<MTLComputeCommandEncoder>)self; }
- (void)setComputePipelineState:(id<MTLComputePipelineState>)pipeline {}
- (void)setBuffer:(id<MTLBuffer>)buffer offset:(NSUInteger)offset atIndex:(NSUInteger)index {
    check(index < 5 && offset == 0, "Unexpected buffer binding");
    buffers_[index] = (id)buffer;
}
- (void)setBytes:(const void*)bytes length:(NSUInteger)length atIndex:(NSUInteger)index {
    check(index == 5 && length == sizeof(shape_), "Unexpected shape binding");
    std::memcpy(shape_.data(), bytes, length);
}
- (void)dispatchThreadgroups:(MTLSize)grid threadsPerThreadgroup:(MTLSize)group {
    check(grid.width * group.width >= shape_[0] && grid.height == shape_[2]
          && group.height == 1 && group.depth == 1, "Wrong batched grid");
    ++dispatches;
}
- (void)endEncoding {}
- (void)commit {
    const auto* targets = static_cast<const uint64_t*>(((id<MTLBuffer>)buffers_[0]).contents);
    const auto* counts = static_cast<const int32_t*>(((id<MTLBuffer>)buffers_[1]).contents);
    const auto* queries = static_cast<const uint64_t*>(((id<MTLBuffer>)buffers_[2]).contents);
    const auto* cutoffs = static_cast<const int32_t*>(((id<MTLBuffer>)buffers_[3]).contents);
    auto* selected = static_cast<int32_t*>(((id<MTLBuffer>)buffers_[4]).contents);
    for (uint32_t query = 0; query < shape_[2]; ++query) {
        for (uint32_t target = 0; target < shape_[0]; ++target) {
            int32_t common = 0;
            for (uint32_t word = 0; word < shape_[1]; ++word)
                common += __builtin_popcountll(targets[size_t(word) * shape_[0] + target]
                                               & queries[size_t(query) * shape_[1] + word]);
            selected[size_t(query) * shape_[0] + target] =
                common >= cutoffs[size_t(query) * shape_[3] + counts[target]] ? common : -1;
        }
    }
}
- (void)waitUntilCompleted {}
- (void*)contents { return failure == "contents" ? nullptr : storage_.mutableBytes; }
- (NSUInteger)maxBufferLength { return NSUInteger(1) << 30; }
- (NSUInteger)maxTotalThreadsPerThreadgroup { return 256; }
- (MTLCommandBufferStatus)status { return failure == "status" ? MTLCommandBufferStatusError : MTLCommandBufferStatusCompleted; }
@end

static id<MTLDevice> fingerprintTestDevice() {
    ++deviceProbes;
    return failure == "device" ? nil : (id<MTLDevice>)[[FingerprintMetalObject alloc] initWithLength:0 buffer:false];
}
#define MTLCreateSystemDefaultDevice fingerprintTestDevice
#include "../metal/metal_fingerprint.mm"
#undef MTLCreateSystemDefaultDevice

namespace smsd {
namespace metal {
bool metalAvailable() noexcept { return !gpu_kern::detail::disabledByEnvironment(); }
std::string metalDeviceInfo() { return "Metal protocol double"; }
}
}

static void checkHostLabelScreens(const smsd_screen_accelerator_v1* api) {
    const int32_t query[]{3, 0, 5, 1};
    const uint64_t offsets[]{0, 0, 2, 4, 6};
    const smsd_label_count entries[]{{0, 2}, {2, 7}, {1, 9}, {3, 1}, {0, 4}, {2, 1}};
    int32_t overlaps[4]{-1, -1, -1, -1};
    char error[256]{};
    check(api->label_overlaps(query, 4, 4, offsets, entries, overlaps, error, sizeof error) == 0,
          "Metal label screens must supply host overlaps");
    const int32_t expected[]{0, 7, 1, 4};
    check(std::equal(std::begin(overlaps), std::end(overlaps), std::begin(expected)),
          "Host label screens must sum the smaller label counts in target order");
    check(api->label_overlaps(nullptr, 0, 0, nullptr, nullptr, nullptr, error, sizeof error) == 0,
          "Empty host label screens must succeed");

    smsd::batch::TargetCorpus corpus;
    corpus.addTargets({smsd::parseSMILES("CCO"), smsd::parseSMILES("N"),
                       smsd::parseSMILES("CCO"), smsd::parseSMILES("CO")});
    const auto molecule = smsd::parseSMILES("CO");
    const smsd::ChemOptions chemistry;
    smsd::screen_accelerator::detach();
    const auto failures = smsd::screen_accelerator::failureCount();
    const auto latest = smsd::screen_accelerator::latestFailure();
    for (double threshold : {0.0, 0.5, 1.0}) {
        const auto cpu = corpus.screenTargets(molecule, threshold, chemistry, 1);
        check(smsd::screen_accelerator::attach(api, api->smsd_version).empty(),
              "Metal fingerprint provider must attach");
        const auto attached = corpus.screenTargets(molecule, threshold, chemistry, 1);
        smsd::screen_accelerator::detach();
        check(attached == cpu, "Attached Metal label screens must preserve CPU hits and target slots");
    }
    check(smsd::screen_accelerator::failureCount() == failures
              && smsd::screen_accelerator::latestFailure() == latest,
          "Supported host label screens must not hide genuine accelerator failures");
}

int main(int argc, char** argv) {
    failure = argc > 1 ? argv[1] : "contract";
    if (failure == "disabled") setenv("SMSD_DISABLE_GPU", "1", 1);
    else unsetenv("SMSD_DISABLE_GPU");
    try {
        @autoreleasepool {
            const auto* api = smsd::metal::fingerprintAccelerator();
            char error[256]{};
            uint64_t minimum = 0;
            check(api->abi_version == SMSD_SCREEN_ACCELERATOR_ABI && api->struct_size == sizeof(*api), "Wrong ABI prefix");
            check(api->minimum_targets(1, &minimum, error, sizeof error) == 0 && minimum == UINT64_MAX, "Unmeasured automatic crossover");
            if (failure == "contract") checkHostLabelScreens(api);
            check(deviceProbes == 0, "Factory or crossover probed the device");
            void* empty = nullptr;
            check(api->create_fingerprints(0, 2, nullptr, nullptr, &empty, error, sizeof error) == 0, "Empty resident creation failed");
            uint64_t emptyOffsets[4]{9, 9, 9, 9};
            check(api->fingerprint_screen(empty, 3, nullptr, nullptr, 0, nullptr, 0, emptyOffsets, error, sizeof error) == 0, "Empty screen failed");
            check(std::all_of(emptyOffsets, emptyOffsets + 4, [](uint64_t value) { return value == 0; }), "Empty offsets were not reset");
            check(api->destroy_resident(empty, error, sizeof error) == 0 && deviceProbes == 0, "Empty screen touched a device");

            uint64_t words[]{0, 3, 1, UINT64_MAX, 5, 0, 0, uint64_t(1) << 63, UINT64_MAX, 10};
            const int32_t counts[]{0, 2, 2, 128, 4};
            void* resident = nullptr;
            const int created = api->create_fingerprints(5, 2, words, counts, &resident, error, sizeof error);
            const bool creationFailure = failure == "queue" || failure == "pipeline" || failure == "function"
                                      || failure == "buffer" || failure == "contents" || failure == "disabled"
                                      || failure == "device" || failure == "exception";
            if (creationFailure) {
                check(created != 0 && resident == nullptr && error[0] != '\0' && liveBuffers == 0, "Creation failure leaked or succeeded");
            } else {
                check(created == 0 && resident != nullptr && liveBuffers == 2, "Resident targets were not retained");
                std::fill(std::begin(words), std::end(words), 0);
                const uint64_t queries[]{0, 0, 3, 10, UINT64_MAX, UINT64_MAX};
                std::array<int32_t, 3 * 129> cutoffs;
                for (int count = 0; count < 129; ++count) {
                    cutoffs[count] = count == 0 ? 0 : 1;
                    cutoffs[129 + count] = 2;
                    cutoffs[258 + count] = count;
                }
                if (failure == "contract") {
                    smsd_fingerprint_hit hit{};
                    uint64_t offsets[2]{};
                    check(api->fingerprint_screen(resident, 1, queries, cutoffs.data(), 129, &hit, 1,
                                                  offsets, error, sizeof error) == 0
                          && offsets[0] == 0 && offsets[1] == 1 && hit.target_index == 0
                          && hit.common_bits == 0, "Single query or zero intersection failed");
                }
                auto screen = [&] {
                    char callError[256]{};
                    smsd_fingerprint_hit hits[15]{};
                    uint64_t offsets[4]{};
                    const int status = api->fingerprint_screen(resident, 3, queries, cutoffs.data(), 129,
                                                               hits, 15, offsets, callError, sizeof callError);
                    if (failure == "command" || failure == "encoder" || failure == "status") {
                        check(status != 0 && callError[0] != '\0', "Command failure was reported as success");
                        return;
                    }
                    check(status == 0, callError);
                    const uint64_t expectedOffsets[]{0, 1, 4, 9};
                    const smsd_fingerprint_hit expected[]{{0, 0}, {1, 2}, {3, 4}, {4, 3},
                                                          {0, 0}, {1, 2}, {2, 2}, {3, 128}, {4, 4}};
                    check(std::equal(offsets, offsets + 4, expectedOffsets), "Wrong ordered query offsets");
                    for (size_t hit = 0; hit < 9; ++hit)
                        check(hits[hit].target_index == expected[hit].target_index
                              && hits[hit].common_bits == expected[hit].common_bits, "Wrong integer hit or target order");
                };
                screen();
                const int allocated = allocations;
                screen();
                check(allocations == allocated, "Warm screen reallocated shared storage");
                if (failure == "concurrent") {
                    std::exception_ptr firstError, secondError;
                    auto repeat = [&](std::exception_ptr& captured) {
                        try { for (int call = 0; call < 10; ++call) screen(); } catch (...) { captured = std::current_exception(); }
                    };
                    std::thread first(repeat, std::ref(firstError)), second(repeat, std::ref(secondError));
                    first.join(); second.join();
                    if (firstError) std::rethrow_exception(firstError);
                    if (secondError) std::rethrow_exception(secondError);
                }
                if (failure == "contract") {
                    uint64_t offsets[4]{77, 77, 77, 77};
                    smsd_fingerprint_hit hits[9]{};
                    check(api->fingerprint_screen(resident, 3, queries, cutoffs.data(), 129, hits, 8,
                                                  offsets, error, sizeof error) != 0 && offsets[0] == 77, "Short output was partly published");
                    check(api->fingerprint_screen(resident, 3, queries, cutoffs.data(), 128, hits, 9,
                                                  offsets, error, sizeof error) != 0, "Short cutoff stride was accepted");
                    check(api->fingerprint_screen(resident, 0, nullptr, nullptr, 0, nullptr, 0,
                                                  offsets, error, sizeof error) == 0 && offsets[0] == 0, "Empty query batch failed");
                }
                check(api->destroy_resident(resident, error, sizeof error) == 0 && liveBuffers == 0, "Resident destruction leaked buffers");
            }
        }
        std::cout << "PASS resident Metal fingerprint " << failure << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
