/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "smsd/metal_gpu.hpp"
#include "smsd/gpu_kernels.hpp"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <vector>

#ifndef SMSD_METAL_SMSD_VERSION
#define SMSD_METAL_SMSD_VERSION "unknown"
#endif

namespace smsd {
namespace metal {
namespace fingerprint_detail {

/** A pair of 32-bit words represents each host uint64_t without 64-bit popcounts. */
static NSString* const shaderSource = @
"#include <metal_stdlib>\n"
"using namespace metal;\n"
"struct Shape { uint targets; uint words; uint queries; uint stride; };\n"
"kernel void fingerprint_screen(\n"
"    device const uint2* targets [[buffer(0)]],\n"
"    device const int* popcounts [[buffer(1)]],\n"
"    device const uint2* queries [[buffer(2)]],\n"
"    device const int* cutoffs [[buffer(3)]],\n"
"    device int* selected [[buffer(4)]],\n"
"    constant Shape& shape [[buffer(5)]],\n"
"    uint2 position [[thread_position_in_grid]]) {\n"
"    uint target = position.x, query = position.y;\n"
"    if (target >= shape.targets || query >= shape.queries) return;\n"
"    int common = 0;\n"
"    for (uint word = 0; word < shape.words; ++word) {\n"
"        uint2 both = targets[ulong(word) * shape.targets + target]\n"
"                   & queries[ulong(query) * shape.words + word];\n"
"        common += int(popcount(both.x) + popcount(both.y));\n"
"    }\n"
"    int cutoff = cutoffs[ulong(query) * shape.stride + uint(popcounts[target])];\n"
"    selected[ulong(query) * shape.targets + target] = common >= cutoff ? common : -1;\n"
"}\n";

struct Shape {
    uint32_t targets, words, queries, stride;
};
static_assert(sizeof(Shape) == 16, "Metal fingerprint dimensions must match the shader");

void require(bool condition, const char* reason) {
    if (!condition) throw std::invalid_argument(reason);
}

size_t product(size_t first, size_t second) {
    require(first == 0 || second <= std::numeric_limits<size_t>::max() / first,
            "Metal fingerprint buffer size overflows");
    return first * second;
}

void message(char* output, size_t size, const char* text) noexcept {
    if (output != nullptr && size != 0) std::snprintf(output, size, "%s", text != nullptr ? text : "Metal screening failed");
}

template<class Body>
int call(char* error, size_t errorSize, Body body) noexcept {
    @autoreleasepool {
        try {
            @try {
                body();
                message(error, errorSize, "");
                return 0;
            } @catch (NSException* exception) {
                message(error, errorSize, exception.reason.UTF8String);
                return 1;
            }
        } catch (const std::exception& exception) {
            message(error, errorSize, exception.what());
        } catch (...) {
            message(error, errorSize, "Metal screening failed");
        }
    }
    return 1;
}

/** The pipeline is created only for a nonempty resident target table. */
struct Context {
    id<MTLDevice> device = nil;
    id<MTLCommandQueue> queue = nil;
    id<MTLComputePipelineState> pipeline = nil;
    std::string failure;

    void init() {
        @autoreleasepool {
            device = MTLCreateSystemDefaultDevice();
            if (!device) { failure = "No Metal device is available"; return; }
            queue = [device newCommandQueue];
            if (!queue) { failure = "Metal command queue creation failed"; return; }
            NSError* error = nil;
            MTLCompileOptions* options = [MTLCompileOptions new];
            options.fastMathEnabled = NO;
            id<MTLLibrary> library = [device newLibraryWithSource:shaderSource options:options error:&error];
            if (!library) { failure = "Metal fingerprint shader compilation failed"; return; }
            id<MTLFunction> function = [library newFunctionWithName:@"fingerprint_screen"];
            if (!function) { failure = "Metal fingerprint function is unavailable"; return; }
            pipeline = [device newComputePipelineStateWithFunction:function error:&error];
            if (!pipeline) failure = "Metal fingerprint pipeline creation failed";
        }
    }

    static Context& get() {
        require(!gpu_kern::detail::disabledByEnvironment() && metalAvailable(), "Metal screening is unavailable or disabled");
        static Context context;
        static std::once_flag once;
        std::call_once(once, [&] { context.init(); });
        require(context.pipeline != nil, context.failure.c_str());
        return context;
    }
};

/** Shared buffers grow only when a later query batch needs more capacity. */
struct Buffer {
    id<MTLBuffer> value = nil;
    size_t capacity = 0;

    void ensure(id<MTLDevice> device, size_t bytes) {
        if (bytes <= capacity) return;
        require(bytes <= device.maxBufferLength, "Metal fingerprint buffer exceeds the device limit");
        id<MTLBuffer> allocated = [device newBufferWithLength:bytes options:MTLResourceStorageModeShared];
        require(allocated != nil && allocated.contents != nullptr, "Metal shared-buffer allocation failed");
        value = allocated;
        capacity = bytes;
    }
};

struct Resident {
    int32_t targets, words, maximumPopcount = 0;
    Buffer fingerprints, popcounts, queries, cutoffs, selected;
    std::mutex mutex;
};

int deviceInfo(char* text, size_t textSize, char* error, size_t errorSize) {
    return call(error, errorSize, [&] {
        require(metalAvailable(), "Metal screening is unavailable or disabled");
        message(text, textSize, metalDeviceInfo().c_str());
    });
}

int labelOverlaps(const int32_t* queryCounts, int32_t, int64_t numTargets, const uint64_t* offsets,
                  const smsd_label_count* entries, int32_t* overlaps, char* error, size_t errorSize) {
    for (int64_t target = 0; target < numTargets; ++target) {
        int32_t overlap = 0;
        for (uint64_t index = offsets[target]; index < offsets[target + 1]; ++index)
            overlap += std::min<int32_t>(entries[index].count, queryCounts[entries[index].label]);
        overlaps[target] = overlap;
    }
    message(error, errorSize, "");
    return 0;
}

int minimumTargets(uint32_t, uint64_t* minimum, char* error, size_t errorSize) {
    return call(error, errorSize, [&] {
        require(minimum != nullptr, "Metal screening needs a crossover output");
        *minimum = UINT64_MAX;
    });
}

int createFingerprints(int64_t numTargets, int32_t numWords, const uint64_t* wordMajor,
                       const int32_t* targetPopcounts, void** output, char* error, size_t errorSize) {
    if (output != nullptr) *output = nullptr;
    return call(error, errorSize, [&] {
        require(output != nullptr && numTargets >= 0 && numTargets <= INT32_MAX
                    && numWords > 0 && numWords <= INT32_MAX / 64, "Invalid Metal fingerprint dimensions");
        require(numTargets == 0 || (wordMajor != nullptr && targetPopcounts != nullptr), "Metal fingerprints need words and popcounts");
        auto resident = std::make_unique<Resident>();
        resident->targets = static_cast<int32_t>(numTargets);
        resident->words = numWords;
        for (int64_t target = 0; target < numTargets; ++target) {
            require(targetPopcounts[target] >= 0 && targetPopcounts[target] <= numWords * 64,
                    "Metal target popcount exceeds its fingerprint capacity");
            resident->maximumPopcount = std::max(resident->maximumPopcount, targetPopcounts[target]);
        }
        if (numTargets > 0) {
            auto& context = Context::get();
            const size_t bytes = product(product(static_cast<size_t>(numTargets), static_cast<size_t>(numWords)), sizeof(uint64_t));
            resident->fingerprints.ensure(context.device, bytes);
            resident->popcounts.ensure(context.device, product(static_cast<size_t>(numTargets), sizeof(int32_t)));
            std::memcpy(resident->fingerprints.value.contents, wordMajor, bytes);
            std::memcpy(resident->popcounts.value.contents, targetPopcounts, static_cast<size_t>(numTargets) * sizeof(int32_t));
        }
        *output = resident.release();
    });
}

int fingerprintScreen(void* handle, int32_t numQueries, const uint64_t* queryWords,
                      const int32_t* minimumCommon, int32_t cutoffStride, smsd_fingerprint_hit* hits,
                      uint64_t hitCapacity, uint64_t* hitOffsets, char* error, size_t errorSize) {
    return call(error, errorSize, [&] {
        require(handle != nullptr && numQueries >= 0 && hitOffsets != nullptr, "Invalid Metal fingerprint query dimensions");
        auto& resident = *static_cast<Resident*>(handle);
        const std::lock_guard<std::mutex> lock(resident.mutex);
        if (numQueries == 0 || resident.targets == 0) {
            std::fill(hitOffsets, hitOffsets + static_cast<size_t>(numQueries) + 1, 0);
            return;
        }
        require(queryWords != nullptr && minimumCommon != nullptr && cutoffStride > resident.maximumPopcount,
                "Metal fingerprint cutoff table is too short");
        auto& context = Context::get();
        const size_t pairs = product(static_cast<size_t>(numQueries), static_cast<size_t>(resident.targets));
        const size_t queryBytes = product(product(static_cast<size_t>(numQueries), static_cast<size_t>(resident.words)), sizeof(uint64_t));
        const size_t cutoffBytes = product(product(static_cast<size_t>(numQueries), static_cast<size_t>(cutoffStride)), sizeof(int32_t));
        resident.queries.ensure(context.device, queryBytes);
        resident.cutoffs.ensure(context.device, cutoffBytes);
        resident.selected.ensure(context.device, product(pairs, sizeof(int32_t)));
        std::memcpy(resident.queries.value.contents, queryWords, queryBytes);
        std::memcpy(resident.cutoffs.value.contents, minimumCommon, cutoffBytes);
        id<MTLCommandBuffer> command = [context.queue commandBuffer];
        require(command != nil, "Metal command buffer creation failed");
        id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
        require(encoder != nil, "Metal compute encoder creation failed");
        [encoder setComputePipelineState:context.pipeline];
        [encoder setBuffer:resident.fingerprints.value offset:0 atIndex:0];
        [encoder setBuffer:resident.popcounts.value offset:0 atIndex:1];
        [encoder setBuffer:resident.queries.value offset:0 atIndex:2];
        [encoder setBuffer:resident.cutoffs.value offset:0 atIndex:3];
        [encoder setBuffer:resident.selected.value offset:0 atIndex:4];
        const Shape shape{static_cast<uint32_t>(resident.targets), static_cast<uint32_t>(resident.words),
                          static_cast<uint32_t>(numQueries), static_cast<uint32_t>(cutoffStride)};
        [encoder setBytes:&shape length:sizeof shape atIndex:5];
        NSUInteger threads = std::max<NSUInteger>(1, std::min<NSUInteger>(256, context.pipeline.maxTotalThreadsPerThreadgroup));
        [encoder dispatchThreadgroups:MTLSizeMake((static_cast<NSUInteger>(resident.targets) + threads - 1) / threads,
                                                 static_cast<NSUInteger>(numQueries), 1)
                 threadsPerThreadgroup:MTLSizeMake(threads, 1, 1)];
        [encoder endEncoding];
        [command commit];
        [command waitUntilCompleted];
        require(command.status == MTLCommandBufferStatusCompleted, "Metal fingerprint command failed");
        const auto* selected = static_cast<const int32_t*>(resident.selected.value.contents);
        require(selected != nullptr, "Metal fingerprint result storage is unavailable");
        if (hits != nullptr && hitCapacity >= pairs
                && pairs <= std::numeric_limits<size_t>::max() / sizeof(smsd_fingerprint_hit)) {
            size_t output = 0;
            hitOffsets[0] = 0;
            for (int32_t query = 0; query < numQueries; ++query) {
                for (int32_t target = 0; target < resident.targets; ++target) {
                    const int32_t common = selected[static_cast<size_t>(query) * resident.targets + target];
                    if (common >= 0) hits[output++] = {target, common};
                }
                hitOffsets[query + 1] = output;
            }
            return;
        }
        std::vector<uint64_t> offsets(static_cast<size_t>(numQueries) + 1, 0);
        for (int32_t query = 0; query < numQueries; ++query) {
            offsets[query + 1] = offsets[query];
            for (int32_t target = 0; target < resident.targets; ++target)
                if (selected[static_cast<size_t>(query) * resident.targets + target] >= 0) ++offsets[query + 1];
        }
        require(offsets.back() <= hitCapacity && offsets.back() <= std::numeric_limits<size_t>::max() / sizeof(smsd_fingerprint_hit)
                    && (offsets.back() == 0 || hits != nullptr), "Metal fingerprint output buffer is too short");
        for (int32_t query = 0; query < numQueries; ++query) {
            size_t output = static_cast<size_t>(offsets[query]);
            for (int32_t target = 0; target < resident.targets; ++target) {
                const int32_t common = selected[static_cast<size_t>(query) * resident.targets + target];
                if (common >= 0) hits[output++] = {target, common};
            }
        }
        std::copy(offsets.begin(), offsets.end(), hitOffsets);
    });
}

int destroyResident(void* handle, char* error, size_t errorSize) {
    return call(error, errorSize, [&] { delete static_cast<Resident*>(handle); });
}

const smsd_screen_accelerator_v1 table = {
    SMSD_SCREEN_ACCELERATOR_ABI, static_cast<uint32_t>(sizeof(smsd_screen_accelerator_v1)),
    SMSD_METAL_SMSD_VERSION, "Metal", &deviceInfo, &labelOverlaps, &minimumTargets,
    &createFingerprints, &fingerprintScreen, nullptr, nullptr, &destroyResident
};

} // namespace fingerprint_detail

const smsd_screen_accelerator_v1* fingerprintAccelerator() noexcept {
    return &fingerprint_detail::table;
}

} // namespace metal
} // namespace smsd
