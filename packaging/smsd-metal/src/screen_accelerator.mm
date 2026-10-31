/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * The screening table of the smsd-metal package, which smsd attaches to screen TargetCorpus
 * libraries on an Apple GPU. Fingerprint screens use the resident Metal screen of the SMSD
 * sources (metal::fingerprintAccelerator); label overlaps run in the kernel below. Both return
 * exact integer intersections and overlaps, from which smsd forms every score with its CPU
 * expressions, so every hit is the hit of the CPU screen.
 */
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "smsd/cuda_device_policy.hpp"
#include "smsd/metal_gpu.hpp"
#include "smsd/screen_accelerator_abi.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>

#ifndef SMSD_METAL_SMSD_VERSION
#error "SMSD_METAL_SMSD_VERSION must name the smsd version the package serves"
#endif

#define SMSD_METAL_EXPORT extern "C" __attribute__((visibility("default")))

namespace {

/** One thread per (target, query) pair; entries are (label, count) pairs of 16-bit integers. */
NSString* const labelShaderSource = @
"#include <metal_stdlib>\n"
"using namespace metal;\n"
"struct Shape { uint targets; uint queries; uint table; uint unused; };\n"
"kernel void label_overlaps(\n"
"    device const ulong* offsets [[buffer(0)]],\n"
"    device const ushort2* entries [[buffer(1)]],\n"
"    device const int* counts [[buffer(2)]],\n"
"    device int* overlaps [[buffer(3)]],\n"
"    constant Shape& shape [[buffer(4)]],\n"
"    uint2 position [[thread_position_in_grid]]) {\n"
"    uint target = position.x, query = position.y;\n"
"    if (target >= shape.targets || query >= shape.queries) return;\n"
"    device const int* row = counts + ulong(query) * shape.table;\n"
"    int overlap = 0;\n"
"    for (ulong entry = offsets[target], end = offsets[target + 1]; entry < end; ++entry) {\n"
"        ushort2 item = entries[entry];\n"
"        if (uint(item.x) < shape.table) overlap += min(int(item.y), row[item.x]);\n"
"    }\n"
"    overlaps[ulong(query) * shape.targets + target] = overlap;\n"
"}\n";

struct Shape {
    uint32_t targets, queries, table, unused;
};
static_assert(sizeof(Shape) == 16, "the label kernel dimensions must match the shader");
static_assert(sizeof(smsd_label_count) == 4 && offsetof(smsd_label_count, label) == 0
                  && offsetof(smsd_label_count, count) == 2,
              "label entries must have the layout the shader reads");

/** The largest label table: labels are 16-bit, so a longer table holds nothing more. */
constexpr int32_t maximumTable = 1 << 16;

/** The automatic fingerprint-screening crossover measured on an Apple M5, and that device's name. */
constexpr uint64_t measuredMinimumTargets = 100000;
constexpr const char* measuredDevice = "Apple M5";

void require(bool condition, const char* reason) {
    if (!condition) throw std::invalid_argument(reason);
}

size_t product(size_t first, size_t second) {
    require(first == 0 || second <= std::numeric_limits<size_t>::max() / first,
            "Metal label buffer size overflows");
    return first * second;
}

void message(char* output, size_t size, const char* text) noexcept {
    if (output != nullptr && size != 0) std::snprintf(output, size, "%s", text != nullptr ? text : "Metal screening failed");
}

/** Runs body, turning an exception into a non-zero status and its reason. */
template<class Body>
int guarded(char* error, size_t errorSize, Body body) noexcept {
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

/** Whether SMSD_DISABLE_GPU keeps every screen on the CPU, by the rule smsd applies. */
bool disabled() noexcept {
    return smsd::cuda_policy::disabled(std::getenv("SMSD_DISABLE_GPU"));
}

/** The system GPU, probed once; nil without one. */
id<MTLDevice> systemDevice() {
    static id<MTLDevice> device = nil;
    static std::once_flag once;
    std::call_once(once, [] {
#if defined(__arm64__) || defined(__aarch64__)
        @autoreleasepool { device = MTLCreateSystemDefaultDevice(); }
#endif
    });
    return device;
}

/** The label pipeline, compiled on first use. */
struct LabelContext {
    id<MTLDevice> device = nil;
    id<MTLCommandQueue> queue = nil;
    id<MTLComputePipelineState> pipeline = nil;
    std::string failure;

    void init() {
        @autoreleasepool {
            device = systemDevice();
            if (!device) { failure = "No Metal device is available"; return; }
            queue = [device newCommandQueue];
            if (!queue) { failure = "Metal command queue creation failed"; return; }
            NSError* error = nil;
            MTLCompileOptions* options = [MTLCompileOptions new];
            options.fastMathEnabled = NO;
            id<MTLLibrary> library = [device newLibraryWithSource:labelShaderSource options:options error:&error];
            if (!library) { failure = "Metal label shader compilation failed"; return; }
            id<MTLFunction> function = [library newFunctionWithName:@"label_overlaps"];
            if (!function) { failure = "Metal label function is unavailable"; return; }
            pipeline = [device newComputePipelineStateWithFunction:function error:&error];
            if (!pipeline) failure = "Metal label pipeline creation failed";
        }
    }

    static LabelContext& get() {
        require(!disabled(), "Metal screening is disabled by SMSD_DISABLE_GPU");
        static LabelContext context;
        static std::once_flag once;
        std::call_once(once, [&] { context.init(); });
        require(context.pipeline != nil, context.failure.c_str());
        return context;
    }
};

/** A shared buffer that grows only when a later batch needs more; never empty. */
struct Buffer {
    id<MTLBuffer> value = nil;
    size_t capacity = 0;

    void ensure(id<MTLDevice> device, size_t bytes) {
        bytes = std::max<size_t>(bytes, sizeof(int32_t));
        if (bytes <= capacity) return;
        require(bytes <= device.maxBufferLength, "Metal label buffer exceeds the device limit");
        id<MTLBuffer> allocated = [device newBufferWithLength:bytes options:MTLResourceStorageModeShared];
        require(allocated != nil && allocated.contents != nullptr, "Metal shared-buffer allocation failed");
        value = allocated;
        capacity = bytes;
    }
};

const smsd_screen_accelerator_v1& fingerprints() noexcept {
    return *smsd::metal::fingerprintAccelerator();
}

/** A resident handle: the SMSD fingerprint table, or the label histograms held here. */
struct Resident {
    enum class Kind { Fingerprints, Labels };

    explicit Resident(Kind held) : kind(held) {}

    const Kind kind;
    void* fingerprintTable = nullptr;
    int32_t targets = 0;
    Buffer offsets, entries, counts, overlaps;
    std::mutex mutex;
};

Resident& residentOf(void* handle, Resident::Kind kind) {
    require(handle != nullptr && static_cast<Resident*>(handle)->kind == kind,
            kind == Resident::Kind::Labels ? "the handle holds no Metal label table"
                                           : "the handle holds no Metal fingerprint table");
    return *static_cast<Resident*>(handle);
}

int deviceInfo(char* text, size_t textSize, char* error, size_t errorSize) {
    return guarded(error, errorSize, [&] {
        id<MTLDevice> device = systemDevice();
        require(device != nil, "No Metal device is available");
        const char* name = device.name.UTF8String;
        if (text != nullptr && textSize != 0)
            std::snprintf(text, textSize, "%s (Metal, %llu MB)", name != nullptr ? name : "Apple GPU",
                          static_cast<unsigned long long>(device.recommendedMaxWorkingSetSize / (1024 * 1024)));
    });
}

/**
 * Automatic fingerprint screening uses an Apple M5 from 100,000 targets. SMSD_METAL_MIN_TARGETS,
 * read as SMSD_CUDA_MIN_TARGETS is, extends it to any device from N targets. Label and substructure
 * screens have no measured crossover and are used only when a caller asks for the device.
 */
int minimumTargets(uint32_t operation, uint64_t* minimum, char* error, size_t errorSize) {
    return guarded(error, errorSize, [&] {
        require(minimum != nullptr && (operation == 1 || operation == 2 || operation == 3), "unknown screening operation");
        *minimum = smsd::cuda_policy::kNoAutomaticUse;
        if (operation != 1 || disabled()) return;
        const uint64_t requested = smsd::cuda_policy::requestedMinimumTargets(std::getenv("SMSD_METAL_MIN_TARGETS"));
        if (requested != 0) {
            *minimum = requested;
            return;
        }
        id<MTLDevice> device = systemDevice();
        if (device != nil && smsd::cuda_policy::sameText(device.name.UTF8String, measuredDevice))
            *minimum = measuredMinimumTargets;
    });
}

int createFingerprints(int64_t numTargets, int32_t numWords, const uint64_t* wordMajor,
                       const int32_t* targetPopcounts, void** output, char* error, size_t errorSize) {
    if (output != nullptr) *output = nullptr;
    std::unique_ptr<Resident> resident;
    if (int status = guarded(error, errorSize, [&] {
            require(output != nullptr, "Metal fingerprints need a handle output");
            resident = std::make_unique<Resident>(Resident::Kind::Fingerprints);
        })) return status;
    if (int status = fingerprints().create_fingerprints(numTargets, numWords, wordMajor, targetPopcounts,
                                                        &resident->fingerprintTable, error, errorSize)) return status;
    *output = resident.release();
    return 0;
}

int fingerprintScreen(void* handle, int32_t numQueries, const uint64_t* queryWords,
                      const int32_t* minimumCommon, int32_t cutoffStride, smsd_fingerprint_hit* hits,
                      uint64_t hitCapacity, uint64_t* hitOffsets, char* error, size_t errorSize) {
    Resident* resident = nullptr;
    if (int status = guarded(error, errorSize, [&] {
            resident = &residentOf(handle, Resident::Kind::Fingerprints);
        })) return status;
    return fingerprints().fingerprint_screen(resident->fingerprintTable, numQueries, queryWords, minimumCommon,
                                             cutoffStride, hits, hitCapacity, hitOffsets, error, errorSize);
}

int createLabels(int64_t numTargets, const uint64_t* offsets, const smsd_label_count* entries,
                 void** output, char* error, size_t errorSize) {
    if (output != nullptr) *output = nullptr;
    return guarded(error, errorSize, [&] {
        require(output != nullptr && numTargets >= 0 && numTargets <= std::numeric_limits<int32_t>::max()
                    && offsets != nullptr, "Invalid Metal label dimensions");
        require(offsets[0] == 0, "Metal label offsets must start at zero");
        for (int64_t target = 0; target < numTargets; ++target)
            require(offsets[target] <= offsets[target + 1], "Metal label offsets are not ordered");
        const uint64_t entryCount = offsets[numTargets];
        require(entryCount == 0 || entries != nullptr, "Metal labels need histogram entries");
        require(entryCount <= std::numeric_limits<size_t>::max() / sizeof(smsd_label_count), "Metal label buffer size overflows");
        auto resident = std::make_unique<Resident>(Resident::Kind::Labels);
        resident->targets = static_cast<int32_t>(numTargets);
        if (numTargets > 0) {
            auto& context = LabelContext::get();
            const size_t offsetBytes = product(static_cast<size_t>(numTargets) + 1, sizeof(uint64_t));
            const size_t entryBytes = product(static_cast<size_t>(entryCount), sizeof(smsd_label_count));
            resident->offsets.ensure(context.device, offsetBytes);
            resident->entries.ensure(context.device, entryBytes);
            std::memcpy(resident->offsets.value.contents, offsets, offsetBytes);
            if (entryBytes != 0) std::memcpy(resident->entries.value.contents, entries, entryBytes);
        }
        *output = resident.release();
    });
}

int residentLabelOverlaps(void* handle, int32_t numQueries, const int32_t* queryCounts,
                          int32_t tableSize, int32_t* overlaps, char* error, size_t errorSize) {
    return guarded(error, errorSize, [&] {
        Resident& resident = residentOf(handle, Resident::Kind::Labels);
        require(numQueries >= 0 && tableSize > 0 && tableSize <= maximumTable, "Invalid Metal label query dimensions");
        const std::lock_guard<std::mutex> lock(resident.mutex);
        if (numQueries == 0 || resident.targets == 0) return;
        require(queryCounts != nullptr && overlaps != nullptr, "Metal label screening needs query and output buffers");
        auto& context = LabelContext::get();
        const size_t countBytes = product(product(static_cast<size_t>(numQueries), static_cast<size_t>(tableSize)), sizeof(int32_t));
        const size_t overlapBytes = product(product(static_cast<size_t>(numQueries), static_cast<size_t>(resident.targets)), sizeof(int32_t));
        resident.counts.ensure(context.device, countBytes);
        resident.overlaps.ensure(context.device, overlapBytes);
        std::memcpy(resident.counts.value.contents, queryCounts, countBytes);
        id<MTLCommandBuffer> command = [context.queue commandBuffer];
        require(command != nil, "Metal command buffer creation failed");
        id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
        require(encoder != nil, "Metal compute encoder creation failed");
        [encoder setComputePipelineState:context.pipeline];
        [encoder setBuffer:resident.offsets.value offset:0 atIndex:0];
        [encoder setBuffer:resident.entries.value offset:0 atIndex:1];
        [encoder setBuffer:resident.counts.value offset:0 atIndex:2];
        [encoder setBuffer:resident.overlaps.value offset:0 atIndex:3];
        const Shape shape{static_cast<uint32_t>(resident.targets), static_cast<uint32_t>(numQueries),
                          static_cast<uint32_t>(tableSize), 0};
        [encoder setBytes:&shape length:sizeof shape atIndex:4];
        const NSUInteger threads = std::max<NSUInteger>(1, std::min<NSUInteger>(256, context.pipeline.maxTotalThreadsPerThreadgroup));
        [encoder dispatchThreadgroups:MTLSizeMake((static_cast<NSUInteger>(resident.targets) + threads - 1) / threads,
                                                 static_cast<NSUInteger>(numQueries), 1)
                 threadsPerThreadgroup:MTLSizeMake(threads, 1, 1)];
        [encoder endEncoding];
        [command commit];
        [command waitUntilCompleted];
        require(command.status == MTLCommandBufferStatusCompleted, "Metal label command failed");
        std::memcpy(overlaps, resident.overlaps.value.contents, overlapBytes);
    });
}

int destroyResident(void* handle, char* error, size_t errorSize) {
    if (handle == nullptr) return guarded(error, errorSize, [] {});
    std::unique_ptr<Resident> resident(static_cast<Resident*>(handle));
    if (resident->kind == Resident::Kind::Fingerprints && resident->fingerprintTable != nullptr)
        return fingerprints().destroy_resident(resident->fingerprintTable, error, errorSize);
    return guarded(error, errorSize, [] {});
}

/** One query through a label table held for the call. */
int labelOverlaps(const int32_t* queryCounts, int32_t tableSize, int64_t numTargets, const uint64_t* offsets,
                  const smsd_label_count* entries, int32_t* overlaps, char* error, size_t errorSize) {
    void* handle = nullptr;
    if (int status = createLabels(numTargets, offsets, entries, &handle, error, errorSize)) return status;
    const int status = residentLabelOverlaps(handle, 1, queryCounts, tableSize, overlaps, error, errorSize);
    char ignored[8] = {};
    destroyResident(handle, ignored, sizeof ignored);
    return status;
}

const smsd_screen_accelerator_v1 table = {
    SMSD_SCREEN_ACCELERATOR_ABI,
    static_cast<uint32_t>(sizeof(smsd_screen_accelerator_v1)),
    SMSD_METAL_SMSD_VERSION,
    "Metal",
    &deviceInfo,
    &labelOverlaps,
    &minimumTargets,
    &createFingerprints,
    &fingerprintScreen,
    &createLabels,
    &residentLabelOverlaps,
    &destroyResident,
};

} // namespace

/** The table smsd attaches; valid while the library stays loaded. */
SMSD_METAL_EXPORT const smsd_screen_accelerator_v1* smsd_metal_screen_accelerator(void) { return &table; }

/** 1 when an Apple GPU can run the screens; 0 with the reason in error. */
SMSD_METAL_EXPORT int smsd_metal_device_available(char* error, size_t errorSize) {
#if defined(__arm64__) || defined(__aarch64__)
    if (systemDevice() == nil) {
        message(error, errorSize, "no Metal device was found");
        return 0;
    }
    message(error, errorSize, "");
    return 1;
#else
    message(error, errorSize, "Metal screening needs an Apple silicon Mac");
    return 0;
#endif
}

/** Writes a description of the device; 0 on success. */
SMSD_METAL_EXPORT int smsd_metal_device_info(char* text, size_t textSize, char* error, size_t errorSize) {
    return deviceInfo(text, textSize, error, errorSize);
}

/** Compiles both screening kernels without screening data. */
SMSD_METAL_EXPORT int smsd_metal_preload(char* error, size_t errorSize) {
    if (int status = guarded(error, errorSize, [] { LabelContext::get(); })) return status;
    // A one-target table makes the fingerprint screen compile its kernel; it is released at once.
    const uint64_t words[1] = {1};
    const int32_t popcounts[1] = {1};
    void* handle = nullptr;
    if (int status = fingerprints().create_fingerprints(1, 1, words, popcounts, &handle, error, errorSize)) return status;
    return fingerprints().destroy_resident(handle, error, errorSize);
}
