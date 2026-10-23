/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Integer screening on immutable device-resident fingerprints and label tables.
 * Each corpus owns its allocations. Each call uses its own stream and scratch
 * buffers; the host computes scores and keeps the CPU's ordering and arithmetic.
 */
#include "smsd/screen_accelerator_abi.h"

#include <cuda_runtime.h>
#include <cub/device/device_select.cuh>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <mutex>
#include <new>
#include <stdexcept>
#include <vector>

#if defined(_WIN32)
#define SMSD_CUDA_EXPORT extern "C" __declspec(dllexport)
#else
#define SMSD_CUDA_EXPORT extern "C" __attribute__((visibility("default")))
#endif

#ifndef SMSD_CUDA_SMSD_VERSION
#error "SMSD_CUDA_SMSD_VERSION must name the smsd version this library serves"
#endif

namespace {

constexpr int kThreads = 256;
constexpr int kMaxBlocks = 65535;
constexpr int32_t kMaxTable = 48 * 1024 / static_cast<int32_t>(sizeof(int32_t));

void describe(char* error, size_t size, const char* what, cudaError_t status) {
    if (error != nullptr && size != 0)
        std::snprintf(error, size, "%s: %s", what, cudaGetErrorString(status));
}

/** The first CUDA failure, including failures while releasing a call's resources. */
struct ErrorState {
    cudaError_t code = cudaSuccess;
    const char* operation = nullptr;

    void note(cudaError_t status, const char* what) noexcept {
        if (code == cudaSuccess && status != cudaSuccess) { code = status; operation = what; }
    }

    struct Failure {};

    void check(cudaError_t status, const char* what) {
        note(status, what);
        if (status != cudaSuccess) throw Failure{};
    }
};

template<class Function>
int guarded(char* error, size_t errorSize, Function&& run) noexcept {
    ErrorState state;
    try {
        run(state);
    } catch (const ErrorState::Failure&) {
    } catch (const std::bad_alloc&) {
        if (error != nullptr && errorSize != 0) std::snprintf(error, errorSize, "host allocation failed");
        return 1;
    } catch (const std::exception& failure) {
        if (error != nullptr && errorSize != 0) std::snprintf(error, errorSize, "%s", failure.what());
        return 1;
    } catch (...) {
        if (error != nullptr && errorSize != 0) std::snprintf(error, errorSize, "screening failed");
        return 1;
    }
    if (state.code != cudaSuccess) { describe(error, errorSize, state.operation, state.code); return 1; }
    return 0;
}

size_t bytesFor(uint64_t count, size_t itemSize) {
    if (count > std::numeric_limits<size_t>::max() / itemSize)
        throw std::invalid_argument("screening buffer size exceeds the address space");
    return static_cast<size_t>(count) * itemSize;
}

void require(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(message);
}

#ifdef SMSD_CUDA_TESTING
/** Failure injection exists only in a separately compiled test library. */
thread_local int injectedFailure = 0;

bool failAllocation() noexcept {
    if (injectedFailure == 1 || injectedFailure == -1) {
        injectedFailure = 0;
        return true;
    }
    if (injectedFailure < -1) ++injectedFailure;
    return false;
}
#endif

void beforeWork(ErrorState& state) {
#ifdef SMSD_CUDA_TESTING
    if (injectedFailure == 2) {
        injectedFailure = 0;
        state.check(cudaErrorUnknown, "injected device loss");
    }
#else
    (void)state;
#endif
}

cudaError_t allocateDevice(void** result, size_t bytes) {
#ifdef SMSD_CUDA_TESTING
    if (failAllocation()) { *result = nullptr; return cudaErrorMemoryAllocation; }
#endif
    return cudaMalloc(result, bytes);
}

cudaError_t allocatePinned(void** result, size_t bytes) {
#ifdef SMSD_CUDA_TESTING
    if (failAllocation()) { *result = nullptr; return cudaErrorMemoryAllocation; }
#endif
    return cudaHostAlloc(result, bytes, cudaHostAllocDefault);
}

/** Select the allocation's device on this thread and restore the caller's device. */
class DeviceScope {
public:
    DeviceScope(ErrorState& state, int device) : state_(state) {
        state_.check(cudaGetDevice(&previous_), "cudaGetDevice");
        if (previous_ != device) {
            state_.check(cudaSetDevice(device), "cudaSetDevice");
            restore_ = true;
        }
    }

    ~DeviceScope() {
        if (restore_) state_.note(cudaSetDevice(previous_), "restore CUDA device");
    }

private:
    ErrorState& state_;
    int previous_ = -1;
    bool restore_ = false;
};

/** Four resident staging slots use at most 64 MiB, including replacement allocation. */
class PinnedStaging {
public:
    std::unique_lock<std::mutex> tryAcquire() { return std::unique_lock<std::mutex>(mutex_, std::try_to_lock); }

    void* allocate(ErrorState& state, size_t index, size_t bytes) {
        if (disabled_ || index >= slots_.size()) return nullptr;
        Slot& slot = slots_[index];
        if (slot.capacity >= bytes) return slot.buffer;
        // Keep the old allocation until replacement succeeds; exceeding the cap uses call-owned staging.
        if (bytes > kMaximumBytes - retainedBytes_) return nullptr;
        void* replacement = nullptr;
        state.check(allocatePinned(&replacement, bytes), "cudaHostAlloc retained staging");
        if (releaseSlot(state, slot, "cudaFreeHost replaced staging") != cudaSuccess) {
            // A free may report an earlier asynchronous error; do not reuse or retry uncertain pointers.
            disabled_ = true;
            quarantined_ = {replacement, bytes, false};
            retainedBytes_ += bytes;
            releaseSlot(state, quarantined_, "cudaFreeHost replacement after failed growth");
            throw ErrorState::Failure{};
        }
        retainedBytes_ += bytes - slot.capacity;
        slot = {replacement, bytes, false};
        return replacement;
    }

    bool empty() const noexcept { return retainedBytes_ == 0; }

    void release(ErrorState& state) noexcept {
        for (Slot& slot : slots_) {
            releaseSlot(state, slot, "cudaFreeHost resident staging");
            slot = {};
        }
        releaseSlot(state, quarantined_, "cudaFreeHost quarantined staging");
        quarantined_ = {};
        retainedBytes_ = 0;
    }

private:
    struct Slot { void* buffer = nullptr; size_t capacity = 0; bool releaseAttempted = false; };

    static cudaError_t releaseSlot(ErrorState& state, Slot& slot, const char* operation) noexcept {
        if (slot.buffer == nullptr || slot.releaseAttempted) return cudaSuccess;
        slot.releaseAttempted = true;
        const cudaError_t status = cudaFreeHost(slot.buffer);
        state.note(status, operation);
        if (status == cudaSuccess) slot.buffer = nullptr;
        return status;
    }

    static constexpr size_t kMaximumBytes = size_t(64) * 1024 * 1024;
    std::mutex mutex_;
    std::array<Slot, 4> slots_{};
    Slot quarantined_;
    size_t retainedBytes_ = 0;
    bool disabled_ = false;
};

/** A call owns its stream and device scratch, and may lease resident pinned staging. */
class Call {
public:
    explicit Call(ErrorState& state, PinnedStaging* retained = nullptr)
        : state_(state), retained_(retained), lease_(retained != nullptr ? retained->tryAcquire() : std::unique_lock<std::mutex>()) {
        device_.reserve(12);
        pinned_.reserve(6);
        state_.check(cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking), "cudaStreamCreateWithFlags");
    }

    ~Call() {
        if (pending_) state_.note(cudaStreamSynchronize(stream_), "finish CUDA call");
        for (void* buffer : device_) state_.note(cudaFree(buffer), "cudaFree");
        for (void* buffer : pinned_) state_.note(cudaFreeHost(buffer), "cudaFreeHost");
        state_.note(cudaStreamDestroy(stream_), "cudaStreamDestroy");
    }

    Call(const Call&) = delete;
    Call& operator=(const Call&) = delete;

    template<class T>
    T* device(uint64_t count) { return static_cast<T*>(allocate(count, sizeof(T), false)); }

    template<class T>
    T* pinned(uint64_t count) { return static_cast<T*>(allocate(count, sizeof(T), true)); }

    template<class T>
    void upload(T* destination, const T* source, uint64_t count) {
        if (count == 0) return;
        T* staging = pinned<T>(count);
        std::memcpy(staging, source, bytesFor(count, sizeof(T)));
        copy(destination, staging, bytesFor(count, sizeof(T)), cudaMemcpyHostToDevice);
    }

    void copy(void* destination, const void* source, size_t bytes, cudaMemcpyKind kind) {
        if (bytes == 0) return;
        pending_ = true;
        state_.check(cudaMemcpyAsync(destination, source, bytes, kind, stream_), "cudaMemcpyAsync");
    }

    void launched() { pending_ = true; state_.check(cudaGetLastError(), "kernel launch"); }
    void pending() noexcept { pending_ = true; }
    cudaStream_t stream() const noexcept { return stream_; }

    void finish() {
        state_.check(cudaStreamSynchronize(stream_), "cudaStreamSynchronize");
        pending_ = false;
    }

private:
    void* allocate(uint64_t count, size_t itemSize, bool host) {
        if (count == 0) return nullptr;
        void* buffer = nullptr;
        const size_t bytes = bytesFor(count, itemSize);
        if (host && lease_.owns_lock()) {
            buffer = retained_->allocate(state_, nextPinnedSlot_++, bytes);
            if (buffer != nullptr) return buffer;
        }
        state_.check(host ? allocatePinned(&buffer, bytes) : allocateDevice(&buffer, bytes),
                     host ? "cudaHostAlloc" : "cudaMalloc");
        try {
            (host ? pinned_ : device_).push_back(buffer);
        } catch (...) {
            state_.note(host ? cudaFreeHost(buffer) : cudaFree(buffer), "release unrecorded allocation");
            throw;
        }
        return buffer;
    }

    ErrorState& state_;
    PinnedStaging* retained_;
    std::unique_lock<std::mutex> lease_;
    size_t nextPinnedSlot_ = 0;
    cudaStream_t stream_ = nullptr;
    bool pending_ = false;
    std::vector<void*> device_, pinned_;
};

enum class ResidentKind { Fingerprints, Labels };

/** Immutable device data and bounded host staging, released by the handle's deleter. */
struct Resident {
    ResidentKind kind;
    int device = -1;
    int32_t targets = 0;
    int32_t words = 0;
    int32_t maximumPopcount = 0;
    uint64_t* fingerprints = nullptr;
    int32_t* popcounts = nullptr;
    uint64_t* offsets = nullptr;
    smsd_label_count* entries = nullptr;
    std::vector<void*> allocations;
    PinnedStaging pinned;

    explicit Resident(ResidentKind value) : kind(value) { allocations.reserve(2); }

    template<class T>
    T* allocate(ErrorState& state, uint64_t count) {
        if (count == 0) return nullptr;
        void* buffer = nullptr;
        state.check(allocateDevice(&buffer, bytesFor(count, sizeof(T))), "cudaMalloc resident corpus");
        try { allocations.push_back(buffer); }
        catch (...) { state.note(cudaFree(buffer), "release unrecorded resident allocation"); throw; }
        return static_cast<T*>(buffer);
    }

    void release(ErrorState& state) noexcept {
        if (allocations.empty() && pinned.empty()) return;
        int previous = -1;
        state.note(cudaGetDevice(&previous), "cudaGetDevice during cleanup");
        const bool change = previous != device;
        if (change) state.note(cudaSetDevice(device), "select resident device during cleanup");
        for (void* buffer : allocations) state.note(cudaFree(buffer), "cudaFree resident corpus");
        allocations.clear();
        pinned.release(state);
        if (change && previous >= 0) state.note(cudaSetDevice(previous), "restore device after cleanup");
    }
};

struct ResidentDeleter {
    ErrorState* state;
    void operator()(Resident* resident) const noexcept {
        if (resident != nullptr) { resident->release(*state); delete resident; }
    }
};

using ResidentOwner = std::unique_ptr<Resident, ResidentDeleter>;

Resident& residentOf(void* handle, ResidentKind kind) {
    require(handle != nullptr, "no resident corpus was given");
    Resident& resident = *static_cast<Resident*>(handle);
    require(resident.kind == kind, "the resident corpus has another screening type");
    return resident;
}

int blocksFor(int32_t targets) {
    return static_cast<int>(std::min<int64_t>(kMaxBlocks, (static_cast<int64_t>(targets) + kThreads - 1) / kThreads));
}

/** Consecutive threads read consecutive targets from each fingerprint word plane. */
__global__ __launch_bounds__(kThreads) void fingerprintKernel(
    int32_t targets, int32_t words, const uint64_t* __restrict__ fingerprints,
    const int32_t* __restrict__ targetPopcounts, int32_t queries,
    const uint64_t* __restrict__ queryWords, const int32_t* __restrict__ cutoffs, int32_t cutoffStride,
    smsd_fingerprint_hit* __restrict__ candidates, unsigned char* __restrict__ selected) {
    const int64_t stride = static_cast<int64_t>(gridDim.x) * blockDim.x;
    for (int64_t queryIndex = blockIdx.y; queryIndex < queries; queryIndex += gridDim.y) {
        for (int64_t targetIndex = static_cast<int64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
             targetIndex < targets; targetIndex += stride) {
            int32_t common = 0;
            for (int32_t wordIndex = 0; wordIndex < words; ++wordIndex)
                common += __popcll(queryWords[queryIndex * words + wordIndex]
                                   & fingerprints[static_cast<int64_t>(wordIndex) * targets + targetIndex]);
            const int64_t pair = queryIndex * targets + targetIndex;
            candidates[pair] = {static_cast<int32_t>(targetIndex), common};
            selected[pair] = common >= cutoffs[queryIndex * cutoffStride + targetPopcounts[targetIndex]];
        }
    }
}

/** Integer label overlaps, with one shared query table per block. */
__global__ __launch_bounds__(kThreads) void labelOverlapKernel(
    const int32_t* __restrict__ queryCounts, int32_t tableSize, int32_t queries,
    int32_t targets, const uint64_t* __restrict__ offsets,
    const smsd_label_count* __restrict__ entries, int32_t* __restrict__ overlaps) {
    extern __shared__ int32_t table[];
    const int64_t stride = static_cast<int64_t>(gridDim.x) * blockDim.x;
    for (int64_t queryIndex = blockIdx.y; queryIndex < queries; queryIndex += gridDim.y) {
        for (int32_t label = threadIdx.x; label < tableSize; label += blockDim.x)
            table[label] = queryCounts[queryIndex * tableSize + label];
        __syncthreads();
        for (int64_t targetIndex = static_cast<int64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
             targetIndex < targets; targetIndex += stride) {
            int32_t overlap = 0;
            const uint64_t end = offsets[targetIndex + 1];
            for (uint64_t entryIndex = offsets[targetIndex]; entryIndex < end; ++entryIndex) {
                const smsd_label_count entry = entries[entryIndex];
                if (entry.label < tableSize) overlap += min(static_cast<int32_t>(entry.count), table[entry.label]);
            }
            overlaps[queryIndex * targets + targetIndex] = overlap;
        }
        __syncthreads();
    }
}

int deviceInfo(char* text, size_t textSize, char* error, size_t errorSize) {
    return guarded(error, errorSize, [&](ErrorState& state) {
        int device = 0;
        state.check(cudaGetDevice(&device), "cudaGetDevice");
        cudaDeviceProp properties;
        state.check(cudaGetDeviceProperties(&properties, device), "cudaGetDeviceProperties");
        if (text != nullptr && textSize != 0)
            std::snprintf(text, textSize, "%s (compute %d.%d, %zu MB)", properties.name, properties.major,
                          properties.minor, static_cast<size_t>(properties.totalGlobalMem / (1024 * 1024)));
    });
}

int minimumTargets(uint32_t operation, uint64_t* minimum, char* error, size_t errorSize) {
    return guarded(error, errorSize, [&](ErrorState& state) {
        require(minimum != nullptr && (operation == 1 || operation == 2 || operation == 3), "unknown screening operation");
        *minimum = std::numeric_limits<uint64_t>::max();
        if (operation != 1) return;
        int device = 0;
        state.check(cudaGetDevice(&device), "cudaGetDevice");
        cudaDeviceProp properties;
        state.check(cudaGetDeviceProperties(&properties, device), "cudaGetDeviceProperties");
        if (properties.major == 7 && properties.minor == 5 && std::strcmp(properties.name, "Tesla T4") == 0)
            *minimum = 100000;
    });
}

int createFingerprints(int64_t numTargets, int32_t numWords, const uint64_t* wordMajor,
                       const int32_t* targetPopcounts, void** handle, char* error, size_t errorSize) {
    if (handle != nullptr) *handle = nullptr;
    return guarded(error, errorSize, [&](ErrorState& state) {
        require(handle != nullptr && numTargets >= 0 && numTargets <= std::numeric_limits<int32_t>::max()
                    && numWords > 0 && numWords <= std::numeric_limits<int32_t>::max() / 64,
                "invalid resident fingerprint dimensions");
        require(numTargets == 0 || (wordMajor != nullptr && targetPopcounts != nullptr),
                "resident fingerprints need words and popcounts");
        beforeWork(state);
        ResidentOwner resident(new Resident(ResidentKind::Fingerprints), ResidentDeleter{&state});
        state.check(cudaGetDevice(&resident->device), "cudaGetDevice");
        resident->targets = static_cast<int32_t>(numTargets);
        resident->words = numWords;
        for (int64_t targetIndex = 0; targetIndex < numTargets; ++targetIndex) {
            require(targetPopcounts[targetIndex] >= 0 && targetPopcounts[targetIndex] <= numWords * 64,
                    "target fingerprint popcount exceeds its capacity");
            resident->maximumPopcount = std::max(resident->maximumPopcount, targetPopcounts[targetIndex]);
        }
        const uint64_t words = static_cast<uint64_t>(numTargets) * numWords;
        resident->fingerprints = resident->allocate<uint64_t>(state, words);
        resident->popcounts = resident->allocate<int32_t>(state, numTargets);
        {
            Call call(state);
            call.upload(resident->fingerprints, wordMajor, words);
            call.upload(resident->popcounts, targetPopcounts, numTargets);
            call.finish();
        }
        if (state.code != cudaSuccess) throw ErrorState::Failure{};
        *handle = resident.release();
    });
}

int fingerprintScreen(void* handle, int32_t numQueries, const uint64_t* queryWords,
                      const int32_t* minimumCommon, int32_t cutoffStride, smsd_fingerprint_hit* hits,
                      uint64_t hitCapacity, uint64_t* hitOffsets, char* error, size_t errorSize) {
    return guarded(error, errorSize, [&](ErrorState& state) {
        Resident& resident = residentOf(handle, ResidentKind::Fingerprints);
        require(numQueries >= 0 && hitOffsets != nullptr, "invalid fingerprint query dimensions");
        std::fill(hitOffsets, hitOffsets + static_cast<size_t>(numQueries) + 1, uint64_t(0));
        if (numQueries == 0 || resident.targets == 0) return;
        require(numQueries == 0 || (queryWords != nullptr && minimumCommon != nullptr
                    && cutoffStride > resident.maximumPopcount), "fingerprint cutoff table is too short");
        beforeWork(state);
        DeviceScope device(state, resident.device);
        Call call(state, &resident.pinned);
        const uint64_t pairs = static_cast<uint64_t>(resident.targets) * numQueries;
        const uint64_t queryWordCount = static_cast<uint64_t>(resident.words) * numQueries;
        const uint64_t cutoffCount = static_cast<uint64_t>(cutoffStride) * numQueries;
        auto* deviceQueries = call.device<uint64_t>(queryWordCount);
        auto* deviceCutoffs = call.device<int32_t>(cutoffCount);
        auto* candidates = call.device<smsd_fingerprint_hit>(pairs);
        auto* flags = call.device<unsigned char>(pairs);
        auto* selected = call.device<smsd_fingerprint_hit>(pairs);
        auto* selectedCounts = call.device<int32_t>(numQueries);
        auto* hostCounts = call.pinned<int32_t>(numQueries);
        size_t temporaryBytes = 0;
        state.check(cub::DeviceSelect::Flagged(nullptr, temporaryBytes, candidates, flags, selected,
                    selectedCounts, resident.targets, call.stream()), "CUB selection storage");
        auto* temporary = call.device<unsigned char>(temporaryBytes);
        call.upload(deviceQueries, queryWords, queryWordCount);
        call.upload(deviceCutoffs, minimumCommon, cutoffCount);
        fingerprintKernel<<<dim3(blocksFor(resident.targets), std::min(numQueries, kMaxBlocks)),
                            kThreads, 0, call.stream()>>>(resident.targets, resident.words,
            resident.fingerprints, resident.popcounts, numQueries, deviceQueries, deviceCutoffs,
            cutoffStride, candidates, flags);
        call.launched();
        for (int32_t queryIndex = 0; queryIndex < numQueries; ++queryIndex) {
            const size_t start = static_cast<size_t>(queryIndex) * resident.targets;
            call.pending();
            state.check(cub::DeviceSelect::Flagged(temporary, temporaryBytes, candidates + start, flags + start,
                        selected + start, selectedCounts + queryIndex, resident.targets, call.stream()),
                        "CUB ordered fingerprint selection");
        }
        call.copy(hostCounts, selectedCounts, bytesFor(numQueries, sizeof(int32_t)), cudaMemcpyDeviceToHost);
        call.finish();
        for (int32_t queryIndex = 0; queryIndex < numQueries; ++queryIndex) {
            require(hostCounts[queryIndex] >= 0 && hostCounts[queryIndex] <= resident.targets,
                    "invalid compacted fingerprint hit count");
            hitOffsets[queryIndex + 1] = hitOffsets[queryIndex] + static_cast<uint64_t>(hostCounts[queryIndex]);
        }
        const uint64_t total = hitOffsets[numQueries];
        require(total <= hitCapacity && (total == 0 || hits != nullptr), "fingerprint output buffer is too short");
        if (total == 0) return;
        auto* hostHits = call.pinned<smsd_fingerprint_hit>(total);
        for (int32_t queryIndex = 0; queryIndex < numQueries; ++queryIndex)
            if (hostCounts[queryIndex] != 0)
                call.copy(hostHits + hitOffsets[queryIndex], selected + static_cast<size_t>(queryIndex) * resident.targets,
                          bytesFor(hostCounts[queryIndex], sizeof(smsd_fingerprint_hit)), cudaMemcpyDeviceToHost);
        call.finish();
        std::memcpy(hits, hostHits, bytesFor(total, sizeof(smsd_fingerprint_hit)));
    });
}

int createLabels(int64_t numTargets, const uint64_t* offsets, const smsd_label_count* entries,
                 void** handle, char* error, size_t errorSize) {
    if (handle != nullptr) *handle = nullptr;
    return guarded(error, errorSize, [&](ErrorState& state) {
        require(handle != nullptr && numTargets >= 0 && numTargets <= std::numeric_limits<int32_t>::max()
                    && offsets != nullptr, "invalid resident label dimensions");
        require(offsets[0] == 0, "resident label offsets must start at zero");
        for (int64_t targetIndex = 0; targetIndex < numTargets; ++targetIndex)
            require(offsets[targetIndex] <= offsets[targetIndex + 1], "resident label offsets are not ordered");
        const uint64_t entryCount = offsets[numTargets];
        require(entryCount == 0 || entries != nullptr, "resident labels need histogram entries");
        beforeWork(state);
        ResidentOwner resident(new Resident(ResidentKind::Labels), ResidentDeleter{&state});
        state.check(cudaGetDevice(&resident->device), "cudaGetDevice");
        resident->targets = static_cast<int32_t>(numTargets);
        resident->offsets = resident->allocate<uint64_t>(state, static_cast<uint64_t>(numTargets) + 1);
        resident->entries = resident->allocate<smsd_label_count>(state, entryCount);
        {
            Call call(state);
            call.upload(resident->offsets, offsets, static_cast<uint64_t>(numTargets) + 1);
            call.upload(resident->entries, entries, entryCount);
            call.finish();
        }
        if (state.code != cudaSuccess) throw ErrorState::Failure{};
        *handle = resident.release();
    });
}

int residentLabelOverlaps(void* handle, int32_t numQueries, const int32_t* queryCounts,
                          int32_t tableSize, int32_t* overlaps, char* error, size_t errorSize) {
    return guarded(error, errorSize, [&](ErrorState& state) {
        Resident& resident = residentOf(handle, ResidentKind::Labels);
        require(numQueries >= 0 && tableSize > 0 && tableSize <= kMaxTable, "invalid label query dimensions");
        if (numQueries == 0 || resident.targets == 0) return;
        require(queryCounts != nullptr && overlaps != nullptr, "label screening needs query and output buffers");
        beforeWork(state);
        DeviceScope device(state, resident.device);
        Call call(state, &resident.pinned);
        const uint64_t querySize = static_cast<uint64_t>(numQueries) * tableSize;
        const uint64_t outputSize = static_cast<uint64_t>(numQueries) * resident.targets;
        auto* deviceQueries = call.device<int32_t>(querySize);
        auto* deviceOverlaps = call.device<int32_t>(outputSize);
        auto* hostOverlaps = call.pinned<int32_t>(outputSize);
        call.upload(deviceQueries, queryCounts, querySize);
        labelOverlapKernel<<<dim3(blocksFor(resident.targets), std::min(numQueries, kMaxBlocks)), kThreads,
                             bytesFor(tableSize, sizeof(int32_t)), call.stream()>>>(deviceQueries,
            tableSize, numQueries, resident.targets, resident.offsets, resident.entries, deviceOverlaps);
        call.launched();
        call.copy(hostOverlaps, deviceOverlaps, bytesFor(outputSize, sizeof(int32_t)), cudaMemcpyDeviceToHost);
        call.finish();
        std::memcpy(overlaps, hostOverlaps, bytesFor(outputSize, sizeof(int32_t)));
    });
}

int destroyResident(void* handle, char* error, size_t errorSize) {
    return guarded(error, errorSize, [&](ErrorState& state) {
        ResidentOwner resident(static_cast<Resident*>(handle), ResidentDeleter{&state});
    });
}

/** The original callback has one call-scoped stream and retains no allocations. */
int labelOverlaps(const int32_t* queryCounts, int32_t tableSize, int64_t numTargets,
                  const uint64_t* offsets, const smsd_label_count* entries, int32_t* overlaps,
                  char* error, size_t errorSize) {
    if (numTargets == 0) return 0;
    return guarded(error, errorSize, [&](ErrorState& state) {
        require(numTargets > 0 && numTargets <= std::numeric_limits<int32_t>::max()
                    && tableSize > 0 && tableSize <= kMaxTable && queryCounts != nullptr
                    && offsets != nullptr && overlaps != nullptr, "invalid label screening dimensions");
        require(offsets[0] == 0, "label offsets must start at zero");
        for (int64_t targetIndex = 0; targetIndex < numTargets; ++targetIndex)
            require(offsets[targetIndex] <= offsets[targetIndex + 1], "label offsets are not ordered");
        const uint64_t entryCount = offsets[numTargets];
        require(entryCount == 0 || entries != nullptr, "label screening needs histogram entries");
        beforeWork(state);
        Call call(state);
        auto* deviceOffsets = call.device<uint64_t>(static_cast<uint64_t>(numTargets) + 1);
        auto* deviceEntries = call.device<smsd_label_count>(entryCount);
        auto* deviceQueries = call.device<int32_t>(tableSize);
        auto* deviceOverlaps = call.device<int32_t>(numTargets);
        auto* hostOverlaps = call.pinned<int32_t>(numTargets);
        call.upload(deviceOffsets, offsets, static_cast<uint64_t>(numTargets) + 1);
        call.upload(deviceEntries, entries, entryCount);
        call.upload(deviceQueries, queryCounts, tableSize);
        labelOverlapKernel<<<dim3(blocksFor(static_cast<int32_t>(numTargets)), 1), kThreads,
                             bytesFor(tableSize, sizeof(int32_t)), call.stream()>>>(deviceQueries,
            tableSize, 1, static_cast<int32_t>(numTargets), deviceOffsets, deviceEntries, deviceOverlaps);
        call.launched();
        call.copy(hostOverlaps, deviceOverlaps, bytesFor(numTargets, sizeof(int32_t)), cudaMemcpyDeviceToHost);
        call.finish();
        std::memcpy(overlaps, hostOverlaps, bytesFor(numTargets, sizeof(int32_t)));
    });
}

const smsd_screen_accelerator_v1 kTable = {
    SMSD_SCREEN_ACCELERATOR_ABI,
    static_cast<uint32_t>(sizeof(smsd_screen_accelerator_v1)),
    SMSD_CUDA_SMSD_VERSION,
    "CUDA",
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
SMSD_CUDA_EXPORT const smsd_screen_accelerator_v1* smsd_cuda_screen_accelerator(void) { return &kTable; }

/** The device count, or -1 with the CUDA failure written to error. */
SMSD_CUDA_EXPORT int smsd_cuda_device_count(char* error, size_t errorSize) {
    int count = 0;
    const cudaError_t status = cudaGetDeviceCount(&count);
    if (status != cudaSuccess) { describe(error, errorSize, "cudaGetDeviceCount", status); return -1; }
    return count;
}

/** Writes a description of the current device; 0 on success. */
SMSD_CUDA_EXPORT int smsd_cuda_device_info(char* text, size_t textSize, char* error, size_t errorSize) {
    return deviceInfo(text, textSize, error, errorSize);
}

#ifdef SMSD_CUDA_TESTING
/** Test builds fail the next allocation (1), Nth allocation (-N), or next device operation (2). */
SMSD_CUDA_EXPORT void smsd_cuda_test_fail_next(int failure) { injectedFailure = failure; }

/** The active blocks per multiprocessor for each test kernel on the current device. */
SMSD_CUDA_EXPORT int smsd_cuda_test_occupancy(int32_t tableSize, int* fingerprintBlocks, int* labelBlocks,
                                            char* error, size_t errorSize) {
    return guarded(error, errorSize, [&](ErrorState& state) {
        require(tableSize > 0 && tableSize <= kMaxTable && fingerprintBlocks != nullptr && labelBlocks != nullptr,
                "invalid occupancy output buffers");
        state.check(cudaOccupancyMaxActiveBlocksPerMultiprocessor(fingerprintBlocks, fingerprintKernel, kThreads, 0),
                    "fingerprint kernel occupancy");
        state.check(cudaOccupancyMaxActiveBlocksPerMultiprocessor(labelBlocks, labelOverlapKernel, kThreads,
                    bytesFor(tableSize, sizeof(int32_t))), "label kernel occupancy");
    });
}
#endif
