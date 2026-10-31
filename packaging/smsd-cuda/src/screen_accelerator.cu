/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Integer screening on immutable device-resident fingerprints and label tables.
 * Each corpus owns its allocations. Calls lease bounded workspace or use their
 * own stream and scratch. Ordered screening leaves scores on the host; optional top-k ranks
 * binary64 scores on the device and returns the original integer intersections.
 */
#include "smsd/screen_accelerator_abi.h"
#include "smsd/cuda_device_policy.hpp"
#include <cuda_runtime.h>
#include "binary_mma_layout.hpp"
#include <cub/device/device_select.cuh>
#include <cub/device/device_radix_sort.cuh>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
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
constexpr int kFingerprintQueryTile = 4;
constexpr int kFingerprintSharedBytes = 32 * 1024;

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

/** Whether the device can run the screens; when it cannot, the reason is written to reason. */
bool deviceSupported(ErrorState& state, int device, char* reason, size_t reasonSize) {
    int major = 0;
    int minor = 0;
    state.check(cudaDeviceGetAttribute(&major, cudaDevAttrComputeCapabilityMajor, device), "compute capability");
    state.check(cudaDeviceGetAttribute(&minor, cudaDevAttrComputeCapabilityMinor, device), "compute capability");
    if (smsd::cuda_policy::supportedCapability(major, minor)) return true;
    cudaDeviceProp properties;
    state.check(cudaGetDeviceProperties(&properties, device), "cudaGetDeviceProperties");
    return smsd::cuda_policy::explainCapability(properties.name, major, minor, reason, reasonSize);
}

/** Refuses work on a current device below the supported compute capability, with the reason. */
void requireSupportedDevice(ErrorState& state) {
    int device = 0;
    state.check(cudaGetDevice(&device), "cudaGetDevice");
    char reason[256] = {};
    if (!deviceSupported(state, device, reason, sizeof reason)) throw std::runtime_error(reason);
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

/** One non-blocking lease retains at most 256 MiB of device scratch and one graph. */
class DeviceWorkspace {
public:
    std::unique_lock<std::mutex> tryAcquire() { return std::unique_lock<std::mutex>(mutex_, std::try_to_lock); }
    bool usable() const noexcept { return !disabled_; }
    bool empty() const noexcept { return retainedBytes_ == 0 && stream_ == nullptr && graph_ == nullptr; }

    cudaStream_t stream(ErrorState& state) {
        if (stream_ == nullptr) {
            cudaStream_t created = nullptr;
            const cudaError_t status = cudaStreamCreateWithFlags(&created, cudaStreamNonBlocking);
            state.note(status, "create retained CUDA stream");
            if (status != cudaSuccess) {
                disabled_ = true;
                if (created != nullptr) state.note(cudaStreamDestroy(created), "discard uncertain retained CUDA stream");
                throw ErrorState::Failure{};
            }
            stream_ = created;
        }
        return stream_;
    }

    void* allocate(ErrorState& state, size_t index, size_t bytes) {
        if (disabled_ || index >= slots_.size()) return nullptr;
        Slot& slot = slots_[index];
        if (slot.capacity >= bytes) return slot.buffer;
        // Include both old and replacement buffers while growing a slot.
        if (bytes > kMaximumBytes - retainedBytes_) return nullptr;
        discardGraph(state);
        void* replacement = nullptr;
        state.check(allocateDevice(&replacement, bytes), "cudaMalloc retained scratch");
        spare_ = {replacement, bytes, false};
        retainedBytes_ += bytes;
        const size_t oldCapacity = slot.capacity;
        if (releaseSlot(state, slot) != cudaSuccess) throw ErrorState::Failure{};
        retainedBytes_ -= oldCapacity;
        slot = spare_;
        spare_ = {};
        return slot.buffer;
    }

    bool graphMatches(int32_t queries, int32_t cutoffStride) const noexcept {
        return graph_ != nullptr && graphQueries_ == queries && graphCutoffStride_ == cutoffStride;
    }

    void beginGraph(ErrorState& state) {
        discardGraph(state);
        state.check(cudaStreamBeginCapture(stream_, cudaStreamCaptureModeThreadLocal), "begin fingerprint graph capture");
        capturing_ = true;
    }

    void finishGraph(ErrorState& state, int32_t queries, int32_t cutoffStride) {
        cudaGraph_t captured = nullptr;
        const cudaError_t ended = cudaStreamEndCapture(stream_, &captured);
        capturing_ = false;
        state.note(ended, "end fingerprint graph capture");
        if (ended != cudaSuccess) {
            if (captured != nullptr) state.note(cudaGraphDestroy(captured), "discard invalid fingerprint graph");
            throw ErrorState::Failure{};
        }
        cudaGraphExec_t executable = nullptr;
        const cudaError_t instantiated = cudaGraphInstantiateWithFlags(&executable, captured, 0);
        state.note(instantiated, "instantiate fingerprint graph");
        state.note(cudaGraphDestroy(captured), "release captured fingerprint graph");
        if (state.code != cudaSuccess) {
            if (executable != nullptr) state.note(cudaGraphExecDestroy(executable), "discard fingerprint graph executable");
            throw ErrorState::Failure{};
        }
        graph_ = executable;
        graphQueries_ = queries;
        graphCutoffStride_ = cutoffStride;
    }

    void launchGraph(ErrorState& state) {
        state.check(cudaGraphLaunch(graph_, stream_), "launch fingerprint graph");
    }

    /** End an abandoned capture before synchronising or releasing its stream. */
    void cancelCapture(ErrorState& state) noexcept {
        if (!capturing_) return;
        cudaGraph_t captured = nullptr;
        state.note(cudaStreamEndCapture(stream_, &captured), "cancel fingerprint graph capture");
        capturing_ = false;
        if (captured != nullptr) state.note(cudaGraphDestroy(captured), "discard cancelled fingerprint graph");
    }

    void release(ErrorState& state) noexcept {
        cancelCapture(state);
        if (graph_ != nullptr) {
            const cudaGraphExec_t graph = graph_;
            graph_ = nullptr;
            state.note(cudaGraphExecDestroy(graph), "release fingerprint graph executable");
        }
        for (Slot& slot : slots_) releaseSlot(state, slot);
        releaseSlot(state, spare_);
        retainedBytes_ = 0;
        if (stream_ != nullptr) {
            const cudaStream_t stream = stream_;
            stream_ = nullptr;
            state.note(cudaStreamDestroy(stream), "release retained CUDA stream");
        }
    }

    /** Uncertain resources are released once and never leased again. */
    void disable(ErrorState& state) noexcept { disabled_ = true; release(state); }

private:
    struct Slot { void* buffer = nullptr; size_t capacity = 0; bool releaseAttempted = false; };
    static constexpr size_t kMaximumBytes = size_t(256) * 1024 * 1024;

    static cudaError_t releaseSlot(ErrorState& state, Slot& slot) noexcept {
        if (slot.buffer == nullptr || slot.releaseAttempted) return cudaSuccess;
        slot.releaseAttempted = true;
        const cudaError_t status = cudaFree(slot.buffer);
        state.note(status, "cudaFree retained scratch");
        if (status == cudaSuccess) slot = {};
        return status;
    }

    void discardGraph(ErrorState& state) {
        if (graph_ == nullptr) return;
        const cudaGraphExec_t graph = graph_;
        graph_ = nullptr;
        state.check(cudaGraphExecDestroy(graph), "replace fingerprint graph executable");
    }

    std::mutex mutex_;
    std::array<Slot, 12> slots_{};
    Slot spare_;
    size_t retainedBytes_ = 0;
    cudaStream_t stream_ = nullptr;
    cudaGraphExec_t graph_ = nullptr;
    int32_t graphQueries_ = 0, graphCutoffStride_ = 0;
    bool capturing_ = false;
    bool disabled_ = false;
};

/** A call owns resources or leases one resident workspace; concurrent calls do not wait. */
class Call {
public:
    explicit Call(ErrorState& state, PinnedStaging* retained = nullptr, DeviceWorkspace* workspace = nullptr)
        : state_(state), retained_(retained), lease_(retained != nullptr ? retained->tryAcquire() : std::unique_lock<std::mutex>()),
          workspace_(workspace), workspaceLease_(workspace != nullptr ? workspace->tryAcquire() : std::unique_lock<std::mutex>()) {
        device_.reserve(12);
        pinned_.reserve(6);
        if (workspaceLease_.owns_lock() && !workspace_->usable()) workspaceLease_.unlock();
        if (workspaceLease_.owns_lock()) stream_ = workspace_->stream(state_);
        else state_.check(cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking), "cudaStreamCreateWithFlags");
    }

    ~Call() {
        if (workspaceLease_.owns_lock()) workspace_->cancelCapture(state_);
        if (pending_) state_.note(cudaStreamSynchronize(stream_), "finish CUDA call");
        for (void* buffer : device_) state_.note(cudaFree(buffer), "cudaFree");
        for (void* buffer : pinned_) state_.note(cudaFreeHost(buffer), "cudaFreeHost");
        if (workspaceLease_.owns_lock()) {
            if (state_.code != cudaSuccess) workspace_->disable(state_);
        } else state_.note(cudaStreamDestroy(stream_), "cudaStreamDestroy");
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
    DeviceWorkspace* graphWorkspace() const noexcept {
        return workspaceLease_.owns_lock() && allDevicesPooled_ ? workspace_ : nullptr;
    }

    void finish() {
        state_.check(cudaStreamSynchronize(stream_), "cudaStreamSynchronize");
        pending_ = false;
    }

private:
    void* allocate(uint64_t count, size_t itemSize, bool host) {
        if (count == 0) return nullptr;
        void* buffer = nullptr;
        const size_t bytes = bytesFor(count, itemSize);
        if (!host && workspaceLease_.owns_lock()) {
            buffer = workspace_->allocate(state_, nextDeviceSlot_++, bytes);
            if (buffer != nullptr) return buffer;
            allDevicesPooled_ = false;
        }
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
    DeviceWorkspace* workspace_;
    std::unique_lock<std::mutex> workspaceLease_;
    size_t nextDeviceSlot_ = 0;
    bool allDevicesPooled_ = true;
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
    DeviceWorkspace workspace;

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
        if (allocations.empty() && pinned.empty() && workspace.empty()) return;
        int previous = -1;
        state.note(cudaGetDevice(&previous), "cudaGetDevice during cleanup");
        const bool change = previous != device;
        if (change) state.note(cudaSetDevice(device), "select resident device during cleanup");
        workspace.release(state);
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

/** Reuse each target word across a tile of queries; output remains query-major. */
template<int QueryTile>
__global__ __launch_bounds__(kThreads) void fingerprintTileKernel(
    int32_t targets, int32_t words, const uint64_t* __restrict__ fingerprints,
    const int32_t* __restrict__ targetPopcounts, int32_t queries,
    const uint64_t* __restrict__ queryWords, const int32_t* __restrict__ cutoffs, int32_t cutoffStride,
    smsd_fingerprint_hit* __restrict__ candidates, unsigned char* __restrict__ selected) {
    extern __shared__ uint64_t tileWords[];
    const int64_t targetStride = static_cast<int64_t>(gridDim.x) * blockDim.x;
    const int64_t queryStride = static_cast<int64_t>(gridDim.y) * QueryTile;
    for (int64_t queryBegin = static_cast<int64_t>(blockIdx.y) * QueryTile;
         queryBegin < queries; queryBegin += queryStride) {
        const int64_t remainingQueries = static_cast<int64_t>(queries) - queryBegin;
        const int queryCount = remainingQueries < QueryTile
            ? static_cast<int>(remainingQueries) : QueryTile;
        for (int64_t word = threadIdx.x; word < static_cast<int64_t>(queryCount) * words; word += blockDim.x)
            tileWords[word] = queryWords[queryBegin * words + word];
        __syncthreads();
        for (int64_t targetIndex = static_cast<int64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
             targetIndex < targets; targetIndex += targetStride) {
            int32_t common[QueryTile] = {};
            for (int32_t wordIndex = 0; wordIndex < words; ++wordIndex) {
                const uint64_t targetWord = fingerprints[static_cast<int64_t>(wordIndex) * targets + targetIndex];
#pragma unroll
                for (int query = 0; query < QueryTile; ++query)
                    if (query < queryCount)
                        common[query] += __popcll(tileWords[static_cast<int64_t>(query) * words + wordIndex] & targetWord);
            }
            const int32_t targetBits = targetPopcounts[targetIndex];
#pragma unroll
            for (int query = 0; query < QueryTile; ++query) {
                if (query < queryCount) {
                    const int64_t queryIndex = queryBegin + query;
                    const int64_t pair = queryIndex * targets + targetIndex;
                    candidates[pair] = {static_cast<int32_t>(targetIndex), common[query]};
                    selected[pair] = common[query] >= cutoffs[queryIndex * cutoffStride + targetBits];
                }
            }
        }
        __syncthreads();
    }
}

#ifdef SMSD_CUDA_BINARY_MMA_PROTOTYPE
/** One warp computes 16 query by 8 target intersections without floating-point arithmetic. */
__global__ __launch_bounds__(kThreads) void fingerprintBinaryMmaKernel(
    int32_t targets, int32_t words, const uint64_t* __restrict__ fingerprints,
    const int32_t* __restrict__ targetPopcounts, int32_t queries,
    const uint64_t* __restrict__ queryWords, const int32_t* __restrict__ cutoffs, int32_t cutoffStride,
    smsd_fingerprint_hit* __restrict__ candidates, unsigned char* __restrict__ selected) {
#if __CUDA_ARCH__ >= 800 && __CUDA_ARCH__ < 1000
    using namespace smsd_cuda_binary_mma;
    const int lane = threadIdx.x & 31;
    const int warp = threadIdx.x >> 5;
    constexpr int warps = kThreads / 32;
    const int64_t targetStride = static_cast<int64_t>(gridDim.x) * warps * 8;
    const int64_t queryStride = static_cast<int64_t>(gridDim.y) * 16;
    for (int64_t queryBegin = static_cast<int64_t>(blockIdx.y) * 16;
         queryBegin < queries; queryBegin += queryStride) {
        for (int64_t targetBegin = (static_cast<int64_t>(blockIdx.x) * warps + warp) * 8;
             targetBegin < targets; targetBegin += targetStride) {
            int32_t common0 = 0, common1 = 0, common2 = 0, common3 = 0;
            for (int wordBegin = 0; wordBegin < words; wordBegin += 4) {
                const uint32_t query0 = queryFragment(queryWords, queries, words, queryBegin, wordBegin, lane, 0);
                const uint32_t query1 = queryFragment(queryWords, queries, words, queryBegin, wordBegin, lane, 1);
                const uint32_t query2 = queryFragment(queryWords, queries, words, queryBegin, wordBegin, lane, 2);
                const uint32_t query3 = queryFragment(queryWords, queries, words, queryBegin, wordBegin, lane, 3);
                const uint32_t target0 = targetFragment(fingerprints, targets, words, targetBegin, wordBegin, lane, 0);
                const uint32_t target1 = targetFragment(fingerprints, targets, words, targetBegin, wordBegin, lane, 1);
                asm volatile(
                    "mma.sync.aligned.m16n8k256.row.col.s32.b1.b1.s32.and.popc "
                    "{%0, %1, %2, %3}, {%4, %5, %6, %7}, {%8, %9}, {%0, %1, %2, %3};"
                    : "+r"(common0), "+r"(common1), "+r"(common2), "+r"(common3)
                    : "r"(query0), "r"(query1), "r"(query2), "r"(query3), "r"(target0), "r"(target1));
            }
            const int32_t common[4] = {common0, common1, common2, common3};
#pragma unroll
            for (int fragment = 0; fragment < 4; ++fragment) {
                const int64_t query = queryBegin + resultRow(lane, fragment);
                const int64_t target = targetBegin + resultColumn(lane, fragment);
                if (query < queries && target < targets) {
                    const int64_t pair = query * targets + target;
                    candidates[pair] = {static_cast<int32_t>(target), common[fragment]};
                    selected[pair] = common[fragment] >= cutoffs[query * cutoffStride + targetPopcounts[target]];
                }
            }
        }
    }
#else
    // Host dispatch excludes architectures whose binary MMA path is not qualified.
    asm volatile("trap;");
#endif
}

bool binaryMmaRequested(ErrorState& state, int device) {
    const char* requested = std::getenv("SMSD_CUDA_BINARY_MMA");
    if (requested == nullptr || std::strcmp(requested, "1") != 0) return false;
    cudaDeviceProp properties;
    state.check(cudaGetDeviceProperties(&properties, device), "binary MMA device properties");
    return properties.major >= 8 && properties.major < 10;
}
#endif

/** Query popcounts are exact integers, independent of host scoring or packing conventions. */
__global__ __launch_bounds__(kThreads) void fingerprintQueryPopcounts(
    int32_t words, int32_t queries, const uint64_t* queryWords, int32_t* popcounts) {
    __shared__ int32_t partial[kThreads];
    for (int64_t query = blockIdx.x; query < queries; query += gridDim.x) {
        int32_t count = 0;
        for (int64_t word = threadIdx.x; word < words; word += blockDim.x)
            count += __popcll(queryWords[query * words + word]);
        partial[threadIdx.x] = count;
        __syncthreads();
        for (int width = kThreads / 2; width != 0; width /= 2) {
            if (threadIdx.x < width) partial[threadIdx.x] += partial[threadIdx.x + width];
            __syncthreads();
        }
        if (threadIdx.x == 0) popcounts[query] = partial[0];
        __syncthreads();
    }
}

/** Stable radix sorting retains ascending target index when binary64 scores tie. */
__global__ __launch_bounds__(kThreads) void fingerprintScoreKeys(
    int32_t count, const smsd_fingerprint_hit* hits, const int32_t* queryPopcount,
    const int32_t* targetPopcounts, double* scores) {
    const int64_t stride = static_cast<int64_t>(gridDim.x) * blockDim.x;
    for (int64_t index = static_cast<int64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
         index < count; index += stride) {
        const smsd_fingerprint_hit hit = hits[index];
        const int64_t unionBits = static_cast<int64_t>(*queryPopcount) + targetPopcounts[hit.target_index] - hit.common_bits;
        scores[index] = unionBits == 0 ? 1.0 : __ddiv_rn(static_cast<double>(hit.common_bits),
                                                        static_cast<double>(unionBits));
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
        *minimum = smsd::cuda_policy::automaticMinimumTargets(operation, properties.name, properties.major,
            properties.minor, std::getenv("SMSD_DISABLE_GPU"), std::getenv("SMSD_CUDA_MIN_TARGETS"));
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
        requireSupportedDevice(state);
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
        const char* poolSetting = std::getenv("SMSD_CUDA_POOL");
        const bool poolEnabled = poolSetting == nullptr || std::strcmp(poolSetting, "0") != 0;
        Call call(state, &resident.pinned, poolEnabled ? &resident.workspace : nullptr);
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
        const size_t wordBytes = bytesFor(resident.words, sizeof(uint64_t));
        int queryTile = kFingerprintQueryTile;
        if (numQueries >= 16 && wordBytes <= kFingerprintSharedBytes / 16) queryTile = 16;
        else if (numQueries >= 8 && wordBytes <= kFingerprintSharedBytes / 8) queryTile = 8;
        const size_t sharedBytes = bytesFor(static_cast<uint64_t>(resident.words) * queryTile,
                                            sizeof(uint64_t));
        bool useBinaryMma = false;
#ifdef SMSD_CUDA_BINARY_MMA_PROTOTYPE
        useBinaryMma = binaryMmaRequested(state, resident.device);
#endif
        const auto screenAndSelect = [&] {
#ifdef SMSD_CUDA_BINARY_MMA_PROTOTYPE
            if (useBinaryMma) {
                constexpr int targetsPerBlock = (kThreads / 32) * 8;
                const int blocks = static_cast<int>(std::min<int64_t>(kMaxBlocks,
                    (static_cast<int64_t>(resident.targets) + targetsPerBlock - 1) / targetsPerBlock));
                const int queryTiles = (numQueries - 1) / 16 + 1;
                fingerprintBinaryMmaKernel<<<dim3(blocks, std::min(queryTiles, kMaxBlocks)),
                                              kThreads, 0, call.stream()>>>(
                    resident.targets, resident.words, resident.fingerprints, resident.popcounts,
                    numQueries, deviceQueries, deviceCutoffs, cutoffStride, candidates, flags);
            } else
#endif
            if (numQueries > 1 && sharedBytes <= kFingerprintSharedBytes) {
                const int queryTiles = (numQueries - 1) / queryTile + 1;
                const dim3 grid(blocksFor(resident.targets), std::min(queryTiles, kMaxBlocks));
                if (queryTile == 16) {
                    fingerprintTileKernel<16><<<grid, kThreads, sharedBytes, call.stream()>>>(
                        resident.targets, resident.words, resident.fingerprints, resident.popcounts,
                        numQueries, deviceQueries, deviceCutoffs, cutoffStride, candidates, flags);
                } else if (queryTile == 8) {
                    fingerprintTileKernel<8><<<grid, kThreads, sharedBytes, call.stream()>>>(
                        resident.targets, resident.words, resident.fingerprints, resident.popcounts,
                        numQueries, deviceQueries, deviceCutoffs, cutoffStride, candidates, flags);
                } else {
                    fingerprintTileKernel<kFingerprintQueryTile><<<grid, kThreads, sharedBytes, call.stream()>>>(
                        resident.targets, resident.words, resident.fingerprints, resident.popcounts,
                        numQueries, deviceQueries, deviceCutoffs, cutoffStride, candidates, flags);
                }
            } else {
                fingerprintKernel<<<dim3(blocksFor(resident.targets), std::min(numQueries, kMaxBlocks)),
                                    kThreads, 0, call.stream()>>>(resident.targets, resident.words,
                    resident.fingerprints, resident.popcounts, numQueries, deviceQueries, deviceCutoffs,
                    cutoffStride, candidates, flags);
            }
            call.launched();
            for (int32_t queryIndex = 0; queryIndex < numQueries; ++queryIndex) {
                const size_t start = static_cast<size_t>(queryIndex) * resident.targets;
                call.pending();
                state.check(cub::DeviceSelect::Flagged(temporary, temporaryBytes, candidates + start, flags + start,
                            selected + start, selectedCounts + queryIndex, resident.targets, call.stream()),
                            "CUB ordered fingerprint selection");
            }
        };
        const char* graphSetting = std::getenv("SMSD_CUDA_GRAPHS");
        DeviceWorkspace* graphWorkspace = call.graphWorkspace();
        const bool useGraph = graphWorkspace != nullptr && !useBinaryMma && numQueries <= 64
            && graphSetting != nullptr && std::strcmp(graphSetting, "1") == 0;
        if (useGraph) {
            if (!graphWorkspace->graphMatches(numQueries, cutoffStride)) {
                graphWorkspace->beginGraph(state);
                screenAndSelect();
                graphWorkspace->finishGraph(state, numQueries, cutoffStride);
            }
            call.pending();
            graphWorkspace->launchGraph(state);
        } else {
            screenAndSelect();
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

/** Top-k scratch is linear in one target row; only selected prefixes return to the host. */
int fingerprintTopK(void* handle, int32_t numQueries, const uint64_t* queryWords,
                    const int32_t* minimumCommon, int32_t cutoffStride, int32_t k,
                    smsd_fingerprint_hit* hits, uint64_t hitCapacity, uint64_t* hitOffsets,
                    char* error, size_t errorSize) {
    return guarded(error, errorSize, [&](ErrorState& state) {
        Resident& resident = residentOf(handle, ResidentKind::Fingerprints);
        require(numQueries >= 0 && k >= 0 && hitOffsets != nullptr, "invalid fingerprint top-k dimensions");
        std::fill(hitOffsets, hitOffsets + static_cast<size_t>(numQueries) + 1, uint64_t(0));
        if (numQueries == 0 || resident.targets == 0 || k == 0) return;
        require(queryWords != nullptr && minimumCommon != nullptr && cutoffStride > resident.maximumPopcount,
                "fingerprint cutoff table is too short");
        beforeWork(state);
        DeviceScope device(state, resident.device);
        Call call(state, &resident.pinned);
        const uint64_t queryWordCount = static_cast<uint64_t>(resident.words) * numQueries;
        const uint64_t cutoffCount = static_cast<uint64_t>(cutoffStride) * numQueries;
        const uint64_t maximumHits = static_cast<uint64_t>(numQueries) * std::min(k, resident.targets);
        auto* deviceQueries = call.device<uint64_t>(queryWordCount);
        auto* deviceCutoffs = call.device<int32_t>(cutoffCount);
        auto* queryPopcounts = call.device<int32_t>(numQueries);
        auto* candidates = call.device<smsd_fingerprint_hit>(resident.targets);
        auto* flags = call.device<unsigned char>(resident.targets);
        auto* selected = call.device<smsd_fingerprint_hit>(resident.targets);
        auto* sorted = call.device<smsd_fingerprint_hit>(resident.targets);
        auto* scores = call.device<double>(resident.targets);
        auto* sortedScores = call.device<double>(resident.targets);
        auto* selectedCount = call.device<int32_t>(1);
        auto* hostCount = call.pinned<int32_t>(1);
        auto* hostHits = call.pinned<smsd_fingerprint_hit>(std::min(maximumHits, hitCapacity));
        size_t selectBytes = 0, sortBytes = 0;
        state.check(cub::DeviceSelect::Flagged(nullptr, selectBytes, candidates, flags, selected,
                    selectedCount, resident.targets, call.stream()), "CUB top-k selection storage");
        state.check(cub::DeviceRadixSort::SortPairsDescending(nullptr, sortBytes, scores, sortedScores,
                    selected, sorted, resident.targets, 0, 64, call.stream()), "CUB top-k sort storage");
        const size_t temporaryCapacity = std::max(selectBytes, sortBytes);
        auto* temporary = call.device<unsigned char>(temporaryCapacity);
        call.upload(deviceQueries, queryWords, queryWordCount);
        call.upload(deviceCutoffs, minimumCommon, cutoffCount);
        fingerprintQueryPopcounts<<<std::min(numQueries, kMaxBlocks), kThreads, 0, call.stream()>>>(
            resident.words, numQueries, deviceQueries, queryPopcounts);
        call.launched();
        for (int32_t query = 0; query < numQueries; ++query) {
            fingerprintKernel<<<blocksFor(resident.targets), kThreads, 0, call.stream()>>>(
                resident.targets, resident.words, resident.fingerprints, resident.popcounts, 1,
                deviceQueries + static_cast<int64_t>(query) * resident.words,
                deviceCutoffs + static_cast<int64_t>(query) * cutoffStride, cutoffStride, candidates, flags);
            call.launched();
            size_t available = temporaryCapacity;
            call.pending();
            state.check(cub::DeviceSelect::Flagged(temporary, available, candidates, flags, selected,
                        selectedCount, resident.targets, call.stream()), "CUB top-k ordered selection");
            call.copy(hostCount, selectedCount, sizeof(int32_t), cudaMemcpyDeviceToHost);
            call.finish();
            require(*hostCount >= 0 && *hostCount <= resident.targets, "invalid fingerprint top-k hit count");
            const int32_t kept = std::min(k, *hostCount);
            hitOffsets[query + 1] = hitOffsets[query] + static_cast<uint64_t>(kept);
            require(hitOffsets[query + 1] <= hitCapacity && (kept == 0 || hits != nullptr),
                    "fingerprint top-k output buffer is too short");
            if (kept == 0) continue;
            // Keep query popcounts on the device, including for long fingerprints.
            fingerprintScoreKeys<<<blocksFor(*hostCount), kThreads, 0, call.stream()>>>(
                *hostCount, selected, queryPopcounts + query, resident.popcounts, scores);
            call.launched();
            available = temporaryCapacity;
            call.pending();
            state.check(cub::DeviceRadixSort::SortPairsDescending(temporary, available, scores, sortedScores,
                        selected, sorted, *hostCount, 0, 64, call.stream()), "CUB stable fingerprint top-k sort");
            call.copy(hostHits + hitOffsets[query], sorted, bytesFor(kept, sizeof(smsd_fingerprint_hit)),
                      cudaMemcpyDeviceToHost);
        }
        call.finish();
        if (hitOffsets[numQueries] != 0)
            std::memcpy(hits, hostHits, bytesFor(hitOffsets[numQueries], sizeof(smsd_fingerprint_hit)));
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
        requireSupportedDevice(state);
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
        requireSupportedDevice(state);
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

/** Optional package entry point; the version-one screening table remains unchanged. */
SMSD_CUDA_EXPORT int smsd_cuda_fingerprint_top_k(
    void* resident, int32_t numQueries, const uint64_t* queryWords, const int32_t* minimumCommon,
    int32_t cutoffStride, int32_t k, smsd_fingerprint_hit* hits, uint64_t hitCapacity,
    uint64_t* hitOffsets, char* error, size_t errorSize) {
    return fingerprintTopK(resident, numQueries, queryWords, minimumCommon, cutoffStride, k,
                           hits, hitCapacity, hitOffsets, error, errorSize);
}

/** Initialise the current device and load screening kernels without running a screen. */
SMSD_CUDA_EXPORT int smsd_cuda_preload(char* error, size_t errorSize) {
    return guarded(error, errorSize, [&](ErrorState& state) {
        beforeWork(state);
        requireSupportedDevice(state);
        state.check(cudaFree(nullptr), "initialise CUDA context");
        cudaFuncAttributes attributes;
        state.check(cudaFuncGetAttributes(&attributes, fingerprintKernel), "load fingerprint kernel");
        state.check(cudaFuncGetAttributes(&attributes, fingerprintTileKernel<4>), "load four-query kernel");
        state.check(cudaFuncGetAttributes(&attributes, fingerprintTileKernel<8>), "load eight-query kernel");
        state.check(cudaFuncGetAttributes(&attributes, fingerprintTileKernel<16>), "load sixteen-query kernel");
        state.check(cudaFuncGetAttributes(&attributes, labelOverlapKernel), "load label kernel");
        state.check(cudaFuncGetAttributes(&attributes, fingerprintQueryPopcounts), "load query popcount kernel");
        state.check(cudaFuncGetAttributes(&attributes, fingerprintScoreKeys), "load top-k score kernel");
#ifdef SMSD_CUDA_BINARY_MMA_PROTOTYPE
        int device = 0;
        state.check(cudaGetDevice(&device), "cudaGetDevice");
        if (binaryMmaRequested(state, device))
            state.check(cudaFuncGetAttributes(&attributes, fingerprintBinaryMmaKernel), "load binary MMA kernel");
#endif
    });
}

/** Whether this build and the current device can exercise the opt-in binary MMA prototype. */
SMSD_CUDA_EXPORT int smsd_cuda_binary_mma_available(char* error, size_t errorSize) {
#ifdef SMSD_CUDA_BINARY_MMA_PROTOTYPE
    int available = 0;
    const int status = guarded(error, errorSize, [&](ErrorState& state) {
        int device = 0;
        state.check(cudaGetDevice(&device), "cudaGetDevice");
        cudaDeviceProp properties;
        state.check(cudaGetDeviceProperties(&properties, device), "binary MMA device properties");
        available = properties.major >= 8 && properties.major < 10;
    });
    return status == 0 ? available : -1;
#else
    (void)error;
    (void)errorSize;
    return 0;
#endif
}

/** The device count, or -1 with the CUDA failure written to error. */
SMSD_CUDA_EXPORT int smsd_cuda_device_count(char* error, size_t errorSize) {
    int count = 0;
    const cudaError_t status = cudaGetDeviceCount(&count);
    if (status != cudaSuccess) { describe(error, errorSize, "cudaGetDeviceCount", status); return -1; }
    return count;
}

/** 1 when the current device can run the screens; 0 with the reason, or -1 with the CUDA failure, in error. */
SMSD_CUDA_EXPORT int smsd_cuda_device_supported(char* error, size_t errorSize) {
    int supported = 0;
    const int status = guarded(error, errorSize, [&](ErrorState& state) {
        int device = 0;
        state.check(cudaGetDevice(&device), "cudaGetDevice");
        supported = deviceSupported(state, device, error, errorSize) ? 1 : 0;
    });
    return status == 0 ? supported : -1;
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
