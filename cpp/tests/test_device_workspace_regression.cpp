/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/gpu_kernels.hpp"
#include <cassert>
#include <cstdlib>
#include <iostream>

struct Runtime {
    static bool failNext;
    static bool leavePointerOnFailure;
    static int allocations, releases;
    static bool allocate(void** pointer, size_t bytes) {
        ++allocations;
        if (failNext) { failNext = false; if (!leavePointerOnFailure) *pointer = nullptr; return false; }
        *pointer = std::malloc(bytes);
        return *pointer != nullptr;
    }
    static void release(void* pointer) { ++releases; std::free(pointer); }
};
bool Runtime::failNext = false;
bool Runtime::leavePointerOnFailure = false;
int Runtime::allocations = 0;
int Runtime::releases = 0;

int main() {
    {
        smsd::gpu_kern::detail::DeviceWorkspace<Runtime> workspace;
        Runtime::failNext = true;
        assert(!workspace.ensure(64) && workspace.ptr == nullptr && workspace.capacity == 0);
        assert(workspace.ensure(16));
        void* original = workspace.ptr;
        Runtime::failNext = true;
        assert(!workspace.ensure(64));
        assert(workspace.ptr == original && workspace.capacity == 16 && Runtime::releases == 0);
        Runtime::leavePointerOnFailure = true;
        Runtime::failNext = true;
        assert(!workspace.ensure(64) && workspace.ptr == original && workspace.capacity == 16);
        const int attempts = Runtime::allocations;
        assert(workspace.ensure(8) && workspace.ptr == original && Runtime::allocations == attempts);
        assert(workspace.ensure(64) && workspace.ptr != original && workspace.capacity == 64);
        assert(Runtime::releases == 1);
    }
    assert(Runtime::releases == 2);
    std::cout << "PASS failed allocation, preserved capacity, reuse, recovery and release\n";
}
