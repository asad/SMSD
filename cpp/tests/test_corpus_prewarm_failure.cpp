/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/batch.hpp"

#include <atomic>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <new>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

static std::atomic<bool> failWorkerAllocation{false};
static std::atomic<int> failedAllocations{0};

void* operator new(std::size_t bytes) {
#ifdef _OPENMP
    if (omp_in_parallel() && failWorkerAllocation.exchange(false)) {
        ++failedAllocations;
        throw std::bad_alloc();
    }
#endif
    if (void* allocation = std::malloc(bytes == 0 ? 1 : bytes)) return allocation;
    throw std::bad_alloc();
}
void* operator new[](std::size_t bytes) { return ::operator new(bytes); }
void operator delete(void* allocation) noexcept { std::free(allocation); }
void operator delete[](void* allocation) noexcept { std::free(allocation); }
void operator delete(void* allocation, std::size_t) noexcept { std::free(allocation); }
void operator delete[](void* allocation, std::size_t) noexcept { std::free(allocation); }
// Sanitiser runtimes supply their own nothrow forms, so route these through the same allocator.
void* operator new(std::size_t size, const std::nothrow_t&) noexcept {
    try { return ::operator new(size); } catch (...) { return nullptr; }
}
void* operator new[](std::size_t size, const std::nothrow_t&) noexcept {
    try { return ::operator new[](size); } catch (...) { return nullptr; }
}
void operator delete(void* allocation, const std::nothrow_t&) noexcept { std::free(allocation); }
void operator delete[](void* allocation, const std::nothrow_t&) noexcept { std::free(allocation); }

int main(int argc, char** argv) {
    const bool warmed = argc > 1 && std::string(argv[1]) == "warm";
    smsd::batch::TargetCorpus corpus;
    for (int target = 0; target < 64; ++target)
        corpus.addTarget(smsd::parseSMILES("CC(=O)Oc1ccccc1C(=O)O"));
    if (warmed) {
        corpus.prewarm(2);
        assert(corpus.isPrewarmed());
    }
#ifdef _OPENMP
    omp_set_dynamic(0);
    // Two workers even when OMP_NUM_THREADS=1, as in the wheel builds.
    omp_set_num_threads(2);
    #pragma omp parallel num_threads(2)
    {
        assert(omp_get_num_threads() == 2);
    }
    failWorkerAllocation = true;
    bool propagated = false;
    try {
        corpus.prewarm(2);
    } catch (const std::bad_alloc&) {
        propagated = true;
    }
    failWorkerAllocation = false;
    assert(propagated && failedAllocations == 1);
    assert(!corpus.isPrewarmed());
#endif
    corpus.prewarm(2);
    assert(corpus.isPrewarmed() && corpus.size() == 64);
    const auto parallel = corpus.screen(smsd::parseSMILES("CC(=O)Oc1ccccc1C(=O)O"), 0.0);
    assert(parallel.size() == 64);
    corpus.prewarm(1);
    const auto serial = corpus.screen(smsd::parseSMILES("CC(=O)Oc1ccccc1C(=O)O"), 0.0);
    assert(parallel.size() == serial.size());
    std::cout << "PASS corpus prewarm exception propagation and reuse " << (warmed ? "warm" : "cold") << '\n';
}
