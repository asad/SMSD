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
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef _MSC_VER
#include <crtdbg.h>
#include <stdlib.h>
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

int main() {
#ifdef _MSC_VER
    // A failed check reports on stderr and ends the run instead of waiting on a dialog.
    _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDERR);
    _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDERR);
    _set_abort_behavior(0, _WRITE_ABORT_MSG | _CALL_REPORTFAULT);
#endif
#ifdef _OPENMP
    omp_set_dynamic(0);
    // Two workers even when OMP_NUM_THREADS=1, as in the wheel builds.
    omp_set_num_threads(2);
    #pragma omp parallel num_threads(2)
    {
        assert(omp_get_num_threads() == 2);
    }
#endif
    // The first index that fails is reported, as on one thread, and every index before it is visited.
    for (int workers : {1, 2}) {
        std::vector<int> visited(4096, 0);
        std::string reported;
        try {
            smsd::batch::detail::forEachIndex(4096, workers, 64, [&](int i) {
                visited[i] = 1;
                if (i == 1000 || i == 3000) throw std::runtime_error(std::to_string(i));
            });
        } catch (const std::runtime_error& failure) {
            reported = failure.what();
        }
        assert(reported == "1000");
        for (int i = 0; i <= 1000; ++i) assert(visited[i] == 1);
    }

    std::vector<smsd::MolGraph> molecules;
    for (int i = 0; i < 64; ++i)
        molecules.push_back(smsd::parseSMILES(i % 2 ? "CC(=O)Oc1ccccc1C(=O)O" : "c1ccc2ccccc2c1"));
    const auto serial = smsd::batch::batchFingerprint(molecules, 7, 1024, 1);
    // Failing the first allocation made inside a worker depends on where the standard library
    // allocates, which differs under MSVC, so this check runs with the other compilers.
#if defined(_OPENMP) && !defined(_MSC_VER)
    // A worker whose allocation fails reports to the caller instead of ending the process.
    failWorkerAllocation = true;
    bool propagated = false;
    try {
        smsd::batch::batchFingerprint(molecules, 7, 1024, 2);
    } catch (const std::bad_alloc&) {
        propagated = true;
    }
    failWorkerAllocation = false;
    assert(propagated && failedAllocations == 1);
#endif
    // Parallel runs return the serial results.
    assert(smsd::batch::batchFingerprint(molecules, 7, 1024, 2) == serial);
    std::vector<std::vector<uint64_t>> targets;
    for (int i = 0; i < 4096; ++i) targets.push_back(serial[i % serial.size()]);
    const auto subsets = smsd::batch::batchFingerprintScreen(serial[1], targets, 2);
    assert(subsets.size() == 2048 && subsets == smsd::batch::batchFingerprintScreen(serial[1], targets, 1));
    smsd::batch::TargetCorpus corpus;
    for (int i = 0; i < 2048; ++i) corpus.addTarget(molecules[i % molecules.size()]);
    corpus.prewarm(2);
    const auto query = smsd::parseSMILES("CC(=O)Oc1ccccc1C(=O)O");
    const auto similar = corpus.screen(query, 0.5, 2);
    assert(similar.size() == 1024 && similar == corpus.screen(query, 0.5, 1));
    std::cout << "PASS batch worker exception propagation and parallel results\n";
}
