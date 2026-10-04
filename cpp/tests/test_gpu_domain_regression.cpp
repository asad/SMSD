/* SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * A deterministic kernel double exercises accelerator domains on every platform.
 */
#ifndef SMSD_TEST_NATIVE_GPU
#define SMSD_ENABLE_METAL
#endif
#include "smsd/vf2pp.hpp"
#include <iostream>
#include <stdexcept>

#ifndef SMSD_TEST_NATIVE_GPU
namespace smsd::gpu_kern {
static bool available = true;
static int calls = 0;
bool graphKernelsAvailable() noexcept { return available; }
bool domainInit(int nq, int nt, int words,
    const int* qz, const int* qc, const int* qr, const int* qa,
    const int* tz, const int* tc, const int* tr, const int* ta,
    bool type, bool charge, bool ring, bool aromatic, uint64_t* out) {
    ++calls;
    for (int i = 0; i < nq; ++i)
        for (int j = 0; j < nt; ++j)
            if ((!type || qz[i] == tz[j]) && (!charge || qc[i] == tc[j])
                && (!ring || !qr[i] || tr[j]) && (!aromatic || qa[i] == ta[j]))
                out[i * words + (j >> 6)] |= uint64_t(1) << (j & 63);
    return true;
}
}
#endif

static smsd::MolGraph atoms(int n, int z) {
    return smsd::MolGraph::Builder().atomCount(n)
        .atomicNumbers(std::vector<int>(n, z))
        .setNeighbors(std::vector<std::vector<int>>(n)).build(false);
}

#ifdef SMSD_TEST_NATIVE_GPU
static void verifyNativeKernel(const smsd::MolGraph& q, const smsd::MolGraph& t,
                               const smsd::ChemOptions& c) {
    const int words = (t.n + 63) / 64;
    std::vector<int> qr(q.n), qa(q.n), tr(t.n), ta(t.n);
    for (int i = 0; i < q.n; ++i) { qr[i] = q.ring[i]; qa[i] = q.aromatic[i]; }
    for (int j = 0; j < t.n; ++j) { tr[j] = t.ring[j]; ta[j] = t.aromatic[j]; }
    std::vector<uint64_t> domain(static_cast<size_t>(q.n) * words, 0);
    const bool type = c.matchAtomType && !c.tautomerAware;
    const bool aromatic = c.aromaticityMode == smsd::ChemOptions::AromaticityMode::STRICT;
    if (!smsd::gpu_kern::domainInit(q.n, t.n, words,
            q.atomicNum.data(), q.formalCharge.data(), qr.data(), qa.data(),
            t.atomicNum.data(), t.formalCharge.data(), tr.data(), ta.data(),
            type, c.matchFormalCharge, c.ringMatchesRingOnly, aromatic, domain.data()))
        throw std::runtime_error("native domain kernel failed; CPU fallback cannot verify GPU");
    for (int i = 0; i < q.n; ++i)
        for (int j = 0; j < words * 64; ++j) {
            const bool expected = j < t.n
                && (!type || q.atomicNum[i] == t.atomicNum[j])
                && (!c.matchFormalCharge || q.formalCharge[i] == t.formalCharge[j])
                && (!c.ringMatchesRingOnly || !qr[i] || tr[j])
                && (!aromatic || qa[i] == ta[j]);
            const bool actual = (domain[static_cast<size_t>(i) * words + (j >> 6)]
                                 & (uint64_t(1) << (j & 63))) != 0;
            if (actual != expected)
                throw std::runtime_error("native domain kernel differs from basic compatibility oracle");
        }
}
#endif

static void check(const smsd::MolGraph& q, const smsd::MolGraph& t,
                  smsd::ChemOptions c, bool expected, const char* name) {
    for (auto engine : {smsd::ChemOptions::MatcherEngine::VF2,
                        smsd::ChemOptions::MatcherEngine::VF2PP}) {
        c.matcherEngine = engine;
#ifndef SMSD_TEST_NATIVE_GPU
        smsd::gpu_kern::available = false;
        const bool cpu = smsd::isSubstructure(q, t, c);
        smsd::gpu_kern::available = true;
        const int before = smsd::gpu_kern::calls;
        const bool gpu = smsd::isSubstructure(q, t, c);
        if (cpu != expected || gpu != cpu || smsd::gpu_kern::calls == before)
            throw std::runtime_error(name);
#else
        verifyNativeKernel(q, t, c);
        if (smsd::isSubstructure(q, t, c) != expected)
            throw std::runtime_error(name);
#endif
    }
}

int main() {
    try {
#ifdef SMSD_TEST_NATIVE_GPU
        if (!smsd::gpu_kern::graphKernelsAvailable()) {
            std::cout << "Native GPU unavailable\n";
            return 77;
        }
#endif
        smsd::ChemOptions c;
        auto q = atoms(33, 6), t = atoms(65, 6);
        q.massNumber[0] = 13;
        std::fill(t.massNumber.begin(), t.massNumber.end(), 12);
        c.matchIsotope = true;
        check(q, t, c, false, "accelerator must preserve isotope rejection");

        q = atoms(33, 6); t = atoms(65, 6); c = {};
        q.tetraChirality[0] = 1;
        std::fill(t.tetraChirality.begin(), t.tetraChirality.end(), 2);
        c.useChirality = true;
        check(q, t, c, false, "accelerator must preserve atom chirality constraints");

        q = atoms(33, 6); t = atoms(65, 6); c = {};
        std::fill(t.ring.begin(), t.ring.end(), 1);
        c.ringMatchesRingOnly = true;
        check(q, t, c, false, "accelerator ring-domain policy must match CPU");

        q = atoms(33, 6); t = atoms(65, 7); c = {};
        q.tautomerClass.assign(q.n, 0); t.tautomerClass.assign(t.n, 0);
        c.tautomerAware = true;
        check(q, t, c, true, "accelerator must retain cross-element tautomer candidates");
        std::cout << "CPU/accelerator domain policy parity passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
