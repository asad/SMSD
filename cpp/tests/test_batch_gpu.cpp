/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * Consolidated batch processing and GPU test suite.
 * OpenMP parallel batch, fingerprints, screen-and-match, and GPU/Metal tests.
 *
 * Merges: test_batch.cpp, test_gpu.cpp
 */

#include "smsd/smsd.hpp"
#include "smsd/batch.hpp"
#include "smsd/smiles_parser.hpp"

#if defined(SMSD_ENABLE_CUDA) || defined(SMSD_ENABLE_METAL)
#include "smsd/gpu.hpp"
#define HAS_GPU 1
#else
#define HAS_GPU 0
#endif

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// ============================================================================
// Test harness (self-registering)
// ============================================================================

static int g_pass = 0, g_fail = 0;

#define TEST(name)                                                            \
    static void test_##name();                                                \
    struct Register_##name {                                                   \
        Register_##name() { tests().push_back({#name, test_##name}); }        \
    };                                                                        \
    static Register_##name reg_##name;                                        \
    static void test_##name()

#define ASSERT_TRUE(expr)                                                     \
    do {                                                                       \
        if (!(expr)) {                                                         \
            std::cerr << "  FAIL: " << #expr                                  \
                      << " at " << __FILE__ << ":" << __LINE__ << "\n";       \
            throw std::runtime_error("assertion failed");                      \
        }                                                                      \
    } while (0)

#define ASSERT_EQ(a, b) ASSERT_TRUE((a) == (b))

#define CHECK(cond, msg)                                           \
    do {                                                           \
        if (cond) {                                                \
            ++g_pass;                                              \
        } else {                                                   \
            ++g_fail;                                              \
            std::fprintf(stderr, "FAIL  %s  (%s:%d)\n",           \
                         (msg), __FILE__, __LINE__);               \
        }                                                          \
    } while (0)

struct TestEntry {
    std::string name;
    void (*func)();
};

static std::vector<TestEntry>& tests() {
    static std::vector<TestEntry> t;
    return t;
}

// ============================================================================
// Molecule builders
// ============================================================================

static void finalize(smsd::MolGraph& g, int bondOrd, bool bondRing, bool bondArom) {
    int n = g.n;
    int words = (n + 63) >> 6;
    g.adjLong.assign(n, std::vector<uint64_t>(words, 0));
    g.useDense = (n <= 200);
    if (g.useDense) {
        g.bondOrdMatrix.assign(n, std::vector<int>(n, 0));
        g.bondRingMatrix.assign(n, std::vector<bool>(n, false));
        g.bondAromMatrix.assign(n, std::vector<bool>(n, false));
    }
    for (int i = 0; i < n; ++i) {
        for (int j : g.neighbors[i]) {
            g.adjLong[i][j >> 6] |= uint64_t(1) << (j & 63);
            if (g.useDense) {
                g.bondOrdMatrix[i][j] = bondOrd;
                g.bondRingMatrix[i][j] = bondRing;
                g.bondAromMatrix[i][j] = bondArom;
            }
        }
    }
    g.morganRank.assign(n, 0);
    g.tautomerClass.assign(n, -1);
    g.orbit.resize(n);
    g.ringCount.assign(n, 0);
    g.tetraChirality.assign(n, 0);
    g.canonicalLabel.resize(n);
    for (int i = 0; i < n; ++i) {
        g.morganRank[i] = g.label[i];
        g.orbit[i] = i;
        g.canonicalLabel[i] = i;
    }
}

static smsd::MolGraph makeBenzene() {
    smsd::MolGraph g;
    g.n = 6;
    g.atomicNum = {6,6,6,6,6,6};
    g.formalCharge = {0,0,0,0,0,0};
    g.aromatic = {true,true,true,true,true,true};
    g.ring = {true,true,true,true,true,true};
    g.massNumber = {0,0,0,0,0,0};
    g.label.resize(6);
    g.degree.resize(6);
    g.neighbors.resize(6);
    for (int i = 0; i < 6; ++i) {
        g.label[i] = (6 << 2) | 2 | 1;
        g.neighbors[i] = {(i+5)%6, (i+1)%6};
        g.degree[i] = 2;
    }
    finalize(g, 4, true, true);
    return g;
}

static smsd::MolGraph makePhenol() {
    smsd::MolGraph g;
    g.n = 7;
    g.atomicNum = {6,6,6,6,6,6,8};
    g.formalCharge = {0,0,0,0,0,0,0};
    g.aromatic = {true,true,true,true,true,true,false};
    g.ring = {true,true,true,true,true,true,false};
    g.massNumber = {0,0,0,0,0,0,0};
    g.label.resize(7);
    g.degree.resize(7);
    g.neighbors.resize(7);
    for (int i = 0; i < 6; ++i) {
        g.label[i] = (6 << 2) | 2 | 1;
        g.neighbors[i] = {(i+5)%6, (i+1)%6};
        g.degree[i] = 2;
    }
    g.label[6] = (8 << 2);
    g.neighbors[0].push_back(6);
    g.neighbors[6] = {0};
    g.degree[0] = 3;
    g.degree[6] = 1;
    finalize(g, 1, false, false);
    for (int i = 0; i < 6; ++i) {
        int j = (i + 1) % 6;
        g.bondOrdMatrix[i][j] = 4;
        g.bondOrdMatrix[j][i] = 4;
        g.bondRingMatrix[i][j] = true;
        g.bondRingMatrix[j][i] = true;
        g.bondAromMatrix[i][j] = true;
        g.bondAromMatrix[j][i] = true;
    }
    return g;
}

static smsd::MolGraph makeChain(int n) {
    smsd::MolGraph g;
    g.n = n;
    g.atomicNum.assign(n, 6);
    g.formalCharge.assign(n, 0);
    g.aromatic.assign(n, false);
    g.ring.assign(n, false);
    g.massNumber.assign(n, 0);
    g.label.resize(n);
    g.degree.resize(n);
    g.neighbors.resize(n);
    for (int i = 0; i < n; ++i) {
        g.label[i] = (6 << 2);
        g.neighbors[i] = {};
        if (i > 0) g.neighbors[i].push_back(i - 1);
        if (i < n - 1) g.neighbors[i].push_back(i + 1);
        g.degree[i] = static_cast<int>(g.neighbors[i].size());
    }
    finalize(g, 1, false, false);
    return g;
}

static smsd::MolGraph makeNaphthalene() {
    smsd::MolGraph g;
    g.n = 10;
    g.atomicNum.assign(10, 6);
    g.formalCharge.assign(10, 0);
    g.aromatic.assign(10, true);
    g.ring.assign(10, true);
    g.massNumber.assign(10, 0);
    g.label.resize(10);
    g.degree.resize(10);
    g.neighbors.resize(10);
    for (int i = 0; i < 10; ++i) g.label[i] = (6 << 2) | 2 | 1;
    int edges[][2] = {
        {0,1},{1,2},{2,3},{3,4},{4,5},{5,0},
        {5,6},{6,7},{7,8},{8,9},{9,4}
    };
    for (auto& e : edges) {
        g.neighbors[e[0]].push_back(e[1]);
        g.neighbors[e[1]].push_back(e[0]);
    }
    for (int i = 0; i < 10; ++i)
        g.degree[i] = static_cast<int>(g.neighbors[i].size());
    finalize(g, 4, true, true);
    return g;
}

static smsd::MolGraph parse(const char* smi) {
    return smsd::parseSMILES(smi);
}

// ============================================================================
// SECTION 1: Batch processing tests (from test_batch.cpp)
// ============================================================================

TEST(batchSubstructure_basic) {
    auto benzene = makeBenzene();
    auto phenol = makePhenol();
    auto chain3 = makeChain(3);
    smsd::ChemOptions opts;
    opts.matchBondOrder = smsd::ChemOptions::BondOrderMode::LOOSE;
    opts.aromaticityMode = smsd::ChemOptions::AromaticityMode::FLEXIBLE;
    opts.ringMatchesRingOnly = false;
    std::vector<smsd::MolGraph> targets = {phenol, chain3};
    auto results = smsd::batch::batchSubstructure(benzene, targets, opts);
    ASSERT_EQ(results.size(), 2u);
    ASSERT_TRUE(results[0]);
    ASSERT_TRUE(!results[1]);
}

TEST(batchSubstructure_100_molecules) {
    auto query = makeChain(3);
    smsd::ChemOptions opts;
    opts.matchBondOrder = smsd::ChemOptions::BondOrderMode::LOOSE;
    opts.ringMatchesRingOnly = false;
    std::vector<smsd::MolGraph> targets;
    targets.reserve(100);
    for (int i = 0; i < 100; ++i)
        targets.push_back(makeChain(i + 2));
    auto results = smsd::batch::batchSubstructure(query, targets, opts);
    ASSERT_EQ(results.size(), 100u);
    ASSERT_TRUE(!results[0]);
    for (int i = 1; i < 100; ++i)
        ASSERT_TRUE(results[i]);
}

TEST(batchSubstructure_empty) {
    auto query = makeBenzene();
    smsd::ChemOptions opts;
    std::vector<smsd::MolGraph> empty;
    auto results = smsd::batch::batchSubstructure(query, empty, opts);
    ASSERT_EQ(results.size(), 0u);
}

TEST(batchMCS_basic) {
    auto benzene = makeBenzene();
    auto phenol = makePhenol();
    smsd::ChemOptions chem;
    chem.matchBondOrder = smsd::ChemOptions::BondOrderMode::LOOSE;
    chem.aromaticityMode = smsd::ChemOptions::AromaticityMode::FLEXIBLE;
    chem.ringMatchesRingOnly = false;
    smsd::MCSOptions mcsOpts;
    mcsOpts.timeoutMs = 5000;
    std::vector<smsd::MolGraph> targets = {phenol};
    auto results = smsd::batch::batchMCS(benzene, targets, chem, mcsOpts);
    ASSERT_EQ(results.size(), 1u);
    ASSERT_TRUE(results[0].size() >= 6u);
}

TEST(batchMCS_50_pairs) {
    auto query = makeChain(4);
    smsd::ChemOptions chem;
    chem.matchBondOrder = smsd::ChemOptions::BondOrderMode::LOOSE;
    chem.ringMatchesRingOnly = false;
    smsd::MCSOptions mcsOpts;
    mcsOpts.timeoutMs = 2000;
    std::vector<smsd::MolGraph> targets;
    targets.reserve(50);
    for (int i = 0; i < 50; ++i)
        targets.push_back(makeChain(i + 3));
    auto results = smsd::batch::batchMCS(query, targets, chem, mcsOpts);
    ASSERT_EQ(results.size(), 50u);
    for (int i = 0; i < 50; ++i) {
        int targetN = i + 3;
        int expectedMCS = std::min(4, targetN);
        ASSERT_TRUE(static_cast<int>(results[i].size()) >= expectedMCS);
    }
}

TEST(screenAndMatch_pipeline) {
    auto query = makeBenzene();
    smsd::ChemOptions chem;
    chem.matchBondOrder = smsd::ChemOptions::BondOrderMode::LOOSE;
    chem.aromaticityMode = smsd::ChemOptions::AromaticityMode::FLEXIBLE;
    chem.ringMatchesRingOnly = false;
    smsd::MCSOptions mcsOpts;
    mcsOpts.timeoutMs = 5000;
    std::vector<smsd::MolGraph> targets;
    targets.push_back(makePhenol());
    targets.push_back(makeNaphthalene());
    targets.push_back(makeChain(3));
    targets.push_back(makeChain(5));
    auto results = smsd::batch::screenAndMatch(query, targets, chem, mcsOpts, 0.1);
    bool foundPhenol = false, foundNaph = false;
    for (auto& [idx, mapping] : results) {
        if (idx == 0) foundPhenol = true;
        if (idx == 1) foundNaph = true;
        ASSERT_TRUE(!mapping.empty());
    }
    ASSERT_TRUE(foundPhenol);
    ASSERT_TRUE(foundNaph);
}

TEST(batchFingerprint_basic) {
    std::vector<smsd::MolGraph> mols;
    mols.push_back(makeBenzene());
    mols.push_back(makePhenol());
    mols.push_back(makeChain(5));
    auto fps = smsd::batch::batchFingerprint(mols, 7, 1024);
    ASSERT_EQ(fps.size(), 3u);
    for (auto& fp : fps) ASSERT_EQ(fp.size(), 16u);
    for (auto& fp : fps) {
        uint64_t orAll = 0;
        for (auto w : fp) orAll |= w;
        ASSERT_TRUE(orAll != 0);
    }
}

TEST(batchFingerprint_tanimoto) {
    auto benzene = makeBenzene();
    auto phenol = makePhenol();
    auto chain10 = makeChain(10);
    auto fp_benz = smsd::batch::detail::computePathFingerprint(benzene, 7, 1024);
    auto fp_phen = smsd::batch::detail::computePathFingerprint(phenol, 7, 1024);
    auto fp_chain = smsd::batch::detail::computePathFingerprint(chain10, 7, 1024);
    double sim_benz_phen = smsd::batch::fingerprintTanimoto(fp_benz, fp_phen);
    double sim_benz_chain = smsd::batch::fingerprintTanimoto(fp_benz, fp_chain);
    ASSERT_TRUE(sim_benz_phen > sim_benz_chain);
    ASSERT_TRUE(sim_benz_phen >= 0.0 && sim_benz_phen <= 1.0);
    ASSERT_TRUE(sim_benz_chain >= 0.0 && sim_benz_chain <= 1.0);
}

TEST(batchFingerprint_subset) {
    auto benzene = makeBenzene();
    auto phenol = makePhenol();
    auto fp_benz = smsd::batch::detail::computePathFingerprint(benzene, 7, 1024);
    auto fp_phen = smsd::batch::detail::computePathFingerprint(phenol, 7, 1024);
    ASSERT_TRUE(smsd::batch::fingerprintSubset(fp_benz, fp_benz));
    ASSERT_TRUE(smsd::batch::fingerprintSubset(fp_phen, fp_phen));
    std::vector<uint64_t> empty;
    ASSERT_TRUE(smsd::batch::fingerprintSubset(empty, fp_benz));
}

TEST(batchFingerprintScreen_cpu) {
    auto query = makeChain(3);
    auto fp_query = smsd::batch::detail::computePathFingerprint(query, 7, 1024);
    std::vector<smsd::MolGraph> targets;
    for (int i = 0; i < 20; ++i)
        targets.push_back(makeChain(i + 2));
    auto targetFPs = smsd::batch::batchFingerprint(targets, 7, 1024);
    auto hits = smsd::batch::batchFingerprintScreen(fp_query, targetFPs);
    ASSERT_TRUE(!hits.empty());
}

TEST(openmp_sequential_consistency) {
    auto query = makeChain(3);
    smsd::ChemOptions opts;
    opts.matchBondOrder = smsd::ChemOptions::BondOrderMode::LOOSE;
    opts.ringMatchesRingOnly = false;
    std::vector<smsd::MolGraph> targets;
    for (int i = 0; i < 30; ++i)
        targets.push_back(makeChain(i + 2));
    auto seq = smsd::batch::batchSubstructure(query, targets, opts, 1);
    auto par = smsd::batch::batchSubstructure(query, targets, opts, 0);
    ASSERT_EQ(seq.size(), par.size());
    for (size_t i = 0; i < seq.size(); ++i)
        ASSERT_EQ(seq[i], par[i]);
}

TEST(performance_batch_1000) {
    auto query = makeChain(4);
    smsd::ChemOptions opts;
    opts.matchBondOrder = smsd::ChemOptions::BondOrderMode::LOOSE;
    opts.ringMatchesRingOnly = false;
    std::vector<smsd::MolGraph> targets;
    targets.reserve(1000);
    for (int i = 0; i < 1000; ++i)
        targets.push_back(makeChain((i % 10) + 3));
    auto t0 = std::chrono::high_resolution_clock::now();
    auto results = smsd::batch::batchSubstructure(query, targets, opts);
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "    batch substruct 1000 molecules: " << ms << " ms\n";
    ASSERT_TRUE(ms < 1000.0);
    ASSERT_EQ(results.size(), 1000u);
}

TEST(performance_batch_fingerprint_1000) {
    std::vector<smsd::MolGraph> mols;
    mols.reserve(1000);
    for (int i = 0; i < 1000; ++i)
        mols.push_back(makeChain((i % 15) + 3));
    auto t0 = std::chrono::high_resolution_clock::now();
    auto fps = smsd::batch::batchFingerprint(mols, 7, 1024);
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "    batch fingerprint 1000 molecules: " << ms << " ms\n";
    ASSERT_TRUE(ms < 1000.0);
    ASSERT_EQ(fps.size(), 1000u);
}

TEST(batchMCS_empty_inputs) {
    smsd::MolGraph empty;
    empty.n = 0;
    smsd::ChemOptions chem;
    smsd::MCSOptions mcsOpts;
    std::vector<smsd::MolGraph> targets = {makeBenzene()};
    auto r1 = smsd::batch::batchMCS(empty, targets, chem, mcsOpts);
    ASSERT_EQ(r1.size(), 1u);
    ASSERT_TRUE(r1[0].empty());
    std::vector<smsd::MolGraph> noTargets;
    auto r2 = smsd::batch::batchMCS(makeBenzene(), noTargets, chem, mcsOpts);
    ASSERT_EQ(r2.size(), 0u);
}

TEST(screenAndMatch_high_threshold) {
    auto query = makeBenzene();
    smsd::ChemOptions chem;
    chem.matchBondOrder = smsd::ChemOptions::BondOrderMode::LOOSE;
    chem.ringMatchesRingOnly = false;
    smsd::MCSOptions mcsOpts;
    std::vector<smsd::MolGraph> targets;
    targets.push_back(makeChain(3));
    targets.push_back(makeChain(5));
    auto results = smsd::batch::screenAndMatch(query, targets, chem, mcsOpts, 0.99);
    ASSERT_TRUE(results.empty());
}

TEST(fingerprint_different_sizes) {
    auto mol = makeBenzene();
    auto fp512 = smsd::batch::detail::computePathFingerprint(mol, 7, 512);
    auto fp2048 = smsd::batch::detail::computePathFingerprint(mol, 7, 2048);
    ASSERT_EQ(fp512.size(), 8u);
    ASSERT_EQ(fp2048.size(), 32u);
    uint64_t or512 = 0, or2048 = 0;
    for (auto w : fp512) or512 |= w;
    for (auto w : fp2048) or2048 |= w;
    ASSERT_TRUE(or512 != 0);
    ASSERT_TRUE(or2048 != 0);
}

// ============================================================================
// SECTION 2: GPU / Metal tests (from test_gpu.cpp)
// Compiled only when SMSD_ENABLE_CUDA or SMSD_ENABLE_METAL is defined.
// ============================================================================

#if HAS_GPU

static std::vector<smsd::gpu::ScreenHit> gpuScreen(
    const smsd::MolGraph&              q,
    const std::vector<smsd::MolGraph>& ts,
    double                             threshold = 0.0)
{
    auto hits = smsd::gpu::batchRascalScreenAuto(q, ts, threshold);
    std::sort(hits.begin(), hits.end(),
              [](const smsd::gpu::ScreenHit& a, const smsd::gpu::ScreenHit& b) {
                  return a.index < b.index;
              });
    return hits;
}

static std::vector<double> cpuRef(
    const smsd::MolGraph&              q,
    const std::vector<smsd::MolGraph>& ts)
{
    std::vector<double> r;
    r.reserve(ts.size());
    for (auto& t : ts)
        r.push_back(smsd::similarityUpperBound(q, t));
    return r;
}

static void runGpuTests() {
    std::printf("\n  -- GPU / Metal tests --\n");

    // Device info
    {
        std::string info = smsd::gpu::deviceInfo();
        CHECK(!info.empty(), "deviceInfo() non-empty");
        bool ok = info.find("GPU:") != std::string::npos
               || info.find("CPU:") != std::string::npos
               || info.find("Metal") != std::string::npos;
        CHECK(ok, "deviceInfo() contains GPU/CPU/Metal prefix");
        std::printf("  deviceInfo: %s\n", info.c_str());
        std::printf("  isAvailable: %s\n", smsd::gpu::isAvailable() ? "true" : "false");
    }

    // Correctness (small)
    {
        const std::vector<std::string> smiles = {
            "c1ccccc1", "Oc1ccccc1", "CC(=O)Oc1ccccc1C(=O)O",
            "C", "CC", "c1ccc2ccccc2c1",
        };
        std::vector<smsd::MolGraph> mols;
        for (auto& s : smiles) mols.push_back(parse(s.c_str()));
        auto query = parse("c1ccccc1");
        auto hits = gpuScreen(query, mols, 0.0);
        auto ref  = cpuRef(query, mols);
        CHECK(hits.size() == mols.size(), "all targets returned at threshold=0");
        bool scores_match = true;
        for (auto& h : hits) {
            double diff = std::abs(h.score - ref[h.index]);
            if (diff > 0.02) {
                scores_match = false;
                std::fprintf(stderr, "  score mismatch at index %d: gpu=%.4f cpu=%.4f\n",
                             h.index, h.score, ref[h.index]);
            }
        }
        CHECK(scores_match, "GPU scores match CPU reference (2% tolerance)");
    }

    // Threshold filtering
    {
        std::vector<smsd::MolGraph> targets = {
            parse("c1ccccc1"), parse("C"), parse("c1ccc2ccccc2c1"),
        };
        auto query = parse("c1ccccc1");
        auto hits = gpuScreen(query, targets, 0.5);
        bool has_benzene = false, has_naphthalene = false, has_methane = false;
        for (auto& h : hits) {
            if (h.index == 0) has_benzene = true;
            if (h.index == 1) has_methane = true;
            if (h.index == 2) has_naphthalene = true;
        }
        CHECK(has_benzene, "threshold=0.5: benzene included");
        CHECK(has_naphthalene, "threshold=0.5: naphthalene included");
        CHECK(!has_methane, "threshold=0.5: methane excluded");
    }

    // Empty targets
    {
        auto query = parse("c1ccccc1");
        std::vector<smsd::MolGraph> empty;
        auto hits = gpuScreen(query, empty, 0.0);
        CHECK(hits.empty(), "empty targets -> empty hits");
    }

    // Self-match
    {
        auto mol = parse("CC(C)Cc1ccc(CC(C)C(=O)O)cc1");
        auto hits = gpuScreen(mol, {mol}, 0.0);
        CHECK(hits.size() == 1, "self-match returns 1 hit");
        if (!hits.empty())
            CHECK(std::abs(hits[0].score - 1.0) < 0.01, "self-match score ~1.0");
    }

    // Large batch
    {
        std::vector<std::string> pool = {
            "c1ccccc1", "Oc1ccccc1", "CC(=O)O", "c1ccncc1", "C1CCCCC1",
            "CC(C)C", "CCCO", "c1ccc(O)cc1", "CCCCC", "c1ccc2ccccc2c1",
        };
        std::vector<smsd::MolGraph> targets;
        targets.reserve(500);
        for (int i = 0; i < 500; ++i)
            targets.push_back(parse(pool[i % pool.size()].c_str()));
        auto query = parse("c1ccccc1");
        auto hits  = gpuScreen(query, targets, 0.0);
        CHECK(hits.size() == 500, "large batch: 500 hits at threshold=0");
    }
}

#endif // HAS_GPU

// ============================================================================
// main
// ============================================================================

int main() {
    std::cout << "=== SMSD Batch & GPU Test Suite ===\n";
    std::cout << "===================================\n";

#ifdef _OPENMP
    std::cout << "OpenMP: enabled (max threads = "
              << omp_get_max_threads() << ")\n";
#else
    std::cout << "OpenMP: not available (sequential fallback)\n";
#endif
    std::cout << "\n";

    // Run self-registered batch tests
    for (auto& t : tests()) {
        std::cout << "  [RUN]  " << t.name << "\n";
        try {
            t.func();
            std::cout << "  [PASS] " << t.name << "\n";
            g_pass++;
        } catch (const std::exception& e) {
            std::cout << "  [FAIL] " << t.name << ": " << e.what() << "\n";
            g_fail++;
        }
    }

#if HAS_GPU
    runGpuTests();
#else
    std::cout << "\n  -- GPU tests skipped (no GPU backend compiled) --\n";
#endif

    std::cout << "\n" << g_pass << " passed, " << g_fail << " failed.\n";
    return g_fail > 0 ? 1 : 0;
}
