/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/smarts_parser.hpp"
#include "smsd/batch.hpp"
#include "smsd/ring_finder.hpp"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
#include <thread>

#if defined(__has_feature)
#if __has_feature(address_sanitizer) || __has_feature(undefined_behavior_sanitizer)
#define SMSD_TEST_SANITISED 1
#endif
#endif
#if defined(__SANITIZE_ADDRESS__) && !defined(SMSD_TEST_SANITISED)
#define SMSD_TEST_SANITISED 1
#endif

/** The time allowed to perceive the large rings of one test, in milliseconds. */
static long long largeRingLimitMillis() {
#if defined(SMSD_TEST_SANITISED) || !(defined(__OPTIMIZE__) || (defined(_MSC_VER) && !defined(_DEBUG)))
    return 90000;
#else
    return 20000;
#endif
}

namespace {
bool limitAllocations = false;
size_t remainingAllocation = 0;
}

void* operator new(std::size_t size) {
    if (limitAllocations) {
        if (size > remainingAllocation) throw std::bad_alloc();
        remainingAllocation -= size;
    }
    if (void* memory = std::malloc(size == 0 ? 1 : size)) return memory;
    throw std::bad_alloc();
}
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete(void* memory) noexcept { std::free(memory); }
void operator delete[](void* memory) noexcept { std::free(memory); }
void operator delete(void* memory, std::size_t) noexcept { std::free(memory); }
void operator delete[](void* memory, std::size_t) noexcept { std::free(memory); }

using smsd_test::require;

template<class Function> static void requireInvalid(Function&& function) {
    try { function(); }
    catch (const std::invalid_argument&) { return; }
    throw std::runtime_error("out-of-range or incomplete input was accepted");
}

static void testDeepBranches() {
    constexpr int depth = 100000;
    const auto nested = smsd::parseSMILES("C" + std::string(depth, '(') + "O" + std::string(depth, ')') + "N");
    require(nested.n == 3 && nested.bondOrder(0, 1) == 1 && nested.bondOrder(0, 2) == 1,
            "nested branches changed their parent atom");
    requireInvalid([&] { smsd::parseSMILES(std::string(depth, '(')); });
    smsd::ParseOptions options;
    options.lenient = true;
    require(smsd::parseSMILES("C" + std::string(depth, '(') + "O", options).n == 2,
            "lenient recovery lost an unfinished branch");
}

static void testAtomNegation() {
    const auto carbon = smsd::parseSMILES("C");
    const auto oxygen = smsd::parseSMILES("O");
    // The longest negation chains that fit MAX_SMARTS_LENGTH.
    for (int count : {99996, 99997}) {
        const auto query = smsd::parseSMARTS("[" + std::string(count, '!') + "C]");
        require(smsd_test::patternMatches(query, carbon) == (count % 2 == 0), "atom negation parity changed");
        require(smsd_test::patternMatches(query, oxygen) == (count % 2 != 0), "negated atom matched incorrectly");
    }
}

static void testBondNegation() {
    const auto single = smsd::parseSMILES("CC");
    const auto multiple = smsd::parseSMILES("C=C");
    for (int count : {99996, 99997}) {
        const auto query = smsd::parseSMARTS("C" + std::string(count, '!') + "-C");
        require(smsd_test::patternMatches(query, single) == (count % 2 == 0), "bond negation parity changed");
        require(smsd_test::patternMatches(query, multiple) == (count % 2 != 0), "negated bond matched incorrectly");
    }
}

static void testNumericBounds() {
    for (const char* text : {"[99999999999C]", "[CH99999999999]", "[C+99999999999]",
            "[C-99999999999]", "[C:99999999999]", "[R99999999999]"})
        requireInvalid([&] { smsd::parseSMILES(text); });
    for (const char* text : {"[#99999999999]", "[H99999999999]", "[X99999999999]",
            "[+99999999999]", "[-99999999999]", "[99999999999C]", "[:99999999999]",
            "[D{1-99999999999}]", "C%99999999999CC%99999999999"})
        requireInvalid([&] { smsd::parseSMARTS(text); });
    // OpenSMILES allows one digit for a hydrogen count, two for a charge within -15 to +15 and
    // three for an isotope; the atom class is unbounded.
    for (const char* text : {"[CH10]", "[CH00]", "[CH2147483647]", "[C+16]", "[C-16]", "[C+128]", "[C+016]",
            "[C++++++++++++++++]", "[1000C]", "[65536C]", "[2147483647C+2147483647:2147483647]", "[O+2147483647]"})
        requireInvalid([&] { smsd::parseSMILES(text); });
    const auto bounds = smsd::parseSMILES("[999CH9+15:2147483647]");
    require(bounds.massNumber[0] == 999 && bounds.hydrogenCount[0] == 9 && bounds.formalCharge[0] == 15
            && bounds.atomClass[0] == std::numeric_limits<int>::max(), "OpenSMILES bracket limits changed");
    require(smsd::parseSMILES("[C-15]").formalCharge[0] == -15 && smsd::parseSMILES("[C+++++++++++++++]").formalCharge[0] == 15
            && smsd::parseSMILES("[C+00]").formalCharge[0] == 0, "charge limits changed");
    require(smsd::toSMILES(bounds) == "[999CH9+15:2147483647]", "bracket atom at the limits did not round trip");
    // SMARTS keeps its own, wider numeric fields, which a graph built directly can hold.
    const auto largeHydrogen = smsd::MolGraph::Builder().atomCount(1).atomicNumbers({6})
        .hydrogenCounts({std::numeric_limits<int>::max()}).setNeighbors({{}}).build();
    require(smsd_test::patternMatches(smsd::parseSMARTS("[H2147483647;X2147483647;v2147483647]"), largeHydrogen),
            "representable boundary properties changed");
    const auto largeCharge = smsd::MolGraph::Builder().atomCount(1).atomicNumbers({8})
        .formalCharges({std::numeric_limits<int>::max()}).massNumbers({std::numeric_limits<int>::max()})
        .atomClasses({std::numeric_limits<int>::max()}).hydrogenCounts({0}).setNeighbors({{}}).build();
    require(smsd_test::patternMatches(smsd::parseSMARTS("[2147483647O+2147483647:2147483647]"), largeCharge),
            "boundary isotope, charge and class query changed");
    require(smsd_test::patternMatches(smsd::parseSMARTS("[O+2147483647;H0;X0;v0]"), largeCharge), "boundary charge changed atom properties");
}

// SMARTS length, atoms (recursive SMARTS included) and expression nodes are bounded, and the
// bounds are inclusive.
static void testSmartsLimits() {
    static_assert(smsd::MAX_SMARTS_LENGTH == 100000 && smsd::MAX_SMARTS_EXPRESSION_NODES == 100000
                  && smsd::MAX_SMARTS_NODES_PER_EXPRESSION == 10000 && smsd::MAX_ATOMS == 10000,
                  "documented input limits changed");
    const auto chain = smsd::parseSMILES(std::string(20, 'C'));
    require(smsd_test::patternMatches(smsd::parseSMARTS("[" + std::string(99997, '!') + "C]"), smsd::parseSMILES("N")),
            "the longest SMARTS was rejected");
    requireInvalid([] { smsd::parseSMARTS("[" + std::string(99998, '!') + "C]"); });
    require(smsd::parseSMARTS(std::string(10000, 'C')).atomCount() == 10000, "the largest SMARTS was rejected");
    requireInvalid([] { smsd::parseSMARTS(std::string(10001, 'C')); });
    require(smsd::parseSMARTS("[$(" + std::string(9999, 'C') + ")]").atomCount() == 1,
            "recursive SMARTS at the atom limit was rejected");
    requireInvalid([] { smsd::parseSMARTS("[$(" + std::string(10000, 'C') + ")]"); });
    // [C R R ...] has one primitive per symbol and one implicit AND between neighbours, so
    // 4999 Rs make 9999 nodes.
    const std::string ringAtom = "[C" + std::string(4999, 'R') + "]";
    require(smsd_test::patternMatches(smsd::parseSMARTS(ringAtom), smsd::parseSMILES("C1CC1")),
            "an expression within the node limit was rejected");
    require(!smsd_test::patternMatches(smsd::parseSMARTS(ringAtom), chain), "a ring expression matched a chain");
    requireInvalid([] { smsd::parseSMARTS("[C" + std::string(5000, 'R') + "]"); });
    // Ten such atoms and their nine bonds make 99999 nodes; three more nodes exceed the pattern limit.
    std::string tenAtoms;
    for (int atom = 0; atom < 10; ++atom) tenAtoms += ringAtom;
    require(smsd::parseSMARTS(tenAtoms).atomCount() == 10, "a pattern at the node limit was rejected");
    requireInvalid([&] { smsd::parseSMARTS(tenAtoms + ".[CR]"); });
    // Recursive SMARTS count their own nodes.
    std::string nested = "[";
    for (int part = 0; part < 11; ++part) nested += (part ? ";$(" : "$(") + ringAtom + ")";
    requireInvalid([&] { smsd::parseSMARTS(nested + "]"); });
}

// Ring perception of very large rings needs no deep call stack and bounded work; it runs on
// a secondary thread, whose stack is 512 KB on macOS.
static void testLargeRings() {
    std::string failure;
    std::thread worker([&] {
        try {
            const auto started = std::chrono::steady_clock::now();
            const auto ring = smsd::parseSMILES("C1" + std::string(smsd::MAX_ATOMS - 2, 'C') + "C1");
            if (ring.n != smsd::MAX_ATOMS || ring.computeRings().size() != 1 || ring.computeRelevantCycles().size() != 1
                || smsd::computeSSSR(ring).size() != 1 || smsd::computeRelevantCycles(ring).size() != 1)
                failure = "a macrocycle at the atom limit lost its ring";
            // A macrocycle through a benzene ring is one ring system of two rings.
            const auto fused = smsd::parseSMILES("C1" + std::string(2994, 'C') + "c2ccccc12");
            if (fused.computeRings().size() != 2 || smsd::computeSSSR(fused).size() != 2
                || smsd::computeRelevantCycles(fused).size() != 2)
                failure = "a fused macrocycle lost a ring";
            const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - started).count();
            // Perception with unbounded work would take minutes on a ring of 4,000 atoms, even in an
            // optimised build; sanitised and unoptimised builds run the bounded work several times
            // more slowly, so their limit is wider.
            if (failure.empty() && elapsed > largeRingLimitMillis())
                failure = "large-ring perception took " + std::to_string(elapsed) + " ms";
        } catch (const std::exception& error) {
            failure = error.what();
        }
    });
    worker.join();
    if (!failure.empty()) throw std::runtime_error(failure);
    // A 40-ring acene, a long fused ladder, keeps a bounded aromatic-cycle walk; an unbounded
    // walk would take tens of seconds for 30 rings.
    constexpr int rings = 40;
    auto label = [](int ring) { return ring < 10 ? std::to_string(ring) : "%" + std::to_string(ring); };
    std::string acene = "c" + label(1) + "ccc" + label(2);
    for (int ring = 3; ring <= rings; ++ring) acene += "cc" + label(ring);
    acene += "cccc";
    for (int ring = rings; ring >= 3; --ring) acene += "c" + label(ring) + "c";
    acene += "c" + label(2) + "c" + label(1);
    const auto started = std::chrono::steady_clock::now();
    const auto ladder = smsd::parseSMILES(acene);
    const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - started).count();
    require(ladder.n == 4 * rings + 2 && ladder.computeRings().size() == rings
            && ladder.computeRelevantCycles().size() == rings, "a long acene lost a ring");
    for (int atom = 0; atom < ladder.n; ++atom) require(ladder.aromatic[atom] != 0, "a long acene lost aromaticity");
    require(elapsed < 5000, "parsing a 40-ring acene took too long");
}

static void testBondPrecedence() {
    struct Example { const char* query; const char* target; bool expected; };
    for (const auto& example : {
            Example{"C=,#@C", "C=C", true}, {"C=,#@C", "C#C", false},
            {"C=,#@C", "C1#CCCCC1", true}, {"C=,#@C", "CC", false},
            {"C-&!@C", "CC", true}, {"C-&!@C", "C1CC1", false},
            {"C-,=;@C", "C=C", false}, {"C-,=;@C", "C1CC1", true},
            {"C-,=;@C", "C1=CCCCC1", true}})
        require(smsd_test::patternMatches(smsd::parseSMARTS(example.query), smsd::parseSMILES(example.target)) == example.expected,
                "bond Boolean expression has wrong precedence");
    for (const char* text : {"C-&C", "C-,C", "C-;C", "C!C"})
        requireInvalid([&] { smsd::parseSMARTS(text); });
}

static void testSparseStereo() {
    const std::string text = "C/C=C/" + std::string(4093, 'C');
    remainingAllocation = 32 * 1024 * 1024;
    limitAllocations = true;
    try {
        const auto graph = smsd::parseSMILES(text);
        limitAllocations = false;
        require(graph.n == 4096 && graph.dbStereo(1, 2) == 2, "large chain lost its alkene configuration");
    } catch (...) {
        limitAllocations = false;
        throw;
    }
    require(smsd::parseSMILES("C/C=C/C").dbStereo(1, 2) == 2, "small graph alkene configuration changed");
}

int main(int argc, char** argv) {
    return smsd_test::runGroups(argc, argv, {
        {"branches", testDeepBranches}, {"atom-negation", testAtomNegation},
        {"bond-negation", testBondNegation}, {"numbers", testNumericBounds},
        {"bond-precedence", testBondPrecedence}, {"sparse-stereo", testSparseStereo},
        {"smarts-limits", testSmartsLimits}, {"large-rings", testLargeRings}});
}
