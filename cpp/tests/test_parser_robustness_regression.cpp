/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "timing_allowance.hpp"  // also defines SMSD_TEST_SANITISED
#include "smsd/smiles_parser.hpp"
#include "smsd/smarts_parser.hpp"
#include "smsd/batch.hpp"
#include "smsd/mol_reader.hpp"
#include "smsd/ring_finder.hpp"
#include "smsd/fingerprint_options.hpp"

#include <chrono>
#include <clocale>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
#include <thread>
#if defined(__unix__) || defined(__APPLE__)
#include <pthread.h>
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
// Sanitiser runtimes supply their own nothrow forms, so route these through the same allocator.
void* operator new(std::size_t size, const std::nothrow_t&) noexcept {
    try { return ::operator new(size); } catch (...) { return nullptr; }
}
void* operator new[](std::size_t size, const std::nothrow_t&) noexcept {
    try { return ::operator new[](size); } catch (...) { return nullptr; }
}
void operator delete(void* memory, const std::nothrow_t&) noexcept { std::free(memory); }
void operator delete[](void* memory, const std::nothrow_t&) noexcept { std::free(memory); }

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
    require(elapsed < smsd_test::timingAllowance(5000), "parsing a 40-ring acene took too long");
}

// The Kekule structure of an aromatic macrocycle near the atom limit, which perception assigns
// to every supplied aromatic system, needs no deep call stack: it parses on a thread of 256 KB,
// twice the smallest default thread stack of common platforms.
static void testAromaticMacrocycleOnASmallStack() {
    std::string failure;
    std::function<void()> parse = [&failure] {
        try {
            const auto ring = smsd::parseSMILES("c1" + std::string(smsd::MAX_ATOMS - 4, 'c') + "c1");
            int aromatic = 0;
            for (int atom = 0; atom < ring.n; ++atom) aromatic += ring.aromatic[atom] != 0;
            if (ring.n != smsd::MAX_ATOMS - 2 || aromatic != ring.n) failure = "an aromatic macrocycle lost aromaticity";
        } catch (const std::exception& error) {
            failure = error.what();
        }
    };
#if defined(__unix__) || defined(__APPLE__)
#if defined(SMSD_TEST_SANITISED)
    const size_t stackBytes = 4096 * 1024;  // sanitiser frames are several times larger
#else
    const size_t stackBytes = 256 * 1024;
#endif
    pthread_attr_t attributes;
    require(pthread_attr_init(&attributes) == 0 && pthread_attr_setstacksize(&attributes, stackBytes) == 0,
            "cannot set a thread stack size");
    pthread_t worker;
    void* (*run)(void*) = [](void* body) -> void* {
        (*static_cast<std::function<void()>*>(body))();
        return nullptr;
    };
    const int created = pthread_create(&worker, &attributes, run, &parse);
    pthread_attr_destroy(&attributes);
    require(created == 0, "cannot start a thread");
    pthread_join(worker, nullptr);
#else
    std::thread(parse).join();
#endif
    if (!failure.empty()) throw std::runtime_error(failure);
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

// Bytes above 0x7F belong to neither language; every position that classifies a character
// rejects them, and lenient SMILES recovery either recovers or rejects them as well.
static void testHighBytes() {
    smsd::ParseOptions lenient;
    lenient.lenient = true;
    for (int byte = 0x80; byte <= 0xFF; ++byte) {
        const std::string high(1, static_cast<char>(byte));
        for (const std::string& text : {"C" + high, high + "C", "[CH" + high + "]", "[C+" + high + "]",
                                        "C%" + high + high + "C", "C1CC" + high, "[" + high + "]"}) {
            requireInvalid([&] { smsd::parseSMILES(text); });
            try { smsd::parseSMILES(text, lenient); } catch (const std::invalid_argument&) {}
        }
        for (const std::string& text : {"C" + high, high + "C", "[" + high + "]", "[C" + high + "]",
                                        "[#6]" + high, "[C;" + high + "]", "C%" + high + "C"})
            requireInvalid([&] { smsd::parseSMARTS(text); });
    }
}

/** A copy of text with its one occurrence of `from` replaced by `to`. */
static std::string replaced(const std::string& text, const std::string& from, const std::string& to) {
    const auto at = text.find(from);
    require(at != std::string::npos && text.find(from, at + 1) == std::string::npos, "fixture text not unique");
    return text.substr(0, at) + to + text.substr(at + from.size());
}

// V3000 integer fields must be integers in range, coordinates decimal numbers, and a block may hold
// at most MAX_ATOMS atoms; anything else is a std::invalid_argument, as for the V2000 reader.
static void testV3000Fields() {
    const std::string block = smsd::writeMolBlockV3000(smsd::parseSMILES("C[NH3+]"));
    const auto graph = smsd::readMolBlock(block);
    require(graph.n == 2 && graph.formalCharge[1] == 1, "V3000 fixture no longer reads");
    const std::string atom = "M  V30 1 C 0 0 0 0\n", bond = "M  V30 1 1 1 2\n";
    for (const std::string& charge : {"99999999999", "-99999999999", "1x", "1.5", "", "+"})
        requireInvalid([&] { smsd::readMolBlock(replaced(block, "CHG=1", "CHG=" + charge)); });
    for (const std::string& field : {"MASS=", "RAD=", "HCOUNT=", "VAL=", "CLASS=", "CFG="})
        requireInvalid([&] { smsd::readMolBlock(replaced(block, "CHG=1", field + "4294967297")); });
    for (const std::string& line : {"M  V30 99999999999 C 0 0 0 0\n", "M  V30 1 C 1e999 0 0 0\n",
                                    "M  V30 1 C nan 0 0 0\n", "M  V30 1 C 0 inf 0 0\n", "M  V30 1 C 0 0 0x 0\n"})
        requireInvalid([&] { smsd::readMolBlock(replaced(block, atom, line)); });
    for (const std::string& line : {"M  V30 1 x 1 2\n", "M  V30 1 1 1 99999999999\n", "M  V30 1 1 -1x 2\n"})
        requireInvalid([&] { smsd::readMolBlock(replaced(block, bond, line)); });
    const auto exact = smsd::readMolBlock(replaced(block, atom, "M  V30 1 C 1.25 -2.5E1 +0 0\n"));
    require(exact.n == 2, "V3000 coordinates in exponent form no longer read");
    // Coordinates are decimal numbers read by one rule on every platform: a magnitude below 1e-307
    // reads as zero, and one beyond the range of a double, or a number in another form, is refused.
    for (const std::string& line : {"M  V30 1 C 1e309 0 0 0\n", "M  V30 1 C -1.8e308 0 0 0\n", "M  V30 1 C 0x10 0 0 0\n",
                                    "M  V30 1 C infinity 0 0 0\n", "M  V30 1 C . 0 0 0\n", "M  V30 1 C 1e 0 0 0\n"})
        requireInvalid([&] { smsd::readMolBlock(replaced(block, atom, line)); });
    require(smsd::readMolBlock(replaced(block, atom, "M  V30 1 C 1e-310 -4.9e-324 .5 0\n")).n == 2,
            "V3000 coordinates below the normal range no longer read");
    using smsd::mol_detail::v3000Coordinate;
    require(v3000Coordinate("1e-310") == 0.0 && std::signbit(v3000Coordinate("-4.9e-324"))
                && v3000Coordinate("1E-307") == 1e-307 && v3000Coordinate("1.7e308") == 1.7e308
                && v3000Coordinate("0.00012e3") == 0.12 && v3000Coordinate("-2.5E1") == -25.0 && v3000Coordinate("5.") == 5.0,
            "V3000 coordinates read differently");
    // Atom properties must lie in the ranges the writers can give them.
    for (const std::string& field : {"CHG=2147483647", "CHG=-2147483648", "CHG=100", "CHG=-100", "MASS=1000",
                                     "MASS=-1", "RAD=4", "HCOUNT=2147483647", "HCOUNT=100", "HCOUNT=-2",
                                     "VAL=16", "VAL=2147483647", "CFG=4"})
        requireInvalid([&] { smsd::readMolBlock(replaced(block, "CHG=1", field)); });
    for (int charge : {-99, 99})
        require(smsd::readMolBlock(replaced(block, "CHG=1", "CHG=" + std::to_string(charge))).formalCharge[1] == charge,
                "a V3000 charge within the SMILES range no longer reads");
    requireInvalid([&] { smsd::readMolBlock(replaced(block, bond, "M  V30 1 1 1 2 CFG=4\n")); });

    std::string atoms;
    for (int index = 1; index <= smsd::MAX_ATOMS; ++index)
        atoms += "M  V30 " + std::to_string(index) + " C 0 0 0 0\n";
    const auto head = block.substr(0, block.find("M  V30 BEGIN ATOM\n") + 18);
    const std::string tail = "M  V30 END ATOM\nM  V30 END CTAB\nM  END\n";
    require(smsd::readMolBlock(head + atoms + tail).n == smsd::MAX_ATOMS, "V3000 block at MAX_ATOMS rejected");
    const std::string extra = "M  V30 " + std::to_string(smsd::MAX_ATOMS + 1) + " C 0 0 0 0\n";
    requireInvalid([&] { smsd::readMolBlock(head + atoms + extra + tail); });
}

// Charges and hydrogen counts beyond what the readers accept, as a builder may set them, are
// written and fingerprinted without overflow, and a MOL valence field that cannot hold the
// hydrogen count is refused rather than overflowed.
static void testExtremeAtomProperties() {
    auto graphOf = [](int charge, int hydrogens) {
        return smsd::MolGraph::Builder().atomCount(2).atomicNumbers({6, 8}).formalCharges({charge, 0})
            .hydrogenCounts({hydrogens, 1}).setNeighbors({{1}, {0}}).setBondOrders({{1}, {1}}).build();
    };
    for (int charge : {std::numeric_limits<int>::min(), std::numeric_limits<int>::max()}) {
        const auto graph = graphOf(charge, 0);
        require(smsd::toSMILES(graph).find(std::to_string(charge < 0 ? -int64_t(charge) : int64_t(charge)))
                    != std::string::npos, "an extreme charge was not written");
        smsd::batch::FingerprintOptions path;
        path.family = smsd::batch::FingerprintFamily::PATH;
        require(!smsd::batch::computeFingerprint(graph, path).empty(), "an extreme charge has no path fingerprint");
    }
    requireInvalid([&] { smsd::writeMolBlock(graphOf(0, std::numeric_limits<int>::max())); });
}

// Coordinates are read in the classic locale, as the format writes them, whatever locale the
// program has set; the test is skipped where no decimal-comma locale is installed.
static void testCoordinatesIgnoreTheLocale() {
    const std::string previous = std::setlocale(LC_NUMERIC, nullptr);
    const char* found = nullptr;
    for (const char* name : {"de_DE.UTF-8", "de_DE.utf8", "fr_FR.UTF-8", "fr_FR.utf8"})
        if ((found = std::setlocale(LC_NUMERIC, name)) != nullptr) break;
    if (found == nullptr) return;
    std::string failure;
    try {
        const std::string v3000 = smsd::writeMolBlockV3000(smsd::parseSMILES("CC"));
        const std::string atom = "M  V30 1 C 0 0 0 0\n";
        const auto at = v3000.find(atom);
        if (at == std::string::npos) failure = "V3000 fixture changed";
        else {
            const auto fractional = smsd::readMolBlock(v3000.substr(0, at) + "M  V30 1 C 1.5 -2.25 0 0\n"
                                                       + v3000.substr(at + atom.size()));
            if (fractional.n != 2) failure = "a fractional V3000 coordinate no longer reads";
        }
    } catch (const std::exception& error) {
        failure = std::string("a fractional V3000 coordinate was rejected: ") + error.what();
    }
    std::setlocale(LC_NUMERIC, previous.c_str());
    if (!failure.empty()) throw std::runtime_error(failure);
}

int main(int argc, char** argv) {
    return smsd_test::runGroups(argc, argv, {
        {"branches", testDeepBranches}, {"atom-negation", testAtomNegation},
        {"bond-negation", testBondNegation}, {"numbers", testNumericBounds},
        {"bond-precedence", testBondPrecedence}, {"sparse-stereo", testSparseStereo},
        {"smarts-limits", testSmartsLimits}, {"large-rings", testLargeRings},
        {"aromatic-macrocycle-stack", testAromaticMacrocycleOnASmallStack},
        {"high-bytes", testHighBytes}, {"v3000-fields", testV3000Fields},
        {"extreme-atom-properties", testExtremeAtomProperties}, {"coordinate-locale", testCoordinatesIgnoreTheLocale}});
}
