#include "search_test_helpers.hpp"
/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#include "smsd/smiles_parser.hpp"
#include "smsd/smarts_parser.hpp"
#include "smsd/batch.hpp"

#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>

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

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

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
    for (int count : {100000, 100001}) {
        const auto query = smsd::parseSMARTS("[" + std::string(count, '!') + "C]");
        require(smsd_test::patternMatches(query, carbon) == (count % 2 == 0), "atom negation parity changed");
        require(smsd_test::patternMatches(query, oxygen) == (count % 2 != 0), "negated atom matched incorrectly");
    }
}

static void testBondNegation() {
    const auto single = smsd::parseSMILES("CC");
    const auto multiple = smsd::parseSMILES("C=C");
    for (int count : {100000, 100001}) {
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
    const auto largeHydrogen = smsd::parseSMILES("[CH2147483647]");
    require(largeHydrogen.hydrogenCount[0] == std::numeric_limits<int>::max(), "representable H count changed");
    require(smsd_test::patternMatches(smsd::parseSMARTS("[H2147483647;X2147483647;v2147483647]"), largeHydrogen),
            "representable boundary properties changed");
    requireInvalid([] { smsd::parseSMILES("[CH2147483647]C"); });
    requireInvalid([] { smsd::parseSMILES("[H][CH2147483647]"); });
    const auto largeClass = smsd::parseSMILES("[2147483647C+2147483647:2147483647]");
    require(largeClass.massNumber[0] == std::numeric_limits<int>::max()
            && largeClass.atomClass[0] == std::numeric_limits<int>::max()
            && largeClass.formalCharge[0] == std::numeric_limits<int>::max(), "representable numeric field changed");
    require(smsd_test::patternMatches(smsd::parseSMARTS("[2147483647C+2147483647:2147483647]"), largeClass),
            "boundary isotope, charge and class query changed");
    const auto largeCharge = smsd::parseSMILES("[O+2147483647]");
    require(smsd_test::patternMatches(smsd::parseSMARTS("[O+2147483647;H0;X0;v0]"), largeCharge), "boundary charge changed atom properties");
    require(!smsd::batch::fingerprintFromSmiles("[O+2147483647]").empty(), "boundary charge broke ECFP");
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
    int failures = 0;
    const std::pair<const char*, void (*)()> tests[] = {
        {"branches", testDeepBranches}, {"atom-negation", testAtomNegation},
        {"bond-negation", testBondNegation}, {"numbers", testNumericBounds},
        {"bond-precedence", testBondPrecedence}, {"sparse-stereo", testSparseStereo}
    };
    for (const auto& test : tests) {
        if (argc > 1 && argv[1] != std::string(test.first)) continue;
        try { test.second(); std::cout << "PASS " << test.first << '\n'; }
        catch (const std::exception& error) {
            ++failures;
            std::cerr << "FAIL " << test.first << ": " << error.what() << '\n';
        }
    }
    return failures == 0 ? 0 : 1;
}
