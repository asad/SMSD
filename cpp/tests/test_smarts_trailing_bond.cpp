/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/smarts_parser.hpp"

#include <cstddef>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

int main() {
    int failures = 0, rejected = 0, valid = 0;
    auto requireInvalid = [&](const std::string& pattern, const char* context, const std::string& bond,
                              const std::string& message = "Trailing bond in SMARTS") {
        try {
            (void)smsd::parseSMARTS(pattern);
            ++failures;
            std::cout << "FAIL " << context << " trailing bond " << bond << ": accepted malformed SMARTS\n";
        } catch (const std::invalid_argument& error) {
            if (error.what() != message) {
                ++failures;
                std::cout << "FAIL " << context << " trailing bond " << bond << ": wrong message: " << error.what() << '\n';
                return;
            }
            ++rejected;
            std::cout << "PASS " << context << " trailing bond " << bond << '\n';
        } catch (const std::exception& error) {
            ++failures;
            std::cout << "FAIL " << context << " trailing bond " << bond << ": wrong exception: " << error.what() << '\n';
        } catch (...) {
            ++failures;
            std::cout << "FAIL " << context << " trailing bond " << bond << ": unknown exception\n";
        }
    };
    auto requireValid = [&](const std::string& pattern, int expectedAtoms, std::size_t expectedBonds, const char* context,
                            const std::string& bond) {
        try {
            const auto query = smsd::parseSMARTS(pattern);
            if (query.atomCount() != expectedAtoms || query.bonds.size() != expectedBonds) {
                ++failures;
                std::cout << "FAIL " << context << " valid bond " << bond << ": changed query shape\n";
            } else {
                ++valid;
                std::cout << "PASS " << context << " valid bond " << bond << '\n';
            }
        } catch (const std::exception& error) {
            ++failures;
            std::cout << "FAIL " << context << " valid bond " << bond << ": " << error.what() << '\n';
        } catch (...) {
            ++failures;
            std::cout << "FAIL " << context << " valid bond " << bond << ": unknown exception\n";
        }
    };

    for (const std::string bond : {"-", "=", "#", ":", "~", "@", "/", "\\", "/?", "\\?"}) {
        const std::string dangling = "C" + bond;
        const std::string recursive = "[$(" + dangling + ")]";
        const std::pair<const char*, std::string> malformed[] = {
            {"literal", dangling}, {"whitespace", dangling + " \t\r\n"},
            {"recursive", recursive}, {"nested recursive", "[$(" + recursive + ")]"}
        };
        for (const auto& input : malformed) requireInvalid(input.second, input.first, bond);

        const std::string complete = bond == ":" ? "c:c" : "C" + bond + "C";
        const std::string completeRecursive = "[$(" + complete + ")]";
        requireValid(complete, 2, 1, "literal", bond);
        requireValid(complete + " \t\r\n", 2, 1, "whitespace", bond);
        requireValid(completeRecursive, 1, 0, "recursive", bond);
        requireValid("[$(" + completeRecursive + ")]", 1, 0, "nested recursive", bond);
    }
    requireValid("C", 1, 0, "ordinary atom", "none");
    requireValid("", 0, 0, "empty SMARTS", "none");

    // Other unfinished patterns keep their earlier messages word for word.
    const std::pair<const char*, const char*> unfinished[] = {
        {"C(", "Unmatched '(' in SMARTS"}, {"C1CC", "Unclosed ring in SMARTS"},
        {"[C", "Expected ']' in bracket atom at pos 2"}, {"C)C", "Unmatched ')' in SMARTS"},
        {"CC-,=", "Trailing bond in SMARTS"}, {"CC!-", "Trailing bond in SMARTS"},
        {"CC-;@", "Trailing bond in SMARTS"}
    };
    for (const auto& input : unfinished) requireInvalid(input.first, "unfinished", input.first, input.second);
    // A bond written at a ring-closure digit is not trailing.
    for (const char* closure : {"C-1CC1", "C1CC-1", "C=1CC1", "C1CC=1"})
        requireValid(closure, 3, 3, "ring-closure", closure);
    std::cout << "Trailing-bond cases rejected=" << rejected << "/47, valid=" << valid
              << "/46, failures=" << failures << '\n';
    return failures == 0 ? 0 : 1;
}
