/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Stereo kept by extractSubgraph and R-group decomposition, from the fixture also read by the
 * Java and Python tests.
 */
#include "smsd/smsd.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef SMSD_SUBGRAPH_STEREO_JSON
#error "SMSD_SUBGRAPH_STEREO_JSON must name the shared fixture"
#endif

namespace {

/** A JSON value; only what the fixture uses. */
struct Json {
    enum class Kind { Null, Bool, Number, String, Array, Object } kind = Kind::Null;
    bool boolean = false;
    double number = 0;
    std::string text;
    std::vector<Json> items;
    std::map<std::string, Json> members;

    const Json& operator[](const std::string& key) const {
        const auto found = members.find(key);
        if (kind != Kind::Object || found == members.end()) throw std::runtime_error("Missing key " + key);
        return found->second;
    }
};

class JsonReader {
    const std::string& source_;
    size_t at_ = 0;

    void space() { while (at_ < source_.size() && std::isspace(static_cast<unsigned char>(source_[at_]))) ++at_; }
    char peek() { space(); if (at_ >= source_.size()) throw std::runtime_error("Unexpected end of JSON"); return source_[at_]; }
    void expect(char c) { if (peek() != c) throw std::runtime_error(std::string("Expected ") + c); ++at_; }
    std::string string() {
        expect('"');
        std::string value;
        while (source_.at(at_) != '"') {
            if (source_[at_] == '\\') {
                const char escaped = source_.at(++at_);
                value += escaped == 'n' ? '\n' : escaped == 't' ? '\t' : escaped;
            } else {
                value += source_[at_];
            }
            ++at_;
        }
        ++at_;
        return value;
    }

public:
    explicit JsonReader(const std::string& source) : source_(source) {}

    Json value() {
        Json result;
        const char c = peek();
        if (c == '{') {
            result.kind = Json::Kind::Object;
            ++at_;
            if (peek() == '}') { ++at_; return result; }
            for (;;) {
                const std::string key = string();
                expect(':');
                result.members[key] = value();
                if (peek() == ',') { ++at_; continue; }
                expect('}');
                return result;
            }
        }
        if (c == '[') {
            result.kind = Json::Kind::Array;
            ++at_;
            if (peek() == ']') { ++at_; return result; }
            for (;;) {
                result.items.push_back(value());
                if (peek() == ',') { ++at_; continue; }
                expect(']');
                return result;
            }
        }
        if (c == '"') { result.kind = Json::Kind::String; result.text = string(); return result; }
        if (source_.compare(at_, 4, "true") == 0) { at_ += 4; result.kind = Json::Kind::Bool; result.boolean = true; return result; }
        if (source_.compare(at_, 5, "false") == 0) { at_ += 5; result.kind = Json::Kind::Bool; return result; }
        if (source_.compare(at_, 4, "null") == 0) { at_ += 4; return result; }
        size_t used = 0;
        result.kind = Json::Kind::Number;
        result.number = std::stod(source_.substr(at_, 32), &used);
        at_ += used;
        return result;
    }
};

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL " << message << '\n';
    }
}

bool hasStereo(const std::string& smiles) {
    return smiles.find('@') != std::string::npos || smiles.find('/') != std::string::npos
        || smiles.find('\\') != std::string::npos;
}

/** The canonical SMILES of a fragment against that of the SMILES it should equal. */
void compare(const smsd::MolGraph& fragment, const std::string& expected, bool stereo, const std::string& name) {
    const std::string actual = smsd::toSMILES(fragment);
    const std::string wanted = smsd::toSMILES(smsd::parseSMILES(expected));
    check(actual == wanted, name + ": expected " + wanted + ", got " + actual);
    check(hasStereo(actual) == stereo, name + ": stereo " + (stereo ? "lost" : "invented") + " in " + actual);
}

void subgraphs(const Json& cases) {
    for (const auto& item : cases.items) {
        const auto parent = smsd::parseSMILES(item["smiles"].text);
        check(hasStereo(smsd::toSMILES(parent)), item["name"].text + ": the parent carries no stereo");
        std::vector<int> atoms;
        for (const auto& atom : item["atoms"].items) atoms.push_back(static_cast<int>(atom.number));
        compare(smsd::extractSubgraph(parent, atoms), item["expected"].text, item["stereo"].boolean, item["name"].text);
    }
}

void rgroups(const Json& cases) {
    for (const auto& item : cases.items) {
        const auto results = smsd::Substructure::decomposeRGroups(
            smsd::parseSMILES(item["core"].text), {smsd::parseSMILES(item["smiles"].text)});
        check(results.size() == 1 && !results.front().search.matches.empty(), item["name"].text + ": the core does not match");
        if (results.size() != 1) continue;
        std::vector<std::string> actual, wanted;
        bool stereo = false;
        for (const auto& group : results.front().rgroups) {
            actual.push_back(smsd::toSMILES(group.second));
            stereo = stereo || hasStereo(actual.back());
        }
        for (const auto& expected : item["expected"].items) wanted.push_back(smsd::toSMILES(smsd::parseSMILES(expected.text)));
        std::sort(actual.begin(), actual.end());
        std::sort(wanted.begin(), wanted.end());
        check(actual == wanted, item["name"].text + ": the R-groups differ");
        check(stereo == item["stereo"].boolean, item["name"].text + ": R-group stereo differs");
    }
}

/** A centre whose implicit hydrogen changes slot with the atom order keeps its configuration. */
void permutedImplicitHydrogen() {
    const auto parent = smsd::parseSMILES("N[C@@H](C)C(=O)O");
    const auto fragment = smsd::extractSubgraph(parent, {3, 4, 1, 0, 2});
    const auto inOrder = smsd::extractSubgraph(parent, {0, 1, 2, 3, 4});
    check(fragment.tetraChirality[2] != 0 && inOrder.tetraChirality[1] != 0, "permuted centre: tag lost");
    check(smsd::toSMILES(fragment) == smsd::toSMILES(inOrder), "permuted centre: configuration changed with atom order");
}

} // namespace

int main() {
    std::ifstream input(SMSD_SUBGRAPH_STEREO_JSON);
    if (!input) {
        std::cerr << "FAIL cannot read " << SMSD_SUBGRAPH_STEREO_JSON << '\n';
        return 1;
    }
    std::stringstream buffer;
    buffer << input.rdbuf();
    const std::string source = buffer.str();
    try {
        const Json fixture = JsonReader(source).value();
        subgraphs(fixture["subgraphs"]);
        rgroups(fixture["rgroups"]);
        permutedImplicitHydrogen();
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
    if (failures != 0) {
        std::cerr << failures << " subgraph stereo check(s) failed\n";
        return 1;
    }
    std::cout << "subgraph stereo checks passed\n";
    return 0;
}
