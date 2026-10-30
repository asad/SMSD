/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * MOL and SDF readers give E/Z only to double bonds that can carry it, so a
 * drawn molecule has the canonical SMILES of the same molecule read as SMILES.
 * The cases are shared with the Java and Python tests.
 */
#include "smsd/mol_reader.hpp"
#include "smsd/smiles_parser.hpp"

#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#ifndef SMSD_MOL_EZ_ELIGIBILITY_JSON
#error "SMSD_MOL_EZ_ELIGIBILITY_JSON must name the shared fixture"
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

/** The number of double bonds that carry a configuration. */
int configuredBonds(const smsd::MolGraph& graph) {
    int count = 0;
    for (int atom = 0; atom < graph.n; ++atom)
        for (int neighbor : graph.neighbors[atom])
            if (atom < neighbor && graph.dbStereo(atom, neighbor) != 0) ++count;
    return count;
}

std::string block(const std::vector<std::string>& lines) {
    std::string text;
    for (const auto& line : lines) text += line + "\n";
    return text;
}

void sharedCases(const Json& fixture) {
    for (const auto& item : fixture["cases"].items) {
        const std::string name = item["name"].text;
        std::vector<std::string> lines;
        for (const auto& line : item["molblock"].items) lines.push_back(line.text);
        const auto drawn = smsd::readMolBlock(block(lines));
        const auto written = smsd::parseSMILES(item["smiles"].text);
        check(configuredBonds(drawn) == static_cast<int>(item["ez_bonds"].number),
              name + ": got " + std::to_string(configuredBonds(drawn)) + " configured double bonds");
        check(drawn.toCanonicalSmiles() == written.toCanonicalSmiles(),
              name + ": MOL gives " + drawn.toCanonicalSmiles() + ", SMILES gives " + written.toCanonicalSmiles());
    }
}

/** A V2000 block from atoms (symbol, x, y) and bonds (first, second, order), 1-based. */
std::string v2000(const std::vector<std::tuple<std::string, double, double>>& atoms,
                  const std::vector<std::array<int, 3>>& bonds, const std::string& tail = "") {
    std::ostringstream text;
    text << "drawn\n  SMSD              2D\n\n";
    char line[96];
    std::snprintf(line, sizeof(line), "%3d%3d  0  0  0  0  0  0  0  0999 V2000\n",
                  static_cast<int>(atoms.size()), static_cast<int>(bonds.size()));
    text << line;
    for (const auto& atom : atoms) {
        std::snprintf(line, sizeof(line), "%10.4f%10.4f%10.4f %-3s 0  0  0  0  0  0  0  0  0  0  0  0\n",
                      std::get<1>(atom), std::get<2>(atom), 0.0, std::get<0>(atom).c_str());
        text << line;
    }
    for (const auto& bond : bonds) {
        std::snprintf(line, sizeof(line), "%3d%3d%3d  0\n", bond[0], bond[1], bond[2]);
        text << line;
    }
    return text.str() + tail + "M  END\n";
}

void inlineCases() {
    // The V3000 reader shares the rule: a drawn cyclohexene has no E/Z.
    const std::string cyclohexene =
        "cyclohexene\n  SMSD              2D\n\n  0  0  0  0  0  0            999 V3000\n"
        "M  V30 BEGIN CTAB\nM  V30 COUNTS 6 6 0 0 0\nM  V30 BEGIN ATOM\n"
        "M  V30 1 C 1.5 0 0 0\nM  V30 2 C 0.75 1.299 0 0\nM  V30 3 C -0.75 1.299 0 0\n"
        "M  V30 4 C -1.5 0 0 0\nM  V30 5 C -0.75 -1.299 0 0\nM  V30 6 C 0.75 -1.299 0 0\n"
        "M  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 2 1 2\nM  V30 2 1 2 3\nM  V30 3 1 3 4\n"
        "M  V30 4 1 4 5\nM  V30 5 1 5 6\nM  V30 6 1 6 1\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
    const auto v3000 = smsd::readMolBlock(cyclohexene);
    check(configuredBonds(v3000) == 0, "V3000 cyclohexene must not carry E/Z");
    check(v3000.toCanonicalSmiles() == smsd::parseSMILES("C1CC=CCC1").toCanonicalSmiles(),
          "V3000 cyclohexene must give the SMILES-derived canonical SMILES");

    // A configuration written on a small-ring bond in SMILES does not survive a MOL round trip,
    // since the reader does not read one from the ring drawing.
    const auto tagged = smsd::parseSMILES("C1/C=C\\CCC1");
    for (const auto& text : {smsd::writeMolBlock(tagged), smsd::writeMolBlockV3000(tagged)}) {
        const auto restored = smsd::readMolBlock(text);
        check(configuredBonds(restored) == 0, "A cyclohexene MOL round trip must not carry E/Z");
        check(restored.toCanonicalSmiles() == smsd::parseSMILES("C1CC=CCC1").toCanonicalSmiles(),
              "A cyclohexene MOL round trip must give the unconfigured canonical SMILES");
    }

    // A seven-membered ring is still below eight atoms.
    std::vector<std::tuple<std::string, double, double>> ring;
    std::vector<std::array<int, 3>> ringBonds;
    for (int atom = 0; atom < 7; ++atom) {
        const double angle = 2 * std::acos(-1.0) * atom / 7;
        ring.emplace_back("C", 1.7 * std::cos(angle), 1.7 * std::sin(angle));
        ringBonds.push_back({atom + 1, (atom + 1) % 7 + 1, atom == 0 ? 2 : 1});
    }
    check(configuredBonds(smsd::readMolBlock(v2000(ring, ringBonds))) == 0, "cycloheptene must not carry E/Z");

    // Two hydrogens are alike only with the same mass number: =CHD keeps E/Z, =CD2 has none.
    const std::vector<std::array<int, 3>> propene = {{{1, 2, 2}}, {{2, 3, 1}}, {{1, 4, 1}}, {{1, 5, 1}}};
    const std::vector<std::tuple<std::string, double, double>> hd = {
        {"C", 0.0, 0.0}, {"C", 1.5, 0.0}, {"C", 2.25, 1.299}, {"H", -0.75, 1.299}, {"H", -0.75, -1.299}};
    check(configuredBonds(smsd::readMolBlock(v2000(hd, propene, "M  ISO  1   5   2\n"))) == 1,
          "1-deuteriopropene must keep its E/Z");
    check(configuredBonds(smsd::readMolBlock(v2000(hd, propene, "M  ISO  2   4   2   5   2\n"))) == 0,
          "1,1-dideuteriopropene must not carry E/Z");
    // One drawn hydrogen and one implicit hydrogen are alike.
    const std::vector<std::tuple<std::string, double, double>> oneDrawn = {
        {"C", 0.0, 0.0}, {"C", 1.5, 0.0}, {"C", 2.25, 1.299}, {"H", -0.75, 1.299}};
    check(configuredBonds(smsd::readMolBlock(v2000(oneDrawn, {{{1, 2, 2}}, {{2, 3, 1}}, {{1, 4, 1}}}))) == 0,
          "A methylene with one drawn hydrogen must not carry E/Z");
}

} // namespace

int main() {
    std::ifstream input(SMSD_MOL_EZ_ELIGIBILITY_JSON);
    if (!input) {
        std::cerr << "FAIL cannot read " << SMSD_MOL_EZ_ELIGIBILITY_JSON << '\n';
        return 1;
    }
    std::stringstream buffer;
    buffer << input.rdbuf();
    const std::string source = buffer.str();
    try {
        sharedCases(JsonReader(source).value());
        inlineCases();
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
    if (failures != 0) return 1;
    std::cout << "PASS MOL E/Z eligibility\n";
    return 0;
}
