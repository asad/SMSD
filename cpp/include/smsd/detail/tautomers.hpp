/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * MolGraph tautomer classes and pKa-informed relevance weights with their solvent
 * and pH corrections, and the tautomer confidence of a selected mapping.
 *
 * Part of mol_graph.hpp, which includes this header; including it first includes
 * the whole of mol_graph.hpp.
 */
#pragma once

#include "smsd/mol_graph.hpp"

namespace smsd {

// ============================================================================
// pH correction of tautomer weights (Henderson-Hasselbalch)
//
// The ratio of neutral tautomers does not depend on pH; the pH sets the fraction of the class
// that is in its neutral form. Each class takes one site model from the ionisable sites of its
// atoms: with pKa_a the lowest acid pKa and pKa_b the highest base pKa (that of the conjugate
// acid BH+) among them, the neutral fraction of the amphoteric species is
// f0 = 1 / (1 + 10^(pH - pKa_a) + 10^(pKa_b - pH)), and a class with sites of one kind only
// drops the other term, which leaves 1 / (1 + 10^(pH - pKa)) for an acid and
// 1 / (1 + 10^(pKa - pH)) for a base.
// ============================================================================

namespace detail {

/**
 * The neutral fraction of a class with an acid site of pKa {@code acidPKa} when
 * {@code hasAcid} and a base site of pKa {@code basePKa} when {@code hasBase}. The fraction is
 * computed in double precision, the acid term added before the base term, and rounded to single
 * precision, steps that the Java library repeats so that both report weights that agree within the
 * shared golden tolerance of 1e-12. The relation takes the dominant site of each kind and leaves
 * out the zwitterion, so it holds for pKa_b below pKa_a.
 */
inline float neutralFraction(bool hasAcid, float acidPKa, bool hasBase, float basePKa, double pH) {
    double denominator = 1.0;
    if (hasAcid) denominator += std::pow(10.0, pH - static_cast<double>(acidPKa));
    if (hasBase) denominator += std::pow(10.0, static_cast<double>(basePKa) - pH);
    return static_cast<float>(1.0 / denominator);
}

/** A tautomer weight scaled by the neutral fraction of its class, clamped to [0.1, 0.99]. */
inline float neutralFormWeight(float weight, float fraction) {
    return std::max(0.1f, std::min(0.99f, weight * fraction));
}

/**
 * The neutral fraction of each atom's class, for the atoms that have an ionisable site; atoms
 * without a site, or outside every class, get NaN and keep their weight at every pH.
 */
inline std::vector<float> tautomerNeutralFractions(const MolGraph& graph, double pH) {
    const int n = graph.n;
    std::vector<float> fractions(static_cast<size_t>(n), std::numeric_limits<float>::quiet_NaN());
    if (static_cast<int>(graph.tautomerClass.size()) < n || static_cast<int>(graph.tautomerSite.size()) < n
        || static_cast<int>(graph.tautomerPKa.size()) < n) return fractions;
    int classes = 0;
    for (int atom = 0; atom < n; ++atom) classes = std::max(classes, graph.tautomerClass[atom] + 1);
    std::vector<char> hasAcid(static_cast<size_t>(classes), 0), hasBase(static_cast<size_t>(classes), 0);
    std::vector<float> acidPKa(static_cast<size_t>(classes), 0.0f), basePKa(static_cast<size_t>(classes), 0.0f);
    for (int atom = 0; atom < n; ++atom) {
        const int cls = graph.tautomerClass[atom];
        if (cls < 0) continue;
        const float pKa = graph.tautomerPKa[atom];
        if (graph.tautomerSite[atom] == MolGraph::TAUTOMER_SITE_ACID) {
            if (!hasAcid[cls] || pKa < acidPKa[cls]) acidPKa[cls] = pKa;
            hasAcid[cls] = 1;
        } else if (graph.tautomerSite[atom] == MolGraph::TAUTOMER_SITE_BASE) {
            if (!hasBase[cls] || pKa > basePKa[cls]) basePKa[cls] = pKa;
            hasBase[cls] = 1;
        }
    }
    for (int atom = 0; atom < n; ++atom) {
        const int cls = graph.tautomerClass[atom];
        if (cls < 0 || graph.tautomerSite[atom] == MolGraph::TAUTOMER_SITE_NONE) continue;
        fractions[atom] = neutralFraction(hasAcid[cls] != 0, acidPKa[cls], hasBase[cls] != 0, basePKa[cls], pH);
    }
    return fractions;
}

} // namespace detail

inline void MolGraph::adjustWeightsForPH(double pH) {
    const std::vector<float> fractions = detail::tautomerNeutralFractions(*this, pH);
    for (int i = 0; i < n; ++i) {
        if (tautomerClass[i] == -1 || tautomerWeight[i] >= 1.0f || std::isnan(fractions[i])) continue;
        tautomerWeight[i] = detail::neutralFormWeight(tautomerWeight[i], fractions[i]);
    }
}

// ============================================================================
// Tautomer classes, motif weights and sites
//
// Every motif pattern is evaluated over all of its matches. The atoms of one match are atoms whose
// hydrogens and bond orders a proton shift of that motif moves together, so they share a class;
// matches that share an atom join their classes, so that atoms which interconvert through any
// chain of tautomeric motifs form one class, in the manner of the mobile-hydrogen groups of InChI
// but without requiring a mobile hydrogen. The class partition is the union-find closure of the
// matches, a function of the graph alone: it depends on neither the order of the atoms nor that of
// their neighbours. tautomerClassNumbers() then numbers the classes by an invariant of the graph.
//
// Each atom takes the weight of the first motif, in the fixed priority order of tautomerMotifs(),
// whose pattern covers it, and the ionisable site and pKa that the matches of that motif covering
// it carry. Every motif claims the atoms of all of its matches before the next motif is tried, so
// this choice does not depend on the atom order either, and the atoms that some motif covers are
// exactly the atoms of the classes. A match of a motif with an acid site carries the site only
// when an atom of the match that would release the proton bears a hydrogen: the N-O oxygen of a
// hydroxamic acid, so that an N-alkoxy amide has none; the nitrogen of an imide, urea, amide,
// thioamide or thioamide/iminothiol; the bridging carbon of a 1,3-diketone; the alpha carbon of a
// ketone or thione; the ring nitrogen or the oxygen of a pyridone/hydroxypyridine or pyridinone;
// and a phenol oxygen with no second heavy-atom neighbour, so that anisoles, aryl esters and
// lactam or lactone carbonyls have none. When every such atom with a hydrogen belongs to an earlier
// motif, the match carries the lowest acid site among them instead, so that each hydrogen keeps
// the pKa of its own motif; an atom covered by several matches of its motif takes the lowest acid
// pKa among them. No motif carries a base site. The order puts the
// most specific motif first: hydroxamic acid, imide and urea before amide; 1,3-diketone,
// enaminone and 1,5-keto-enol before keto/enol and amide; the aromatic lactam motifs
// (pyridone/hydroxypyridine, pyridinone, hydroxypyrimidine) before amide and phenol; quinone
// methide before phenol; tetrazole and triazole before pyrazole; benzimidazole before imidazole;
// guanidine before amidine. The imide motif comes before the barbituric motif, whose ring nitrogen
// between two ring carbonyls is a cyclic imide, so that the nitrogen and the carbonyl carbons keep
// the imide N-H site; the barbituric motif claims the carbonyl oxygens. Of the motifs with the
// same pattern, the pyridone/hydroxypyridine motif, which carries the acid site, comes before the
// hydroxypyrimidine motif, and the furanose/pyranose motif before the lactol motif. The Java
// library uses the same patterns, sites, order and numbering.
// ============================================================================

namespace detail {

/**
 * The aqueous N-H pKa of pyrimidin-4(3H)-one, which a pyridone/hydroxypyridine or pyridinone match
 * takes in place of the 2-pyridone value when the other ring neighbour of its nitrogen is a carbon
 * bonded to a second aromatic ring nitrogen.
 */
inline constexpr float TPKA_PYRIMIDINONE = 8.6f;

/** The weight, ionisable site and pKa that the motif priority gives each atom, and the classes. */
struct TautomerMotifs {
    std::vector<float> weight;   // 1 when no motif covers the atom
    std::vector<int8_t> site;    // MolGraph::TAUTOMER_SITE_*
    std::vector<float> pKa;      // 0 when no match of the atom's motif carries a site
    std::vector<char> covered;   // whether a motif covers the atom, that is whether it has a class
    /**
     * A union-find over the atoms in which the atoms of every motif match share a set, so the sets
     * of the covered atoms are the tautomer classes; read them with unionFindRoot().
     */
    std::vector<int> joined;
};

/** The motif and the class set of every atom of a graph; computeTautomerClasses() numbers the sets. */
inline TautomerMotifs tautomerMotifs(const MolGraph& g) {
    using G = MolGraph;
    const int n = g.n;
    const auto& atomicNum = g.atomicNum;
    const auto& aromatic = g.aromatic;
    const auto& ring = g.ring;
    const auto& neighbors = g.neighbors;
    TautomerMotifs motifs{std::vector<float>(static_cast<size_t>(n), 1.0f),
                          std::vector<int8_t>(static_cast<size_t>(n), G::TAUTOMER_SITE_NONE),
                          std::vector<float>(static_cast<size_t>(n), 0.0f),
                          std::vector<char>(static_cast<size_t>(n), 0),
                          std::vector<int>(static_cast<size_t>(n))};
    std::iota(motifs.joined.begin(), motifs.joined.end(), 0);
    const auto ACID = G::TAUTOMER_SITE_ACID, NONE = G::TAUTOMER_SITE_NONE;
    // The motif being tried and the motif that claimed each atom, the first whose pattern covers it.
    int motifIndex = -1;
    std::vector<int> owner(static_cast<size_t>(n), -1);
    float motifWeight = 1.0f;
    // The site that the match being claimed carries, set by carry() before its atoms are claimed.
    int8_t matchSite = NONE;
    float matchPKa = 0.0f;
    auto claim = [&](int atom) {
        if (!motifs.covered[atom]) {
            motifs.covered[atom] = 1;
            owner[atom] = motifIndex;
            motifs.weight[atom] = motifWeight;
            motifs.site[atom] = matchSite;
            motifs.pKa[atom] = matchPKa;
        } else if (owner[atom] == motifIndex && matchSite == ACID
                   && (motifs.site[atom] != ACID || matchPKa < motifs.pKa[atom])) {
            // Of the matches of its motif that cover an atom, the one with the lowest acid pKa
            // sets its site, whatever the order in which they are found.
            motifs.site[atom] = ACID;
            motifs.pKa[atom] = matchPKa;
        }
    };
    // beginSite(), siteAtom() for each atom of a match that would release the proton, then
    // carry(pKa): the match carries the acid site of pKa `pKa` when such an atom with a hydrogen
    // belongs to this motif (or to none yet, so that the match claims it), else the lowest acid
    // site of the earlier motifs of those atoms, else no site.
    bool ownSite = false, earlierSite = false;
    float earlierPKa = 0.0f;
    auto beginSite = [&]() { ownSite = earlierSite = false; earlierPKa = 0.0f; };
    auto siteAtom = [&](int atom) {
        if (g.attachedHydrogenCount(atom) == 0) return;
        if (!motifs.covered[atom] || owner[atom] == motifIndex) {
            ownSite = true;
        } else if (motifs.site[atom] == ACID && (!earlierSite || motifs.pKa[atom] < earlierPKa)) {
            earlierSite = true;
            earlierPKa = motifs.pKa[atom];
        }
    };
    auto carry = [&](float pKa) {
        matchSite = ownSite || earlierSite ? ACID : NONE;
        matchPKa = ownSite ? pKa : earlierSite ? earlierPKa : 0.0f;
    };
    auto carries = [&](float pKa, std::initializer_list<int> atoms) {
        beginSite();
        for (int atom : atoms) siteAtom(atom);
        carry(pKa);
    };
    // join(atom, with) adds an atom to the match of `with`; group() claims the atoms of a match.
    auto join = [&](int atom, int with) { claim(atom); unionFindJoin(motifs.joined, atom, with); };
    auto group = [&](std::initializer_list<int> atoms) {
        const int first = *atoms.begin();
        claim(first);
        for (int atom : atoms) join(atom, first);
    };
    // A motif without a site claims its matches with none; a motif with one calls carry() or
    // carries() before it claims each match.
    auto motif = [&](float weight, const auto& match) {
        ++motifIndex;
        motifWeight = weight;
        matchSite = NONE;
        matchPKa = 0.0f;
        for (int i = 0; i < n; ++i) match(i);
    };
    auto bond = [&](int a, int b) { return g.bondOrder(a, b); };
    auto hasNitrogenBesides = [&](int carbon, int except) {
        for (int k : neighbors[carbon]) if (k != except && atomicNum[k] == 7) return true;
        return false;
    };
    auto isHydroxyl = [&](int nitrogen, int o) {   // the N-O of a hydroxamic acid
        return atomicNum[o] == 8 && !aromatic[o] && bond(nitrogen, o) == 1;
    };
    auto isCarbonyl = [&](int carbon, int except) {   // a C=O carbon, whose oxygen is not `except`
        for (int k : neighbors[carbon]) if (k != except && atomicNum[k] == 8 && bond(carbon, k) == 2) return true;
        return false;
    };
    // The para carbons c of an aromatic ring carbon j (path j-a-b-c), as T4 and T24 walk them.
    auto forParaCarbons = [&](int i, int j, const auto& visit) {
        for (int a : neighbors[j]) {
            if (a == i || atomicNum[a] != 6 || !aromatic[a] || !ring[a]) continue;
            for (int b : neighbors[a]) {
                if (b == j || atomicNum[b] != 6 || !aromatic[b] || !ring[b]) continue;
                for (int c : neighbors[b])
                    if (c != a && c != j && atomicNum[c] == 6 && aromatic[c] && ring[c]) visit(c);
            }
        }
    };
    auto aromaticRingN = [&](int atom) { return atomicNum[atom] == 7 && aromatic[atom] && ring[atom]; };
    // The acid pKa of a pyridone/hydroxypyridine or pyridinone match of ring nitrogen `nitrogen`
    // and oxygen-bearing carbon `carbon`: that of pyrimidin-4(3H)-one when the other ring neighbour
    // of the nitrogen is a carbon bonded to a second aromatic ring nitrogen, else `pKa`.
    auto lactamPKa = [&](int nitrogen, int carbon, float pKa) {
        for (int m : neighbors[nitrogen]) {
            if (m == carbon || atomicNum[m] != 6 || !aromatic[m] || !ring[m]) continue;
            for (int p : neighbors[m]) if (p != nitrogen && aromaticRingN(p)) return TPKA_PYRIMIDINONE;
        }
        return pKa;
    };
    auto heavyNeighbours = [&](int atom) {
        int count = 0;
        for (int k : neighbors[atom]) count += atomicNum[k] > 1;
        return count;
    };

    // T13 hydroxamic acid N(-OH)-C=O: the N, its N-O oxygens and every C=O on the N; the site is
    // the O-H of an N-O oxygen
    motif(G::TW_HYDROXAMIC, [&](int i) {
        if (atomicNum[i] != 7) return;
        bool hydroxyl = false, carbonyl = false;
        for (int o : neighbors[i]) hydroxyl = hydroxyl || isHydroxyl(i, o);
        for (int j : neighbors[i]) carbonyl = carbonyl || (atomicNum[j] == 6 && isCarbonyl(j, i));
        if (!hydroxyl || !carbonyl) return;
        beginSite();
        for (int o : neighbors[i]) if (isHydroxyl(i, o)) siteAtom(o);
        carry(G::TPKA_HYDROXAMIC);
        claim(i);
        for (int o : neighbors[i]) if (isHydroxyl(i, o)) join(o, i);
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6) continue;
            for (int k : neighbors[j])
                if (k != i && atomicNum[k] == 8 && bond(j, k) == 2) { join(j, i); join(k, i); }
        }
    });
    // T12 imide C(=O)-N-C(=O): the N and its carbonyl carbons; the site is the N-H
    motif(G::TW_IMIDE, [&](int i) {
        if (atomicNum[i] != 7) return;
        int carbonyls = 0;
        for (int j : neighbors[i]) if (atomicNum[j] == 6 && isCarbonyl(j, i)) ++carbonyls;
        if (carbonyls < 2) return;
        carries(G::TPKA_IMIDE, {i});
        claim(i);
        for (int j : neighbors[i]) if (atomicNum[j] == 6 && isCarbonyl(j, i)) join(j, i);
    });
    // T28 barbituric: ring N between ring carbonyls, with their oxygens
    motif(G::TW_BARBITURIC, [&](int i) {
        if (atomicNum[i] != 7 || !ring[i]) return;
        int carbonyls = 0;
        for (int j : neighbors[i]) if (atomicNum[j] == 6 && ring[j] && isCarbonyl(j, i)) ++carbonyls;
        if (carbonyls < 2) return;
        claim(i);
        for (int c : neighbors[i]) {
            if (atomicNum[c] != 6 || !ring[c] || !isCarbonyl(c, i)) continue;
            join(c, i);
            for (int k : neighbors[c]) if (atomicNum[k] == 8 && bond(c, k) == 2) join(k, i);
        }
    });
    // T10 urea O=C(-N)(-N): the O, its carbon and that carbon's nitrogens; the site is an N-H
    motif(G::TW_UREA, [&](int i) {
        if (atomicNum[i] != 8 || aromatic[i]) return;
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6 || bond(i, j) != 2) continue;
            int nitrogens = 0;
            for (int k : neighbors[j]) if (k != i && atomicNum[k] == 7) ++nitrogens;
            if (nitrogens < 2) continue;
            beginSite();
            for (int k : neighbors[j]) if (k != i && atomicNum[k] == 7) siteAtom(k);
            carry(G::TPKA_UREA);
            group({i, j});
            for (int k : neighbors[j]) if (k != i && atomicNum[k] == 7) join(k, i);
        }
    });
    // T5 1,3-diketone O=C-C-C=O: both oxygens and the bridging carbon, whose C-H is the site
    motif(G::TW_DIKETONE_ENOL, [&](int i) {
        if (atomicNum[i] != 8 || aromatic[i]) return;
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6 || bond(i, j) != 2) continue;
            for (int bridge : neighbors[j]) {
                if (bridge == i || atomicNum[bridge] != 6 || aromatic[bridge]) continue;
                for (int m : neighbors[bridge]) {
                    if (m == j || m == i || atomicNum[m] != 6) continue;
                    if (std::find(neighbors[j].begin(), neighbors[j].end(), m) != neighbors[j].end()) continue;
                    for (int o2 : neighbors[m])
                        if (o2 != bridge && atomicNum[o2] == 8 && !aromatic[o2] && bond(m, o2) == 2) {
                            carries(G::TPKA_DIKETONE_ENOL, {bridge});
                            group({i, bridge, o2});
                        }
                }
            }
        }
    });
    // T11 beta-enaminone N-C=C-C=O
    motif(G::TW_ENAMINONE, [&](int i) {
        if (atomicNum[i] != 7 || aromatic[i]) return;
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6 || bond(i, j) != 1) continue;
            for (int k : neighbors[j]) {
                if (k == i || atomicNum[k] != 6 || bond(j, k) != 2) continue;
                for (int m : neighbors[k]) {
                    if (m == j || atomicNum[m] != 6) continue;
                    for (int o : neighbors[m])
                        if (o != k && atomicNum[o] == 8 && bond(m, o) == 2) group({i, j, k, m, o});
                }
            }
        }
    });
    // T21 1,5-keto-enol O=C-C=C-C
    motif(G::TW_15_KETO_ENOL, [&](int i) {
        if (atomicNum[i] != 8 || aromatic[i]) return;
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6 || bond(i, j) != 2) continue;
            for (int k : neighbors[j]) {
                if (k == i || atomicNum[k] != 6 || bond(j, k) != 1) continue;
                for (int m : neighbors[k]) {
                    if (m == j || atomicNum[m] != 6 || bond(k, m) != 2) continue;
                    for (int p : neighbors[m])
                        if (p != k && atomicNum[p] == 6 && !aromatic[p]) group({i, j, k, m, p});
                }
            }
        }
    });
    // T27 pyridone/hydroxypyridine: aromatic ring N, aromatic ring C, single-bonded O; the site is
    // the N-H or O-H
    motif(G::TW_PYRIDONE_HYDROXYPYRIDINE, [&](int i) {
        if (!aromaticRingN(i)) return;
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6 || !aromatic[j] || !ring[j]) continue;
            for (int k : neighbors[j])
                if (k != i && atomicNum[k] == 8 && !aromatic[k] && bond(j, k) == 1) {
                    carries(lactamPKa(i, j, G::TPKA_PYRIDONE_HYDROXYPYRIDINE), {i, k});
                    group({i, j, k});
                }
        }
    });
    // T6 pyridinone: aromatic ring N next to an aromatic C bearing an exocyclic O; the site is the
    // N-H or O-H
    motif(G::TW_PYRIDINONE, [&](int i) {
        if (!aromaticRingN(i)) return;
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6 || !aromatic[j]) continue;
            for (int k : neighbors[j])
                if (k != i && atomicNum[k] == 8 && !aromatic[k]) {
                    carries(lactamPKa(i, j, G::TPKA_PYRIDINONE), {i, k});
                    group({i, k});
                }
        }
    });
    // T14 hydroxypyrimidine: exocyclic O on an aromatic ring C next to an aromatic ring N
    motif(G::TW_HYDROXYPYRIMIDINE, [&](int i) {
        if (atomicNum[i] != 8 || aromatic[i]) return;
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6 || !aromatic[j] || !ring[j] || bond(i, j) != 1) continue;
            for (int k : neighbors[j])
                if (k != i && aromaticRingN(k)) group({i, j, k});
        }
    });
    // T1 keto/enol C=O with its alpha carbons, and its nitrogens for an amide; the site is the N-H
    // of an amide and the alpha C-H of a ketone
    auto carbonyl = [&](bool amide, float pKa) {
        return [&, amide, pKa](int i) {
            if (atomicNum[i] != 8 || aromatic[i]) return;
            for (int j : neighbors[i]) {
                if (atomicNum[j] != 6 || bond(i, j) != 2 || hasNitrogenBesides(j, i) != amide) continue;
                for (int k : neighbors[j])
                    if (k != i && ((atomicNum[k] == 6 && !aromatic[k]) || atomicNum[k] == 7)) {
                        beginSite();
                        if (!amide || atomicNum[k] == 7) siteAtom(k);
                        carry(pKa);
                        group({i, k});
                    }
            }
        };
    };
    motif(G::TW_AMIDE_IMIDIC, carbonyl(true, G::TPKA_AMIDE_IMIDIC));
    motif(G::TW_KETO_ENOL, carbonyl(false, G::TPKA_KETO_ENOL));
    // T2 thione/thiol C=S with the carbon and nitrogen neighbours of its carbon; the site is the N-H
    // of a thioamide and the non-aromatic alpha C-H of a thione
    auto thiocarbonyl = [&](bool amide, float pKa) {
        return [&, amide, pKa](int i) {
            if (atomicNum[i] != 16 || aromatic[i]) return;
            for (int j : neighbors[i]) {
                if (atomicNum[j] != 6 || bond(i, j) != 2 || hasNitrogenBesides(j, i) != amide) continue;
                for (int k : neighbors[j])
                    if (k != i && (atomicNum[k] == 6 || atomicNum[k] == 7)) {
                        beginSite();
                        if (amide ? atomicNum[k] == 7 : !aromatic[k]) siteAtom(k);
                        carry(pKa);
                        group({i, k});
                    }
            }
        };
    };
    motif(G::TW_THIOAMIDE, thiocarbonyl(true, G::TPKA_THIOAMIDE));
    motif(G::TW_THIONE_THIOL, thiocarbonyl(false, G::TPKA_THIONE_THIOL));
    // T16 thioamide/iminothiol S=C-N; the site is the N-H
    motif(G::TW_THIOAMIDE_IMINOTHIOL, [&](int i) {
        if (atomicNum[i] != 16 || aromatic[i]) return;
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6 || bond(i, j) != 2) continue;
            for (int k : neighbors[j])
                if (k != i && atomicNum[k] == 7) {
                    carries(G::TPKA_THIOAMIDE_IMINOTHIOL, {k});
                    group({i, j, k});
                }
        }
    });
    // T3 nitroso/oxime N=C-O
    motif(G::TW_NITROSO_OXIME, [&](int i) {
        if (atomicNum[i] != 7 || ring[i]) return;
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6 || bond(i, j) != 2) continue;
            for (int k : neighbors[j])
                if (k != i && atomicNum[k] == 8) group({i, k});
        }
    });
    // T24 phenol/quinone methide: para carbon with an exocyclic C=C; no site
    motif(G::TW_PHENOL_QUINONE_METHIDE, [&](int i) {
        if (atomicNum[i] != 8 || aromatic[i]) return;
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6 || !aromatic[j] || !ring[j] || bond(i, j) != 1) continue;
            forParaCarbons(i, j, [&](int c) {
                for (int d : neighbors[c])
                    if (atomicNum[d] == 6 && !ring[d] && bond(c, d) == 2) group({i, j, c, d});
            });
        }
    });
    // T4 phenol/quinone: O on an aromatic ring C, its para carbon and that carbon's oxygens; the site
    // is the O-H of an oxygen with no second heavy-atom neighbour
    motif(G::TW_PHENOL_QUINONE, [&](int i) {
        if (atomicNum[i] != 8 || aromatic[i]) return;
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6 || !aromatic[j] || !ring[j]) continue;
            forParaCarbons(i, j, [&](int c) {
                beginSite();
                if (heavyNeighbours(i) == 1) siteAtom(i);
                carry(G::TPKA_PHENOL_QUINONE);
                group({i, j, c});
                for (int d : neighbors[c]) if (atomicNum[d] == 8 && !aromatic[d]) join(d, i);
            });
        }
    });
    // T15 tetrazole/triazole: aromatic ring N with two or more aromatic ring N neighbours
    motif(G::TW_TETRAZOLE_NH, [&](int i) {
        if (!aromaticRingN(i)) return;
        int count = 0;
        for (int j : neighbors[i]) if (aromaticRingN(j)) ++count;
        if (count < 2) return;
        claim(i);
        for (int j : neighbors[i]) if (aromaticRingN(j)) join(j, i);
    });
    // T25 1,2,3-triazole NH shift N-N-N
    motif(G::TW_TRIAZOLE_NH_SHIFT, [&](int i) {
        if (!aromaticRingN(i)) return;
        int only = -1, count = 0;
        for (int j : neighbors[i]) if (aromaticRingN(j)) { only = j; ++count; }
        if (count != 1) return;
        for (int k : neighbors[only])
            if (k != i && aromaticRingN(k)) group({i, only, k});
    });
    // T26 benzimidazole N-C-N through a fusion carbon
    motif(G::TW_BENZIMIDAZOLE_NH, [&](int i) {
        if (!aromaticRingN(i)) return;
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6 || !aromatic[j] || !ring[j]) continue;
            int ringNeighbours = 0;
            for (int k : neighbors[j]) if (ring[k] && aromatic[k]) ++ringNeighbours;
            if (ringNeighbours < 3) continue;
            for (int k : neighbors[j])
                if (k != i && aromaticRingN(k)) group({i, j, k});
        }
    });
    // T7 pyrazole N-N, then imidazole N-C-N in a five-membered ring
    motif(G::TW_IMIDAZOLE_NH, [&](int i) {
        if (!aromaticRingN(i)) return;
        for (int j : neighbors[i]) if (aromaticRingN(j)) group({i, j});
    });
    motif(G::TW_IMIDAZOLE_NH, [&](int i) {
        if (!aromaticRingN(i)) return;
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6 || !aromatic[j] || !ring[j]) continue;
            for (int k : neighbors[j]) {
                if (k == i || !aromaticRingN(k)) continue;
                bool in5 = false;
                for (const auto& r : g.computeRings()) {
                    if (r.size() != 5) continue;
                    if (std::find(r.begin(), r.end(), i) != r.end() && std::find(r.begin(), r.end(), k) != r.end()) {
                        in5 = true;
                        break;
                    }
                }
                if (in5) group({i, k});
            }
        }
    });
    // T9 guanidine: C with three or more non-aromatic N
    motif(G::TW_GUANIDINE, [&](int i) {
        if (atomicNum[i] != 6) return;
        int count = 0;
        for (int j : neighbors[i]) if (atomicNum[j] == 7 && !aromatic[j]) ++count;
        if (count < 3) return;
        claim(i);
        for (int j : neighbors[i]) if (atomicNum[j] == 7 && !aromatic[j]) join(j, i);
    });
    // T8 amidine N=C-N
    motif(G::TW_AMIDINE, [&](int i) {
        if (atomicNum[i] != 7 || aromatic[i]) return;
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6 || bond(i, j) != 2) continue;
            for (int k : neighbors[j])
                if (k != i && atomicNum[k] == 7 && !aromatic[k]) group({i, j, k});
        }
    });
    // T17 phosphonate P with two or more O
    motif(G::TW_PHOSPHONATE, [&](int i) {
        if (atomicNum[i] != 15) return;
        int count = 0;
        for (int j : neighbors[i]) if (atomicNum[j] == 8) ++count;
        if (count < 2) return;
        claim(i);
        for (int j : neighbors[i]) if (atomicNum[j] == 8) join(j, i);
    });
    // T19 nitro/aci-nitro: non-aromatic N with two or more O and a C
    motif(G::TW_NITRO_ACI, [&](int i) {
        if (atomicNum[i] != 7 || aromatic[i]) return;
        int oxygens = 0, carbons = 0;
        for (int j : neighbors[i]) {
            if (atomicNum[j] == 8) ++oxygens;
            else if (atomicNum[j] == 6) ++carbons;
        }
        if (oxygens < 2 || carbons == 0) return;
        claim(i);
        for (int j : neighbors[i]) if (atomicNum[j] == 8 || atomicNum[j] == 6) join(j, i);
    });
    // T22 furanose/pyranose: ring O, ring C, exocyclic single-bonded O
    motif(G::TW_FURANOSE_PYRANOSE, [&](int i) {
        if (atomicNum[i] != 8 || !ring[i] || aromatic[i]) return;
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6 || !ring[j]) continue;
            for (int k : neighbors[j])
                if (k != i && atomicNum[k] == 8 && !ring[k] && bond(j, k) == 1) group({i, j, k});
        }
    });
    // T23 lactol: ring C bonded to a ring O and an exocyclic single-bonded O
    motif(G::TW_LACTOL, [&](int i) {
        if (atomicNum[i] != 6 || !ring[i]) return;
        auto ringOxygen = [&](int j) { return atomicNum[j] == 8 && ring[j] && !aromatic[j]; };
        auto exocyclicOxygen = [&](int j) {
            return atomicNum[j] == 8 && !ring[j] && !aromatic[j] && bond(i, j) == 1;
        };
        bool hasRingO = false, hasExoO = false;
        for (int j : neighbors[i]) { hasRingO = hasRingO || ringOxygen(j); hasExoO = hasExoO || exocyclicOxygen(j); }
        if (!hasRingO || !hasExoO) return;
        claim(i);
        for (int j : neighbors[i]) if (ringOxygen(j) || exocyclicOxygen(j)) join(j, i);
    });
    // T29 allyl shift X(-H)-C=C-C, X = O, S or N
    motif(G::TW_ALLYL_SHIFT, [&](int i) {
        if ((atomicNum[i] != 8 && atomicNum[i] != 16 && atomicNum[i] != 7) || aromatic[i] || ring[i]) return;
        int sumBO = 0;
        for (int nb : neighbors[i]) sumBO += bond(i, nb);
        const int maxValence = (atomicNum[i] == 7) ? 3 : 2;
        const int64_t charge = g.formalCharge[i];
        if (!(sumBO < int64_t(maxValence) + (charge < 0 ? -charge : charge))) return;
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6 || aromatic[j] || bond(i, j) != 1) continue;
            for (int k : neighbors[j]) {
                if (k == i || atomicNum[k] != 6 || bond(j, k) != 2) continue;
                for (int m : neighbors[k])
                    if (m != j && atomicNum[m] == 6 && !aromatic[m] && bond(k, m) == 1) group({i, j, k, m});
            }
        }
    });
    // T30 selenol/selenoketone C=Se with the carbon and nitrogen neighbours of its carbon
    motif(G::TW_SELENOL_SELENOKETONE, [&](int i) {
        if (atomicNum[i] != 34 || aromatic[i]) return;
        for (int j : neighbors[i]) {
            if (atomicNum[j] != 6 || bond(i, j) != 2) continue;
            group({i, j});
            for (int k : neighbors[j]) if (k != i && (atomicNum[k] == 6 || atomicNum[k] == 7)) join(k, i);
        }
    });
    return motifs;
}

// ----------------------------------------------------------------------------
// Class numbering
//
// The classes are numbered in increasing order of a key, the sorted colours of their atoms. The
// colour of an atom starts from its element, charge, mass number, hydrogen count, aromaticity,
// ring membership, degree and whether a motif covers it, and each round of refinement hashes into
// it the bond orders and colours of its neighbours, as in Weisfeiler-Leman refinement; rounds stop
// once the keys differ or a round splits no colour. When classes still share the smallest key
// the refinement cannot tell them apart (in practice a symmetry of the graph exchanges them): the
// one with the lowest atom index takes the next number, its atoms are given colours of their own
// and the refinement resumes, so that the classes tied to it are told apart consistently. The
// numbering is therefore the same for every order of the atoms, except among classes that the
// refinement cannot tell apart.
// Colours are 64-bit hashes cut to 62 bits, so that the signed comparison of the Java library
// gives the same order; the Java library repeats every step.
// ----------------------------------------------------------------------------

/** The 62-bit hash of `value` combined into `hash`, shared bit for bit with the Java library. */
inline uint64_t tautomerColourMix(uint64_t hash, int64_t value) {
    uint64_t x = (hash ^ (static_cast<uint64_t>(value) + 0x9E3779B97F4A7C15ULL)) * 0xFF51AFD7ED558CCDULL;
    x ^= x >> 33;
    x *= 0xC4CEB9FE1A85EC53ULL;
    x ^= x >> 33;
    return x & 0x3FFFFFFFFFFFFFFFULL;
}

/** The number of distinct values in `colours`. */
inline size_t tautomerDistinctColours(std::vector<uint64_t> colours) {
    std::sort(colours.begin(), colours.end());
    return static_cast<size_t>(std::unique(colours.begin(), colours.end()) - colours.begin());
}

/**
 * The tautomer class of every atom: -1 for an atom that no motif covers, otherwise the number of
 * its set in `motifs.joined`, numbered from 0 by the rule above.
 */
inline std::vector<int> tautomerClassNumbers(const MolGraph& g, TautomerMotifs& motifs) {
    const int n = g.n;
    std::vector<int> numbers(static_cast<size_t>(n), -1);
    // The atoms of each class, in increasing atom order.
    std::vector<int> classOfRoot(static_cast<size_t>(n), -1);
    std::vector<std::vector<int>> members;
    for (int atom = 0; atom < n; ++atom) {
        if (!motifs.covered[atom]) continue;
        const int root = unionFindRoot(motifs.joined, atom);
        if (classOfRoot[root] < 0) {
            classOfRoot[root] = static_cast<int>(members.size());
            members.emplace_back();
        }
        members[static_cast<size_t>(classOfRoot[root])].push_back(atom);
    }
    if (members.size() == 1)
        for (int atom : members[0]) numbers[atom] = 0;
    if (members.size() <= 1) return numbers;

    auto property = [](const auto& values, int atom) -> int64_t {
        return atom < static_cast<int>(values.size()) ? static_cast<int64_t>(values[atom]) : 0;
    };
    std::vector<uint64_t> colours(static_cast<size_t>(n));
    for (int atom = 0; atom < n; ++atom) {
        uint64_t colour = tautomerColourMix(0x5441555443ULL, g.atomicNum[atom]);
        colour = tautomerColourMix(colour, property(g.formalCharge, atom));
        colour = tautomerColourMix(colour, property(g.massNumber, atom));
        colour = tautomerColourMix(colour, g.implicitHydrogenCount(atom));
        colour = tautomerColourMix(colour, g.aromatic[atom] ? 1 : 0);
        colour = tautomerColourMix(colour, g.ring[atom] ? 1 : 0);
        colour = tautomerColourMix(colour, static_cast<int64_t>(g.neighbors[atom].size()));
        colours[atom] = tautomerColourMix(colour, motifs.covered[atom] ? 1 : 0);
    }
    size_t distinct = tautomerDistinctColours(colours);
    std::vector<uint64_t> refined(static_cast<size_t>(n)), around;
    // One round of refinement; false when it splits no colour, which leaves the colours as they were.
    auto refine = [&]() {
        for (int atom = 0; atom < n; ++atom) {
            around.clear();
            for (int neighbor : g.neighbors[atom])
                around.push_back(tautomerColourMix(tautomerColourMix(0x424F4E44ULL, g.bondOrder(atom, neighbor)),
                                                   static_cast<int64_t>(colours[neighbor])));
            std::sort(around.begin(), around.end());
            uint64_t colour = tautomerColourMix(0x524F554E44ULL, static_cast<int64_t>(colours[atom]));
            for (uint64_t value : around) colour = tautomerColourMix(colour, static_cast<int64_t>(value));
            refined[atom] = colour;
        }
        const size_t split = tautomerDistinctColours(refined);
        if (split <= distinct) return false;
        colours.swap(refined);
        distinct = split;
        return true;
    };

    std::vector<int> remaining(members.size());
    std::iota(remaining.begin(), remaining.end(), 0);
    std::vector<std::vector<uint64_t>> keys(members.size());
    auto sortRemaining = [&]() {
        bool tied = false;
        for (int cls : remaining) {
            auto& key = keys[static_cast<size_t>(cls)];
            key.clear();
            for (int atom : members[static_cast<size_t>(cls)]) key.push_back(colours[atom]);
            std::sort(key.begin(), key.end());
        }
        std::sort(remaining.begin(), remaining.end(), [&](int first, int second) {
            if (keys[first] != keys[second]) return keys[first] < keys[second];
            return members[first][0] < members[second][0];
        });
        for (size_t index = 1; index < remaining.size(); ++index)
            tied = tied || keys[remaining[index - 1]] == keys[remaining[index]];
        return tied;
    };
    int next = 0;
    while (!remaining.empty()) {
        while (sortRemaining() && refine()) {}
        size_t position = 0;
        bool individualised = false;
        for (; position < remaining.size() && !individualised; ++position) {
            const int cls = remaining[position];
            const int number = next++;
            for (int atom : members[static_cast<size_t>(cls)]) numbers[atom] = number;
            if (position + 1 < remaining.size() && keys[cls] == keys[remaining[position + 1]]) {
                for (int atom : members[static_cast<size_t>(cls)])
                    colours[atom] = tautomerColourMix(tautomerColourMix(0x494E444956ULL, number),
                                                      static_cast<int64_t>(colours[atom]));
                distinct = tautomerDistinctColours(colours);
                individualised = true;
            }
        }
        remaining.erase(remaining.begin(), remaining.begin() + static_cast<std::ptrdiff_t>(position));
    }
    return numbers;
}

} // namespace detail

// ============================================================================
// Tautomer class detection: pKa-informed transforms (Tier 1)
// ============================================================================

inline void MolGraph::computeTautomerClasses() {
    detail::TautomerMotifs motifs = detail::tautomerMotifs(*this);
    tautomerClass = detail::tautomerClassNumbers(*this, motifs);
    // Every atom a motif covers is in a class, so the weight, site and pKa of each atom are those
    // of its motif; atoms outside every class keep weight 1 and no site.
    tautomerWeight = std::move(motifs.weight);
    tautomerSite = std::move(motifs.site);
    tautomerPKa = std::move(motifs.pKa);
}

// ============================================================================
// Tier 2: Solvent-dependent tautomer weight corrections
//
// Aprotic/non-aqueous solvents shift keto-enol and related equilibria.
// Literature pKa shifts (Bordwell 1988; Olmstead 1977):
//   DMSO:       keto-enol +5 pKa units  -> enol more prevalent (shift -0.15)
//               amide/lactam +3          -> imidic slightly more prevalent
//   Methanol:   keto-enol +2             -> small shift (-0.05)
//   Chloroform: keto-enol +8 (aprotic)   -> much more enol (-0.25)
// ============================================================================

inline void MolGraph::applySolventCorrection(ChemOptions::Solvent solvent) {
    if (solvent == ChemOptions::Solvent::AQUEOUS) return;

    float ketoShift  = 0.0f;
    float amideShift = 0.0f;
    float generalShift = 0.0f;

    if (solvent == ChemOptions::Solvent::DMSO) {
        ketoShift    = -0.15f;
        amideShift   = -0.08f;
        generalShift = -0.05f;
    } else if (solvent == ChemOptions::Solvent::CHLOROFORM) {
        ketoShift    = -0.25f;
        amideShift   = -0.12f;
        generalShift = -0.08f;
    } else if (solvent == ChemOptions::Solvent::METHANOL) {
        ketoShift    = -0.05f;
        amideShift   = -0.02f;
        generalShift = -0.01f;
    } else if (solvent == ChemOptions::Solvent::ACETONITRILE) {
        ketoShift    = -0.10f;
        amideShift   = -0.05f;
        generalShift = -0.03f;
    } else if (solvent == ChemOptions::Solvent::DIETHYL_ETHER) {
        ketoShift    = -0.20f;
        amideShift   = -0.10f;
        generalShift = -0.06f;
    }

    for (int i = 0; i < n; i++) {
        if (tautomerClass[i] == -1) continue;
        float w = tautomerWeight[i];
        if (w >= 1.0f) continue;

        float shift;
        if (w >= 0.96f) {
            shift = (w >= 0.97f) ? amideShift : ketoShift;
            if (std::abs(w - TW_KETO_ENOL) < 0.005f) shift = ketoShift;
        } else if (w >= 0.93f) {
            shift = amideShift;
        } else if (w >= 0.84f) {
            shift = generalShift;
        } else {
            shift = generalShift;
        }

        tautomerWeight[i] = std::max(0.05f, std::min(0.99f, w + shift));
    }
}

// ===========================================================================
// Tautomer confidence -- heuristic quality of selected tautomeric relaxations.
// ===========================================================================

namespace detail {

/**
 * The rule shared by both tautomer-confidence forms.  A query atom is relaxed when its element
 * or aromaticity differs from its target atom, or when it ends a selected bond whose order or
 * aromaticity differs.  The confidence is the geometric mean, over relaxed atoms with a
 * tautomer class in either molecule, of the smaller baseline weight after the solvent and pH
 * context; it is 1 when nothing is relaxed or when the chemistry is not tautomer-aware.
 * Out-of-range atoms throw std::out_of_range and bonds outside the mapping throw
 * std::invalid_argument.
 */
inline double selectedTautomerConfidence(const MolGraph& query, const MolGraph& target,
                                         const std::map<int, int>& atoms,
                                         const std::vector<std::array<int, 4>>& bonds,
                                         const ChemOptions& chemistry) {
    if (!std::isfinite(chemistry.pH) || chemistry.pH < 0 || chemistry.pH > 14)
        throw std::invalid_argument("pH must be finite and between 0 and 14");
    switch (chemistry.solvent) {
        case ChemOptions::Solvent::AQUEOUS: case ChemOptions::Solvent::DMSO:
        case ChemOptions::Solvent::METHANOL: case ChemOptions::Solvent::CHLOROFORM:
        case ChemOptions::Solvent::ACETONITRILE: case ChemOptions::Solvent::DIETHYL_ETHER: break;
        default: throw std::invalid_argument("Unknown solvent");
    }
    if (!chemistry.tautomerAware) return 1.0;
    std::vector<bool> relaxed(query.n, false);
    for (const auto& atom : atoms) {
        if (atom.first < 0 || atom.first >= query.n || atom.second < 0 || atom.second >= target.n)
            throw std::out_of_range("Tautomer confidence atom index out of range");
        relaxed[atom.first] = query.atomicNum[atom.first] != target.atomicNum[atom.second]
            || query.aromatic[atom.first] != target.aromatic[atom.second];
    }
    for (const auto& bond : bonds) {
        if (bond[0] < 0 || bond[1] < 0 || bond[0] >= query.n || bond[1] >= query.n
            || bond[2] < 0 || bond[3] < 0 || bond[2] >= target.n || bond[3] >= target.n)
            throw std::out_of_range("Tautomer confidence bond atom index out of range");
        const auto first = atoms.find(bond[0]), second = atoms.find(bond[1]);
        if (first == atoms.end() || second == atoms.end()
            || first->second != bond[2] || second->second != bond[3]
            || !query.hasBond(bond[0], bond[1]) || !target.hasBond(bond[2], bond[3]))
            throw std::invalid_argument("Tautomer confidence bond is outside the selected mapping");
        if (query.bondOrder(bond[0], bond[1]) != target.bondOrder(bond[2], bond[3])
            || query.bondAromatic(bond[0], bond[1]) != target.bondAromatic(bond[2], bond[3]))
            relaxed[bond[0]] = relaxed[bond[1]] = true;
    }
    if (std::none_of(relaxed.begin(), relaxed.end(), [](bool value) { return value; })) return 1.0;
    auto baseline = [](const MolGraph& source) {
        MolGraph copy;
        {
            const std::lock_guard<std::recursive_mutex> cacheLock(source.lazyCacheMutex());
            copy = source;
        }
        copy.computeTautomerClasses();
        return copy;
    };
    // Recompute local baseline weights because Builder graphs may carry a solvent context.
    MolGraph queryBaseline = baseline(query), targetBaseline = baseline(target);
    queryBaseline.applySolventCorrection(chemistry.solvent);
    targetBaseline.applySolventCorrection(chemistry.solvent);
    queryBaseline.adjustWeightsForPH(chemistry.pH);
    targetBaseline.adjustWeightsForPH(chemistry.pH);
    double logSum = 0.0;
    size_t count = 0;
    for (const auto& atom : atoms) {
        if (relaxed[atom.first] && (queryBaseline.tautomerClass[atom.first] >= 0
                                   || targetBaseline.tautomerClass[atom.second] >= 0)) {
            const double weight = std::min(queryBaseline.tautomerWeight[atom.first],
                                           targetBaseline.tautomerWeight[atom.second]);
            logSum += std::log(std::max(weight, 1e-10));
            ++count;
        }
    }
    return count == 0 ? 1.0 : std::exp(logSum / static_cast<double>(count));
}

/** The mapped query bonds whose target atoms are bonded, as {query, query, target, target}. */
inline std::vector<std::array<int, 4>> inducedBondMapping(const MolGraph& query, const MolGraph& target,
                                                          const std::map<int, int>& atoms) {
    std::vector<std::array<int, 4>> bonds;
    for (const auto& atom : atoms)
        for (int neighbor : query.neighbors[atom.first]) {
            if (neighbor <= atom.first) continue;
            const auto other = atoms.find(neighbor);
            if (other != atoms.end() && target.hasBond(atom.second, other->second))
                bonds.push_back({atom.first, neighbor, atom.second, other->second});
        }
    return bonds;
}

} // namespace detail

/**
 * Tautomer confidence of an atom mapping alone.
 *
 * @deprecated Use MCS::computeTautomerConfidence(query, target, match, chemistry), which takes
 * the selected bonds and the solvent and pH context.  This form selects every mapped bond
 * present in both molecules and applies the tautomer-aware rule at pH 7.4 in water, so it
 * returns the same value as the four-argument form for that selection.  Out-of-range atoms
 * throw std::out_of_range.
 */
[[deprecated("Use MCS::computeTautomerConfidence(query, target, match, chemistry).")]]
inline double computeTautomerConfidence(
        const MolGraph& gq,
        const MolGraph& gt,
        const std::unordered_map<int,int>& mcs) {
    const std::map<int, int> atoms(mcs.begin(), mcs.end());
    for (const auto& atom : atoms)
        if (atom.first < 0 || atom.first >= gq.n || atom.second < 0 || atom.second >= gt.n)
            throw std::out_of_range("Tautomer confidence atom index out of range");
    ChemOptions chemistry;
    chemistry.tautomerAware = true;
    return detail::selectedTautomerConfidence(gq, gt, atoms, detail::inducedBondMapping(gq, gt, atoms),
                                              chemistry);
}

} // namespace smsd
