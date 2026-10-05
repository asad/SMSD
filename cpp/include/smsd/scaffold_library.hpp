/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. */
#pragma once
#ifndef SMSD_SCAFFOLD_LIBRARY_HPP
#define SMSD_SCAFFOLD_LIBRARY_HPP

#include <array>
#include <cstddef>
#include <cstring>

namespace smsd {

// ============================================================================
// Pharmaceutical scaffold library — constexpr reference data for drug design.
//
// Contents:
//   1.  Top 50 ring systems by frequency in FDA-approved drugs
//   2.  Top 30 Murcko scaffolds from approved drugs
//   3.  Pharmacophore SMARTS patterns (8 standard + 10 extended families)
//   4.  Drug-like property boundaries (Lipinski, Veber, average drug stats)
//   5.  Natural product scaffold classes (terpenoids, alkaloids, polyketides,
//       phenylpropanoids, steroids)
//   6.  Privileged structures in medicinal chemistry
//   7.  Amino acid side chains (20 standard, with pKa/hydrophobicity)
//   8.  Common functional groups with SMARTS patterns
//   9.  Pharmaceutical statistics constants
//
// All data is constexpr / header-only.  Thread-safe by construction.
// ============================================================================

// ---- Forward helpers -------------------------------------------------------

namespace detail {

/// Simple constexpr string equality (null-terminated C-strings).
constexpr bool streq(const char* a, const char* b) {
    while (*a && *b) {
        if (*a != *b) return false;
        ++a;
        ++b;
    }
    return *a == *b;
}

} // namespace detail

// ============================================================================
//  1. Ring Systems — Top 50 by frequency in FDA-approved drugs
// ============================================================================

struct RingSystemEntry {
    const char* smiles;       ///< Canonical SMILES of the ring system
    const char* name;         ///< Common chemical name
    int         ringSize;     ///< Number of atoms in the ring
    int         heteroCount;  ///< Number of heteroatoms (non-carbon ring atoms)
    int         rank;         ///< Frequency rank (1 = most common)
};

inline constexpr std::array<RingSystemEntry, 50> RING_SYSTEMS = {{
    { "c1ccccc1",                          "benzene",            6,  0,  1 },
    { "c1ccncc1",                          "pyridine",           6,  1,  2 },
    { "C1CCNCC1",                          "piperidine",         6,  1,  3 },
    { "C1CNCCN1",                          "piperazine",         6,  2,  4 },
    { "C1CCCCC1",                          "cyclohexane",        6,  0,  5 },
    { "c1ncccn1",                          "pyrimidine",         6,  2,  6 },
    { "C1CCOCC1",                          "tetrahydropyran",    6,  1,  7 },
    { "C1CC1",                             "cyclopropane",       3,  0,  8 },
    { "C1CCNC1",                           "pyrrolidine",        5,  1,  9 },
    { "c1c[nH]cn1",                        "imidazole",          5,  2, 10 },
    { "c1cc[nH]n1",                        "pyrazole",           5,  2, 11 },
    { "c1ccsc1",                           "thiophene",          5,  1, 12 },
    { "c1ccoc1",                           "furan",              5,  1, 13 },
    { "c1cn[nH]n1",                        "1,2,4-triazole",     5,  3, 14 },
    { "C1COCCN1",                          "morpholine",         6,  2, 15 },
    { "c1cscn1",                           "thiazole",           5,  2, 16 },
    { "c1cocn1",                           "oxazole",            5,  2, 17 },
    { "c1ccc2[nH]ccc2c1",                  "indole",             9,  1, 18 },
    { "c1cc[nH]c1",                        "pyrrole",            5,  1, 19 },
    { "c1ccc2ccccc2c1",                    "naphthalene",       10,  0, 20 },
    { "c1ccc2[nH]cnc2c1",                  "benzimidazole",      9,  2, 21 },
    { "c1ccc2ncccc2c1",                    "quinoline",         10,  1, 22 },
    { "c1ccc2ccncc2c1",                    "isoquinoline",      10,  1, 23 },
    { "c1cnccn1",                          "pyrazine",           6,  2, 24 },
    { "c1nnn[nH]1",                        "tetrazole",          5,  4, 25 },
    { "C1CCCC1",                           "cyclopentane",       5,  0, 26 },
    { "c1conc1",                           "isoxazole",          5,  2, 27 },
    { "c1nonn1",                           "1,2,4-oxadiazole",   5,  3, 28 },
    { "c1nsnn1",                           "1,2,4-thiadiazole",  5,  3, 29 },
    { "C1CCCO1",                           "tetrahydrofuran",    5,  1, 30 },
    { "c1nc2[nH]cnc2n1",                   "purine",             9,  4, 31 },
    { "C1=CNC=CC1",                        "dihydropyridine",    6,  1, 32 },
    { "c1ccc2ncncc2c1",                    "quinazoline",       10,  2, 33 },
    { "C1CCN1",                            "azetidine",          4,  1, 34 },
    { "c1ccnnc1",                          "pyridazine",         6,  2, 35 },
    { "c1ccc2occc2c1",                     "benzofuran",         9,  1, 36 },
    { "c1ccc2sccc2c1",                     "benzothiophene",     9,  1, 37 },
    { "c1ccc2c(c1)[nH]c3ccccc32",          "carbazole",         13,  1, 38 },
    { "c1ccc2nc3ccccc3cc2c1",              "acridine",          14,  1, 39 },
    { "O=C1CCN1",                          "beta-lactam",        4,  1, 40 },
    { "O=c1ccoc2ccccc12",                  "chromone",          10,  1, 41 },
    { "O=c1ccc2ccccc2o1",                  "coumarin",          10,  1, 42 },
    { "C1CSCN1",                           "thiazolidine",       5,  2, 43 },
    { "C1COC1",                            "oxetane",            4,  1, 44 },
    { "C1CNCN1",                           "imidazolidine",      5,  2, 45 },
    { "C1CCC1",                            "cyclobutane",        4,  0, 46 },
    { "O=C1CCCN1",                         "2-pyrrolidinone",    5,  1, 47 },
    { "O=C1CCCCN1",                        "2-piperidinone",     6,  1, 48 },
    { "C1CCCCCC1",                         "cycloheptane",       7,  0, 49 },
    { "c1cnc2ccccc2n1",                    "quinoxaline",       10,  2, 50 },
}};

// ============================================================================
//  2. Murcko Scaffolds — Top 30 from approved drugs
// ============================================================================

struct MurckoScaffoldEntry {
    const char* smiles;   ///< Canonical SMILES of the scaffold
    const char* name;     ///< Common name
    int         rank;     ///< Frequency rank (1 = most common)
};

inline constexpr std::array<MurckoScaffoldEntry, 30> MURCKO_SCAFFOLDS = {{
    { "c1ccccc1",                                    "benzene",                   1 },
    { "c1ccc(-c2ccccc2)cc1",                         "biphenyl",                  2 },
    { "c1ccc2ccccc2c1",                              "naphthalene",               3 },
    { "C1CCN(CC1)c1ccccc1",                          "phenylpiperidine",          4 },
    { "c1ccc(cc1)N1CCNCC1",                          "phenylpiperazine",          5 },
    { "C(c1ccccc1)c1ccccc1",                         "diphenylmethane",           6 },
    { "c1ccncc1",                                    "pyridine",                  7 },
    { "c1ccc(-c2ccccn2)cc1",                         "phenyl-pyridine",           8 },
    { "c1ccc2[nH]ccc2c1",                            "indole",                    9 },
    { "c1ccc2ncccc2c1",                              "quinoline",                10 },
    { "c1ccc2[nH]cnc2c1",                            "benzimidazole",            11 },
    { "c1ccc(-c2ccncn2)cc1",                         "phenyl-pyrimidine",        12 },
    { "C1CCC2C(C1)CCC1C2CCC2CCCCC12",               "steroid-gonane",           13 },
    { "C1=CNC=CC1",                                  "dihydropyridine",          14 },
    { "O=c1[nH]c(=O)c2[nH]cn(c2n1)",                "xanthine",                 15 },
    // -- scaffolds 16-30 (extended set) --
    { "O=C1CC2SC(C)(C)C(N2C1)C(=O)O",               "beta-lactam-thiazolidine", 16 },
    { "O=C(O)c1cn2c3ccccc3c(=O)c2cc1",              "fluoroquinolone",          17 },
    { "NC(=O)c1ccccc1",                              "benzamide",                18 },
    { "NCCOc1ccccc1",                                "phenoxyethylamine",        19 },
    { "NS(=O)(=O)c1ccccc1",                          "sulfonamide-benzene",      20 },
    { "OC(=O)Cc1ccccc1",                             "phenylpropanoic-acid",     21 },
    { "c1ccc(-c2cn[nH]n2)cc1",                       "phenyl-triazole",          22 },
    { "c1ccc(-c2nnn[nH]2)cc1",                       "phenyl-tetrazole",         23 },
    { "c1ccc2ccncc2c1",                              "isoquinoline",             24 },
    { "c1nc2[nH]cnc2c(=O)[nH]1",                    "purine",                   25 },
    { "c1cnc2ncncc2n1",                              "pteridine",                26 },
    { "O=C1CCOc2ccccc21",                            "chromanone",               27 },
    { "c1ccc2nc3ccccc3cc2c1",                        "acridine",                 28 },
    { "c1ccc2c(c1)[nH]c1ccccc12",                    "carbazole",                29 },
    { "O=c1cc(-c2ccccc2)oc2ccccc12",                 "flavone",                  30 },
}};

// ============================================================================
//  3. Pharmacophore SMARTS — 8 standard feature families
// ============================================================================

struct PharmacophorePattern {
    const char* featureName;   ///< Feature family name
    const char* smarts;        ///< SMARTS pattern for feature detection
};

inline constexpr std::array<PharmacophorePattern, 18> PHARMACOPHORE_PATTERNS = {{
    // -- 8 standard feature families --
    { "Donor",
      "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]" },
    { "Acceptor",
      "[$([O;H1;v2]),$([O;H0;v2;!$(O=N-*)]),$([O;-;!$(*-N=O)]),"
      "$([o;+0;!$(o~[Nv5])]),$([n;+0;!$(n(cc)cc)])]" },
    { "NegIonizable",
      "[C,S](=[O,S,P])-[O;H1,H0&-1]" },
    { "PosIonizable",
      "[#7;+;!$([N+]-[O-])]" },
    { "Aromatic5",
      "a1aaaa1" },
    { "Aromatic6",
      "a1aaaaa1" },
    { "Hydrophobe",
      "[D3,D4;$([#6;+0;!$([#6]~[#7,#8,F]);!$([F,Cl,Br,I])])]" },
    { "ZnBinder",
      "[S;D1]-[#6]" },
    // -- 10 extended pharmacophore patterns --
    { "Amide",
      "[NX3][CX3](=[OX1])[#6]" },
    { "Sulfonamide",
      "[$([#16X4](=[OX1])(=[OX1])([#6])[NX3])]" },
    { "CarboxylicAcid",
      "[CX3](=O)[OX2H1]" },
    { "PrimaryAmine",
      "[NX3;H2;!$(NC=O)]" },
    { "SecondaryAmine",
      "[NX3;H1;!$(NC=O)]" },
    { "TertiaryAmine",
      "[NX3;H0;!$(NC=O);!$(N=O)]" },
    { "Hydroxyl",
      "[OX2H]" },
    { "Thiol",
      "[SX2H]" },
    { "Ester",
      "[#6][CX3](=O)[OX2H0][#6]" },
    { "Ether",
      "[OD2]([#6])[#6]" },
}};

// ============================================================================
//  4. Drug-like property boundaries
// ============================================================================

/// Lipinski Rule-of-Five thresholds.
struct LipinskiLimits {
    double maxMW;    ///< Molecular weight upper bound (Da)
    double maxLogP;  ///< Calculated LogP upper bound
    int    maxHBD;   ///< Hydrogen-bond donor count upper bound
    int    maxHBA;   ///< Hydrogen-bond acceptor count upper bound
};

inline constexpr LipinskiLimits LIPINSKI_LIMITS = { 500.0, 5.0, 5, 10 };

/// Veber oral bioavailability thresholds.
struct VeberLimits {
    int    maxRotBonds;  ///< Rotatable bond count upper bound
    double maxTPSA;      ///< Topological polar surface area upper bound (A^2)
};

inline constexpr VeberLimits VEBER_LIMITS = { 10, 140.0 };

/// Average properties of FDA-approved small-molecule drugs.
struct AverageDrugProperties {
    double meanMW;     ///< Mean molecular weight (Da)
    double meanLogP;   ///< Mean LogP
    double meanHBD;    ///< Mean hydrogen-bond donor count
    double meanHBA;    ///< Mean hydrogen-bond acceptor count
};

inline constexpr AverageDrugProperties AVG_DRUG_PROPS = { 396.0, 2.5, 1.8, 4.2 };

// ============================================================================
//  5. Natural Product Scaffold Classes
// ============================================================================
//
//  Organised by biosynthetic origin.  Each entry carries the class, subclass,
//  a representative SMILES, the carbon-backbone length, and the biological
//  source kingdom most commonly associated with the scaffold.
// ============================================================================

struct NaturalProductScaffold {
    const char* scaffoldClass;  ///< Biosynthetic class (e.g. "Terpenoid")
    const char* subclass;       ///< Subclass (e.g. "monoterpene")
    const char* smiles;         ///< Representative SMILES
    int         carbonBackbone; ///< Carbon backbone length (isoprene units x 5)
    const char* sourceKingdom;  ///< Dominant source kingdom
};

inline constexpr std::array<NaturalProductScaffold, 23> NATURAL_PRODUCT_SCAFFOLDS = {{
    // -- Terpenoids (isoprene-derived, mevalonate / MEP pathway) --
    { "Terpenoid", "monoterpene",    "CC1=CCC(CC1)C(C)=C",       10, "Plantae"    },
    { "Terpenoid", "sesquiterpene",  "CC1=CC2C(CC1)C(C)=CCC2C",  15, "Plantae"    },
    { "Terpenoid", "diterpene",      "CC(=CCCC(C)=CCCC(C)=CCCC(C)=C)C", 20, "Plantae" },
    { "Terpenoid", "triterpene",     "CC1(C)CCC2(C)C(C1)CCC1C3(C)CCC(O)C(C)(C)C3CCC12C", 30, "Plantae" },
    { "Terpenoid", "tetraterpene",   "CC(=CC=CC(C)=CC=CC(C)=CC=CC=C(C)C=CC=C(C)C=CC=C(C)C)C", 40, "Plantae" },
    // -- Alkaloids (nitrogen-containing, diverse biosynthetic origins) --
    { "Alkaloid",  "indole",         "c1ccc2[nH]ccc2c1",         9,  "Plantae"    },
    { "Alkaloid",  "isoquinoline",   "c1ccc2ccncc2c1",          10,  "Plantae"    },
    { "Alkaloid",  "pyridine",       "c1ccncc1",                 6,  "Fungi"      },
    { "Alkaloid",  "tropane",        "C1CC2CCC(C1)N2",           8,  "Plantae"    },
    { "Alkaloid",  "quinoline",      "c1ccc2ncccc2c1",          10,  "Plantae"    },
    { "Alkaloid",  "purine",         "c1nc2[nH]cnc2c(=O)[nH]1",  9,  "universal"  },
    // -- Polyketides (acetate / malonate pathway) --
    { "Polyketide", "macrolide",     "O=C1CCCCCCCCCCCCOC1",     14,  "Bacteria"   },
    { "Polyketide", "tetracycline",  "C1CC2C(CC1O)C(=O)c1cccc(O)c1C2(O)O", 18, "Bacteria" },
    { "Polyketide", "anthracycline", "O=C1c2cccc(O)c2C(=O)c2cccc(O)c21", 14, "Bacteria" },
    // -- Phenylpropanoids (shikimate / phenylpropanoid pathway) --
    { "Phenylpropanoid", "flavonoid",  "O=C1CC(c2ccccc2)Oc2ccccc21",     15, "Plantae" },
    { "Phenylpropanoid", "coumarin",   "O=c1ccc2ccccc2o1",               10, "Plantae" },
    { "Phenylpropanoid", "lignan",     "C(c1ccc(O)cc1)c1ccc(O)cc1",      18, "Plantae" },
    { "Phenylpropanoid", "stilbene",   "C(=Cc1ccccc1)c1ccccc1",          14, "Plantae" },
    // -- Steroids (modified triterpene, cholesterol pathway) --
    { "Steroid", "androstane",  "C1CCC2C(C1)CCC1C2CCC2CCCCC12",    19, "Animalia" },
    { "Steroid", "estrane",     "C1CCC2C(C1)CCC1C2CCc2ccccc21",    18, "Animalia" },
    { "Steroid", "pregnane",    "CC(=O)C1CCC2C1CCC1C2CCC2CCCCC12", 21, "Animalia" },
    { "Steroid", "cholestane",  "CC(CCCC(C)C)C1CCC2C1CCC1C2CCC2CCCCC12", 27, "Animalia" },
    { "Steroid", "ergostane",   "CC(C)C(C)CCC(C)C1CCC2C1CCC1C2CCC2CCCCC12", 28, "Fungi" },
}};

// ============================================================================
//  6. Privileged Structures in Medicinal Chemistry
// ============================================================================
//
//  Privileged scaffolds recur across diverse target classes (GPCRs, kinases,
//  nuclear receptors, ion channels).  They provide starting points that are
//  biased toward bioactivity while retaining drug-like ADME properties.
// ============================================================================

struct PrivilegedStructure {
    const char* smiles;        ///< Canonical SMILES
    const char* name;          ///< Common name
    const char* targetClass;   ///< Primary target class association
};

inline constexpr std::array<PrivilegedStructure, 16> PRIVILEGED_STRUCTURES = {{
    // -- Acyclic privileged motifs --
    { "NC(=O)c1ccccc1",                       "benzamide",               "Kinase / HDAC"        },
    { "NS(=O)(=O)c1ccccc1",                   "sulfonamide",             "Carbonic anhydrase"   },
    { "NNC(=O)c1ccccc1",                      "hydrazide",               "MAO / Tuberculosis"   },
    { "NC(=O)Nc1ccccc1",                      "urea",                    "Kinase / sEH"         },
    { "NC(=O)Oc1ccccc1",                      "carbamate",               "AChE / BChE"          },
    { "O=C1CN(c2ccccc2)C(=O)O1",              "oxazolidinone",           "Ribosome (antibacterial)" },
    // -- Biaryl / conjugated privileged motifs --
    { "c1ccc(-c2ccccc2)cc1",                   "biphenyl",                "ARB / GPCR"           },
    { "C(=Cc1ccccc1)c1ccccc1",                "stilbene",                "Estrogen receptor"    },
    { "O=C(/C=C/c1ccccc1)c1ccccc1",           "chalcone",                "NF-kB / tubulin"      },
    // -- Fused bicyclic privileged motifs --
    { "c1ccc2[nH]cnc2c1",                      "benzimidazole",           "PPI / antiparasitic"  },
    { "c1ccc2ocnc2c1",                         "benzoxazole",             "CB2 / anticonvulsant" },
    { "c1ccc2scnc2c1",                         "benzothiazole",           "Amyloid / kinase"     },
    // -- Partially saturated privileged motifs --
    { "C1=CNC=CC1",                            "dihydropyridine",         "L-type Ca channel"    },
    { "c1ccc2c(c1)CCNC2",                      "tetrahydroisoquinoline",  "Opioid / dopamine"    },
    // -- Additional privileged heterocycles --
    { "O=c1ccoc2ccccc12",                      "chromone",                "Kinase / antifungal"  },
    { "c1ccnc2[nH]ccc12",                      "pyrrolo-pyridine",        "Kinase (CDK/JAK)"     },
}};

// ============================================================================
//  7. Amino Acid Side Chains — 20 Standard Amino Acids
// ============================================================================
//
//  Each entry: full name, three-letter code, one-letter code, side-chain
//  SMILES (R-group only, not including backbone), molecular weight of the
//  full amino acid, side-chain pKa (0.0 if not ionisable at pH 7), and
//  Kyte-Doolittle hydrophobicity index.
// ============================================================================

struct AminoAcidEntry {
    const char* name;        ///< Full name
    const char* code3;       ///< Three-letter code
    char        code1;       ///< One-letter code
    const char* sideChain;   ///< Side-chain SMILES (R-group)
    double      mw;          ///< Molecular weight of full amino acid (Da)
    double      pKaSide;     ///< Side-chain pKa (0.0 = not titratable)
    double      hydropathy;  ///< Kyte-Doolittle hydropathy index
};

inline constexpr std::array<AminoAcidEntry, 20> AMINO_ACIDS = {{
    { "Glycine",       "Gly", 'G', "[H]",                  75.03,  0.0, -0.4 },
    { "Alanine",       "Ala", 'A', "C",                    89.09,  0.0,  1.8 },
    { "Valine",        "Val", 'V', "CC(C)",               117.15,  0.0,  4.2 },
    { "Leucine",       "Leu", 'L', "CC(C)C",              131.17,  0.0,  3.8 },
    { "Isoleucine",    "Ile", 'I', "C(CC)C",              131.17,  0.0,  4.5 },
    { "Proline",       "Pro", 'P', "C1CNC1",              115.13,  0.0, -1.6 },
    { "Phenylalanine", "Phe", 'F', "Cc1ccccc1",           165.19,  0.0,  2.8 },
    { "Tryptophan",    "Trp", 'W', "Cc1c[nH]c2ccccc12",  204.23,  0.0, -0.9 },
    { "Methionine",    "Met", 'M', "CCSC",                149.21,  0.0,  1.9 },
    { "Serine",        "Ser", 'S', "CO",                  105.09,  0.0, -0.8 },
    { "Threonine",     "Thr", 'T', "C(C)O",               119.12,  0.0, -0.7 },
    { "Cysteine",      "Cys", 'C', "CS",                  121.16,  8.3,  2.5 },
    { "Tyrosine",      "Tyr", 'Y', "Cc1ccc(O)cc1",        181.19, 10.1, -1.3 },
    { "Asparagine",    "Asn", 'N', "CC(=O)N",             132.12,  0.0, -3.5 },
    { "Glutamine",     "Gln", 'Q', "CCC(=O)N",            146.15,  0.0, -3.5 },
    { "Aspartate",     "Asp", 'D', "CC(=O)O",             133.10,  3.65,-3.5 },
    { "Glutamate",     "Glu", 'E', "CCC(=O)O",            147.13,  4.25,-3.5 },
    { "Lysine",        "Lys", 'K', "CCCCN",               146.19, 10.5, -3.9 },
    { "Arginine",      "Arg", 'R', "CCCNC(=N)N",          174.20, 12.5, -4.5 },
    { "Histidine",     "His", 'H', "Cc1c[nH]cn1",         155.16,  6.0, -3.2 },
}};

// ============================================================================
//  8. Common Functional Groups with SMARTS Patterns
// ============================================================================
//
//  Each entry: name, SMARTS pattern for substructure matching, and the
//  typical aqueous pKa range (pKaLow..pKaHigh).  A range of 0.0/0.0
//  indicates the group is not commonly ionisable in aqueous solution.
// ============================================================================

struct FunctionalGroupEntry {
    const char* name;       ///< Functional group name
    const char* smarts;     ///< SMARTS substructure pattern
    double      pKaLow;    ///< Typical pKa lower bound (0.0 = N/A)
    double      pKaHigh;   ///< Typical pKa upper bound (0.0 = N/A)
};

inline constexpr std::array<FunctionalGroupEntry, 22> FUNCTIONAL_GROUPS = {{
    { "hydroxyl",     "[OX2H]",                            15.0, 19.0 },
    { "carboxyl",     "[CX3](=O)[OX2H1]",                  2.0,  5.0 },
    { "amino",        "[NX3;H2;!$(NC=O)]",                  9.0, 11.0 },
    { "amide",        "[NX3][CX3](=[OX1])[#6]",             0.0,  0.0 },
    { "ester",        "[#6][CX3](=O)[OX2H0][#6]",           0.0,  0.0 },
    { "ether",        "[OD2]([#6])[#6]",                     0.0,  0.0 },
    { "thiol",        "[SX2H]",                              8.0, 11.0 },
    { "sulfide",      "[#16X2H0]",                           0.0,  0.0 },
    { "disulfide",    "[#16X2][#16X2]",                      0.0,  0.0 },
    { "nitro",        "[NX3+](=O)[O-]",                      0.0,  0.0 },
    { "nitrile",      "[NX1]#[CX2]",                         0.0,  0.0 },
    { "aldehyde",     "[CX3H1](=O)[#6]",                     0.0,  0.0 },
    { "ketone",       "[#6][CX3](=O)[#6]",                   0.0,  0.0 },
    { "imine",        "[CX3;$([C]([#6])[#6]),$([CH][#6])]=[NX2][#6]", 0.0, 0.0 },
    { "phosphate",    "[PX4](=[OX1])([OX2H])([OX2H])[OX2H]", 1.0, 13.0 },
    { "sulfonate",    "[#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1-]", -1.5, 1.5 },
    { "sulfonamide",  "[#16X4](=[OX1])(=[OX1])([#6])[NX3]",  0.0,  0.0 },
    { "fluoride",     "[FX1][#6]",                            0.0,  0.0 },
    { "chloride",     "[ClX1][#6]",                           0.0,  0.0 },
    { "bromide",      "[BrX1][#6]",                           0.0,  0.0 },
    { "iodide",       "[IX1][#6]",                            0.0,  0.0 },
    { "epoxide",      "[OX2r3]1[#6r3][#6r3]1",               0.0,  0.0 },
}};

// ============================================================================
//  9. Pharmaceutical Statistics Constants
// ============================================================================
//
//  Curated from published analyses of FDA-approved small-molecule drugs.
//  Useful for benchmarking compound libraries and chemical space coverage.
// ============================================================================

struct PharmaceuticalStats {
    /// Fraction of FDA-approved drugs containing at least one N-heterocycle.
    double heterocyclicPrevalence;
    /// Fraction of FDA-approved drugs containing aza-heterocycles.
    double azaHeterocyclePrevalence;
    /// Average number of ring systems per approved drug.
    double avgRingCountPerDrug;
    /// Average Fsp3 (fraction sp3 carbons) for natural products.
    double naturalProductAvgFsp3;
    /// Average Fsp3 for synthetic drugs.
    double syntheticDrugAvgFsp3;
    /// Fraction of ring-system occurrences covered by the top-10 systems.
    double top10RingCoverage;
    /// Average molecular weight of FDA oral drugs (Da).
    double avgMWOralDrugs;
    /// Average heavy atom count of FDA oral drugs.
    double avgHeavyAtomCount;
    /// Average number of chiral centres per approved drug.
    double avgChiralCentres;
    /// Fraction of approved drugs containing at least one stereogenic centre.
    double fractionWithStereoCentre;
};

inline constexpr PharmaceuticalStats PHARMA_STATS = {
    0.59,   // 59% of FDA drugs contain N-heterocycles
    0.80,   // 80% contain aza-heterocycles
    2.1,    // average ring count per drug
    0.47,   // natural product average Fsp3
    0.32,   // synthetic drug average Fsp3
    0.90,   // top 10 ring systems cover 90% of occurrences
    371.0,  // average MW of FDA oral drugs (Da)
    26.0,   // average heavy atom count
    1.6,    // average chiral centres per drug
    0.56,   // 56% of approved drugs have at least one stereocentre
};

// ============================================================================
//  Lookup / query functions
// ============================================================================

/// Return the total number of ring systems in the library.
constexpr std::size_t ringSystemCount() noexcept {
    return RING_SYSTEMS.size();
}

/// Return the total number of Murcko scaffolds in the library.
constexpr std::size_t murckoScaffoldCount() noexcept {
    return MURCKO_SCAFFOLDS.size();
}

/// Return the total number of natural product scaffolds in the library.
constexpr std::size_t naturalProductScaffoldCount() noexcept {
    return NATURAL_PRODUCT_SCAFFOLDS.size();
}

/// Return the total number of privileged structures in the library.
constexpr std::size_t privilegedStructureCount() noexcept {
    return PRIVILEGED_STRUCTURES.size();
}

/// Return the total number of amino acid entries.
constexpr std::size_t aminoAcidCount() noexcept {
    return AMINO_ACIDS.size();
}

/// Return the total number of functional group entries.
constexpr std::size_t functionalGroupCount() noexcept {
    return FUNCTIONAL_GROUPS.size();
}

/// Look up a ring system by common name.
/// Returns nullptr if not found.
constexpr const RingSystemEntry* getRingSystemByName(const char* name) noexcept {
    for (const auto& e : RING_SYSTEMS) {
        if (detail::streq(e.name, name)) return &e;
    }
    return nullptr;
}

/// Look up a ring system by frequency rank (1-based).
/// Returns nullptr if rank is out of range.
constexpr const RingSystemEntry* getRingSystemByRank(int rank) noexcept {
    if (rank < 1 || rank > static_cast<int>(RING_SYSTEMS.size())) return nullptr;
    // The array is sorted by rank, so index = rank - 1.
    return &RING_SYSTEMS[static_cast<std::size_t>(rank - 1)];
}

/// Look up a Murcko scaffold by common name.
/// Returns nullptr if not found.
constexpr const MurckoScaffoldEntry* getMurckoScaffold(const char* name) noexcept {
    for (const auto& e : MURCKO_SCAFFOLDS) {
        if (detail::streq(e.name, name)) return &e;
    }
    return nullptr;
}

/// Look up a natural product scaffold by subclass name.
/// Returns nullptr if not found.
constexpr const NaturalProductScaffold*
getNaturalProductBySubclass(const char* subclass) noexcept {
    for (const auto& e : NATURAL_PRODUCT_SCAFFOLDS) {
        if (detail::streq(e.subclass, subclass)) return &e;
    }
    return nullptr;
}

/// Look up a privileged structure by common name.
/// Returns nullptr if not found.
constexpr const PrivilegedStructure*
getPrivilegedStructure(const char* name) noexcept {
    for (const auto& e : PRIVILEGED_STRUCTURES) {
        if (detail::streq(e.name, name)) return &e;
    }
    return nullptr;
}

/// Look up an amino acid by one-letter code.
/// Returns nullptr if the code is not recognised.
constexpr const AminoAcidEntry* getAminoAcidByCode1(char code) noexcept {
    for (const auto& e : AMINO_ACIDS) {
        if (e.code1 == code) return &e;
    }
    return nullptr;
}

/// Look up an amino acid by three-letter code.
/// Returns nullptr if the code is not recognised.
constexpr const AminoAcidEntry* getAminoAcidByCode3(const char* code3) noexcept {
    for (const auto& e : AMINO_ACIDS) {
        if (detail::streq(e.code3, code3)) return &e;
    }
    return nullptr;
}

/// Look up a functional group by name.
/// Returns nullptr if not found.
constexpr const FunctionalGroupEntry*
getFunctionalGroup(const char* name) noexcept {
    for (const auto& e : FUNCTIONAL_GROUPS) {
        if (detail::streq(e.name, name)) return &e;
    }
    return nullptr;
}

/// Check whether a SMILES string matches any of the top-50 ring systems.
constexpr bool isCommonRingSystem(const char* smiles) noexcept {
    for (const auto& e : RING_SYSTEMS) {
        if (detail::streq(e.smiles, smiles)) return true;
    }
    return false;
}

/// Get the SMARTS pattern for a pharmacophore feature family by name.
/// Returns nullptr if the feature name is not found.
constexpr const char* getPharmacophorePattern(const char* featureName) noexcept {
    for (const auto& p : PHARMACOPHORE_PATTERNS) {
        if (detail::streq(p.featureName, featureName)) return p.smarts;
    }
    return nullptr;
}

// ============================================================================
//  Drug-likeness compliance checks
// ============================================================================

/// Check Lipinski Rule-of-Five compliance.
/// Returns true if the molecule satisfies all four conditions.
constexpr bool isLipinskiCompliant(double mw, double logp,
                                   int hbd, int hba) noexcept {
    return mw  <= LIPINSKI_LIMITS.maxMW
        && logp <= LIPINSKI_LIMITS.maxLogP
        && hbd  <= LIPINSKI_LIMITS.maxHBD
        && hba  <= LIPINSKI_LIMITS.maxHBA;
}

/// Check Veber oral bioavailability compliance.
/// Returns true if rotatable bonds and TPSA are within thresholds.
constexpr bool isVeberCompliant(int rotBonds, double tpsa) noexcept {
    return rotBonds <= VEBER_LIMITS.maxRotBonds
        && tpsa     <= VEBER_LIMITS.maxTPSA;
}

} // namespace smsd

#endif // SMSD_SCAFFOLD_LIBRARY_HPP
