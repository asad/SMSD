/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. */
#pragma once
#ifndef SMSD_PERIODIC_TABLE_HPP
#define SMSD_PERIODIC_TABLE_HPP

// ============================================================================
// Authoritative periodic-table element property lookup for SMSD.
//
// Data sources:
//   Atomic weights      — IUPAC 2021 standard atomic weights
//   Electronegativity   — Pauling scale (Allred & Rochow revised values)
//   Covalent radii      — Cordero et al. 2008, Dalton Trans. 2832-2838
//   Van der Waals radii — Bondi 1964 / Alvarez 2013 (main group extension)
//   Default valences    — OpenSMILES / Daylight convention
//
// Header-only, constexpr where possible, C++17.
// ============================================================================

#include <array>
#include <cmath>
#include <string>
#include <unordered_map>

namespace smsd {

// ============================================================================
// Element data structure
// ============================================================================

struct ElementData {
    int         atomicNum;
    const char* symbol;
    const char* name;
    double      atomicWeight;       // standard atomic weight (IUPAC 2021)
    double      electronegativity;  // Pauling scale (0.0 if unknown)
    double      covalentRadius;     // Angstroms, single bond (Cordero 2008)
    double      vanDerWaalsRadius;  // Angstroms (Bondi 1964 / Alvarez 2013)
    int         maxValence;         // maximum common valence
    int         defaultValences[4]; // common valences, terminated by -1
    int         group;              // periodic table group 1-18 (0 = lanthanide/actinide)
    int         period;             // periodic table period 1-7
    const char* block;              // "s", "p", "d", "f"
    const char* electronConfig;     // abbreviated electron configuration
};

// ============================================================================
// Full periodic table: elements 0 (dummy) through 118 (Oganesson)
//
// Index 0 is a wildcard/dummy atom used for query atoms ("*" in SMILES).
// Elements 1-86 carry experimentally verified data.
// Elements 87-118 use best available data or calculated estimates.
// ============================================================================

constexpr ElementData PERIODIC_TABLE[119] = {
    // Z=0: Dummy/wildcard atom
    {  0, "*",  "Dummy",        0.0,   0.0,  0.00, 0.00,  0, {-1,-1,-1,-1},  0, 0, "",  ""},

    // ---- Period 1 ----
    {  1, "H",  "Hydrogen",     1.008,  2.20, 0.31, 1.20,  1, { 1,-1,-1,-1},  1, 1, "s", "1s1"},
    {  2, "He", "Helium",       4.0026, 0.00, 0.28, 1.40,  0, {-1,-1,-1,-1}, 18, 1, "s", "1s2"},

    // ---- Period 2 ----
    {  3, "Li", "Lithium",      6.94,   0.98, 1.28, 2.12,  1, { 1,-1,-1,-1},  1, 2, "s", "[He]2s1"},
    {  4, "Be", "Beryllium",    9.0122, 1.57, 0.96, 1.85,  2, { 2,-1,-1,-1},  2, 2, "s", "[He]2s2"},
    {  5, "B",  "Boron",       10.81,   2.04, 0.84, 1.92,  3, { 3,-1,-1,-1}, 13, 2, "p", "[He]2s2 2p1"},
    {  6, "C",  "Carbon",      12.011,  2.55, 0.76, 1.70,  4, { 4,-1,-1,-1}, 14, 2, "p", "[He]2s2 2p2"},
    {  7, "N",  "Nitrogen",    14.007,  3.04, 0.71, 1.55,  5, { 3, 5,-1,-1}, 15, 2, "p", "[He]2s2 2p3"},
    {  8, "O",  "Oxygen",      15.999,  3.44, 0.66, 1.52,  2, { 2,-1,-1,-1}, 16, 2, "p", "[He]2s2 2p4"},
    {  9, "F",  "Fluorine",    18.998,  3.98, 0.57, 1.47,  1, { 1,-1,-1,-1}, 17, 2, "p", "[He]2s2 2p5"},
    { 10, "Ne", "Neon",        20.180,  0.00, 0.58, 1.54,  0, {-1,-1,-1,-1}, 18, 2, "p", "[He]2s2 2p6"},

    // ---- Period 3 ----
    { 11, "Na", "Sodium",      22.990,  0.93, 1.66, 2.27,  1, { 1,-1,-1,-1},  1, 3, "s", "[Ne]3s1"},
    { 12, "Mg", "Magnesium",   24.305,  1.31, 1.41, 1.73,  2, { 2,-1,-1,-1},  2, 3, "s", "[Ne]3s2"},
    { 13, "Al", "Aluminium",   26.982,  1.61, 1.21, 2.05,  3, { 3,-1,-1,-1}, 13, 3, "p", "[Ne]3s2 3p1"},
    { 14, "Si", "Silicon",     28.086,  1.90, 1.11, 2.10,  4, { 4,-1,-1,-1}, 14, 3, "p", "[Ne]3s2 3p2"},
    { 15, "P",  "Phosphorus",  30.974,  2.19, 1.07, 1.80,  5, { 3, 5,-1,-1}, 15, 3, "p", "[Ne]3s2 3p3"},
    { 16, "S",  "Sulfur",      32.06,   2.58, 1.05, 1.80,  6, { 2, 4, 6,-1}, 16, 3, "p", "[Ne]3s2 3p4"},
    { 17, "Cl", "Chlorine",    35.45,   3.16, 1.02, 1.75,  7, { 1,-1,-1,-1}, 17, 3, "p", "[Ne]3s2 3p5"},
    { 18, "Ar", "Argon",       39.95,   0.00, 1.06, 1.88,  0, {-1,-1,-1,-1}, 18, 3, "p", "[Ne]3s2 3p6"},

    // ---- Period 4 ----
    { 19, "K",  "Potassium",   39.098,  0.82, 2.03, 2.75,  1, { 1,-1,-1,-1},  1, 4, "s", "[Ar]4s1"},
    { 20, "Ca", "Calcium",     40.078,  1.00, 1.76, 2.45,  2, { 2,-1,-1,-1},  2, 4, "s", "[Ar]4s2"},
    { 21, "Sc", "Scandium",    44.956,  1.36, 1.70, 2.25,  3, { 3,-1,-1,-1},  3, 4, "d", "[Ar]3d1 4s2"},
    { 22, "Ti", "Titanium",    47.867,  1.54, 1.60, 2.10,  4, { 4,-1,-1,-1},  4, 4, "d", "[Ar]3d2 4s2"},
    { 23, "V",  "Vanadium",    50.942,  1.63, 1.53, 2.05,  5, { 5,-1,-1,-1},  5, 4, "d", "[Ar]3d3 4s2"},
    { 24, "Cr", "Chromium",    51.996,  1.66, 1.39, 2.00,  6, { 6,-1,-1,-1},  6, 4, "d", "[Ar]3d5 4s1"},
    { 25, "Mn", "Manganese",   54.938,  1.55, 1.50, 2.00,  7, { 2, 4, 7,-1},  7, 4, "d", "[Ar]3d5 4s2"},
    { 26, "Fe", "Iron",        55.845,  1.83, 1.42, 2.00,  6, { 2, 3, 6,-1},  8, 4, "d", "[Ar]3d6 4s2"},
    { 27, "Co", "Cobalt",      58.933,  1.88, 1.38, 1.95,  4, { 2, 3,-1,-1},  9, 4, "d", "[Ar]3d7 4s2"},
    { 28, "Ni", "Nickel",      58.693,  1.91, 1.24, 1.63,  4, { 2,-1,-1,-1}, 10, 4, "d", "[Ar]3d8 4s2"},
    { 29, "Cu", "Copper",      63.546,  1.90, 1.32, 1.40,  4, { 1, 2,-1,-1}, 11, 4, "d", "[Ar]3d10 4s1"},
    { 30, "Zn", "Zinc",        65.38,   1.65, 1.22, 1.39,  2, { 2,-1,-1,-1}, 12, 4, "d", "[Ar]3d10 4s2"},
    { 31, "Ga", "Gallium",     69.723,  1.81, 1.22, 1.87,  3, { 3,-1,-1,-1}, 13, 4, "p", "[Ar]3d10 4s2 4p1"},
    { 32, "Ge", "Germanium",   72.630,  2.01, 1.20, 2.05,  4, { 4,-1,-1,-1}, 14, 4, "p", "[Ar]3d10 4s2 4p2"},
    { 33, "As", "Arsenic",     74.922,  2.18, 1.19, 1.85,  5, { 3, 5,-1,-1}, 15, 4, "p", "[Ar]3d10 4s2 4p3"},
    { 34, "Se", "Selenium",    78.971,  2.55, 1.20, 1.90,  6, { 2, 4, 6,-1}, 16, 4, "p", "[Ar]3d10 4s2 4p4"},
    { 35, "Br", "Bromine",     79.904,  2.96, 1.20, 1.85,  7, { 1,-1,-1,-1}, 17, 4, "p", "[Ar]3d10 4s2 4p5"},
    { 36, "Kr", "Krypton",     83.798,  3.00, 1.16, 2.02,  2, { 2,-1,-1,-1}, 18, 4, "p", "[Ar]3d10 4s2 4p6"},

    // ---- Period 5 ----
    { 37, "Rb", "Rubidium",    85.468,  0.82, 2.20, 3.00,  1, { 1,-1,-1,-1},  1, 5, "s", "[Kr]5s1"},
    { 38, "Sr", "Strontium",   87.62,   0.95, 1.95, 2.60,  2, { 2,-1,-1,-1},  2, 5, "s", "[Kr]5s2"},
    { 39, "Y",  "Yttrium",     88.906,  1.22, 1.90, 2.40,  3, { 3,-1,-1,-1},  3, 5, "d", "[Kr]4d1 5s2"},
    { 40, "Zr", "Zirconium",   91.224,  1.33, 1.75, 2.30,  4, { 4,-1,-1,-1},  4, 5, "d", "[Kr]4d2 5s2"},
    { 41, "Nb", "Niobium",     92.906,  1.60, 1.64, 2.15,  5, { 5,-1,-1,-1},  5, 5, "d", "[Kr]4d4 5s1"},
    { 42, "Mo", "Molybdenum",  95.95,   2.16, 1.54, 2.10,  6, { 6,-1,-1,-1},  6, 5, "d", "[Kr]4d5 5s1"},
    { 43, "Tc", "Technetium",  97.0,    1.90, 1.47, 2.10,  7, { 7,-1,-1,-1},  7, 5, "d", "[Kr]4d5 5s2"},
    { 44, "Ru", "Ruthenium",  101.07,   2.20, 1.46, 2.05,  8, { 3, 4, 8,-1},  8, 5, "d", "[Kr]4d7 5s1"},
    { 45, "Rh", "Rhodium",    102.91,   2.28, 1.42, 2.00,  6, { 3,-1,-1,-1},  9, 5, "d", "[Kr]4d8 5s1"},
    { 46, "Pd", "Palladium",  106.42,   2.20, 1.39, 1.63,  4, { 2, 4,-1,-1}, 10, 5, "d", "[Kr]4d10"},
    { 47, "Ag", "Silver",     107.87,   1.93, 1.45, 1.72,  3, { 1,-1,-1,-1}, 11, 5, "d", "[Kr]4d10 5s1"},
    { 48, "Cd", "Cadmium",    112.41,   1.69, 1.44, 1.62,  2, { 2,-1,-1,-1}, 12, 5, "d", "[Kr]4d10 5s2"},
    { 49, "In", "Indium",     114.82,   1.78, 1.42, 1.93,  3, { 3,-1,-1,-1}, 13, 5, "p", "[Kr]4d10 5s2 5p1"},
    { 50, "Sn", "Tin",        118.71,   1.96, 1.39, 2.17,  4, { 2, 4,-1,-1}, 14, 5, "p", "[Kr]4d10 5s2 5p2"},
    { 51, "Sb", "Antimony",   121.76,   2.05, 1.39, 2.25,  5, { 3, 5,-1,-1}, 15, 5, "p", "[Kr]4d10 5s2 5p3"},
    { 52, "Te", "Tellurium",  127.60,   2.10, 1.38, 2.06,  6, { 2, 4, 6,-1}, 16, 5, "p", "[Kr]4d10 5s2 5p4"},
    { 53, "I",  "Iodine",     126.904,  2.66, 1.39, 1.98,  7, { 1, 3, 5, 7}, 17, 5, "p", "[Kr]4d10 5s2 5p5"},
    { 54, "Xe", "Xenon",      131.29,   2.60, 1.40, 2.16,  8, { 2, 4, 6, 8}, 18, 5, "p", "[Kr]4d10 5s2 5p6"},

    // ---- Period 6: Cs, Ba ----
    { 55, "Cs", "Caesium",     132.91,  0.79, 2.44, 3.15,  1, { 1,-1,-1,-1},  1, 6, "s", "[Xe]6s1"},
    { 56, "Ba", "Barium",     137.33,   0.89, 2.15, 2.70,  2, { 2,-1,-1,-1},  2, 6, "s", "[Xe]6s2"},

    // ---- Lanthanides (Z=57-71) ----
    { 57, "La", "Lanthanum",   138.91,  1.10, 2.07, 2.50,  3, { 3,-1,-1,-1},  3, 6, "f", "[Xe]5d1 6s2"},
    { 58, "Ce", "Cerium",      140.12,  1.12, 2.04, 2.48,  4, { 3, 4,-1,-1},  0, 6, "f", "[Xe]4f1 5d1 6s2"},
    { 59, "Pr", "Praseodymium",140.91,  1.13, 2.03, 2.47,  4, { 3, 4,-1,-1},  0, 6, "f", "[Xe]4f3 6s2"},
    { 60, "Nd", "Neodymium",   144.24,  1.14, 2.01, 2.45,  3, { 3,-1,-1,-1},  0, 6, "f", "[Xe]4f4 6s2"},
    { 61, "Pm", "Promethium",  145.0,   1.13, 1.99, 2.43,  3, { 3,-1,-1,-1},  0, 6, "f", "[Xe]4f5 6s2"},
    { 62, "Sm", "Samarium",    150.36,  1.17, 1.98, 2.42,  3, { 2, 3,-1,-1},  0, 6, "f", "[Xe]4f6 6s2"},
    { 63, "Eu", "Europium",    151.96,  1.12, 1.98, 2.40,  3, { 2, 3,-1,-1},  0, 6, "f", "[Xe]4f7 6s2"},
    { 64, "Gd", "Gadolinium",  157.25,  1.20, 1.96, 2.38,  3, { 3,-1,-1,-1},  0, 6, "f", "[Xe]4f7 5d1 6s2"},
    { 65, "Tb", "Terbium",     158.93,  1.10, 1.94, 2.37,  4, { 3, 4,-1,-1},  0, 6, "f", "[Xe]4f9 6s2"},
    { 66, "Dy", "Dysprosium",  162.50,  1.22, 1.92, 2.35,  3, { 3,-1,-1,-1},  0, 6, "f", "[Xe]4f10 6s2"},
    { 67, "Ho", "Holmium",     164.93,  1.23, 1.92, 2.33,  3, { 3,-1,-1,-1},  0, 6, "f", "[Xe]4f11 6s2"},
    { 68, "Er", "Erbium",      167.26,  1.24, 1.89, 2.32,  3, { 3,-1,-1,-1},  0, 6, "f", "[Xe]4f12 6s2"},
    { 69, "Tm", "Thulium",     168.93,  1.25, 1.90, 2.30,  3, { 2, 3,-1,-1},  0, 6, "f", "[Xe]4f13 6s2"},
    { 70, "Yb", "Ytterbium",   173.05,  1.10, 1.87, 2.28,  3, { 2, 3,-1,-1},  0, 6, "f", "[Xe]4f14 6s2"},
    { 71, "Lu", "Lutetium",    174.97,  1.27, 1.87, 2.25,  3, { 3,-1,-1,-1},  3, 6, "d", "[Xe]4f14 5d1 6s2"},

    // ---- Period 6 d-block (Z=72-80) ----
    { 72, "Hf", "Hafnium",     178.49,  1.30, 1.75, 2.25,  4, { 4,-1,-1,-1},  4, 6, "d", "[Xe]4f14 5d2 6s2"},
    { 73, "Ta", "Tantalum",    180.95,  1.50, 1.70, 2.20,  5, { 5,-1,-1,-1},  5, 6, "d", "[Xe]4f14 5d3 6s2"},
    { 74, "W",  "Tungsten",    183.84,  2.36, 1.62, 2.15,  6, { 6,-1,-1,-1},  6, 6, "d", "[Xe]4f14 5d4 6s2"},
    { 75, "Re", "Rhenium",     186.21,  1.90, 1.51, 2.10,  7, { 4, 7,-1,-1},  7, 6, "d", "[Xe]4f14 5d5 6s2"},
    { 76, "Os", "Osmium",      190.23,  2.20, 1.44, 2.00,  8, { 4, 8,-1,-1},  8, 6, "d", "[Xe]4f14 5d6 6s2"},
    { 77, "Ir", "Iridium",     192.22,  2.20, 1.41, 2.00,  6, { 3, 4, 6,-1},  9, 6, "d", "[Xe]4f14 5d7 6s2"},
    { 78, "Pt", "Platinum",    195.08,  2.28, 1.36, 1.72,  6, { 2, 4,-1,-1}, 10, 6, "d", "[Xe]4f14 5d9 6s1"},
    { 79, "Au", "Gold",        196.97,  2.54, 1.36, 1.66,  5, { 1, 3,-1,-1}, 11, 6, "d", "[Xe]4f14 5d10 6s1"},
    { 80, "Hg", "Mercury",     200.59,  2.00, 1.32, 1.70,  4, { 1, 2,-1,-1}, 12, 6, "d", "[Xe]4f14 5d10 6s2"},

    // ---- Period 6 p-block (Z=81-86) ----
    { 81, "Tl", "Thallium",    204.38,  1.62, 1.45, 1.96,  3, { 1, 3,-1,-1}, 13, 6, "p", "[Xe]4f14 5d10 6s2 6p1"},
    { 82, "Pb", "Lead",        207.2,   2.33, 1.46, 2.02,  4, { 2, 4,-1,-1}, 14, 6, "p", "[Xe]4f14 5d10 6s2 6p2"},
    { 83, "Bi", "Bismuth",     208.98,  2.02, 1.48, 2.35,  5, { 3, 5,-1,-1}, 15, 6, "p", "[Xe]4f14 5d10 6s2 6p3"},
    { 84, "Po", "Polonium",    209.0,   2.00, 1.40, 1.90,  6, { 2, 4,-1,-1}, 16, 6, "p", "[Xe]4f14 5d10 6s2 6p4"},
    { 85, "At", "Astatine",    210.0,   2.20, 1.50, 2.00,  7, { 1, 3, 5, 7}, 17, 6, "p", "[Xe]4f14 5d10 6s2 6p5"},
    { 86, "Rn", "Radon",       222.0,   0.00, 1.50, 2.20,  2, { 2,-1,-1,-1}, 18, 6, "p", "[Xe]4f14 5d10 6s2 6p6"},

    // ---- Period 7: Fr, Ra ----
    { 87, "Fr", "Francium",    223.0,   0.70, 2.60, 3.48,  1, { 1,-1,-1,-1},  1, 7, "s", "[Rn]7s1"},
    { 88, "Ra", "Radium",      226.0,   0.90, 2.21, 2.83,  2, { 2,-1,-1,-1},  2, 7, "s", "[Rn]7s2"},

    // ---- Actinides (Z=89-103) ----
    { 89, "Ac", "Actinium",    227.0,   1.10, 2.15, 2.60,  3, { 3,-1,-1,-1},  3, 7, "f", "[Rn]6d1 7s2"},
    { 90, "Th", "Thorium",     232.04,  1.30, 2.06, 2.45,  4, { 4,-1,-1,-1},  0, 7, "f", "[Rn]6d2 7s2"},
    { 91, "Pa", "Protactinium",231.04,  1.50, 2.00, 2.43,  5, { 4, 5,-1,-1},  0, 7, "f", "[Rn]5f2 6d1 7s2"},
    { 92, "U",  "Uranium",     238.03,  1.38, 1.96, 1.86,  6, { 4, 6,-1,-1},  0, 7, "f", "[Rn]5f3 6d1 7s2"},
    { 93, "Np", "Neptunium",   237.0,   1.36, 1.90, 2.40,  7, { 4, 5, 6,-1},  0, 7, "f", "[Rn]5f4 6d1 7s2"},
    { 94, "Pu", "Plutonium",   244.0,   1.28, 1.87, 2.40,  7, { 3, 4, 5, 6},  0, 7, "f", "[Rn]5f6 7s2"},
    { 95, "Am", "Americium",   243.0,   1.13, 1.80, 2.40,  6, { 3, 4,-1,-1},  0, 7, "f", "[Rn]5f7 7s2"},
    { 96, "Cm", "Curium",      247.0,   1.28, 1.69, 2.40,  4, { 3, 4,-1,-1},  0, 7, "f", "[Rn]5f7 6d1 7s2"},
    { 97, "Bk", "Berkelium",   247.0,   1.30, 1.70, 2.40,  4, { 3, 4,-1,-1},  0, 7, "f", "[Rn]5f9 7s2"},
    { 98, "Cf", "Californium", 251.0,   1.30, 1.70, 2.40,  4, { 3, 4,-1,-1},  0, 7, "f", "[Rn]5f10 7s2"},
    { 99, "Es", "Einsteinium", 252.0,   1.30, 1.70, 2.40,  3, { 3,-1,-1,-1},  0, 7, "f", "[Rn]5f11 7s2"},
    {100, "Fm", "Fermium",     257.0,   1.30, 1.70, 2.40,  3, { 2, 3,-1,-1},  0, 7, "f", "[Rn]5f12 7s2"},
    {101, "Md", "Mendelevium", 258.0,   1.30, 1.70, 2.40,  3, { 2, 3,-1,-1},  0, 7, "f", "[Rn]5f13 7s2"},
    {102, "No", "Nobelium",    259.0,   1.30, 1.70, 2.40,  3, { 2, 3,-1,-1},  0, 7, "f", "[Rn]5f14 7s2"},
    {103, "Lr", "Lawrencium",  266.0,   1.30, 1.70, 2.40,  3, { 3,-1,-1,-1},  3, 7, "d", "[Rn]5f14 7s2 7p1"},

    // ---- Period 7 d-block (Z=104-112) ----
    // Superheavy elements: data from relativistic calculations and limited experiment
    {104, "Rf", "Rutherfordium",267.0,  0.00, 1.57, 2.40,  4, { 4,-1,-1,-1},  4, 7, "d", "[Rn]5f14 6d2 7s2"},
    {105, "Db", "Dubnium",     268.0,   0.00, 1.49, 2.40,  5, { 5,-1,-1,-1},  5, 7, "d", "[Rn]5f14 6d3 7s2"},
    {106, "Sg", "Seaborgium",  269.0,   0.00, 1.43, 2.40,  6, { 6,-1,-1,-1},  6, 7, "d", "[Rn]5f14 6d4 7s2"},
    {107, "Bh", "Bohrium",     270.0,   0.00, 1.41, 2.40,  7, { 7,-1,-1,-1},  7, 7, "d", "[Rn]5f14 6d5 7s2"},
    {108, "Hs", "Hassium",     277.0,   0.00, 1.34, 2.40,  8, { 8,-1,-1,-1},  8, 7, "d", "[Rn]5f14 6d6 7s2"},
    {109, "Mt", "Meitnerium",  278.0,   0.00, 1.29, 2.40,  6, { 3, 6,-1,-1},  9, 7, "d", "[Rn]5f14 6d7 7s2"},
    {110, "Ds", "Darmstadtium",281.0,   0.00, 1.28, 2.40,  6, { 6,-1,-1,-1}, 10, 7, "d", "[Rn]5f14 6d8 7s2"},
    {111, "Rg", "Roentgenium", 282.0,   0.00, 1.21, 2.40,  5, { 3, 5,-1,-1}, 11, 7, "d", "[Rn]5f14 6d9 7s2"},
    {112, "Cn", "Copernicium", 285.0,   0.00, 1.22, 2.40,  4, { 2, 4,-1,-1}, 12, 7, "d", "[Rn]5f14 6d10 7s2"},

    // ---- Period 7 p-block (Z=113-118) ----
    {113, "Nh", "Nihonium",    286.0,   0.00, 1.36, 2.40,  3, { 1, 3,-1,-1}, 13, 7, "p", "[Rn]5f14 6d10 7s2 7p1"},
    {114, "Fl", "Flerovium",   289.0,   0.00, 1.43, 2.40,  4, { 2, 4,-1,-1}, 14, 7, "p", "[Rn]5f14 6d10 7s2 7p2"},
    {115, "Mc", "Moscovium",   290.0,   0.00, 1.62, 2.40,  3, { 1, 3,-1,-1}, 15, 7, "p", "[Rn]5f14 6d10 7s2 7p3"},
    {116, "Lv", "Livermorium", 293.0,   0.00, 1.75, 2.40,  4, { 2, 4,-1,-1}, 16, 7, "p", "[Rn]5f14 6d10 7s2 7p4"},
    {117, "Ts", "Tennessine",  294.0,   0.00, 1.65, 2.40,  3, { 1, 3,-1,-1}, 17, 7, "p", "[Rn]5f14 6d10 7s2 7p5"},
    {118, "Og", "Oganesson",   294.0,   0.00, 1.57, 2.40,  2, { 2,-1,-1,-1}, 18, 7, "p", "[Rn]5f14 6d10 7s2 7p6"},
};

// ============================================================================
// Element property lookup by atomic number
// ============================================================================

inline const ElementData& getElement(int atomicNum) {
    if (atomicNum < 0 || atomicNum > 118) return PERIODIC_TABLE[0];
    return PERIODIC_TABLE[atomicNum];
}

inline const char* symbolFromAtomicNum(int z) {
    return getElement(z).symbol;
}

inline int atomicNumFromSymbol(const std::string& symbol) {
    // Static lookup table built on first call (includes aromatic lowercase forms)
    static const std::unordered_map<std::string, int> table = []() {
        std::unordered_map<std::string, int> m;
        m.reserve(256);
        for (int z = 0; z <= 118; ++z) {
            m[PERIODIC_TABLE[z].symbol] = z;
        }
        // Aromatic lowercase forms used in SMILES notation
        m["b"]  =  5;  m["c"]  =  6;  m["n"]  =  7;  m["o"]  =  8;
        m["p"]  = 15;  m["s"]  = 16;  m["se"] = 34;  m["te"] = 52;
        // Deuterium / Tritium aliases (both map to hydrogen, Z=1)
        m["D"] = 1;  m["T"] = 1;
        return m;
    }();
    auto it = table.find(symbol);
    return (it != table.end()) ? it->second : -1;
}

// ============================================================================
// Physical properties
// ============================================================================

inline double paulingElectronegativity(int z) {
    return getElement(z).electronegativity;
}

inline double covalentRadiusSingle(int z) {
    return getElement(z).covalentRadius;
}

inline double vanDerWaalsRadius(int z) {
    return getElement(z).vanDerWaalsRadius;
}

inline double atomicWeight(int z) {
    return getElement(z).atomicWeight;
}

// ============================================================================
// Valence and bonding
// ============================================================================

inline int maxValence(int z) {
    return getElement(z).maxValence;
}

/// Returns pointer to the defaultValences array (terminated by -1).
inline const int* allowedValences(int z) {
    return getElement(z).defaultValences;
}

/// Compute number of implicit hydrogens using the OpenSMILES valence model.
///
/// For organic-subset atoms at charge 0:
///   implicitH = smallestFittingValence - bondOrderSum
///
/// For charged atoms, the target valence is adjusted:
///   - Cations of N, P, O, S, As, Se: target = v + charge  (valence expansion)
///   - Anions of B, Al: target = v + |charge|  (isoelectronic with Group 14)
///   - All others: target = v - |charge|
inline int implicitHydrogens(int z, int bondOrderSum, int charge) {
    const int* vals = allowedValences(z);
    if (vals[0] == -1) return 0;  // no default valences known

    for (int i = 0; i < 4 && vals[i] != -1; ++i) {
        int v = vals[i];
        int target = v;

        if (charge != 0) {
            // Cations of Group 15/16 elements expand into next valence shell
            if (charge > 0 && (z == 7 || z == 15 || z == 8 || z == 16 ||
                               z == 33 || z == 34)) {
                target = v + charge;
            }
            // Group 13 anions gain valence (isoelectronic with Group 14)
            else if (charge < 0 && (z == 5 || z == 13)) {
                target = v + (-charge);
            }
            // General case: charge reduces available valence
            else {
                target = v - std::abs(charge);
            }
        }

        if (target >= bondOrderSum) {
            int h = target - bondOrderSum;
            return (h >= 0) ? h : 0;
        }
    }
    return 0;
}

// ============================================================================
// Orbital configuration
// ============================================================================

inline const char* electronConfiguration(int z) {
    return getElement(z).electronConfig;
}

inline const char* orbitalBlock(int z) {
    return getElement(z).block;
}

/// Return the number of valence electrons for element Z.
///
/// Valence shell assignment:
///   s-block: group number (1 or 2)
///   p-block: group - 10
///   d-block: varies by filling anomalies; typically group number for groups 3-12
///   f-block (lanthanides/actinides): 3 for most, special cases handled
inline int valenceElectrons(int z) {
    if (z <= 0 || z > 118) return 0;

    const auto& e = PERIODIC_TABLE[z];

    // Noble gases: full octet
    if (e.group == 18) {
        return (z == 2) ? 2 : 8;
    }

    // s-block
    if (e.group == 1) return 1;
    if (e.group == 2) return 2;

    // p-block (groups 13-17)
    if (e.group >= 13 && e.group <= 17) {
        return e.group - 10;
    }

    // d-block (groups 3-12): valence electrons = s + d electrons
    // Approximate: for VSEPR purposes, use the group number
    // but cap at the number of electrons beyond the last noble gas core
    if (e.group >= 3 && e.group <= 12) {
        // Standard d-block: valence = group number for groups 3-7,
        // groups 8-10 have varying oxidation states,
        // groups 11-12: 1 and 2 s-electrons dominate chemistry
        switch (e.group) {
            case  3: return 3;
            case  4: return 4;
            case  5: return 5;
            case  6: return 6;
            case  7: return 7;
            case  8: return 8;
            case  9: return 9;
            case 10: return 10;
            case 11: return 1;   // coinage metals: Cu, Ag, Au (d10 s1)
            case 12: return 2;   // Zn, Cd, Hg (d10 s2 — d electrons inert)
            default: return 2;
        }
    }

    // f-block (lanthanides and actinides, group == 0)
    // Most lanthanides: [Xe]4f^n 5d^0-1 6s^2 => ~3 valence electrons
    // Most actinides: [Rn]5f^n 6d^0-1 7s^2 => varies, approximate as 3
    if (e.group == 0) {
        // Cerium and a few others have 4f participation
        if (z == 58 || z == 59 || z == 65) return 4;  // Ce, Pr, Tb (known +4 states)
        if (z >= 57 && z <= 71) return 3;  // lanthanides: predominantly +3
        if (z >= 89 && z <= 103) {
            // Actinides show wider variation
            if (z == 92) return 6;  // U: prominent +6 state
            if (z == 93) return 5;  // Np
            if (z == 94) return 6;  // Pu
            return 3;               // default for actinides
        }
    }

    return 0;
}

/// Return the number of core (non-valence) electrons for element Z.
/// Core electrons = total electrons - valence electrons.
inline int coreElectrons(int z) {
    if (z <= 0 || z > 118) return 0;
    return z - valenceElectrons(z);
}

// ============================================================================
// Hybridization
// ============================================================================

enum class Hybridization {
    S,        // atomic s orbital (hydrogen)
    SP,       // linear (e.g., CO2 carbon, acetylene)
    SP2,      // trigonal planar (e.g., ethylene, BF3)
    SP3,      // tetrahedral (e.g., methane, water)
    SP3D,     // trigonal bipyramidal (e.g., PF5)
    SP3D2,    // octahedral (e.g., SF6)
    UNKNOWN
};

/// Compute hybridization from atom connectivity.
///
/// Uses the steric number model (VSEPR):
///   effectiveElectrons = valenceElectrons(Z) - charge
///   lonePairs = (effectiveElectrons - bondOrderSum) / 2
///   stericNumber = totalDegree + lonePairs
///
/// Aromatic atoms are assigned SP2 (pi system participation).
inline Hybridization computeHybridization(int atomicNum, int totalDegree,
                                          int bondOrderSum, int charge,
                                          bool isAromatic) {
    if (isAromatic) return Hybridization::SP2;

    int ve = valenceElectrons(atomicNum) - charge;
    if (ve < 0) ve = 0;

    int lp = (ve - bondOrderSum) / 2;
    if (lp < 0) lp = 0;

    int steric = totalDegree + lp;

    // Hydrogen and single-coordinate atoms with no lone pairs
    if (atomicNum == 1) return Hybridization::S;

    switch (steric) {
        case 1:  return Hybridization::SP;    // unusual, e.g., [CH]+ fragment
        case 2:  return Hybridization::SP;    // linear
        case 3:  return Hybridization::SP2;   // trigonal planar
        case 4:  return Hybridization::SP3;   // tetrahedral
        case 5:  return Hybridization::SP3D;  // trigonal bipyramidal
        case 6:  return Hybridization::SP3D2; // octahedral
        default: return Hybridization::UNKNOWN;
    }
}

// ============================================================================
// Classification functions
// ============================================================================

/// True if the element is a metal (alkali, alkaline earth, transition metal,
/// post-transition metal, lanthanide, or actinide).
inline bool isMetal(int z) {
    if (z <= 0 || z > 118) return false;
    // Non-metals: H, He, C, N, O, F, Ne, P, S, Cl, Ar, Se, Br, Kr, I, Xe, At, Rn,
    //             noble gases, and halogens are handled; metalloids are debatable.
    // Metalloids (B, Si, Ge, As, Sb, Te) are NOT classified as metals here.
    // Po (Z=84) is a post-transition metal per IUPAC and IS classified as metal.
    static constexpr bool nonmetal[119] = {
        // 0   1    2    3    4    5    6    7    8    9
        true,true,true,false,false,true,true,true,true,true,  //  0- 9
        true,false,false,false,true,true,true,true,true,false, // 10-19
        false,false,false,false,false,false,false,false,false,false, // 20-29
        false,false,true,true,true,true,true,false,false,false, // 30-39
        false,false,false,false,false,false,false,false,false,false, // 40-49
        false,true,true,true,true,false,false,false,false,false, // 50-59
        false,false,false,false,false,false,false,false,false,false, // 60-69
        false,false,false,false,false,false,false,false,false,false, // 70-79
        false,false,false,false,false,true,true,false,false,false, // 80-89
        false,false,false,false,false,false,false,false,false,false, // 90-99
        false,false,false,false,false,false,false,false,false,false, //100-109
        false,false,false,false,false,false,false,true,true         //110-118
    };
    return !nonmetal[z];
}

/// True if the element is a transition metal (groups 3-12, periods 4-7).
inline bool isTransitionMetal(int z) {
    if (z <= 0 || z > 118) return false;
    const auto& e = PERIODIC_TABLE[z];
    return (e.group >= 3 && e.group <= 12 && e.period >= 4);
}

/// True for halogens: F, Cl, Br, I, At, Ts.
inline bool isHalogen(int z) {
    return z == 9 || z == 17 || z == 35 || z == 53 || z == 85 || z == 117;
}

/// True for noble gases: He, Ne, Ar, Kr, Xe, Rn, Og.
inline bool isNobleGas(int z) {
    return z == 2 || z == 10 || z == 18 || z == 36 || z == 54 || z == 86 || z == 118;
}

/// True for main-group elements (s-block and p-block).
inline bool isMainGroup(int z) {
    if (z <= 0 || z > 118) return false;
    const auto& e = PERIODIC_TABLE[z];
    return (e.block[0] == 's' || e.block[0] == 'p');
}

/// True for lanthanides (Z = 57-71).
inline bool isLanthanide(int z) {
    return z >= 57 && z <= 71;
}

/// True for actinides (Z = 89-103).
inline bool isActinide(int z) {
    return z >= 89 && z <= 103;
}

/// True for the organic subset: B, C, N, O, P, S, F, Cl, Br, I.
/// These atoms can appear without brackets in SMILES notation.
inline bool isOrganicSubset(int z) {
    return z == 5 || z == 6 || z == 7 || z == 8 || z == 15 || z == 16 ||
           z == 9 || z == 17 || z == 35 || z == 53;
}

} // namespace smsd

#endif // SMSD_PERIODIC_TABLE_HPP
