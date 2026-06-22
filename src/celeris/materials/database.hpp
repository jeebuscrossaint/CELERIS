#pragma once
// Built-in material library + a named registry. The analytic models below
// return references to long-lived Material objects (constructed once on first
// use), so callers can hold the reference safely.
//
// Real metalens materials WITH absorption (TiO2, Si, Au, Ag, Al) are NOT
// hardcoded as fabricated numbers — they are loaded from tabulated n,k data
// shipped in data/ (sourced from refractiveindex.info, CC0). The transparent
// dielectrics (glass, sapphire, GaN, Si3N4) use canonical published Sellmeier
// models. Everything is reachable by a short canonical name via by_name().

#include "celeris/materials/material.hpp"

#include <string>
#include <vector>

namespace celeris::materials {

// --- Analytic models (lossless Sellmeier / constant) -----------------------
const Material& air();             // n = 1 exactly
const Material& bk7();             // Schott N-BK7 borosilicate crown glass (Sellmeier)
const Material& fused_silica();    // Fused silica SiO2 (Malitson 1965 Sellmeier)
const Material& silicon_nitride(); // Si3N4 LPCVD (Luke et al. 2015 Sellmeier)
const Material& sapphire();        // Al2O3 ordinary ray (Malitson-Dodge 1972 Sellmeier)
const Material& gallium_nitride(); // alpha-GaN ordinary ray (Barker-Ilegems 1973)

// --- Named registry --------------------------------------------------------
// Metadata for one catalog entry (for `celeris materials` listing + the GUI).
struct MaterialInfo {
    std::string name;         // canonical lookup key, e.g. "au", "c-si", "sio2"
    std::string description;  // human label + literature source
    bool tabulated;           // true = real n,k from a data file (may be lossy)
    double lambda_min_um;     // validity range; 0,0 = broadband analytic
    double lambda_max_um;
    std::string data_file;    // basename of the CSV (tabulated only), else ""
    bool available;           // tabulated: the data file was located on disk
};

// Every material the registry knows about, with metadata. Aliases are NOT
// listed (see by_name) — this is the canonical, de-duplicated set.
std::vector<MaterialInfo> catalog();

// Resolve a material by canonical name or alias (case-insensitive; '-'/'_'
// equivalent). Analytic models are returned by value (Material is value-
// semantic); tabulated ones are loaded from data/ on each call. Throws
// std::runtime_error on an unknown name or a missing tabulated data file.
Material by_name(const std::string& name);

// Try to locate a data file shipped in data/ by basename, searching a handful
// of candidate roots (CELERIS_DATA_DIR, ./data, ../data, ../../data, ...).
// Returns the first existing path, or "" if not found.
std::string find_data_file(const std::string& basename);

} // namespace celeris::materials
