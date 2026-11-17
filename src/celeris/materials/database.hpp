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

/// @brief Built-in materials and the named material registry.
namespace celeris::materials {

// --- Analytic models (lossless Sellmeier / constant) -----------------------
/// @brief Air, n = 1 exactly.
const Material& air();
/// @brief Schott N-BK7 borosilicate crown glass (Sellmeier).
const Material& bk7();
/// @brief Fused silica SiO₂ (Malitson 1965 Sellmeier).
const Material& fused_silica();
/// @brief Si₃N₄, LPCVD (Luke et al. 2015 Sellmeier).
const Material& silicon_nitride();
/// @brief Sapphire Al₂O₃, ordinary ray (Malitson–Dodge 1972 Sellmeier).
const Material& sapphire();
/// @brief α-GaN, ordinary ray (Barker–Ilegems 1973).
const Material& gallium_nitride();

// --- Named registry --------------------------------------------------------
// Metadata for one catalog entry (for `celeris materials` listing + the GUI).
/// @brief Metadata for one material registry entry (for `celeris materials` and the GUI).
struct MaterialInfo {
    std::string name;         ///< Canonical lookup key, e.g. "au", "c-si", "sio2".
    std::string description;  ///< Human label + literature source.
    bool tabulated;           ///< True = real n,k from a data file (may be lossy).
    double lambda_min_um;     ///< Validity range start (µm); 0 with lambda_max_um 0 = broadband analytic.
    double lambda_max_um;     ///< Validity range end (µm).
    std::string data_file;    ///< Basename of the CSV (tabulated only), else "".
    bool available;           ///< Tabulated: the data file was located on disk.
};

/// @brief Every material the registry knows, with metadata.
///
/// Aliases are not listed (see `by_name`); this is the canonical, de-duplicated set.
std::vector<MaterialInfo> catalog();

/// @brief Resolve a material by canonical name or alias.
///
/// Case-insensitive; '-' and '_' are equivalent. Analytic models are returned
/// by value; tabulated ones are loaded from data/ on each call (see
/// `find_data_file`).
/// @param name Canonical name or alias, e.g. "tio2", "a-si", "gan".
/// @throws std::runtime_error on an unknown name or a missing tabulated data file.
Material by_name(const std::string& name);

/// @brief Locate a data file shipped in data/ by basename.
///
/// Searches CELERIS_DATA_DIR, ./data, ../data, ../../data and a few more roots.
/// @param basename File name, e.g. "TiO2_Siefke.csv".
/// @return The first existing path, or "" if not found.
std::string find_data_file(const std::string& basename);

} // namespace celeris::materials
