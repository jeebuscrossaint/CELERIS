#pragma once
/// @file
/// @brief Tabulated materials from a refractive-index CSV file.
///
/// Load a material's dispersion from a CSV of refractive-index data, e.g. the
/// exports from refractiveindex.info. Each data row is:
///     wavelength_um , n [, k]
/// Header lines and anything not starting with a number are skipped; k defaults
/// to 0 (lossless) if absent. Wavelengths must be ascending. This is how real
/// materials (TiO2, Si, GaN, ...) enter CELERIS instead of a constant index.

#include "celeris/materials/material.hpp"

#include <string>

namespace celeris {

/// @brief Load a tabulated material from a CSV of wavelength (µm), n and optional k.
///
/// Rows that do not start with a number are skipped; k defaults to 0. The
/// wavelengths must be ascending.
/// @param path CSV file, e.g. an export from refractiveindex.info.
/// @param name Label for the material.
/// @return The material, interpolated between the tabulated points.
/// @throws std::runtime_error If the file is missing or unreadable, or has
///         fewer than two valid data rows.
Material load_material_csv(const std::string& path, const std::string& name);

} // namespace celeris
