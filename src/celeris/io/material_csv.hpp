#pragma once
// Load a material's dispersion from a CSV of refractive-index data, e.g. the
// exports from refractiveindex.info. Each data row is:
//     wavelength_um , n [, k]
// Header lines and anything not starting with a number are skipped; k defaults
// to 0 (lossless) if absent. Wavelengths must be ascending. This is how real
// materials (TiO2, Si, GaN, ...) enter CELERIS instead of a constant index.

#include "celeris/materials/material.hpp"

#include <string>

namespace celeris {

// Throws std::runtime_error on a missing/unreadable file or fewer than two
// valid data rows.
Material load_material_csv(const std::string& path, const std::string& name);

} // namespace celeris
