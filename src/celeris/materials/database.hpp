#pragma once
// Built-in material library. Returns references to long-lived Material objects
// (constructed once on first use), so callers can hold the reference safely.
//
// These three are exact, well-published dispersion models. Real metalens
// materials with absorption (TiO2, Si, GaN) should be loaded from tabulated
// n,k data (refractiveindex.info) via Material::tabulated — we do NOT hardcode
// fabricated numbers for those.

#include "celeris/materials/material.hpp"

namespace celeris::materials {

const Material& air();           // n = 1 exactly
const Material& bk7();           // Schott N-BK7 borosilicate crown glass (Sellmeier)
const Material& fused_silica();  // Fused silica SiO2 (Malitson 1965 Sellmeier)

} // namespace celeris::materials
