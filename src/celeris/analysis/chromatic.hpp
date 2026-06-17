#pragma once
// Chromatic analysis — how the focus moves with wavelength. A standard
// phase-profile metalens is strongly chromatic: its focal length scales as
// f(λ) ≈ f₀·λ₀/λ. This is the #1 limitation of metalenses (and why achromatic
// designs are a hot research area), so quantifying it is essential.
//
// Model: the pillars' imparted phase is held fixed at the design value
// (geometric-phase approximation) and only the propagation wavenumber k(λ)
// varies. This captures the dominant chromatic focal shift. A fully rigorous
// treatment would rebuild the RCWA library at each wavelength (heavier; future).

#include "celeris/design/metalens.hpp"

#include <vector>

namespace celeris {

struct ChromaticPoint {
    double wavelength_um;
    double focal_length_um;   // where on-axis intensity actually peaks
    double rel_peak;          // peak intensity relative to the design wavelength
};

// Sweep wavelength from lambda_min..lambda_max (n points) and locate the focus
// at each by scanning on-axis intensity. design_focal/design_wavelength are the
// nominal design point; diameter sets the circular aperture.
std::vector<ChromaticPoint> analyze_chromatic(const MetalensDesign& lens,
                                              const UnitCellLibrary& lib,
                                              double design_focal,
                                              double design_wavelength,
                                              double diameter, double lambda_min,
                                              double lambda_max, int n);

} // namespace celeris
