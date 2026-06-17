#pragma once
// Focal-plane analysis — does the designed lens actually focus, and how well?
//
// Takes a finished metalens design, builds its complex aperture field (each
// pillar's transmission from the library), propagates it to the focal plane by
// Rayleigh-Sommerfeld summation, and reports the standard optical metrics:
//   - Strehl ratio: peak intensity vs an ideal (perfect-phase) lens (1.0 = perfect)
//   - focal-spot FWHM, compared to the diffraction limit ~ lambda*f/D
//   - encircled energy within the first Airy null (focusing quality)
// This is the end-to-end proof that the whole pipeline produces a working lens.

#include "celeris/design/metalens.hpp"

namespace celeris {

struct FocalAnalysis {
    double strehl;              // peak(design) / peak(ideal phase), 0..1
    double fwhm_um;            // focal-spot full width at half maximum
    double diffraction_limit_um;  // lambda * f / D for reference
    double encircled_energy;   // fraction of focal-window power within first null
};

// Analyze the focus of `lens` at the given focal length / wavelength / aperture
// diameter. The aperture is masked to a circle of the given diameter.
FocalAnalysis analyze_focus(const MetalensDesign& lens,
                            const UnitCellLibrary& lib, double focal_length_um,
                            double wavelength_um, double diameter_um);

// Full 2D focal-plane intensity map (the point-spread function), sampled on an
// n x n grid spanning +/- half_window_um about the axis at the focal plane.
struct PsfMap {
    int n;
    double half_window_um;
    std::vector<double> intensity;  // row-major, n*n, |E|^2
};
PsfMap compute_psf(const MetalensDesign& lens, const UnitCellLibrary& lib,
                   double focal_length_um, double wavelength_um,
                   double diameter_um, int n, double half_window_um);

} // namespace celeris
