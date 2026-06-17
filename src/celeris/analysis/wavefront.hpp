#pragma once
// Wavefront (OPD) analysis — the core of real optical-design software. The
// metalens imparts a realized phase that deviates from the ideal focusing
// phase; that deviation IS the wavefront error. From it we get the metrics an
// optical engineer actually reads:
//   - OPD map (waves) across the aperture
//   - RMS and peak-to-valley wavefront error (piston removed)
//   - Strehl via the Maréchal approximation  exp(-(2π·RMS)²)  (independent
//     cross-check of the diffraction PSF Strehl)
//   - Zernike decomposition: how much defocus / astigmatism / coma / spherical
//     / ... the residual contains (named aberrations)

#include "celeris/design/metalens.hpp"

#include <string>
#include <vector>

namespace celeris {

struct ZernikeTerm {
    int noll;            // Noll index
    std::string name;    // e.g. "defocus", "coma (x)"
    double coeff_waves;  // fitted coefficient, waves RMS
};

struct WavefrontAnalysis {
    double rms_waves;        // RMS wavefront error (piston removed)
    double pv_waves;         // peak-to-valley
    double strehl_marechal;  // exp(-(2π·RMS)²)
    std::vector<ZernikeTerm> zernike;
    int n;                   // OPD map is n×n at lens-cell resolution
    std::vector<double> opd; // waves; entries outside the aperture set to 0
    std::vector<char> mask;  // 1 inside aperture, 0 outside (same n×n)
};

// Compute the wavefront error of `lens` relative to the ideal focusing phase.
WavefrontAnalysis analyze_wavefront(const MetalensDesign& lens,
                                    const UnitCellLibrary& lib,
                                    double focal_length_um, double wavelength_um,
                                    double diameter_um);

} // namespace celeris
