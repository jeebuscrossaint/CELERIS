#pragma once
/// @file
/// @brief Wavefront (OPD) analysis: realized vs ideal focusing phase, RMS/PV, Zernikes.
///
/// The core of real optical-design software. The metalens imparts a realized
/// phase that deviates from the ideal focusing phase; that deviation IS the
/// wavefront error. From it we get the metrics an optical engineer actually reads:
///   - OPD map (waves) across the aperture
///   - RMS and peak-to-valley wavefront error (piston removed)
///   - Strehl via the Maréchal approximation  exp(-(2π·RMS)²)  (independent
///     cross-check of the diffraction PSF Strehl)
///   - Zernike decomposition: how much defocus / astigmatism / coma / spherical
///     / ... the residual contains (named aberrations)

#include "celeris/design/metalens.hpp"

#include <string>
#include <vector>

namespace celeris {

/// @brief One fitted Zernike aberration term.
struct ZernikeTerm {
    int noll;            ///< Noll index (2..11: tip/tilt through primary spherical).
    std::string name;    ///< Aberration name, e.g. "defocus", "coma (x)".
    double coeff_waves;  ///< Fitted coefficient (waves RMS).
};

/// @brief Wavefront error of a lens and its decomposition.
struct WavefrontAnalysis {
    double rms_waves;        ///< RMS wavefront error, piston removed (waves).
    double pv_waves;         ///< Peak-to-valley wavefront error (waves).
    double strehl_marechal;  ///< Maréchal Strehl estimate exp(-(2π·RMS)²).
    std::vector<ZernikeTerm> zernike;  ///< Least-squares Zernike fit, Noll 2..11.
    int n;                   ///< OPD map is n×n at lens-cell resolution.
    std::vector<double> opd; ///< Row-major OPD in waves, wrapped to (-0.5, 0.5]; 0 outside the aperture.
    std::vector<char> mask;  ///< 1 inside the aperture, 0 outside (same n×n).
};

/// @brief Compute the wavefront error of `lens` relative to the ideal focusing phase.
///
/// The ideal phase is the hyperbolic -k(sqrt(r²+f²) - f). Per cell, the OPD is
/// the realized library phase minus the ideal phase, wrapped to one wave.
/// @param lens            The finished design.
/// @param lib             The library the design was built from.
/// @param focal_length_um Focal length of the ideal reference (µm).
/// @param wavelength_um   Vacuum wavelength (µm).
/// @param diameter_um     Circular aperture diameter (µm).
/// @return The analysis; an empty aperture gives zero error and Strehl 1.
WavefrontAnalysis analyze_wavefront(const MetalensDesign& lens,
                                    const UnitCellLibrary& lib,
                                    double focal_length_um, double wavelength_um,
                                    double diameter_um);

} // namespace celeris
