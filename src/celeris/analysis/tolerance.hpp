#pragma once
/// @file
/// @brief Fabrication tolerance: Monte-Carlo Strehl under pillar-size error.
///
/// How does the focus survive real-world manufacturing error? Lithography + etch
/// leave each pillar a few nanometers off its design size. This Monte-Carlo
/// perturbs every pillar's width by a Gaussian of a given sigma and measures the
/// resulting Strehl, so you can answer "what CD (critical-dimension) control does
/// this design need to hit spec?" — a hard requirement for defense/aerospace
/// procurement.

#include "celeris/design/metalens.hpp"

#include <vector>

namespace celeris {

/// @brief Strehl statistics at one fabrication-error level.
struct ToleranceResult {
    double sigma_nm;      ///< 1-sigma fabrication error on pillar size (nm).
    double mean_strehl;   ///< Mean Strehl over the Monte-Carlo trials.
    double std_strehl;    ///< Standard deviation of Strehl.
    double worst_strehl;  ///< Worst-case trial.
};

/// @brief Monte-Carlo Strehl statistics for each CD-error level.
///
/// For each sigma, runs `trials` realizations: every pillar's fill gets an
/// independent Gaussian error of sigma_nm / period (clamped to 0.05..0.95), and
/// the transmission-weighted Strehl from analyze_focus() is recorded. sigma_nm = 0
/// reproduces the nominal design. Trials run in parallel.
/// @param lens            The finished design.
/// @param lib             The library the design was built from.
/// @param focal_length_um Focal length (µm).
/// @param wavelength_um   Vacuum wavelength (µm).
/// @param diameter_um     Circular aperture diameter (µm).
/// @param sigmas_nm       1-sigma CD errors to evaluate (nm).
/// @param trials          Realizations per sigma (at least 1 is run).
/// @param seed            RNG seed; the same seed gives the same result.
/// @return One ToleranceResult per sigma, in input order.
std::vector<ToleranceResult> analyze_tolerance(
    const MetalensDesign& lens, const UnitCellLibrary& lib, double focal_length_um,
    double wavelength_um, double diameter_um, const std::vector<double>& sigmas_nm,
    int trials, unsigned seed);

} // namespace celeris
