#pragma once
/// @file
/// @brief System-level design optimizer: search period × pillar height by a merit function.
///
/// The "merit function" idea from optical design tools, applied to a metalens.
/// Instead of the user guessing the period and pillar height, this searches
/// that design space and returns the combination that maximizes a merit (focal
/// Strehl, optionally weighted by transmission). Each trial builds a full
/// unit-cell library, assembles the lens, and measures the actual focal
/// performance — so it optimizes the real metric, not a proxy.

#include "celeris/materials/material.hpp"

#include <functional>

namespace celeris {

/// @brief The best design found by optimize_system().
struct SystemOptResult {
    double period_um;       ///< Winning lattice period (µm).
    double thickness_um;    ///< Winning pillar height (µm).
    double strehl;          ///< Focal Strehl ratio of the winning design.
    double mean_amplitude;  ///< Mean |t| of the winning design's pillars.
    double merit;           ///< strehl + efficiency_weight · mean transmission.
};

/// @brief Grid-search period × pillar height and return the best design by merit.
///
/// merit = strehl + efficiency_weight · mean transmission. Uses a coarse library
/// (small M, few samples) for speed; re-run the winning parameters at full
/// quality.
/// @param pillar,background  Pillar and surrounding materials.
/// @param incident,substrate Half-spaces above and below.
/// @param focal_um           Focal length (µm).
/// @param diameter_um        Aperture width (µm).
/// @param wavelength_um      Design wavelength (µm).
/// @param period_lo,period_hi Period range (µm).
/// @param thick_lo,thick_hi  Height range (µm).
/// @param grid               Grid points per axis.
/// @param M                  Harmonic half-count for the library solves.
/// @param fill_samples       Fill samples per library.
/// @param efficiency_weight  Weight on mean transmission in the merit.
/// @param progress           Called with the completed fraction, 0..1 (may be empty).
SystemOptResult optimize_system(const Material& pillar, const Material& background,
                                const Material& incident, const Material& substrate,
                                double focal_um, double diameter_um,
                                double wavelength_um, double period_lo,
                                double period_hi, double thick_lo, double thick_hi,
                                int grid, int M, int fill_samples,
                                double efficiency_weight,
                                const std::function<void(float)>& progress);

} // namespace celeris
