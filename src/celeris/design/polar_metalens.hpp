#pragma once
/// @file
/// @brief Polarization-multiplexed metalens design.
///
/// A rectangular nanopillar imparts independent phase to x- and y-polarized
/// light (form birefringence), so one metasurface can implement TWO different
/// phase profiles at once -- e.g. focus x-polarized light to one point and
/// y-polarized light to another. This is the basis of polarization-splitting /
/// polarization-imaging metalenses.
///
/// We build a 2D library over (fill_x, fill_y) -- two RCWA solves per cell, one
/// per incident polarization -- then for each lattice site pick the pillar whose
/// (phase_x, phase_y) best matches the two target profiles simultaneously.

#include "celeris/core.hpp"
#include "celeris/materials/material.hpp"

#include <vector>

namespace celeris {

/// @brief One sampled rectangular meta-atom and its dual-polarization response.
struct PolarCell {
    double fill_x;   ///< Pillar width in x (fraction of the period).
    double fill_y;   ///< Pillar width in y (fraction of the period).
    double phase_x;  ///< Phase imparted to x-polarized light (rad).
    double phase_y;  ///< Phase imparted to y-polarized light (rad).
    double amp_x;    ///< |t| for x-polarized light.
    double amp_y;    ///< |t| for y-polarized light.
};

/// @brief A (fill_x, fill_y) library of rectangular pillars.
struct PolarizationLibrary {
    double period_um = 0;      ///< Lattice period (µm).
    double wavelength_um = 0;  ///< Design wavelength (µm).
    double thickness_um = 0;   ///< Pillar height (µm).
    std::vector<PolarCell> cells;  ///< n_samples² entries.

    /// @brief Index of the cell best matching both target phases.
    ///
    /// Minimizes the combined phase error to (target_x, target_y), with an
    /// optional amplitude penalty that favors high transmission.
    /// @param target_x,target_y  Target phases for x- and y-polarization (rad).
    /// @param amplitude_weight   Transmission penalty weight.
    int lookup(double target_x, double target_y, double amplitude_weight = 0.0) const;
};

/// @brief Build the (fill_x, fill_y) library.
///
/// n_samples × n_samples rectangular pillars, each solved twice (x- and
/// y-polarized). Solves run in parallel.
/// @param pillar,background  Pillar and surrounding materials.
/// @param incident,substrate Half-spaces above and below.
/// @param period_um          Lattice period (µm).
/// @param wavelength_um      Design wavelength (µm).
/// @param thickness_um       Pillar height (µm).
/// @param fill_min,fill_max  Fill range, applied to both axes.
/// @param n_samples          Samples per axis.
/// @param M                  Harmonic half-count in x and y.
PolarizationLibrary build_polarization_library(
    const Material& pillar, const Material& background, const Material& incident,
    const Material& substrate, double period_um, double wavelength_um,
    double thickness_um, double fill_min, double fill_max, int n_samples, int M);

/// @brief A polarization-multiplexed metalens: one rectangular pillar per site.
struct PolarMetalensDesign {
    int n_cells = 0;             ///< Cells across the (square) aperture.
    double period_um = 0;        ///< Lattice period (µm).
    std::vector<double> fill_x;  ///< n_cells² pillar widths in x, row-major.
    std::vector<double> fill_y;  ///< n_cells² pillar widths in y, row-major.
    std::vector<cdouble> t_x;    ///< Realized per-pillar transmission, x-polarized.
    std::vector<cdouble> t_y;    ///< Realized per-pillar transmission, y-polarized.
    double rms_phase_error_x_deg = 0;  ///< RMS phase error, x-polarized (deg).
    double rms_phase_error_y_deg = 0;  ///< RMS phase error, y-polarized (deg).
    double mean_amp_x = 0;       ///< Mean |t|, x-polarized.
    double mean_amp_y = 0;       ///< Mean |t|, y-polarized.
};

/// @brief Design a lens that focuses x- and y-polarized light at different distances.
///
/// Maps the two ideal hyperbolic phase profiles onto rectangular pillars. Set
/// the focal lengths equal for a polarization-insensitive lens.
/// @param lib               The (fill_x, fill_y) library.
/// @param focal_x_um        Focal length for x-polarized light (µm).
/// @param focal_y_um        Focal length for y-polarized light (µm).
/// @param diameter_um       Aperture width (µm).
/// @param amplitude_weight  Transmission penalty weight.
PolarMetalensDesign design_polarization_metalens(const PolarizationLibrary& lib,
                                                 double focal_x_um,
                                                 double focal_y_um,
                                                 double diameter_um,
                                                 double amplitude_weight = 0.25);

} // namespace celeris
