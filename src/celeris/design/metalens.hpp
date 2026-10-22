#pragma once
/// @file
/// @brief Metalens design: the layer that turns the RCWA solver into a usable tool.
///
/// Two pieces:
///   1. UnitCellLibrary — sweep a nanopillar's size and tabulate the complex
///      transmission (phase + amplitude) it produces. This is the "phase
///      library": the palette of phase shifts available from one fab process.
///   2. design_metalens() — given a target phase profile and aperture, pick at
///      each lattice site the pillar whose library phase best matches. The
///      output is a fabricable pillar map plus the residual phase error
///      (design fidelity).

#include "celeris/design/phase_profile.hpp"
#include "celeris/materials/material.hpp"
#include "celeris/rcwa/rcwa2d.hpp"

#include <vector>

namespace celeris {

/// @brief A meta-atom phase library: transmission vs pillar fill at one height.
struct UnitCellLibrary {
    double period_um;               ///< Lattice period (µm).
    double wavelength_um;           ///< Design wavelength (µm).
    double thickness_um;            ///< Pillar height (µm).
    std::vector<double> fill;       ///< Sampled pillar size fractions.
    std::vector<double> phase;      ///< Transmission phase per sample (rad), in (-π, π].
    std::vector<double> amplitude;  ///< Zeroth-order |t| per sample.

    /// @brief Phase span (max - min of the wrapped phases), in radians.
    ///
    /// Only meaningful when the phase-vs-fill curve does not wrap across ±π.
    /// For the honest "can I hit any target phase" measure use coverage().
    double phase_span() const;
    /// @brief Effective phase coverage on the circle (rad).
    ///
    /// 2π minus the largest angular gap between consecutive available phases.
    /// This is the metric that governs design fidelity: a small largest gap
    /// means any target phase has a near match. Equals phase_span() for a
    /// monotonic non-wrapping sweep, but stays correct when the phase wraps.
    double coverage() const;
    /// @brief Index of the entry whose phase best matches the target (circularly).
    /// @param target_phase_rad  Target phase (rad).
    int lookup(double target_phase_rad) const;
    /// @brief Like lookup(), but trades phase error against transmission.
    ///
    /// Minimizes phase_error² + amplitude_weight·(1 - |t|)². With weight 0 this
    /// is lookup().
    /// @param target_phase_rad  Target phase (rad).
    /// @param amplitude_weight  Weight on the transmission penalty.
    int lookup_weighted(double target_phase_rad, double amplitude_weight) const;
    /// @brief Complex transmission t = |t|·e^{iφ} of the sampled pillar nearest `fill`.
    cdouble transmission_for_fill(double fill) const;
};

/// @brief Build a library by sweeping pillar fill from fill_min to fill_max.
///
/// Each sample is a full 2D RCWA solve at normal incidence, x-polarized.
/// @param pillar,background  Pillar and surrounding materials.
/// @param incident,substrate Half-spaces above and below the pillar layer.
/// @param period_um          Lattice period (µm).
/// @param wavelength_um      Design wavelength (µm).
/// @param thickness_um       Pillar height (µm).
/// @param fill_min,fill_max  Fill range to sweep (fractions of the period).
/// @param n_samples          Number of fill samples.
/// @param M                  Harmonic half-count in x and y (Mx = My = M).
/// @param shape              Meta-atom cross-section.
/// @param shape_param        Shape parameter (see RectCell2D::shape_param).
UnitCellLibrary build_unit_cell_library(const Material& pillar,
                                        const Material& background,
                                        const Material& incident,
                                        const Material& substrate,
                                        double period_um, double wavelength_um,
                                        double thickness_um, double fill_min,
                                        double fill_max, int n_samples, int M,
                                        MetaShape shape = MetaShape::Rectangle,
                                        double shape_param = 0.5);

/// @brief Build a library for a multilayer unit cell.
///
/// The unit cell is a STACK; the library sweeps the fill of layer
/// `active_layer` through the whole stack (caps, AR coatings, and spacers
/// included). Returns the same UnitCellLibrary the rest of the pipeline consumes.
/// @param stack              The unit-cell stack (its periods set the lattice).
/// @param active_layer       Index of the layer whose fill is swept.
/// @param incident,substrate Half-spaces above and below the stack.
/// @param wavelength_um      Design wavelength (µm).
/// @param fill_min,fill_max  Fill range to sweep.
/// @param n_samples          Number of fill samples.
/// @param M                  Harmonic half-count in x and y.
UnitCellLibrary build_unit_cell_library_stack(Rcwa2DStack stack, int active_layer,
                                              const Material& incident,
                                              const Material& substrate,
                                              double wavelength_um, double fill_min,
                                              double fill_max, int n_samples, int M);

/// @brief One pillar height tried by optimize_height_for_2pi().
struct HeightSweepEntry {
    double thickness_um;        ///< Pillar height (µm).
    double coverage_deg;        ///< Effective circular coverage, 360° minus the largest gap.
    double mean_transmittance;  ///< Mean |t|² across the fill sweep.
};

/// @brief Result of searching pillar height for full-2π phase coverage.
///
/// Selection rule: among heights whose coverage clears `coverage_target_deg`,
/// pick the highest transmittance — the lever that lifts the transmission-
/// weighted Strehl. If no height clears the target, the height with the largest
/// coverage is returned.
struct HeightOptResult {
    double best_thickness_um;   ///< Winning pillar height (µm).
    double coverage_deg;        ///< Coverage at the winning height (deg).
    double mean_transmittance;  ///< Mean |t|² at the winning height.
    bool reached_target;        ///< True if the winner cleared coverage_target_deg.
    double coverage_target_deg; ///< The coverage target that was requested (deg).
    std::vector<HeightSweepEntry> sweep;  ///< Every height tried, for reporting.
    UnitCellLibrary best_library;         ///< Library at the winning height, ready for design_metalens().
};

/// @brief Full-2π library builder: sweep pillar height, pick the best single etch.
///
/// A single-height fill sweep can cap phase coverage short of 2π, or reach it
/// only at a low-transmittance height. This sweeps the height over
/// [thick_lo, thick_hi] (one fill library per height, still a single etch),
/// measures each height's coverage and mean transmittance, and applies the
/// selection rule in HeightOptResult.
/// @param pillar,background  Pillar and surrounding materials.
/// @param incident,substrate Half-spaces above and below.
/// @param period_um          Lattice period (µm).
/// @param wavelength_um      Design wavelength (µm).
/// @param thick_lo,thick_hi  Height range to sweep (µm).
/// @param n_heights          Number of heights tried.
/// @param fill_min,fill_max  Fill range for each library.
/// @param fill_samples       Fill samples per library.
/// @param M                  Harmonic half-count in x and y.
/// @param coverage_target_deg Coverage a height must reach to qualify (deg).
/// @param shape              Meta-atom cross-section.
/// @param shape_param        Shape parameter (see RectCell2D::shape_param).
/// @return The full sweep table plus the winning library.
HeightOptResult optimize_height_for_2pi(const Material& pillar,
                                        const Material& background,
                                        const Material& incident,
                                        const Material& substrate,
                                        double period_um, double wavelength_um,
                                        double thick_lo, double thick_hi,
                                        int n_heights, double fill_min,
                                        double fill_max, int fill_samples, int M,
                                        double coverage_target_deg = 330.0,
                                        MetaShape shape = MetaShape::Rectangle,
                                        double shape_param = 0.5);

/// @brief A propagation-phase metalens: one pillar fill per lattice site.
struct MetalensDesign {
    int n_cells;                      ///< Cells across the (square) aperture, round(D / period).
    double period_um;                 ///< Lattice period (µm).
    std::vector<double> fill_map;     ///< n_cells × n_cells pillar fills, row-major.
    double rms_phase_error_deg;       ///< RMS realized-vs-target phase error (deg).
    double mean_amplitude;            ///< Mean |t| of the chosen pillars (efficiency proxy).
};

/// @brief Design a metalens for an arbitrary phase profile (propagation-phase path).
///
/// At each lattice site, picks the library pillar whose phase best matches the
/// target φ(x,y). Unlike the geometric-phase path this carries a finite residual
/// phase error (the library is discrete) and a non-uniform amplitude, both
/// reported in MetalensDesign.
/// @param lib               Phase library to draw pillars from.
/// @param profile           Target phase profile.
/// @param diameter_um       Aperture width (µm).
/// @param amplitude_weight  Transmission penalty weight (see UnitCellLibrary::lookup_weighted()).
MetalensDesign design_metalens(const UnitCellLibrary& lib,
                               const PhaseProfile& profile, double diameter_um,
                               double amplitude_weight = 0.25);

/// @brief Design a focusing metalens of the given focal length and aperture.
///
/// Convenience overload for the hyperbolic profile
/// φ(r) = -(2π/λ)(√(r² + f²) - f).
/// @param lib               Phase library to draw pillars from.
/// @param focal_length_um   Focal length (µm).
/// @param diameter_um       Aperture width (µm).
/// @param amplitude_weight  Transmission penalty weight.
MetalensDesign design_metalens(const UnitCellLibrary& lib,
                               double focal_length_um, double diameter_um,
                               double amplitude_weight = 0.25);

} // namespace celeris
