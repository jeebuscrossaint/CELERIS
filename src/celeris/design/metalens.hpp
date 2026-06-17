#pragma once
// Metalens design: the layer that turns the RCWA solver into a usable tool.
//
// Two pieces:
//   1. UnitCellLibrary — sweep a nanopillar's size and tabulate the complex
//      transmission (phase + amplitude) it produces. This is the "phase
//      library": the palette of phase shifts available from one fab process.
//   2. design_metalens — given a focal length and aperture, compute the ideal
//      focusing phase profile and, at each lattice site, pick the pillar whose
//      library phase best matches. The output is a fabricable pillar map plus
//      the residual phase error (design fidelity).

#include "celeris/materials/material.hpp"

#include <vector>

namespace celeris {

struct UnitCellLibrary {
    double period_um;
    double wavelength_um;
    double thickness_um;
    std::vector<double> fill;       // pillar size fraction sampled
    std::vector<double> phase;      // transmission phase (rad), in (-pi, pi]
    std::vector<double> amplitude;  // |t| of the zeroth order

    // Phase coverage actually achievable (max - min over the sweep), radians.
    double phase_span() const;
    // Index of the library entry whose phase best matches target (circularly).
    int lookup(double target_phase_rad) const;
    // Like lookup, but trades off phase error against transmission: minimizes
    // phase_error^2 + amplitude_weight*(1-|t|)^2. With weight 0 this is lookup.
    int lookup_weighted(double target_phase_rad, double amplitude_weight) const;
    // Complex transmission t = |t|·e^{iφ} of the sampled pillar nearest `fill`.
    cdouble transmission_for_fill(double fill) const;
};

// Build a library by sweeping square-pillar fill from fill_min..fill_max.
// Each sample is a full 2D RCWA solve at normal incidence (x-polarized).
UnitCellLibrary build_unit_cell_library(const Material& pillar,
                                        const Material& background,
                                        const Material& incident,
                                        const Material& substrate,
                                        double period_um, double wavelength_um,
                                        double thickness_um, double fill_min,
                                        double fill_max, int n_samples, int M);

struct MetalensDesign {
    int n_cells;             // cells across the (square) aperture
    double period_um;
    std::vector<double> fill_map;     // n_cells*n_cells row-major pillar fills
    double rms_phase_error_deg;       // realized vs target phase fidelity
    double mean_amplitude;            // avg |t| of chosen pillars (efficiency proxy)
};

// Design a focusing metalens of the given focal length and aperture diameter
// using the supplied library. Target phase is the hyperbolic lens profile
//   phi(r) = -(2*pi/lambda) * (sqrt(r^2 + f^2) - f).
MetalensDesign design_metalens(const UnitCellLibrary& lib,
                               double focal_length_um, double diameter_um,
                               double amplitude_weight = 0.25);

} // namespace celeris
