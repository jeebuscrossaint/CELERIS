#pragma once
// Polarization-multiplexed metalens design. A rectangular nanopillar imparts
// independent phase to x- and y-polarized light (form birefringence), so one
// metasurface can implement TWO different phase profiles at once -- e.g. focus
// x-polarized light to one point and y-polarized light to another. This is the
// basis of polarization-splitting / polarization-imaging metalenses.
//
// We build a 2D library over (fill_x, fill_y) -- two RCWA solves per cell, one
// per incident polarization -- then for each lattice site pick the pillar whose
// (phase_x, phase_y) best matches the two target profiles simultaneously.

#include "celeris/core.hpp"
#include "celeris/materials/material.hpp"

#include <vector>

namespace celeris {

// One sampled rectangular meta-atom and its dual-polarization response.
struct PolarCell {
    double fill_x, fill_y;     // pillar widths (fractions of period)
    double phase_x, phase_y;   // imparted phase per polarization (rad)
    double amp_x, amp_y;       // |t| per polarization
};

struct PolarizationLibrary {
    double period_um = 0, wavelength_um = 0, thickness_um = 0;
    std::vector<PolarCell> cells;  // n_samples^2 entries

    // Index of the cell minimizing the combined phase error to (target_x,
    // target_y), with an optional amplitude penalty (favor high transmission).
    int lookup(double target_x, double target_y, double amplitude_weight = 0.0) const;
};

// Build the (fill_x, fill_y) library: n_samples x n_samples rectangular pillars,
// each solved twice (x- and y-polarized). Solves run in parallel.
PolarizationLibrary build_polarization_library(
    const Material& pillar, const Material& background, const Material& incident,
    const Material& substrate, double period_um, double wavelength_um,
    double thickness_um, double fill_min, double fill_max, int n_samples, int M);

struct PolarMetalensDesign {
    int n_cells = 0;
    double period_um = 0;
    std::vector<double> fill_x;  // n_cells^2, row-major
    std::vector<double> fill_y;
    std::vector<cdouble> t_x;    // realized per-pillar transmission, X-pol
    std::vector<cdouble> t_y;    // realized per-pillar transmission, Y-pol
    double rms_phase_error_x_deg = 0, rms_phase_error_y_deg = 0;
    double mean_amp_x = 0, mean_amp_y = 0;
};

// Design a polarization-multiplexed lens: x-polarized light focuses at focal_x_um,
// y-polarized light at focal_y_um (set them equal for a polarization-insensitive
// lens). Maps the two ideal hyperbolic phase profiles onto rectangular pillars.
PolarMetalensDesign design_polarization_metalens(const PolarizationLibrary& lib,
                                                 double focal_x_um,
                                                 double focal_y_um,
                                                 double diameter_um,
                                                 double amplitude_weight = 0.25);

} // namespace celeris
