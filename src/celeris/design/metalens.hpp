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
#include "celeris/rcwa/rcwa2d.hpp"

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
    // NOTE: this is max-min of the wrapped phases; it is only meaningful when the
    // phase-vs-fill curve does not wrap across +/-pi. For the honest "can I hit
    // any target phase" measure use coverage() below.
    double phase_span() const;
    // Effective phase coverage on the circle = 2*pi minus the largest angular gap
    // between consecutive available phases (radians). This is the metric that
    // governs design fidelity: a small largest-gap means any target phase has a
    // near match, regardless of where the wrapped phases sit. Equals the span for
    // a monotonic non-wrapping sweep, but stays correct when the phase wraps.
    double coverage() const;
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

// Multi-layer "system" path: the unit cell is a STACK; the phase library is
// built by sweeping the fill of layer `active_layer` through the whole stack
// (caps / AR coatings / spacers included). Returns the same UnitCellLibrary
// the rest of the pipeline consumes.
UnitCellLibrary build_unit_cell_library_stack(Rcwa2DStack stack, int active_layer,
                                              const Material& incident,
                                              const Material& substrate,
                                              double wavelength_um, double fill_min,
                                              double fill_max, int n_samples, int M);

// One height tried by the full-2pi sweep: its achievable phase coverage and the
// mean power transmittance (mean |t|^2) over the fill sweep at that height.
struct HeightSweepEntry {
    double thickness_um;
    double coverage_deg;        // effective circular coverage (360 - largest gap)
    double mean_transmittance;  // mean |t|^2 across the fill sweep
};

// Result of searching pillar height for full-2pi phase coverage. Selection rule:
// among heights whose effective coverage clears `coverage_target_deg` (so any
// target phase has a near match), pick the HIGHEST transmittance -- that is the
// lever that lifts the transmission-weighted Strehl. If no height clears the
// target, the one with the largest coverage is returned.
struct HeightOptResult {
    double best_thickness_um;
    double coverage_deg;
    double mean_transmittance;
    bool reached_target;        // did the winner clear coverage_target_deg?
    double coverage_target_deg;
    std::vector<HeightSweepEntry> sweep;  // every height tried (for reporting)
    UnitCellLibrary best_library;         // ready to feed design_metalens
};

// Full-2pi library builder. A single-height fill sweep can cap phase coverage
// short of 2pi (or hit it only at a low-transmittance height), which limits the
// transmission-weighted Strehl. Sweep the pillar HEIGHT over [thick_lo,thick_hi]
// (one fill-library per height, single-etch so still fabricable), measure each
// height's effective coverage and mean transmittance, and pick the height that
// clears the coverage target with the highest transmittance. Returns the full
// table plus the winning library, ready for design_metalens.
HeightOptResult optimize_height_for_2pi(const Material& pillar,
                                        const Material& background,
                                        const Material& incident,
                                        const Material& substrate,
                                        double period_um, double wavelength_um,
                                        double thick_lo, double thick_hi,
                                        int n_heights, double fill_min,
                                        double fill_max, int fill_samples, int M,
                                        double coverage_target_deg = 330.0);

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
