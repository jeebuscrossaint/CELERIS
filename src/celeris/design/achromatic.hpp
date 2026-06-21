#pragma once
// Achromatic (broadband) metalens design -- dispersion engineering.
//
// A standard phase-profile metalens fixes the phase only at the design
// wavelength, so its focal length drifts as f(lambda) ~ f0*lambda0/lambda (the
// #1 limitation of metalenses). To focus a whole BAND at one plane, the phase a
// meta-atom imparts must satisfy, at every frequency omega in the band,
//
//     phi(r, omega) = -(omega/c) * (sqrt(r^2 + f^2) - f) + phi_ref(omega).
//
// Expanding around the center frequency omega0, two terms must be controlled per
// site:
//   * phi(r, omega0)         -- the BASE phase (matched mod 2*pi, as usual);
//   * d phi / d omega |omega0 -- the GROUP DELAY, which MUST vary with radius
//                                (maximum at center, decreasing to the edge).
// The required group-delay SPAN, (1/c)(sqrt(R^2+f^2)-f), grows with aperture and
// bandwidth. The catch: to hit an ARBITRARY (phase, group-delay) pair per site,
// the atom library must cover the 2-D (phase, GD) plane -- a single geometric
// degree of freedom (e.g. fill only) traces just a 1-D curve and cannot set both
// independently. So achromatic design needs a richer library (here: fill x
// height) and a TWO-objective atom selection.
//
// References: Wang et al., Nat. Nanotechnol. 13, 227 (2018); Chen et al.,
// Nat. Nanotechnol. 13, 220 (2018); Shrestha et al., Light Sci. Appl. 7, 85 (2018).

#include "celeris/design/metalens.hpp"
#include "celeris/materials/material.hpp"
#include "celeris/rcwa/rcwa2d.hpp"

#include <vector>

namespace celeris {

// One meta-atom characterized across a wavelength band: its phase/amplitude at
// every band sample plus the derived center phase and group delay.
struct DispersiveAtom {
    double fill = 0;
    double thickness_um = 0;
    double phase0_rad = 0;        // transmission phase at the center wavelength (wrapped)
    double group_delay_fs = 0;    // d(phi)/d(omega) at center, least-squares over the band [fs]
    double mean_amplitude = 0;    // mean |t| over the band (efficiency proxy)
    std::vector<double> phase_rad;   // phi(lambda) at each band sample
    std::vector<double> amplitude;   // |t|(lambda) at each band sample
};

// A meta-atom library characterized over a wavelength band, ready for the
// achromatic two-objective selection. Atoms span a fill x height grid so the
// (phase, group-delay) plane is covered, not just a 1-D curve.
struct DispersiveLibrary {
    std::vector<double> wavelengths_um;  // band samples (ascending)
    double center_wavelength_um = 0;
    int center_index = 0;                // index of the center wavelength
    double period_um = 0;
    int n_fill = 0, n_height = 0;        // grid dimensions (atoms = n_fill*n_height)
    std::vector<DispersiveAtom> atoms;   // one per (fill, height) sample
    double gd_min_fs = 0, gd_max_fs = 0; // group-delay range the library supplies
};

// Build a dispersive library over a fill x height grid: for each (fill, height)
// atom, solve the 2D RCWA at EVERY band wavelength so phi(lambda)/|t|(lambda)
// capture the meta-atom's true material + waveguide dispersion. The group delay
// is the least-squares slope of the (band-unwrapped) phase vs angular frequency.
// `band_wavelengths_um` must be ascending; `center_wavelength_um` selects the
// base-phase reference sample. n_heights==1 gives a single-etch (1-DOF) library.
DispersiveLibrary build_dispersive_library(
    const Material& pillar, const Material& background, const Material& incident,
    const Material& substrate, double period_um,
    const std::vector<double>& band_wavelengths_um, double center_wavelength_um,
    double fill_min, double fill_max, int n_fills, double thick_lo, double thick_hi,
    int n_heights, int M);

struct AchromaticDesign {
    int n_cells = 0;
    double period_um = 0;
    std::vector<int> atom_index;         // per site -> index into DispersiveLibrary::atoms
    std::vector<double> fill_map;         // per site chosen pillar fill (row-major)
    std::vector<double> thickness_map;    // per site chosen pillar height (row-major)
    double center_wavelength_um = 0;
    double focal_length_um = 0;
    double diameter_um = 0;
    double rms_phase_error_deg = 0;      // base-phase residual at the center wavelength
    double rms_group_delay_error_fs = 0; // group-delay residual (the achromatic lever)
    double mean_amplitude = 0;
    double required_gd_span_fs = 0;      // group-delay span the design demanded
    double available_gd_span_fs = 0;     // span the library could supply
    double gd_coverage = 0;              // available/required (>=1 => library is sufficient)
    bool single_height = true;           // did every chosen atom share one height?
    double min_height_um = 0, max_height_um = 0;
};

// Design a focusing metalens by two-objective atom selection: at each site choose
// the atom that best matches BOTH the base focusing phase (mod 2*pi) and the
// radius-dependent group delay. `gd_weight` scales the group-delay objective --
// the group-delay error is converted to the phase error it causes at the band
// EDGE (its physical effect on focus there), so gd_weight ~ 1 balances center-
// wavelength focus against band-edge focus. gd_weight = 0 reproduces a standard
// single-wavelength (dispersion-blind) design from the same library -- the
// baseline for comparison. `amplitude_weight` biases toward higher transmission.
AchromaticDesign design_achromatic_metalens(const DispersiveLibrary& lib,
                                            double focal_length_um,
                                            double diameter_um,
                                            double gd_weight = 1.0,
                                            double amplitude_weight = 0.25);

// Focal length vs wavelength for a design, using the library's STORED per-atom
// band response (no new RCWA solves): at each band wavelength, propagate the
// aperture with each site's atom transmission at that wavelength and locate the
// on-axis intensity peak. This is the honest rigorous chromatic response (each
// atom's true dispersion) and works for per-site varying heights.
struct AchromaticFocalPoint {
    double wavelength_um;
    double focal_length_um;
    double rel_peak;  // on-axis peak relative to the center-wavelength design point
};
std::vector<AchromaticFocalPoint> verify_achromatic_focus(
    const DispersiveLibrary& lib, const AchromaticDesign& design);

// Adapt to the plain MetalensDesign the GDS writer / analysis battery consume
// (in-plane footprints only; a multi-height design's etch depths are not encoded
// in a single GDS layer -- see AchromaticDesign::single_height).
MetalensDesign to_metalens_design(const AchromaticDesign& d);

} // namespace celeris
