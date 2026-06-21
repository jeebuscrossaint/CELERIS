#pragma once
// Achromatic Pancharatnam-Berry metalens -- geometric phase + dispersion engineering.
//
// This is the modern recipe for a broadband achromatic metalens. It combines the
// two design paths CELERIS already has:
//
//   * the GEOMETRIC (PB) phase -- rotating a birefringent atom by theta imprints
//     a phase -handedness*2*theta on the spin-flipped (cross-circular) output that
//     is EXACT (set by an angle, not by hitting a library target) and, crucially,
//     WAVELENGTH-INDEPENDENT (purely topological);
//   * the PROPAGATION phase -- the atom's own dispersive transmission carries a
//     radius-tunable GROUP DELAY d(phi)/d(omega).
//
// The total phase a focusing achromat needs at frequency omega is
//     phi(r, omega) = -(omega/c)*(sqrt(r^2+f^2) - f) + phi_ref(omega)
//                   = phi(r, omega0) + GD_target(r) * (omega - omega0) + ...
// where GD_target(r) = -(1/c)(sqrt(r^2+f^2) - f) varies with radius.
//
// In a PROPAGATION-only achromat (see achromatic.hpp) each atom must hit BOTH the
// base phase phi(r,omega0) AND the group delay -- a two-objective fit that leaves
// a base-phase RESIDUAL (the library is discrete). PB decouples them: the rotation
// sets the base phase EXACTLY for whatever atom is placed, so the atom is chosen
// PURELY for its group delay (one objective). The realized cross transmission is
//     t_cross(r, omega) = a_cross(atom, omega) * exp(i * phi_geo(r)),
//     a_cross(atom, omega) = (t_x(omega) - t_y(omega)) / 2,   phi_geo = -2*hand*theta,
// and theta(r) is set so arg(t_cross(r, omega0)) = phi(r, omega0) exactly. Hence
// the base-phase RMS is ~0 by construction; the achromatic limit is ONLY the
// library's group-delay coverage. Every atom shares one etch depth (rotated
// rectangles of varying footprint) -> a single fabrication step, no grayscale.
//
// References: Wang et al., Nat. Nanotechnol. 13, 227 (2018); Chen et al.,
// Nat. Nanotechnol. 13, 220 (2018) (the geometric-phase achromat).

#include "celeris/core.hpp"
#include "celeris/design/achromatic.hpp"   // AchromaticFocalPoint
#include "celeris/materials/material.hpp"
#include "celeris/rcwa/rcwa2d.hpp"

#include <vector>

namespace celeris {

// One BIREFRINGENT meta-atom characterized across the band by its spin-flip
// (cross-circular) transmission a_cross(omega) = (t_x(omega) - t_y(omega))/2 --
// the amplitude a PB site contributes once rotated. Birefringence (the HWP
// requirement) and dispersion both come from the rectangle's aspect ratio + size.
struct DispersivePbAtom {
    double fill_x = 0, fill_y = 0;
    double thickness_um = 0;
    double phase0_rad = 0;            // arg(a_cross) at the center wavelength (wrapped)
    double group_delay_fs = 0;        // d(arg a_cross)/d(omega), least-squares over band [fs]
    double mean_amplitude = 0;        // mean |a_cross| over the band (apodization proxy)
    double mean_conversion = 0;       // mean |a_cross|^2 (spin-flip efficiency)
    double retardance_center_deg = 0; // wrap(arg t_x - arg t_y) at center (ideal 180)
    std::vector<cdouble> cross;       // a_cross(lambda): amplitude + phase, for verify
};

// A dispersive birefringent library at ONE etch depth, ready for the achromatic
// PB selection. Atoms span a (fill_x, fill_y) grid so the group delay varies
// while the height is fixed -- a fabricable single-etch achromat.
struct DispersivePbLibrary {
    std::vector<double> wavelengths_um;  // band samples (ascending)
    double center_wavelength_um = 0;
    int center_index = 0;
    double period_um = 0;
    double thickness_um = 0;             // the single shared etch depth
    int n_fill = 0;                      // grid side (atoms = n_fill*n_fill)
    std::vector<DispersivePbAtom> atoms;
    double gd_min_fs = 0, gd_max_fs = 0; // group-delay range supplied
};

// Build a dispersive birefringent library over a (fill_x, fill_y) grid at a single
// height: for each rectangle, two RCWA solves (x- and y-polarized) at every band
// wavelength give t_x(omega), t_y(omega); the spin-flip amplitude a_cross and its
// group delay (least-squares slope of the unwrapped phase vs angular frequency)
// follow. The grid spans a range of group delays at ONE etch depth.
// A non-birefringent atom (fill_x == fill_y -> t_x == t_y) has a_cross ~ 0, so its
// phase -- and therefore its group delay -- is meaningless noise; such atoms would
// pollute the selection (a garbage group delay can match the target while the spin-
// flip amplitude is ~0). `min_amplitude` drops atoms whose mean |a_cross| over the
// band is below it -- only genuinely birefringent (usable PB) atoms are kept.
DispersivePbLibrary build_dispersive_pb_library(
    const Material& pillar, const Material& background, const Material& incident,
    const Material& substrate, double period_um,
    const std::vector<double>& band_wavelengths_um, double center_wavelength_um,
    double fill_min, double fill_max, int n_fills, double thickness_um, int M,
    double min_amplitude = 0.30);

struct PbAchromaticDesign {
    int n_cells = 0;
    double period_um = 0;
    double thickness_um = 0;             // single etch
    int handedness = +1;
    std::vector<int> atom_index;         // per site -> DispersivePbLibrary::atoms
    std::vector<double> rotation_rad;    // per-site geometric-phase rotation theta
    std::vector<double> fill_x_map;      // per-site chosen footprint (row-major)
    std::vector<double> fill_y_map;
    double center_wavelength_um = 0;
    double focal_length_um = 0;
    double diameter_um = 0;
    double rms_phase_error_deg = 0;      // base-phase residual at center (~0: geometric is exact)
    double rms_group_delay_error_fs = 0; // group-delay residual (the achromatic lever)
    double mean_amplitude = 0;           // mean |a_cross| over the aperture
    double mean_conversion = 0;          // mean |a_cross|^2 (efficiency cap)
    double required_gd_span_fs = 0;      // group-delay span the design demanded
    double available_gd_span_fs = 0;     // span the library could supply
    double gd_coverage = 0;              // available/required (>=1 => sufficient)
};

// Design an achromatic PB focusing metalens. At each site: (1) pick the atom whose
// GROUP DELAY best matches the radius-dependent target (the ONLY library
// constraint -- geometric phase handles the base phase); (2) set the rotation so
// the base phase is hit EXACTLY. gd_weight = 0 ignores group delay -> the best-
// conversion atom everywhere = a STANDARD (chromatic) PB lens, the baseline;
// gd_weight > 0 engages the dispersion engineering. amplitude_weight biases toward
// higher spin-flip conversion.
PbAchromaticDesign design_pb_achromatic_metalens(
    const DispersivePbLibrary& lib, double focal_length_um, double diameter_um,
    int handedness = +1, double gd_weight = 1.0, double amplitude_weight = 0.25);

// Focal length vs wavelength using the library's STORED per-atom band response
// (no new RCWA): at each band wavelength propagate the aperture with each site's
// realized cross transmission a_cross(omega)*exp(i*phi_geo). The rigorous chromatic
// response; flat = achromatic. Reuses AchromaticFocalPoint (achromatic.hpp).
std::vector<AchromaticFocalPoint> verify_pb_achromatic_focus(
    const DispersivePbLibrary& lib, const PbAchromaticDesign& design);

} // namespace celeris
