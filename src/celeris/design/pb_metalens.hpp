#pragma once
// Pancharatnam-Berry (geometric-phase) metalens design.
//
// A propagation-phase metalens (the `design`/`polardesign` path) varies the
// pillar's SIZE to dial the optical phase. A Pancharatnam-Berry lens instead
// uses ONE fixed birefringent meta-atom -- a half-wave plate (HWP) -- and
// ROTATES it from site to site. Under circularly polarized illumination, an
// HWP rotated by angle theta converts the light to the opposite handedness and
// stamps a purely GEOMETRIC phase of -/+ 2*theta on it (sign = input spin).
//
// Two properties make this powerful: the imprinted phase is EXACT (set by an
// angle, not by hitting a target with a finite library) and the transmission
// amplitude is IDENTICAL at every site (same atom), so a PB lens is inherently
// phase-error-free and apodization-free -- diffraction-limited, capped only by
// the atom's polarization-conversion efficiency. It is the standard route to
// circular-polarization optics: vortex/OAM plates, spin-Hall devices, and the
// achromatic-friendly geometric-phase metalenses (Khorasaninejad 2016).
//
// References: Pancharatnam (1956); Berry (1987); Kang et al. (2012);
// Khorasaninejad & Capasso, Science 352, 1190 (2016).

#include "celeris/core.hpp"
#include "celeris/materials/material.hpp"
#include "celeris/rcwa/rcwa2d.hpp"

#include <vector>

namespace celeris {

// Full 2x2 linear-basis Jones transmission matrix (zeroth transmitted order):
//   [Eout_x; Eout_y] = J [Ein_x; Ein_y].
struct JonesMatrix {
    cdouble xx, xy, yx, yy;
};

// Zeroth-order Jones matrix of a single-layer 2D cell at normal incidence, from
// two RCWA solves (x- and y-polarized illumination).
JonesMatrix solve_jones(const Material& incident, const Rcwa2DStack& stack,
                        const Material& substrate, double wavelength_um, int M);

// The half-wave-plate meta-atom that PB optics is built from: the rectangular
// pillar whose form birefringence best flips incident circular polarization to
// the opposite handedness. Ideal HWP: |t_x| = |t_y| = 1, retardance = 180 deg,
// conversion efficiency = 1.
struct HwpAtom {
    double fill_x = 0, fill_y = 0;
    double thickness_um = 0;
    cdouble t_x = 0, t_y = 0;          // co-pol linear transmissions at rotation 0
    double retardance_deg = 0;         // wrap(arg(t_x) - arg(t_y))
    double conversion_efficiency = 0;  // |t_x - t_y|^2 / 4  (spin-flip, in [0,1])
};

// Search a (fill_x, fill_y) grid of rectangular pillars for the best HWP atom
// (max spin-flip conversion efficiency). Reuses the polarization-library sweep.
HwpAtom find_hwp_atom(const Material& pillar, const Material& background,
                      const Material& incident, const Material& substrate,
                      double period_um, double wavelength_um, double thickness_um,
                      double fill_min, double fill_max, int n_samples, int M);

// One RCWA-measured sample of the geometric-phase relation (RCP illumination).
struct PbVerifyPoint {
    double rotation_deg;     // atom rotation theta
    double cross_phase_deg;  // phase of the spin-flipped output (tracks -2*theta)
    double conversion_eff;   // |cross amplitude|^2 (spin-flipped, focused)
    double copol_leakage;    // |co amplitude|^2 (same handedness, unfocused)
};

// RCWA-verify the PB relation: solve the ROTATED atom at each angle and report
// the spin-flip phase/efficiency. The proof that rotating the atom imprints
// 2*theta (incident handedness = RCP).
std::vector<PbVerifyPoint> verify_pb_phase(
    const Material& pillar, const Material& background, const Material& incident,
    const Material& substrate, double period_um, double wavelength_um,
    const HwpAtom& atom, const std::vector<double>& rotations_rad, int M);

// A PB metasurface stamps ANY target phase profile phi(x,y) -- the rotation map
// theta = -handedness * phi / 2 is profile-agnostic. These are the canonical
// geometric-phase profiles; each maps to a different optical element.
enum class PbProfileKind {
    Focusing,   // hyperbolic lens: phi = -k(sqrt(r^2+f^2) - f)            -> focal spot
    Vortex,     // OAM plate:       phi = l*atan2(y,x) [+ focusing if f>0]  -> donut / OAM beam
    Deflector,  // blazed grating:  phi = k*sin(a)*(x*cos(az)+y*sin(az))    -> tilted beam
    Axicon,     // conical phase:   phi = -k*sin(b)*r                       -> Bessel / line focus
};

// Parameters for a PB phase profile. Only the fields relevant to `kind` are used.
struct PbProfile {
    PbProfileKind kind = PbProfileKind::Focusing;
    double focal_length_um = 50.0;     // Focusing; Vortex (focused vortex; <=0 => pure OAM)
    int topological_charge = 1;        // Vortex: OAM charge l (winds 2*pi*l around the axis)
    double deflect_deg = 10.0;         // Deflector: beam deflection angle from normal
    double deflect_azimuth_deg = 0.0;  // Deflector: in-plane direction of the deflection
    double axicon_deg = 5.0;           // Axicon: cone half-angle (ray bend toward the axis)
};

// Target geometric phase phi(x,y) [rad] for a PB profile at lattice point (x,y),
// using the engine's exp(+i*k*r) propagator sign convention (so Focusing matches
// the hyperbolic-lens phase the propagation check focuses at z=f).
double pb_profile_phase(const PbProfile& profile, double x, double y,
                        double wavelength_um);

struct PbMetalensDesign {
    int n_cells = 0;
    double period_um = 0;
    HwpAtom atom;
    std::vector<double> rotation_rad;  // n_cells^2 row-major, per-site atom rotation
    std::vector<cdouble> t_cross;      // realized spin-flip transmission per site
    double rms_phase_error_deg = 0;    // ~0: geometric phase is exact
    double conversion_efficiency = 0;  // spin-flip efficiency (focusing-efficiency cap)
};

// Design a PB metasurface for an ARBITRARY phase profile: a fixed HWP atom
// rotated per site so the spin-flipped (cross-circular) output carries phi(x,y).
// handedness = +1 for RCP illumination (cross phase = -2*theta), -1 for LCP.
PbMetalensDesign design_pb_metalens(const HwpAtom& atom, double period_um,
                                    double wavelength_um, const PbProfile& profile,
                                    double diameter_um, int handedness = +1);

// Convenience overload: a PB focusing lens (the original signature).
PbMetalensDesign design_pb_metalens(const HwpAtom& atom, double period_um,
                                    double wavelength_um, double focal_length_um,
                                    double diameter_um, int handedness = +1);

} // namespace celeris
