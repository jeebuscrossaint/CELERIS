#pragma once
// Phase profiles phi(x,y): the target wavefront a metasurface imprints.
//
// A profile is just a rule for the target phase at a lattice point (x,y). The
// SAME profile drives both design paths:
//   - the propagation-phase lens (`design_metalens`): vary pillar SIZE to hit
//     phi via a library lookup (finite library -> a small residual phase error);
//   - the geometric-phase / Pancharatnam-Berry lens (`design_pb_metalens`):
//     rotate a fixed half-wave-plate atom by theta = -handedness*phi/2 (exact).
//
// The canonical profiles are analytic (focusing / vortex / deflector / axicon);
// `Freeform` carries an ARBITRARY loaded phi(x,y) map -- a computer-generated
// hologram (CGH) or any freeform wavefront -- sampled bilinearly.

#include "celeris/core.hpp"

#include <string>
#include <vector>

namespace celeris {

enum class PhaseProfileKind {
    Focusing,   // hyperbolic lens: phi = -k(sqrt(r^2+f^2)-f)              -> focal spot
    Vortex,     // OAM plate:       phi = l*atan2(y,x) [+focusing if f>0]   -> donut / OAM beam
    Deflector,  // blazed grating:  phi = k*sin(a)*(x*cos az + y*sin az)    -> tilted beam
    Axicon,     // conical phase:   phi = -k*sin(b)*r                       -> Bessel / line focus
    Freeform,   // arbitrary loaded phi(x,y) map (hologram / CGH), bilinearly sampled
};

// Parameters for a phase profile. Only the fields relevant to `kind` are used.
struct PhaseProfile {
    PhaseProfileKind kind = PhaseProfileKind::Focusing;
    double focal_length_um = 50.0;     // Focusing; Vortex (focused vortex; <=0 => pure OAM)
    int topological_charge = 1;        // Vortex: OAM charge l (winds 2*pi*l around the axis)
    double deflect_deg = 10.0;         // Deflector: beam deflection angle from normal
    double deflect_azimuth_deg = 0.0;  // Deflector: in-plane direction of the deflection
    double axicon_deg = 5.0;           // Axicon: cone half-angle (ray bend toward the axis)

    // Freeform: a square n x n grid of target phase (radians), row-major, that
    // physically spans [-extent/2, +extent/2] in BOTH x and y. Sampled bilinearly,
    // edges clamped. Populated by load_freeform_phase or built in memory.
    std::vector<double> freeform_phase_rad;
    int freeform_n = 0;
    double freeform_extent_um = 0.0;
};

// Target geometric/optical phase phi(x,y) [rad] at lattice point (x,y), using the
// engine's exp(+i*k*r) propagator sign convention (so Focusing converges on-axis
// at z = f, matching the propagation check).
double phase_profile_value(const PhaseProfile& profile, double x, double y,
                           double wavelength_um);

// Load a freeform phase map from a whitespace-separated text grid of phase values
// in RADIANS (the count must be a perfect square -> an n x n grid). `extent_um` is
// the full physical width the map spans in x and y (it is centered on the origin).
// '#' begins a comment to end-of-line. Throws std::runtime_error on a bad file.
PhaseProfile load_freeform_phase(const std::string& path, double extent_um);

} // namespace celeris
