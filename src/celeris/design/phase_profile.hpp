#pragma once
/// @file
/// @brief Phase profiles φ(x,y): the target wavefront a metasurface imprints.
///
/// A profile is just a rule for the target phase at a lattice point (x,y). The
/// SAME profile drives both design paths:
///   - the propagation-phase lens (design_metalens()): vary pillar SIZE to hit
///     φ via a library lookup (finite library -> a small residual phase error);
///   - the geometric-phase / Pancharatnam–Berry lens (design_pb_metalens()):
///     rotate a fixed half-wave-plate atom by θ = -handedness·φ/2 (exact).
///
/// The canonical profiles are analytic (focusing / vortex / deflector / axicon);
/// `Freeform` carries an ARBITRARY loaded φ(x,y) map -- a computer-generated
/// hologram (CGH) or any freeform wavefront -- sampled bilinearly.

#include "celeris/core.hpp"

#include <string>
#include <vector>

namespace celeris {

/// @brief Which analytic (or loaded) phase profile to imprint.
enum class PhaseProfileKind {
    Focusing,   ///< Hyperbolic lens, φ = -k(√(r²+f²) - f): a diffraction-limited focal spot.
    Quadratic,  ///< Parabolic lens, φ = -k r²/(2f): a wide-field-of-view focus.
    Vortex,     ///< OAM plate, φ = l·atan2(y,x), plus focusing if f > 0: a donut / OAM beam.
    Deflector,  ///< Blazed grating, φ = k sin(a)(x cos az + y sin az): a tilted beam.
    Axicon,     ///< Conical phase, φ = -k sin(b) r: a Bessel beam / line focus.
    Freeform,   ///< Arbitrary loaded φ(x,y) map (hologram / CGH), sampled bilinearly.
};

/// @brief Parameters for a phase profile; only the fields relevant to `kind` are used.
struct PhaseProfile {
    PhaseProfileKind kind = PhaseProfileKind::Focusing;  ///< Profile type.
    /// Focal length (µm) for Focusing/Quadratic, and for Vortex (a focused
    /// vortex; ≤ 0 gives a pure OAM plate).
    double focal_length_um = 50.0;
    int topological_charge = 1;        ///< Vortex: OAM charge l (phase winds 2πl around the axis).
    double deflect_deg = 10.0;         ///< Deflector: beam deflection angle from the normal (deg).
    double deflect_azimuth_deg = 0.0;  ///< Deflector: in-plane direction of the deflection (deg).
    double axicon_deg = 5.0;           ///< Axicon: cone half-angle, the ray bend toward the axis (deg).

    /// @brief Freeform: square n×n grid of target phase (rad), row-major.
    ///
    /// The grid spans [-extent/2, +extent/2] in both x and y and is sampled
    /// bilinearly with edges clamped. Populated by load_freeform_phase() or
    /// built in memory.
    std::vector<double> freeform_phase_rad;
    int freeform_n = 0;               ///< Freeform: grid size n.
    double freeform_extent_um = 0.0;  ///< Freeform: full physical width of the grid (µm).
};

/// @brief Target phase φ(x,y) at a lattice point.
///
/// Uses the engine's exp(+ikr) propagator sign convention, so Focusing converges
/// on-axis at z = f, matching the propagation check.
/// @param profile        The profile to evaluate.
/// @param x,y            Position in the lens plane (µm), origin at the lens center.
/// @param wavelength_um  Design wavelength (µm).
/// @return Target phase (rad).
double phase_profile_value(const PhaseProfile& profile, double x, double y,
                           double wavelength_um);

/// @brief Load a freeform phase map from a text file.
///
/// The file is a whitespace-separated grid of phase values in radians; the count
/// must be a perfect square (an n×n grid). `#` begins a comment to end of line.
/// @param path       File to read.
/// @param extent_um  Full physical width the map spans in x and y, centered on
///                   the origin (µm).
/// @return A PhaseProfile of kind Freeform.
/// @throws std::runtime_error on a missing or malformed file.
PhaseProfile load_freeform_phase(const std::string& path, double extent_um);

} // namespace celeris
