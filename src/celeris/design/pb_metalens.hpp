#pragma once
/// @file
/// @brief Pancharatnam-Berry (geometric-phase) metalens design.
///
/// A propagation-phase metalens (the `design`/`polardesign` path) varies the
/// pillar's SIZE to dial the optical phase. A Pancharatnam-Berry lens instead
/// uses ONE fixed birefringent meta-atom -- a half-wave plate (HWP) -- and
/// ROTATES it from site to site. Under circularly polarized illumination, an
/// HWP rotated by angle theta converts the light to the opposite handedness and
/// stamps a purely GEOMETRIC phase of -/+ 2*theta on it (sign = input spin).
///
/// Two properties make this powerful: the imprinted phase is EXACT (set by an
/// angle, not by hitting a target with a finite library) and the transmission
/// amplitude is IDENTICAL at every site (same atom), so a PB lens is inherently
/// phase-error-free and apodization-free -- diffraction-limited, capped only by
/// the atom's polarization-conversion efficiency. It is the standard route to
/// circular-polarization optics: vortex/OAM plates, spin-Hall devices, and the
/// achromatic-friendly geometric-phase metalenses (Khorasaninejad 2016).
///
/// References: Pancharatnam (1956); Berry (1987); Kang et al. (2012);
/// Khorasaninejad & Capasso, Science 352, 1190 (2016).

#include "celeris/core.hpp"
#include "celeris/design/phase_profile.hpp"
#include "celeris/materials/material.hpp"
#include "celeris/rcwa/rcwa2d.hpp"

#include <vector>

namespace celeris {

/// @brief Full 2×2 linear-basis Jones transmission matrix (zeroth transmitted order).
///
/// [Eout_x; Eout_y] = J [Ein_x; Ein_y].
struct JonesMatrix {
    cdouble xx;  ///< x-in → x-out.
    cdouble xy;  ///< y-in → x-out.
    cdouble yx;  ///< x-in → y-out.
    cdouble yy;  ///< y-in → y-out.
};

/// @brief Zeroth-order Jones matrix of a single-layer 2D cell at normal incidence.
///
/// Two RCWA solves, one x-polarized and one y-polarized.
/// @param incident,substrate Half-spaces above and below.
/// @param stack              The 2D cell.
/// @param wavelength_um      Vacuum wavelength (µm).
/// @param M                  Harmonic half-count in x and y.
JonesMatrix solve_jones(const Material& incident, const Rcwa2DStack& stack,
                        const Material& substrate, double wavelength_um, int M);

/// @brief The half-wave-plate meta-atom that PB optics is built from.
///
/// The rectangular pillar whose form birefringence best flips incident circular
/// polarization to the opposite handedness. Ideal HWP: |t_x| = |t_y| = 1,
/// retardance = 180°, conversion efficiency = 1.
struct HwpAtom {
    double fill_x = 0;                 ///< Pillar width in x (fraction of the period).
    double fill_y = 0;                 ///< Pillar width in y (fraction of the period).
    double thickness_um = 0;           ///< Pillar height (µm).
    cdouble t_x = 0;                   ///< Co-pol transmission for x-polarized light, at rotation 0.
    cdouble t_y = 0;                   ///< Co-pol transmission for y-polarized light, at rotation 0.
    double retardance_deg = 0;         ///< wrap(arg(t_x) − arg(t_y)) (deg); ideal 180.
    double conversion_efficiency = 0;  ///< Spin-flip efficiency |t_x − t_y|² / 4, in [0, 1].
};

/// @brief Find the best HWP atom on a (fill_x, fill_y) grid of rectangular pillars.
///
/// Maximizes spin-flip conversion efficiency. Reuses the polarization-library
/// sweep (`build_polarization_library`).
/// @param pillar,background  Pillar and surrounding materials.
/// @param incident,substrate Half-spaces above and below.
/// @param period_um          Lattice period (µm).
/// @param wavelength_um      Design wavelength (µm).
/// @param thickness_um       Pillar height (µm).
/// @param fill_min,fill_max  Fill range, applied to both axes.
/// @param n_samples          Samples per axis.
/// @param M                  Harmonic half-count in x and y.
HwpAtom find_hwp_atom(const Material& pillar, const Material& background,
                      const Material& incident, const Material& substrate,
                      double period_um, double wavelength_um, double thickness_um,
                      double fill_min, double fill_max, int n_samples, int M);

/// @brief One RCWA-measured sample of the geometric-phase relation (RCP illumination).
struct PbVerifyPoint {
    double rotation_deg;     ///< Atom rotation theta (deg).
    double cross_phase_deg;  ///< Phase of the spin-flipped output (deg); tracks −2·theta.
    double conversion_eff;   ///< |cross amplitude|² (spin-flipped, focused).
    double copol_leakage;    ///< |co amplitude|² (same handedness, unfocused).
};

/// @brief RCWA-verify the PB relation by solving the ROTATED atom at each angle.
///
/// Reports the spin-flip phase and efficiency per angle: the check that
/// rotating the atom by theta imprints 2·theta. Incident handedness is RCP.
/// @param pillar,background  Pillar and surrounding materials.
/// @param incident,substrate Half-spaces above and below.
/// @param period_um          Lattice period (µm).
/// @param wavelength_um      Design wavelength (µm).
/// @param atom               The HWP atom to rotate (from `find_hwp_atom`).
/// @param rotations_rad      Rotation angles to solve (rad).
/// @param M                  Harmonic half-count in x and y.
std::vector<PbVerifyPoint> verify_pb_phase(
    const Material& pillar, const Material& background, const Material& incident,
    const Material& substrate, double period_um, double wavelength_um,
    const HwpAtom& atom, const std::vector<double>& rotations_rad, int M);

// A PB metasurface stamps ANY target phase profile phi(x,y) -- the rotation map
// theta = -handedness * phi / 2 is profile-agnostic, so the canonical profiles
// (focusing / vortex / deflector / axicon / freeform) all flow through the same
// design path. The profile type lives in design/phase_profile.hpp (shared with
// the propagation-phase `design_metalens` path).

/// @brief A PB metasurface: one fixed HWP atom, rotated per site.
struct PbMetalensDesign {
    int n_cells = 0;                   ///< Cells across the (square) aperture.
    double period_um = 0;              ///< Lattice period (µm).
    HwpAtom atom;                      ///< The atom placed at every site.
    std::vector<double> rotation_rad;  ///< n_cells² per-site atom rotations (rad), row-major.
    std::vector<cdouble> t_cross;      ///< Realized spin-flip transmission per site.
    double rms_phase_error_deg = 0;    ///< RMS phase error (deg); ~0, geometric phase is exact.
    double conversion_efficiency = 0;  ///< Spin-flip efficiency (the focusing-efficiency cap).
};

/// @brief Design a PB metasurface for an arbitrary phase profile.
///
/// A fixed HWP atom is rotated per site, theta = −handedness·phi/2, so the
/// spin-flipped (cross-circular) output carries phi(x, y).
/// @param atom           The HWP atom (from `find_hwp_atom`).
/// @param period_um      Lattice period (µm).
/// @param wavelength_um  Design wavelength (µm).
/// @param profile        Target phase profile.
/// @param diameter_um    Aperture width (µm).
/// @param handedness     +1 for RCP illumination (cross phase = −2·theta), −1 for LCP.
PbMetalensDesign design_pb_metalens(const HwpAtom& atom, double period_um,
                                    double wavelength_um, const PhaseProfile& profile,
                                    double diameter_um, int handedness = +1);

/// @brief Convenience overload: a PB focusing lens with focal length `focal_length_um` (µm).
///
/// Equivalent to the profile overload with a focusing profile.
PbMetalensDesign design_pb_metalens(const HwpAtom& atom, double period_um,
                                    double wavelength_um, double focal_length_um,
                                    double diameter_um, int handedness = +1);

} // namespace celeris
