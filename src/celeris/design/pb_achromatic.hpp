#pragma once
/// @file
/// @brief Achromatic Pancharatnam-Berry metalens: geometric phase + dispersion engineering.
///
/// This is the modern recipe for a broadband achromatic metalens. It combines the
/// two design paths CELERIS already has:
///
///   - the GEOMETRIC (PB) phase -- rotating a birefringent atom by theta imprints
///     a phase -handedness*2*theta on the spin-flipped (cross-circular) output that
///     is EXACT (set by an angle, not by hitting a library target) and, crucially,
///     WAVELENGTH-INDEPENDENT (purely topological);
///   - the PROPAGATION phase -- the atom's own dispersive transmission carries a
///     radius-tunable GROUP DELAY d(phi)/d(omega).
///
/// The total phase a focusing achromat needs at frequency omega is
///
///     phi(r, omega) = -(omega/c)*(sqrt(r^2+f^2) - f) + phi_ref(omega)
///                   = phi(r, omega0) + GD_target(r) * (omega - omega0) + ...
///
/// where GD_target(r) = -(1/c)(sqrt(r^2+f^2) - f) varies with radius.
///
/// In a PROPAGATION-only achromat (see achromatic.hpp) each atom must hit BOTH the
/// base phase phi(r,omega0) AND the group delay -- a two-objective fit that leaves
/// a base-phase RESIDUAL (the library is discrete). PB decouples them: the rotation
/// sets the base phase EXACTLY for whatever atom is placed, so the atom is chosen
/// PURELY for its group delay (one objective). The realized cross transmission is
///
///     t_cross(r, omega) = a_cross(atom, omega) * exp(i * phi_geo(r)),
///     a_cross(atom, omega) = (t_x(omega) - t_y(omega)) / 2,   phi_geo = -2*hand*theta,
///
/// and theta(r) is set so arg(t_cross(r, omega0)) = phi(r, omega0) exactly. Hence
/// the base-phase RMS is ~0 by construction; the achromatic limit is ONLY the
/// library's group-delay coverage. Every atom shares one etch depth (rotated
/// rectangles of varying footprint) -> a single fabrication step, no grayscale.
///
/// References: Wang et al., Nat. Nanotechnol. 13, 227 (2018); Chen et al.,
/// Nat. Nanotechnol. 13, 220 (2018) (the geometric-phase achromat).

#include "celeris/core.hpp"
#include "celeris/design/achromatic.hpp"   // AchromaticFocalPoint
#include "celeris/materials/material.hpp"
#include "celeris/rcwa/rcwa2d.hpp"

#include <vector>

namespace celeris {

/// @brief One birefringent meta-atom characterized across the band by its spin-flip transmission.
///
/// a_cross(omega) = (t_x(omega) − t_y(omega)) / 2 is the amplitude a PB site
/// contributes once rotated. Birefringence (the HWP requirement) and dispersion
/// both come from the rectangle's aspect ratio and size.
struct DispersivePbAtom {
    double fill_x = 0;                 ///< Pillar width in x (fraction of the period).
    double fill_y = 0;                 ///< Pillar width in y (fraction of the period).
    double thickness_um = 0;           ///< Pillar height (µm).
    double phase0_rad = 0;             ///< arg(a_cross) at the center wavelength (rad, wrapped).
    double group_delay_fs = 0;         ///< d(arg a_cross)/d(omega), least-squares over the band (fs).
    double mean_amplitude = 0;         ///< Mean |a_cross| over the band (apodization proxy).
    double mean_conversion = 0;        ///< Mean |a_cross|² over the band (spin-flip efficiency).
    double retardance_center_deg = 0;  ///< wrap(arg t_x − arg t_y) at center (deg); ideal 180.
    std::vector<cdouble> cross;        ///< a_cross at each band wavelength (amplitude + phase).
};

/// @brief A dispersive birefringent library at ONE etch depth.
///
/// Atoms span a (fill_x, fill_y) grid so the group delay varies while the height
/// is fixed: a fabricable single-etch achromat. Input to
/// `design_pb_achromatic_metalens`.
struct DispersivePbLibrary {
    std::vector<double> wavelengths_um;  ///< Band samples (µm, ascending).
    double center_wavelength_um = 0;     ///< Base-phase reference wavelength (µm).
    int center_index = 0;                ///< Index of the center wavelength in `wavelengths_um`.
    double period_um = 0;                ///< Lattice period (µm).
    double thickness_um = 0;             ///< The single shared etch depth (µm).
    int n_fill = 0;                      ///< Grid side (up to n_fill² atoms before filtering).
    std::vector<DispersivePbAtom> atoms; ///< The usable (birefringent) atoms.
    double gd_min_fs = 0;                ///< Smallest group delay supplied (fs).
    double gd_max_fs = 0;                ///< Largest group delay supplied (fs).
};

/// @brief Build a dispersive birefringent library over a (fill_x, fill_y) grid at a single height.
///
/// For each rectangle, two RCWA solves (x- and y-polarized) at every band
/// wavelength give t_x(omega) and t_y(omega); the spin-flip amplitude a_cross
/// and its group delay (least-squares slope of the unwrapped phase vs angular
/// frequency) follow.
///
/// A non-birefringent atom (fill_x == fill_y, so t_x == t_y) has a_cross ≈ 0,
/// and its phase -- and therefore its group delay -- is noise that could match a
/// target while contributing no light. `min_amplitude` drops atoms whose mean
/// |a_cross| over the band is below it, so only usable PB atoms are kept.
/// @param pillar,background    Pillar and surrounding materials.
/// @param incident,substrate   Half-spaces above and below.
/// @param period_um            Lattice period (µm).
/// @param band_wavelengths_um  Band samples (µm); must be ascending.
/// @param center_wavelength_um Base-phase reference wavelength (µm).
/// @param fill_min,fill_max    Fill range, applied to both axes.
/// @param n_fills              Samples per axis.
/// @param thickness_um         The single pillar height (µm).
/// @param M                    Harmonic half-count in x and y.
/// @param min_amplitude        Minimum mean |a_cross| for an atom to be kept.
DispersivePbLibrary build_dispersive_pb_library(
    const Material& pillar, const Material& background, const Material& incident,
    const Material& substrate, double period_um,
    const std::vector<double>& band_wavelengths_um, double center_wavelength_um,
    double fill_min, double fill_max, int n_fills, double thickness_um, int M,
    double min_amplitude = 0.30);

/// @brief An achromatic PB focusing metalens: a rotated, group-delay-matched atom per site.
struct PbAchromaticDesign {
    int n_cells = 0;                      ///< Cells across the (square) aperture.
    double period_um = 0;                 ///< Lattice period (µm).
    double thickness_um = 0;              ///< The single etch depth (µm).
    int handedness = +1;                  ///< +1 for RCP illumination, −1 for LCP.
    std::vector<int> atom_index;          ///< Per site, index into `DispersivePbLibrary::atoms`.
    std::vector<double> rotation_rad;     ///< Per-site geometric-phase rotation theta (rad).
    std::vector<double> fill_x_map;       ///< Per-site chosen width in x, row-major.
    std::vector<double> fill_y_map;       ///< Per-site chosen width in y, row-major.
    double center_wavelength_um = 0;      ///< Base-phase reference wavelength (µm).
    double focal_length_um = 0;           ///< Target focal length (µm).
    double diameter_um = 0;               ///< Aperture width (µm).
    double rms_phase_error_deg = 0;       ///< Base-phase residual at center (deg); ~0, geometric is exact.
    double rms_group_delay_error_fs = 0;  ///< Group-delay residual (fs); the achromatic lever.
    double mean_amplitude = 0;            ///< Mean |a_cross| over the aperture.
    double mean_conversion = 0;           ///< Mean |a_cross|² over the aperture (efficiency cap).
    double required_gd_span_fs = 0;       ///< Group-delay span the design demanded (fs).
    double available_gd_span_fs = 0;      ///< Group-delay span the library supplies (fs).
    double gd_coverage = 0;               ///< available / required; ≥ 1 means the library suffices.
};

/// @brief Design an achromatic PB focusing metalens.
///
/// At each site: (1) pick the atom whose GROUP DELAY best matches the
/// radius-dependent target -- the only library constraint, since the geometric
/// phase handles the base phase; (2) set the rotation so the base phase is hit
/// exactly.
/// @param lib               Library from `build_dispersive_pb_library`.
/// @param focal_length_um   Target focal length (µm).
/// @param diameter_um       Aperture width (µm).
/// @param handedness        +1 for RCP illumination, −1 for LCP.
/// @param gd_weight         Group-delay objective weight. 0 ignores group delay
///                          and picks the best-amplitude atom everywhere -- a
///                          standard (chromatic) PB lens, the baseline.
/// @param amplitude_weight  Bias toward atoms with higher mean |a_cross|.
PbAchromaticDesign design_pb_achromatic_metalens(
    const DispersivePbLibrary& lib, double focal_length_um, double diameter_um,
    int handedness = +1, double gd_weight = 1.0, double amplitude_weight = 0.25);

/// @brief Focal length vs wavelength from the library's stored per-atom band response.
///
/// No new RCWA: at each band wavelength the aperture is propagated with each
/// site's realized cross transmission a_cross(omega)·exp(i·phi_geo) and the
/// on-axis peak is located. A flat curve means achromatic.
/// @param lib     The library the design was built from.
/// @param design  The design to verify.
/// @return One `AchromaticFocalPoint` (achromatic.hpp) per band wavelength.
std::vector<AchromaticFocalPoint> verify_pb_achromatic_focus(
    const DispersivePbLibrary& lib, const PbAchromaticDesign& design);

} // namespace celeris
