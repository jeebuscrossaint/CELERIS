#pragma once
/// @file
/// @brief Achromatic (broadband) metalens design: dispersion engineering.
///
/// A standard phase-profile metalens fixes the phase only at the design
/// wavelength, so its focal length drifts as f(lambda) ~ f0*lambda0/lambda (the
/// #1 limitation of metalenses). To focus a whole BAND at one plane, the phase a
/// meta-atom imparts must satisfy, at every frequency omega in the band,
///
///     phi(r, omega) = -(omega/c) * (sqrt(r^2 + f^2) - f) + phi_ref(omega).
///
/// Expanding around the center frequency omega0, two terms must be controlled per
/// site:
///   - phi(r, omega0)          -- the BASE phase (matched mod 2*pi, as usual);
///   - d phi / d omega |omega0 -- the GROUP DELAY, which MUST vary with radius
///                                (maximum at center, decreasing to the edge).
///
/// The required group-delay SPAN, (1/c)(sqrt(R^2+f^2)-f), grows with aperture and
/// bandwidth. The catch: to hit an ARBITRARY (phase, group-delay) pair per site,
/// the atom library must cover the 2-D (phase, GD) plane -- a single geometric
/// degree of freedom (e.g. fill only) traces just a 1-D curve and cannot set both
/// independently. So achromatic design needs a richer library (here: fill x
/// height) and a TWO-objective atom selection.
///
/// References: Wang et al., Nat. Nanotechnol. 13, 227 (2018); Chen et al.,
/// Nat. Nanotechnol. 13, 220 (2018); Shrestha et al., Light Sci. Appl. 7, 85 (2018).

#include "celeris/design/metalens.hpp"
#include "celeris/materials/material.hpp"
#include "celeris/rcwa/rcwa2d.hpp"

#include <vector>

namespace celeris {

/// @brief One meta-atom characterized across a wavelength band.
///
/// Its phase and amplitude at every band sample, plus the derived center phase
/// and group delay. The full in-plane geometry is carried so a single-etch
/// (shape-diverse) library can report the mix of shapes it spans the
/// (phase, group-delay) plane with.
struct DispersiveAtom {
    double fill = 0;                 ///< Representative fill (== fill_x for a square atom).
    double thickness_um = 0;         ///< Pillar height (µm).
    MetaShape shape = MetaShape::Rectangle;  ///< In-plane cross-section.
    double fill_x = 0;               ///< Width in x (fraction of the period).
    double fill_y = 0;               ///< Width in y (fraction of the period).
    double shape_param = 0.5;        ///< Cross: arm width / fill; Ring: inner / outer radius.
    double phase0_rad = 0;           ///< Transmission phase at the center wavelength (rad, wrapped).
    double group_delay_fs = 0;       ///< d(phi)/d(omega) at center, least-squares over the band (fs).
    double mean_amplitude = 0;       ///< Mean |t| over the band (efficiency proxy).
    std::vector<double> phase_rad;   ///< phi at each band sample (rad).
    std::vector<double> amplitude;   ///< |t| at each band sample.
};

/// @brief One meta-atom geometry to characterize.
///
/// The dispersive library is a list of these solved across the band: a
/// fill × height grid (the default), or a single height with varied shapes (the
/// single-etch path).
struct MetaAtomSpec {
    MetaShape shape = MetaShape::Rectangle;  ///< In-plane cross-section.
    double fill_x = 0.5;       ///< Width in x (fraction of the period).
    double fill_y = 0.5;       ///< Width in y (fraction of the period).
    double thickness_um = 0.6; ///< Pillar height (µm).
    double shape_param = 0.5;  ///< Cross: arm width / fill; Ring: inner / outer radius.
};

/// @brief A meta-atom library characterized over a wavelength band.
///
/// Input to the two-objective achromatic selection (`design_achromatic_metalens`).
/// Atoms span a fill × height grid so the (phase, group-delay) plane is covered,
/// not just a 1-D curve.
struct DispersiveLibrary {
    std::vector<double> wavelengths_um;  ///< Band samples (µm, ascending).
    double center_wavelength_um = 0;     ///< Base-phase reference wavelength (µm).
    int center_index = 0;                ///< Index of the center wavelength in `wavelengths_um`.
    double period_um = 0;                ///< Lattice period (µm).
    int n_fill = 0;                      ///< Fill samples in the grid.
    int n_height = 0;                    ///< Height samples in the grid (atoms = n_fill·n_height).
    std::vector<DispersiveAtom> atoms;   ///< One per characterized geometry.
    double gd_min_fs = 0;                ///< Smallest group delay supplied (fs).
    double gd_max_fs = 0;                ///< Largest group delay supplied (fs).
};

/// @brief Build a dispersive library over a fill × height grid.
///
/// For each (fill, height) atom the 2D RCWA is solved at EVERY band wavelength,
/// so phi(lambda) and |t|(lambda) capture the meta-atom's true material and
/// waveguide dispersion. The group delay is the least-squares slope of the
/// band-unwrapped phase vs angular frequency. `n_heights == 1` gives a
/// single-etch (one-degree-of-freedom) library.
/// @param pillar,background    Pillar and surrounding materials.
/// @param incident,substrate   Half-spaces above and below.
/// @param period_um            Lattice period (µm).
/// @param band_wavelengths_um  Band samples (µm); must be ascending.
/// @param center_wavelength_um Base-phase reference wavelength (µm).
/// @param fill_min,fill_max    Fill range.
/// @param n_fills              Fill samples.
/// @param thick_lo,thick_hi    Height range (µm).
/// @param n_heights            Height samples.
/// @param M                    Harmonic half-count in x and y.
DispersiveLibrary build_dispersive_library(
    const Material& pillar, const Material& background, const Material& incident,
    const Material& substrate, double period_um,
    const std::vector<double>& band_wavelengths_um, double center_wavelength_um,
    double fill_min, double fill_max, int n_fills, double thick_lo, double thick_hi,
    int n_heights, int M);

/// @brief Characterize an arbitrary list of meta-atom geometries across the band.
///
/// The general builder: `build_dispersive_library` and `build_single_etch_library`
/// are thin wrappers over it. `n_fill`/`n_height` on the returned library are
/// informational only here (the grid wrapper sets them).
/// @param pillar,background    Pillar and surrounding materials.
/// @param incident,substrate   Half-spaces above and below.
/// @param period_um            Lattice period (µm).
/// @param band_wavelengths_um  Band samples (µm); must be ascending.
/// @param center_wavelength_um Base-phase reference wavelength (µm).
/// @param specs                The geometries to solve.
/// @param M                    Harmonic half-count in x and y.
DispersiveLibrary build_dispersive_library_from_specs(
    const Material& pillar, const Material& background, const Material& incident,
    const Material& substrate, double period_um,
    const std::vector<double>& band_wavelengths_um, double center_wavelength_um,
    const std::vector<MetaAtomSpec>& specs, int M);

/// @brief Build a SINGLE-ETCH dispersive library: every atom shares one pillar height.
///
/// The (phase, group-delay) plane is spanned by varying the in-plane SHAPE
/// instead of the depth: a square + circle fill sweep, plus cross
/// (fill × arm-width) and ring (fill × inner-radius) families. One lithography
/// step, one etch depth, no grayscale. The tradeoff vs the fill × height library
/// is a SMALLER group-delay span (taller pillars accumulate more delay), so it is
/// achromatic over a smaller aperture × bandwidth -- reported via
/// `AchromaticDesign::gd_coverage`.
/// @param pillar,background    Pillar and surrounding materials.
/// @param incident,substrate   Half-spaces above and below.
/// @param period_um            Lattice period (µm).
/// @param band_wavelengths_um  Band samples (µm); must be ascending.
/// @param center_wavelength_um Base-phase reference wavelength (µm).
/// @param fill_min,fill_max    Fill range.
/// @param n_fills              Fill samples per family.
/// @param thickness_um         The single pillar height (µm).
/// @param M                    Harmonic half-count in x and y.
DispersiveLibrary build_single_etch_library(
    const Material& pillar, const Material& background, const Material& incident,
    const Material& substrate, double period_um,
    const std::vector<double>& band_wavelengths_um, double center_wavelength_um,
    double fill_min, double fill_max, int n_fills, double thickness_um, int M);

/// @brief An achromatic focusing metalens: one band-characterized atom per site.
struct AchromaticDesign {
    int n_cells = 0;                      ///< Cells across the (square) aperture.
    double period_um = 0;                 ///< Lattice period (µm).
    std::vector<int> atom_index;          ///< Per site, index into `DispersiveLibrary::atoms`.
    std::vector<double> fill_map;         ///< Per-site chosen pillar fill, row-major.
    std::vector<double> thickness_map;    ///< Per-site chosen pillar height (µm), row-major.
    double center_wavelength_um = 0;      ///< Base-phase reference wavelength (µm).
    double focal_length_um = 0;           ///< Target focal length (µm).
    double diameter_um = 0;               ///< Aperture width (µm).
    double rms_phase_error_deg = 0;       ///< Base-phase residual at the center wavelength (deg).
    double rms_group_delay_error_fs = 0;  ///< Group-delay residual (fs); the achromatic lever.
    double mean_amplitude = 0;            ///< Mean |t| over the aperture.
    double required_gd_span_fs = 0;       ///< Group-delay span the design demanded (fs).
    double available_gd_span_fs = 0;      ///< Group-delay span the library supplies (fs).
    double gd_coverage = 0;               ///< available / required; ≥ 1 means the library suffices.
    bool single_height = true;            ///< True if every chosen atom shares one height.
    double min_height_um = 0;             ///< Lowest chosen pillar height (µm).
    double max_height_um = 0;             ///< Tallest chosen pillar height (µm).
};

/// @brief Design a focusing metalens by two-objective atom selection.
///
/// At each site, choose the atom that best matches BOTH the base focusing phase
/// (mod 2π) and the radius-dependent group delay. The group-delay error is
/// converted to the phase error it causes at the band EDGE (its physical effect
/// on focus there), so `gd_weight` ≈ 1 balances center-wavelength focus against
/// band-edge focus.
/// @param lib               Library from one of the builders above.
/// @param focal_length_um   Target focal length (µm).
/// @param diameter_um       Aperture width (µm).
/// @param gd_weight         Group-delay objective weight. 0 reproduces a standard
///                          single-wavelength (dispersion-blind) design from the
///                          same library -- the baseline for comparison.
/// @param amplitude_weight  Bias toward higher transmission.
AchromaticDesign design_achromatic_metalens(const DispersiveLibrary& lib,
                                            double focal_length_um,
                                            double diameter_um,
                                            double gd_weight = 1.0,
                                            double amplitude_weight = 0.25);

/// @brief One sample of a design's focal length vs wavelength.
struct AchromaticFocalPoint {
    double wavelength_um;    ///< Band wavelength (µm).
    double focal_length_um;  ///< Location of the on-axis intensity peak (µm).
    double rel_peak;         ///< On-axis peak relative to the center-wavelength design point.
};

/// @brief Focal length vs wavelength from the library's stored per-atom band response.
///
/// No new RCWA solves: at each band wavelength the aperture is propagated with
/// each site's atom transmission at that wavelength and the on-axis intensity
/// peak is located. This uses each atom's true dispersion and works for per-site
/// varying heights.
/// @param lib     The library the design was built from.
/// @param design  The design to verify.
/// @return One point per band wavelength.
std::vector<AchromaticFocalPoint> verify_achromatic_focus(
    const DispersiveLibrary& lib, const AchromaticDesign& design);

/// @brief Adapt to the plain `MetalensDesign` the GDS writer and analysis battery consume.
///
/// Only in-plane footprints carry over: a multi-height design's etch depths are
/// not encoded in a single GDS layer (see `AchromaticDesign::single_height`).
MetalensDesign to_metalens_design(const AchromaticDesign& d);

} // namespace celeris
