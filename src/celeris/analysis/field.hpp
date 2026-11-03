#pragma once
/// @file
/// @brief Field-of-view analysis: how the focus degrades for off-axis light.
///
/// A single hyperbolic-phase metalens is optimized for on-axis
/// focusing; at a field angle the spot shifts and aberrates (coma), so the
/// usable field of view is limited. This sweeps incidence angle and reports the
/// relative Strehl and spot shift, then the FOV where quality stays acceptable.

#include "celeris/design/metalens.hpp"

#include <vector>

namespace celeris {

/// @brief Off-axis focus quality at one field angle (x–z plane tilt).
struct FieldPoint {
    double angle_deg;      ///< Incidence (field) angle (degrees).
    double rel_strehl;     ///< Peak intensity relative to the on-axis focus.
    double spot_shift_um;  ///< Lateral focal-spot displacement (µm).
};

/// @brief Sweep field angles with the aperture stop at the lens.
///
/// Incident tilt is in the x–z plane. The first entry (0°) is the on-axis
/// reference (rel_strehl = 1).
/// Models the aperture STOP at the lens (full-aperture tilted-plane-wave
/// illumination): this is the right model for a normal lens's usable field, but
/// it does NOT reveal the wide-FOV benefit of a quadratic-phase lens — that
/// requires an offset stop (see analyze_wide_fov()).
/// @param lens            The finished design.
/// @param lib             The library the design was built from.
/// @param focal_length_um Focal length (µm).
/// @param wavelength_um   Vacuum wavelength (µm).
/// @param diameter_um     Circular aperture diameter (µm).
/// @param angles_deg      Field angles to evaluate (degrees).
/// @return One FieldPoint per angle; rel_strehl is normalized to the 0° peak.
std::vector<FieldPoint> analyze_field_of_view(const MetalensDesign& lens,
                                              const UnitCellLibrary& lib,
                                              double focal_length_um,
                                              double wavelength_um,
                                              double diameter_um,
                                              const std::vector<double>& angles_deg);

/// @brief Wide-FOV analysis with an aperture stop placed in front of the lens.
///
/// The stop has diameter stop_diameter_um and sits a distance stop_distance_um IN
/// FRONT of the metasurface (incidence side). A field angle theta then illuminates
/// a DECENTERED patch of the (larger) lens — a disk of radius stop/2 centered at
/// x = stop_distance*tan(theta) — carrying the tilt phase exp(i k sin(theta) x).
/// This is the configuration in which a quadratic-phase lens is wide-FOV: over the
/// small, low-NA stopped patch the tilted quadratic phase is a perfectly
/// RECENTERED parabola (a sharp focus that just translates with angle), whereas a
/// hyperbolic lens shows coma. The stop sets the resolution (~lambda*f/stop); the
/// lens diameter must be large enough to catch the walk-off.
/// @param lens              The finished design.
/// @param lib               The library the design was built from.
/// @param focal_length_um   Focal length (µm).
/// @param wavelength_um     Vacuum wavelength (µm).
/// @param lens_diameter_um  Full lens diameter (µm).
/// @param stop_diameter_um  Aperture-stop diameter (µm).
/// @param stop_distance_um  Stop-to-metasurface distance (µm).
/// @param angles_deg        Field angles to evaluate (degrees).
/// @return One FieldPoint per angle; rel_strehl is normalized to this lens's own
///         on-axis (theta = 0) peak.
std::vector<FieldPoint> analyze_wide_fov(const MetalensDesign& lens,
                                         const UnitCellLibrary& lib,
                                         double focal_length_um,
                                         double wavelength_um,
                                         double lens_diameter_um,
                                         double stop_diameter_um,
                                         double stop_distance_um,
                                         const std::vector<double>& angles_deg);

/// @brief One node of the full 2D field-resolved grid.
///
/// The off-axis focus for a plane wave tilted by (theta_x, theta_y). Reports the
/// relative Strehl (vs the on-axis peak) and the tangential/sagittal spot FWHM,
/// so a grid of these is the full-field PSF/Strehl map (not just the center-row
/// spot-vs-field of the FieldPoint sweep).
struct FieldGridPoint {
    double theta_x_deg;  ///< Field angle about the y axis (tilt in x–z), degrees.
    double theta_y_deg;  ///< Field angle about the x axis (tilt in y–z), degrees.
    double rel_strehl;   ///< Peak intensity relative to the on-axis (0,0) peak.
    double fwhm_x_um;    ///< Spot FWHM along x (µm).
    double fwhm_y_um;    ///< Spot FWHM along y (µm).
    double cx_um;        ///< Chief-ray landing on the focal plane, x (µm).
    double cy_um;        ///< Chief-ray landing on the focal plane, y (µm).
};

/// @brief A square grid of FieldGridPoint covering ±max_angle_deg on both axes.
struct FieldGrid {
    int n;                  ///< Side length = 2*n_half + 1 angles per axis.
    double max_angle_deg;   ///< Extent of each axis (degrees).
    std::vector<FieldGridPoint> points;  ///< Row-major n x n: theta_y outer, theta_x inner.
};

/// @brief Full-field map: off-axis focus quality on a 2D grid of field angles.
///
/// Sweeps an n_half-symmetric grid of field angles (theta_x, theta_y) over
/// ±max_angle_deg. Each point illuminates the full aperture with the tilted plane
/// wave exp(i k (sin θx·x + sin θy·y)) and propagates to a window centered on the
/// chief ray; the PSF is sampled on psf_n x psf_n. The on-axis (0,0) node is the
/// rel_strehl = 1 reference.
/// @param lens            The finished design.
/// @param lib             The library the design was built from.
/// @param focal_length_um Focal length (µm).
/// @param wavelength_um   Vacuum wavelength (µm).
/// @param diameter_um     Circular aperture diameter (µm).
/// @param max_angle_deg   Largest field angle on each axis (degrees).
/// @param n_half          Angles per half-axis; the side length is 2*n_half + 1.
/// @param psf_n           PSF samples per side at each field point.
/// @return The grid, row-major as documented in FieldGrid.
FieldGrid analyze_field_grid(const MetalensDesign& lens,
                             const UnitCellLibrary& lib,
                             double focal_length_um, double wavelength_um,
                             double diameter_um, double max_angle_deg,
                             int n_half, int psf_n);

} // namespace celeris
