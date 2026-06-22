#pragma once
// Field-of-view analysis — how the focus degrades for off-axis (obliquely
// incident) light. A single hyperbolic-phase metalens is optimized for on-axis
// focusing; at a field angle the spot shifts and aberrates (coma), so the
// usable field of view is limited. This sweeps incidence angle and reports the
// relative Strehl and spot shift, then the FOV where quality stays acceptable.

#include "celeris/design/metalens.hpp"

#include <vector>

namespace celeris {

struct FieldPoint {
    double angle_deg;      // incidence (field) angle
    double rel_strehl;     // peak intensity relative to the on-axis focus
    double spot_shift_um;  // lateral focal-spot displacement
};

// Sweep field angles and report off-axis focus quality. Incident tilt is in the
// x–z plane. The first entry (0°) is the on-axis reference (rel_strehl = 1).
// Models the aperture STOP at the lens (full-aperture tilted-plane-wave
// illumination): this is the right model for a normal lens's usable field, but
// it does NOT reveal the wide-FOV benefit of a quadratic-phase lens — that
// requires an offset stop (see analyze_wide_fov).
std::vector<FieldPoint> analyze_field_of_view(const MetalensDesign& lens,
                                              const UnitCellLibrary& lib,
                                              double focal_length_um,
                                              double wavelength_um,
                                              double diameter_um,
                                              const std::vector<double>& angles_deg);

// Wide-FOV analysis with an aperture STOP of diameter stop_diameter_um placed a
// distance stop_distance_um IN FRONT of the metasurface (incidence side). A field
// angle theta then illuminates a DECENTERED patch of the (larger) lens — a disk of
// radius stop/2 centered at x = stop_distance*tan(theta) — carrying the tilt phase
// exp(i k sin(theta) x). This is the configuration in which a quadratic-phase lens
// is wide-FOV: over the small, low-NA stopped patch the tilted quadratic phase is a
// perfectly RECENTERED parabola (a sharp focus that just translates with angle),
// whereas a hyperbolic lens shows coma. The stop sets the resolution (~lambda*f/
// stop), the lens diameter must be large enough to catch the walk-off. rel_strehl
// is normalized to each lens's own on-axis (theta=0) peak.
std::vector<FieldPoint> analyze_wide_fov(const MetalensDesign& lens,
                                         const UnitCellLibrary& lib,
                                         double focal_length_um,
                                         double wavelength_um,
                                         double lens_diameter_um,
                                         double stop_diameter_um,
                                         double stop_distance_um,
                                         const std::vector<double>& angles_deg);

// One node of the full 2D field-resolved grid: the off-axis focus for a plane
// wave tilted by (theta_x, theta_y). Reports the relative Strehl (vs the on-axis
// peak) and the tangential/sagittal spot FWHM, so a grid of these is the full-
// field PSF/Strehl map (not just the center-row spot-vs-field of the FieldPoint
// sweep above). cx/cy are the chief-ray landing on the focal plane.
struct FieldGridPoint {
    double theta_x_deg, theta_y_deg;
    double rel_strehl;
    double fwhm_x_um, fwhm_y_um;
    double cx_um, cy_um;
};

struct FieldGrid {
    int n;                  // (2*n_half+1) angles per axis; n = side length
    double max_angle_deg;   // extent of each axis
    std::vector<FieldGridPoint> points;  // row-major n x n, theta_y outer / theta_x inner
};

// Sweep an n_half-symmetric grid of field angles (theta_x, theta_y) over
// +/- max_angle_deg and report the off-axis focus quality at each. Each point
// illuminates the full aperture with the tilted plane wave exp(i k (sin θx·x +
// sin θy·y)) and propagates to a window centered on the chief ray; the PSF is
// sampled on psf_n x psf_n. This is the full-field map: the on-axis (0,0) node
// is the rel_strehl=1 reference. side length n = 2*n_half + 1.
FieldGrid analyze_field_grid(const MetalensDesign& lens,
                             const UnitCellLibrary& lib,
                             double focal_length_um, double wavelength_um,
                             double diameter_um, double max_angle_deg,
                             int n_half, int psf_n);

} // namespace celeris
