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

} // namespace celeris
