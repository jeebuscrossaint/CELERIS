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
std::vector<FieldPoint> analyze_field_of_view(const MetalensDesign& lens,
                                              const UnitCellLibrary& lib,
                                              double focal_length_um,
                                              double wavelength_um,
                                              double diameter_um,
                                              const std::vector<double>& angles_deg);

} // namespace celeris
