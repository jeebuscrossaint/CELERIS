#pragma once
// Through-focus analysis — the "focus" axis of optical analysis. Scans the
// on-axis intensity along z around the nominal focus to show how sharply the
// lens focuses and over what axial range (depth of focus). Standard tab in
// real optical-design tools.

#include "celeris/design/metalens.hpp"

#include <vector>

namespace celeris {

struct ThroughFocus {
    std::vector<float> z_um;        // axial position
    std::vector<float> intensity;  // on-axis intensity, normalized to its peak
    double z_peak_um;              // where intensity actually peaks
    double dof_um;                 // depth of focus (axial FWHM of the peak)
};

ThroughFocus analyze_through_focus(const MetalensDesign& lens,
                                   const UnitCellLibrary& lib,
                                   double focal_length_um, double wavelength_um,
                                   double diameter_um);

} // namespace celeris
