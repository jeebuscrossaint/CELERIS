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

    // Longitudinal caustic: |E|^2 on an (nz x nx) grid in the y=0 plane (row =
    // z, col = x), normalized to its peak. Empty if it could not be computed.
    // Shows the focusing cone -- the standard "x-z" focus view.
    int caustic_nx = 0, caustic_nz = 0;
    double caustic_xmin = 0, caustic_xmax = 0;  // lateral extent (um)
    double caustic_zmin = 0, caustic_zmax = 0;  // axial extent (um)
    std::vector<float> caustic;
};

ThroughFocus analyze_through_focus(const MetalensDesign& lens,
                                   const UnitCellLibrary& lib,
                                   double focal_length_um, double wavelength_um,
                                   double diameter_um);

} // namespace celeris
