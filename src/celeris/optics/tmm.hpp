#pragma once
// Transfer Matrix Method (TMM) for planar multilayer stacks.
//
// Physics: each layer is represented by a 2x2 "characteristic matrix" relating
// the tangential E and H fields at its two faces. Stacking layers = multiplying
// their matrices. The stack's reflectance/transmittance falls out of the
// product. This is the exact 1D analogue of the S-matrix layer recursion used
// later in RCWA — so getting it right here is direct practice for the engine.
//
// Reference: H. A. Macleod, "Thin-Film Optical Filters."

#include "celeris/core.hpp"
#include "celeris/materials/material.hpp"

#include <vector>

namespace celeris {

// One homogeneous film: a material and a physical thickness (micrometers).
struct Layer {
    Material material;
    double thickness_um;
};

struct TmmResult {
    double R;   // power reflectance  (0..1)
    double T;   // power transmittance (0..1); R + T == 1 for lossless stacks
    cdouble r;  // complex reflection amplitude (carries phase)
    cdouble t;  // complex transmission amplitude
};

// Solve a stack illuminated from `incident` (assumed lossless) at angle
// `theta0_rad` from normal, through `layers` in order, into semi-infinite
// `substrate`, at vacuum wavelength `wavelength_um`, for one polarization.
TmmResult solve_stack(const Material& incident,
                      const std::vector<Layer>& layers,
                      const Material& substrate,
                      double wavelength_um,
                      double theta0_rad,
                      Pol pol);

} // namespace celeris
