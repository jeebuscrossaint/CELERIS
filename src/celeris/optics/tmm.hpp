#pragma once
/// @file
/// @brief Transfer Matrix Method (TMM) for planar multilayer stacks.
///
/// Physics: each layer is represented by a 2x2 "characteristic matrix" relating
/// the tangential E and H fields at its two faces. Stacking layers = multiplying
/// their matrices. The stack's reflectance/transmittance falls out of the
/// product. This is the exact 1D analogue of the S-matrix layer recursion used
/// in RCWA, and the closed-form reference the RCWA solvers are validated against.
///
/// Reference: H. A. Macleod, "Thin-Film Optical Filters."

#include "celeris/core.hpp"
#include "celeris/materials/material.hpp"

#include <vector>

namespace celeris {

/// @brief One homogeneous film.
struct Layer {
    Material material;    ///< Film material.
    double thickness_um;  ///< Physical thickness (µm).
};

/// @brief Power and amplitude response of a multilayer stack.
struct TmmResult {
    double R;   ///< Power reflectance (0..1).
    double T;   ///< Power transmittance (0..1); R + T == 1 for lossless stacks.
    cdouble r;  ///< Complex reflection amplitude (carries phase).
    cdouble t;  ///< Complex transmission amplitude.
};

/// @brief Solve a planar multilayer stack for one polarization.
/// @param incident       Incidence medium (assumed lossless).
/// @param layers         Films in order, from the incident side down.
/// @param substrate      Semi-infinite exit medium.
/// @param wavelength_um  Vacuum wavelength (µm).
/// @param theta0_rad     Incidence angle from the normal.
/// @param pol            TE or TM.
/// @return Reflectance, transmittance, and complex amplitudes.
TmmResult solve_stack(const Material& incident,
                      const std::vector<Layer>& layers,
                      const Material& substrate,
                      double wavelength_um,
                      double theta0_rad,
                      Pol pol);

} // namespace celeris
