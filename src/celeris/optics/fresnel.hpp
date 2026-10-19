#pragma once
/// @file
/// @brief Fresnel equations for a single planar interface between two media.
///
/// The atomic building block of thin-film and grating physics. The Transfer
/// Matrix Method (tmm.hpp) generalizes this to arbitrary multilayer stacks.
///
/// Convention: light travels from medium 1 into medium 2; theta1 is the
/// incidence angle in medium 1 measured from the surface normal (radians).
/// Indices are complex (n + i*k) so absorbing media work without special cases.

#include "celeris/core.hpp"

#include <cmath>

namespace celeris {

/// @brief Complex Fresnel amplitude coefficients for both polarizations.
struct FresnelCoeffs {
    cdouble rs;  ///< Reflection amplitude, TE (s).
    cdouble rp;  ///< Reflection amplitude, TM (p).
    cdouble ts;  ///< Transmission amplitude, TE (s).
    cdouble tp;  ///< Transmission amplitude, TM (p).
};

/// @brief cos(θ2) from Snell's law, kept complex.
///
/// Staying complex means total internal reflection and absorbing media need no
/// special-casing.
/// @param n1      Complex index of the incidence medium.
/// @param n2      Complex index of the transmission medium.
/// @param theta1  Incidence angle in medium 1 (rad).
inline cdouble snell_cos_theta2(cdouble n1, cdouble n2, double theta1) {
    cdouble sin_t1 = std::sin(theta1);
    cdouble sin_t2 = (n1 / n2) * sin_t1;
    return std::sqrt(cdouble{1.0, 0.0} - sin_t2 * sin_t2);
}

/// @brief Fresnel coefficients of the interface from medium 1 into medium 2.
/// @param n1      Complex index of the incidence medium.
/// @param n2      Complex index of the transmission medium.
/// @param theta1  Incidence angle in medium 1 (rad).
/// @return Reflection and transmission amplitudes for TE and TM.
inline FresnelCoeffs fresnel(cdouble n1, cdouble n2, double theta1) {
    cdouble cos_t1 = std::cos(theta1);
    cdouble cos_t2 = snell_cos_theta2(n1, n2, theta1);

    cdouble n1c1 = n1 * cos_t1;
    cdouble n2c2 = n2 * cos_t2;
    cdouble n2c1 = n2 * cos_t1;
    cdouble n1c2 = n1 * cos_t2;

    FresnelCoeffs f;
    f.rs = (n1c1 - n2c2) / (n1c1 + n2c2);
    f.rp = (n2c1 - n1c2) / (n2c1 + n1c2);
    f.ts = (2.0 * n1c1) / (n1c1 + n2c2);
    f.tp = (2.0 * n1c1) / (n2c1 + n1c2);
    return f;
}

/// @brief Power reflectance R = |r|² from a reflection amplitude.
inline double reflectance(cdouble r) { return std::norm(r); }

} // namespace celeris
