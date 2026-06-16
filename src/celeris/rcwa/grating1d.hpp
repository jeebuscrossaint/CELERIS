#pragma once
// A 1D binary (lamellar) grating: a layer periodic in x, invariant in y.
// One period consists of a "ridge" of width fill*period and a "groove" filling
// the rest. This is the simplest periodic structure and the standard vehicle
// for validating an RCWA implementation before moving to 2D metalens pillars.

#include "celeris/core.hpp"
#include "celeris/materials/material.hpp"

namespace celeris {

struct BinaryGrating1D {
    Material ridge;       // high-index region material
    Material groove;      // the rest of the period
    double period_um;     // grating period Λ
    double fill;          // ridge width fraction, 0..1
    double thickness_um;  // layer thickness along propagation (z)

    // h-th Fourier coefficient ε_h of the periodic permittivity ε(x) at the
    // given vacuum wavelength. The ridge is taken centered at x=0 so the
    // coefficients are real multiples of (ε_ridge − ε_groove):
    //   ε_0 = fill·ε_r + (1−fill)·ε_g
    //   ε_h = (ε_r − ε_g) · sin(π h fill) / (π h),  h ≠ 0
    cdouble eps_fourier(int h, double wavelength_um) const;

    // h-th Fourier coefficient of the RECIPROCAL permittivity 1/ε(x). TM
    // polarization needs this (Li's inverse rule): the Toeplitz matrix built
    // from these coefficients, not the matrix inverse of [ε], is what makes
    // TM converge. Same binary form as eps_fourier with ε → 1/ε.
    cdouble eps_inv_fourier(int h, double wavelength_um) const;
};

// One layer of a multilayer stack: the same binary cross-section as
// BinaryGrating1D but WITHOUT a period — every layer in a 1D RCWA stack must
// share a single period (set on the stack) so their Fourier orders align.
// A homogeneous (unpatterned) layer is just ridge == groove.
struct GratingLayer1D {
    Material ridge;
    Material groove;
    double fill;          // ridge width fraction, 0..1
    double thickness_um;

    cdouble eps_fourier(int h, double wavelength_um) const;
    cdouble eps_inv_fourier(int h, double wavelength_um) const;

    // Convenience: a uniform layer of a single material.
    static GratingLayer1D homogeneous(Material m, double thickness_um) {
        return GratingLayer1D{m, m, 1.0, thickness_um};
    }
};

} // namespace celeris
