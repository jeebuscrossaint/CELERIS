#pragma once
/// @file
/// @brief 1D binary (lamellar) gratings: periodic in x, invariant in y.
///
/// One period consists of a "ridge" of width fill*period and a "groove" filling
/// the rest. This is the simplest periodic structure and the standard vehicle
/// for validating an RCWA implementation before moving to 2D metalens pillars.

#include "celeris/core.hpp"
#include "celeris/materials/material.hpp"

namespace celeris {

/// @brief A single binary grating layer with its own period.
struct BinaryGrating1D {
    Material ridge;       ///< High-index region material.
    Material groove;      ///< Material filling the rest of the period.
    double period_um;     ///< Grating period Λ (µm).
    double fill;          ///< Ridge width fraction, 0..1.
    double thickness_um;  ///< Layer thickness along propagation, z (µm).

    /// @brief h-th Fourier coefficient ε_h of the periodic permittivity ε(x).
    ///
    /// The ridge is centered at x=0, so the coefficients are real multiples of
    /// (ε_ridge − ε_groove):
    ///   ε_0 = fill·ε_r + (1−fill)·ε_g,
    ///   ε_h = (ε_r − ε_g) · sin(π h fill) / (π h) for h ≠ 0.
    /// @param h              Fourier order.
    /// @param wavelength_um  Vacuum wavelength (µm), for the material dispersion.
    cdouble eps_fourier(int h, double wavelength_um) const;

    /// @brief h-th Fourier coefficient of the reciprocal permittivity 1/ε(x).
    ///
    /// TM polarization needs this (Li's inverse rule): the Toeplitz matrix built
    /// from these coefficients, not the matrix inverse of [ε], is what makes TM
    /// converge. Same binary form as eps_fourier() with ε → 1/ε.
    /// @param h              Fourier order.
    /// @param wavelength_um  Vacuum wavelength (µm).
    cdouble eps_inv_fourier(int h, double wavelength_um) const;
};

/// @brief One layer of a 1D multilayer stack (no period of its own).
///
/// The same binary cross-section as BinaryGrating1D but without a period:
/// every layer in a 1D RCWA stack must share a single period (set on
/// Rcwa1DStack) so their Fourier orders align. A homogeneous (unpatterned)
/// layer is just ridge == groove.
struct GratingLayer1D {
    Material ridge;       ///< High-index region material.
    Material groove;      ///< Material filling the rest of the period.
    double fill;          ///< Ridge width fraction, 0..1.
    double thickness_um;  ///< Layer thickness (µm).

    /// @brief h-th Fourier coefficient of ε(x); see BinaryGrating1D::eps_fourier().
    cdouble eps_fourier(int h, double wavelength_um) const;
    /// @brief h-th Fourier coefficient of 1/ε(x); see BinaryGrating1D::eps_inv_fourier().
    cdouble eps_inv_fourier(int h, double wavelength_um) const;

    /// @brief A uniform layer of a single material.
    /// @param m             The layer material.
    /// @param thickness_um  Layer thickness (µm).
    static GratingLayer1D homogeneous(Material m, double thickness_um) {
        return GratingLayer1D{m, m, 1.0, thickness_um};
    }
};

} // namespace celeris
