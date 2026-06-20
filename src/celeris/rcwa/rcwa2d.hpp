#pragma once
// 2D RCWA (Fourier Modal Method) for biperiodic structures — the real metalens
// case: nanopillars on a grid, periodic in BOTH x and y.
//
// Unlike 1D, the polarizations couple, so the per-layer eigenproblem is fully
// vectorial: a 2N×2N operator P·Q acting on the tangential E-field harmonics
// [Ex;Ey], where N = (2Mx+1)(2My+1) is the number of 2D Fourier orders. Layer
// stacking reuses the same scattering-matrix machinery as the 1D solver.
//
// This is the "basic" Fourier factorization (matrix inverse of the convolution
// matrix). It is correct and converges; the improved 2D factorization (Li's
// normal-vector method, for faster convergence on high-contrast/metal cells) is
// a later milestone. References: Moharam & Gaylord (1995); Rumpf, FMM notes.

#include "celeris/core.hpp"
#include "celeris/materials/material.hpp"

#include <vector>

namespace celeris {

// One layer: a rectangular pillar of `pillar` material (fractional size
// fill_x × fill_y, centered in the cell) embedded in `background`. A
// homogeneous layer is pillar == background (or fill = 1). The separable
// rectangle has a closed-form 2D Fourier series, and setting fill_y = 1 makes
// the cell invariant in y — which is how we cross-check 2D against the 1D
// solver.
struct RectCell2D {
    Material pillar;
    Material background;
    double fill_x;
    double fill_y;
    double thickness_um;

    // 2D Fourier coefficient ε_{p,q} of the permittivity.
    cdouble eps_fourier(int p, int q, double wavelength_um) const;

    // 2D Fourier coefficient of the RECIPROCAL permittivity (1/ε)_{p,q}, used by
    // Li's inverse-rule factorization for the field component normal to a
    // material interface (what makes high-contrast cells converge).
    cdouble inv_eps_fourier(int p, int q, double wavelength_um) const;

    static RectCell2D homogeneous(Material m, double thickness_um) {
        return RectCell2D{m, m, 1.0, 1.0, thickness_um};
    }
};

struct Rcwa2DStack {
    double period_x_um;
    double period_y_um;
    std::vector<RectCell2D> layers;  // top (incident side) to bottom
};

struct Rcwa2DResult {
    double R;        // total reflected efficiency (all orders)
    double T;        // total transmitted efficiency (all orders)
    double sum_de;   // R + T; == 1 for a lossless stack
    double de_t0;    // zeroth-order transmitted efficiency
    double de_r0;    // zeroth-order reflected efficiency
    cdouble tx0;     // complex zeroth-order transmitted Ex amplitude
    cdouble ty0;     // complex zeroth-order transmitted Ey amplitude
                     // (arg() of these is the phase delay a metalens pillar
                     //  imparts — the core quantity for unit-cell libraries)
};

// Solve a biperiodic stack. Incidence is set by polar/azimuth angles and the
// tangential incident E-field at order 0 (Ex0, Ey0); e.g. (0,1) = E along y,
// (1,0) = E along x. Mx/My are the harmonic half-counts in each direction.
Rcwa2DResult solve_rcwa_2d(const Material& incident,
                           const Rcwa2DStack& stack,
                           const Material& substrate,
                           double wavelength_um,
                           double theta_rad,
                           double phi_rad,
                           cdouble Ex0,
                           cdouble Ey0,
                           int Mx,
                           int My);

} // namespace celeris
