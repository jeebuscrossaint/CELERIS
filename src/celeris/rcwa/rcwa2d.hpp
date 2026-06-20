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

// Meta-atom cross-section shape, centered in the cell and bounded by the
// fill_x × fill_y box. Rectangle has a closed-form separable Fourier series (the
// fast, analytic path); the others are non-separable and go through the
// numerical Li/Fast-Fourier factorization (a sampled-grid solve).
enum class MetaShape {
    Rectangle,  // |x|<fill_x/2 and |y|<fill_y/2  (the original analytic cell)
    Ellipse,    // (x/(fill_x/2))^2 + (y/(fill_y/2))^2 < 1  (circle if fill_x==fill_y)
    Cross,      // plus/cross: two centered bars; arm width = shape_param * fill
    Ring        // annulus of pillar material; inner radius = shape_param * outer
};

// One layer: a pillar of `pillar` material (shape inside the fractional box
// fill_x × fill_y, centered in the cell) embedded in `background`. A homogeneous
// layer is pillar == background (or fill = 1). Setting fill_y = 1 on a Rectangle
// makes the cell invariant in y — how we cross-check 2D against the 1D solver.
struct RectCell2D {
    Material pillar;
    Material background;
    double fill_x;
    double fill_y;
    double thickness_um;
    MetaShape shape = MetaShape::Rectangle;
    double shape_param = 0.5;  // Cross: arm width / fill; Ring: inner / outer radius
    double rotation_rad = 0.0; // in-plane rotation of the meta-atom about the cell
                               // center (CCW). Nonzero breaks the separable analytic
                               // path -> the cell goes through the sampled-grid
                               // factorization. The enabler for Pancharatnam-Berry
                               // (geometric-phase) optics, where a fixed birefringent
                               // atom is ROTATED per site to imprint phase 2*rotation.

    // 2D Fourier coefficient ε_{p,q} of the permittivity (Rectangle, analytic).
    cdouble eps_fourier(int p, int q, double wavelength_um) const;

    // 2D Fourier coefficient of the RECIPROCAL permittivity (1/ε)_{p,q}, used by
    // Li's inverse-rule factorization for the field component normal to a
    // material interface (what makes high-contrast cells converge). Analytic.
    cdouble inv_eps_fourier(int p, int q, double wavelength_um) const;

    // A rotated cell is no longer separable, so even a Rectangle must take the
    // sampled-grid (Laurent) path rather than the analytic inverse-rule path.
    bool is_plain_rect() const {
        return shape == MetaShape::Rectangle && rotation_rad == 0.0;
    }
    // True if the sampled grid point (cell-fraction coords in [-0.5,0.5]) is
    // inside the pillar for this shape.
    bool inside(double xf, double yf) const;
    // Rasterize ε(x,y) onto an ngx*ngy row-major grid (sample centers), for the
    // numerical factorization of non-rectangular shapes.
    void rasterize_eps(std::vector<cdouble>& out, int ngx, int ngy,
                       double wavelength_um) const;

    static RectCell2D homogeneous(Material m, double thickness_um) {
        return RectCell2D{m, m, 1.0, 1.0, thickness_um, MetaShape::Rectangle, 0.5};
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
