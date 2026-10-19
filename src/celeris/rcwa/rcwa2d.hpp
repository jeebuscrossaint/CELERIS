#pragma once
/// @file
/// @brief 2D RCWA (Fourier modal method) for biperiodic structures.
///
/// The real metalens case: nanopillars on a grid, periodic in BOTH x and y.
/// Unlike 1D, the polarizations couple, so the per-layer eigenproblem is fully
/// vectorial: a 2N×2N operator on the tangential H-field harmonics [Hx;Hy]
/// (Liu–Fan / grcwa convention), where N = (2Mx+1)(2My+1) is the number of 2D
/// Fourier orders. Layer stacking reuses the same scattering-matrix machinery as
/// the 1D solver.
///
/// Fourier factorization: plain (unrotated) rectangles use Li's analytic
/// inverse rule, which converges correctly on high-contrast cells; every other
/// shape (and any rotated cell) is rasterized and uses a stable Laurent
/// factorization, which converges more slowly on curved interfaces.
/// References: Moharam & Gaylord (1995); L. Li, JOSA A 14 (1997);
/// Liu & Fan, Comput. Phys. Commun. 183 (2012).

#include "celeris/core.hpp"
#include "celeris/materials/material.hpp"

#include <vector>

namespace celeris {

/// @brief Meta-atom cross-section shape, centered in the cell.
///
/// Bounded by the fill_x × fill_y box. Rectangle has a closed-form separable
/// Fourier series (the fast, analytic path); the others are non-separable and go
/// through the sampled-grid factorization.
enum class MetaShape {
    Rectangle,  ///< |x| < fill_x/2 and |y| < fill_y/2 (the analytic cell).
    Ellipse,    ///< (x/(fill_x/2))² + (y/(fill_y/2))² < 1; a circle if fill_x == fill_y.
    Cross,      ///< Two centered bars; arm width = shape_param · fill.
    Ring        ///< Annulus of pillar material; inner radius = shape_param · outer.
};

/// @brief One layer of a 2D stack: a shaped pillar embedded in a background.
///
/// The pillar (shape inside the fractional box fill_x × fill_y, centered in the
/// cell) is embedded in `background`. A homogeneous layer is pillar ==
/// background (or fill = 1). Setting fill_y = 1 on a Rectangle makes the cell
/// invariant in y — how the 2D solver is cross-checked against the 1D solver.
struct RectCell2D {
    Material pillar;       ///< Pillar (meta-atom) material.
    Material background;   ///< Surrounding material.
    double fill_x;         ///< Pillar extent in x as a fraction of period_x, 0..1.
    double fill_y;         ///< Pillar extent in y as a fraction of period_y, 0..1.
    double thickness_um;   ///< Layer thickness (µm).
    MetaShape shape = MetaShape::Rectangle;  ///< Cross-section shape.
    double shape_param = 0.5;  ///< Cross: arm width / fill; Ring: inner / outer radius.
    /// @brief In-plane CCW rotation of the meta-atom about the cell center (rad).
    ///
    /// Nonzero breaks the separable analytic path, so the cell goes through the
    /// sampled-grid factorization. This is the enabler for Pancharatnam–Berry
    /// (geometric-phase) optics, where a fixed birefringent atom is rotated per
    /// site to imprint a phase of 2 × rotation.
    double rotation_rad = 0.0;

    /// @brief 2D Fourier coefficient ε_{p,q} of the permittivity (Rectangle, analytic).
    /// @param p,q            Fourier orders in x and y.
    /// @param wavelength_um  Vacuum wavelength (µm).
    cdouble eps_fourier(int p, int q, double wavelength_um) const;

    /// @brief 2D Fourier coefficient of the reciprocal permittivity (1/ε)_{p,q}.
    ///
    /// Used by Li's inverse-rule factorization for the field component normal to
    /// a material interface (what makes high-contrast cells converge). Analytic.
    /// @param p,q            Fourier orders in x and y.
    /// @param wavelength_um  Vacuum wavelength (µm).
    cdouble inv_eps_fourier(int p, int q, double wavelength_um) const;

    /// @brief True if the cell can use the analytic inverse-rule path.
    ///
    /// A rotated cell is no longer separable, so even a Rectangle must take the
    /// sampled-grid (Laurent) path.
    bool is_plain_rect() const {
        return shape == MetaShape::Rectangle && rotation_rad == 0.0;
    }
    /// @brief True if a point is inside the pillar.
    /// @param xf,yf  Position in cell-fraction coordinates, each in [-0.5, 0.5].
    bool inside(double xf, double yf) const;
    /// @brief Rasterize ε(x,y) for the numerical factorization of non-rectangular shapes.
    /// @param out            Output: ngx*ngy samples, row-major, at sample centers.
    /// @param ngx,ngy        Grid size in x and y.
    /// @param wavelength_um  Vacuum wavelength (µm).
    void rasterize_eps(std::vector<cdouble>& out, int ngx, int ngy,
                       double wavelength_um) const;

    /// @brief A uniform (unpatterned) layer of a single material.
    /// @param m             The layer material.
    /// @param thickness_um  Layer thickness (µm).
    static RectCell2D homogeneous(Material m, double thickness_um) {
        return RectCell2D{m, m, 1.0, 1.0, thickness_um, MetaShape::Rectangle, 0.5};
    }
};

/// @brief A biperiodic stack: shared x/y periods and an ordered list of layers.
struct Rcwa2DStack {
    double period_x_um;               ///< Period in x (µm).
    double period_y_um;               ///< Period in y (µm).
    std::vector<RectCell2D> layers;   ///< Layers, top (incident side) to bottom.
};

/// @brief Totals and zeroth-order amplitudes of a 2D RCWA solve.
struct Rcwa2DResult {
    double R;        ///< Total reflected efficiency (all orders).
    double T;        ///< Total transmitted efficiency (all orders).
    double sum_de;   ///< R + T; equals 1 for a lossless stack.
    double de_t0;    ///< Zeroth-order transmitted efficiency.
    double de_r0;    ///< Zeroth-order reflected efficiency.
    /// @brief Complex zeroth-order transmitted Ex amplitude.
    ///
    /// arg() of tx0/ty0 is the phase delay a metalens pillar imparts, the core
    /// quantity for unit-cell libraries.
    cdouble tx0;
    cdouble ty0;     ///< Complex zeroth-order transmitted Ey amplitude.
};

/// @brief Energy split for one diffraction order (p,q).
///
/// Optionally filled by solve_rcwa_2d() for an efficiency budget ("where does
/// the light go?"). Evanescent orders carry ~0 real flux and are bookkeeping.
struct OrderEfficiency {
    int p;        ///< Diffraction order index in x.
    int q;        ///< Diffraction order index in y.
    double de_t;  ///< Fraction of incident power transmitted into this order.
    double de_r;  ///< Fraction of incident power reflected into this order.
    bool prop_t;  ///< True if the order propagates (radiates) in the substrate.
    bool prop_r;  ///< True if the order propagates in the incident region.
};

/// @brief Solve a biperiodic stack between two half-spaces.
/// @param incident       Half-space above the stack (light comes from here).
/// @param stack          The layers and their periods.
/// @param substrate      Half-space below the stack.
/// @param wavelength_um  Vacuum wavelength (µm).
/// @param theta_rad      Polar incidence angle from the normal.
/// @param phi_rad        Azimuthal incidence angle, measured from the x axis.
/// @param Ex0,Ey0        Tangential incident E-field at order 0; e.g. (1,0) is
///                       E along x and (0,1) is E along y.
/// @param Mx,My          Harmonic half-counts: orders -Mx..Mx and -My..My.
/// @param orders_out     If non-null, filled with one OrderEfficiency per (p,q);
///                       pass nullptr to skip that bookkeeping on the hot path.
/// @return Total efficiencies and the complex zeroth-order transmission.
Rcwa2DResult solve_rcwa_2d(const Material& incident,
                           const Rcwa2DStack& stack,
                           const Material& substrate,
                           double wavelength_um,
                           double theta_rad,
                           double phi_rad,
                           cdouble Ex0,
                           cdouble Ey0,
                           int Mx,
                           int My,
                           std::vector<OrderEfficiency>* orders_out = nullptr);

} // namespace celeris
