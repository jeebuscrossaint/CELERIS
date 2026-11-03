#pragma once
/// @file
/// @brief Efficiency budget for one meta-atom: where does the light go?
///
/// RCWA conserves energy exactly for a lossless stack (R + T = 1); with a
/// lossy material (metals, absorbing dielectrics) the deficit A = 1 - R - T is
/// real absorption. This breaks the transmitted/reflected power down PER
/// DIFFRACTION ORDER, so you can see how much lands in the useful zeroth order
/// (the wave a metalens phase-profiles) versus higher orders (stray light /
/// diffraction loss — only present when the pitch is not subwavelength). This is
/// the per-element analog of a Zemax surface-efficiency / ghost budget.

#include "celeris/rcwa/rcwa2d.hpp"
#include "celeris/materials/material.hpp"

#include <vector>

namespace celeris {

/// @brief Power budget of one unit cell, totals and per diffraction order.
///
/// All powers are fractions of the incident power.
struct EfficiencyBudget {
    double reflection;     ///< Total R (all orders).
    double transmission;   ///< Total T (all orders).
    double absorption;     ///< 1 - R - T (material loss; ~0 for lossless).
    double t_zero;         ///< Transmitted into the zeroth order (the useful wave).
    double r_zero;         ///< Reflected zeroth order.
    double t_stray;        ///< Transmitted into all higher orders (stray light).
    double r_stray;        ///< Reflected into all higher orders.
    int n_prop_t;          ///< Number of propagating transmitted orders (radiating channels).
    int n_prop_r;          ///< Number of propagating reflected orders.
    std::vector<OrderEfficiency> orders;  ///< Every (p,q) order, sorted by transmitted power.
};

/// @brief Solve one meta-atom at normal incidence and return its energy budget.
///
/// A pillar of `pillar` in `background`, of size fill_x × fill_y and height
/// `thickness_um`, sits on `substrate` and is lit from `incident` with the linear
/// input polarization (Ex0, Ey0).
/// @param incident,substrate   Half-space materials above and below.
/// @param pillar,background    Meta-atom materials.
/// @param period_x_um,period_y_um  Unit-cell periods (µm).
/// @param fill_x,fill_y        Pillar size as fractions of the periods.
/// @param thickness_um         Pillar height (µm).
/// @param wavelength_um        Vacuum wavelength (µm).
/// @param Ex0,Ey0              Incident field components (linear polarization).
/// @param M                    Fourier half-count, used for both axes (Mx = My = M).
/// @param shape                Pillar cross-section.
/// @param shape_param          Shape-specific parameter (see MetaShape).
/// @return The budget; `absorption` is the energy deficit 1 - R - T.
EfficiencyBudget analyze_efficiency(const Material& incident,
                                    const Material& pillar,
                                    const Material& background,
                                    const Material& substrate,
                                    double period_x_um, double period_y_um,
                                    double fill_x, double fill_y,
                                    double thickness_um, double wavelength_um,
                                    cdouble Ex0, cdouble Ey0, int M,
                                    MetaShape shape = MetaShape::Rectangle,
                                    double shape_param = 0.5);

} // namespace celeris
