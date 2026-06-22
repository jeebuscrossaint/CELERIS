#pragma once
// Efficiency budget — "where does the light go?" for a single meta-atom unit
// cell. RCWA conserves energy exactly for a lossless stack (R + T = 1); with a
// lossy material (metals, absorbing dielectrics) the deficit A = 1 - R - T is
// real absorption. This breaks the transmitted/reflected power down PER
// DIFFRACTION ORDER, so you can see how much lands in the useful zeroth order
// (the wave a metalens phase-profiles) versus higher orders (stray light /
// diffraction loss — only present when the pitch is not subwavelength). This is
// the per-element analog of a Zemax surface-efficiency / ghost budget.

#include "celeris/rcwa/rcwa2d.hpp"
#include "celeris/materials/material.hpp"

#include <vector>

namespace celeris {

struct EfficiencyBudget {
    double reflection;     // total R (all orders)
    double transmission;   // total T (all orders)
    double absorption;     // 1 - R - T  (material loss; ~0 for lossless)
    double t_zero;         // transmitted into the zeroth order (the useful wave)
    double r_zero;         // reflected zeroth order
    double t_stray;        // transmitted into all higher orders (stray light)
    double r_stray;        // reflected into all higher orders
    int n_prop_t;          // # of propagating transmitted orders (radiating channels)
    int n_prop_r;          // # of propagating reflected orders
    std::vector<OrderEfficiency> orders;  // every (p,q), sorted by transmitted power
};

// Solve one meta-atom (pillar of `pillar` in `background`, fill_x × fill_y,
// height `thickness_um`, on `substrate`, in `incident`) at normal incidence with
// the given linear input polarization (Ex0,Ey0) and return its energy budget.
// `shape`/`shape_param` select the cross-section. Mx=My=M harmonics.
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
