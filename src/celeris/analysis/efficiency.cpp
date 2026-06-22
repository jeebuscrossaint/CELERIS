#include "celeris/analysis/efficiency.hpp"

#include <algorithm>

namespace celeris {

EfficiencyBudget analyze_efficiency(const Material& incident,
                                    const Material& pillar,
                                    const Material& background,
                                    const Material& substrate,
                                    double period_x_um, double period_y_um,
                                    double fill_x, double fill_y,
                                    double thickness_um, double wavelength_um,
                                    cdouble Ex0, cdouble Ey0, int M,
                                    MetaShape shape, double shape_param) {
    RectCell2D cell{pillar,     background, fill_x, fill_y, thickness_um,
                    shape, shape_param};
    Rcwa2DStack stack{period_x_um, period_y_um, {cell}};
    std::vector<OrderEfficiency> orders;
    Rcwa2DResult r = solve_rcwa_2d(incident, stack, substrate, wavelength_um,
                                   /*theta=*/0.0, /*phi=*/0.0, Ex0, Ey0, M, M,
                                   &orders);

    EfficiencyBudget b{};
    b.reflection = r.R;
    b.transmission = r.T;
    b.absorption = std::max(0.0, 1.0 - r.R - r.T);
    b.t_zero = r.de_t0;
    b.r_zero = r.de_r0;
    b.t_stray = std::max(0.0, r.T - r.de_t0);
    b.r_stray = std::max(0.0, r.R - r.de_r0);
    b.n_prop_t = 0;
    b.n_prop_r = 0;
    for (const OrderEfficiency& o : orders) {
        if (o.prop_t) ++b.n_prop_t;
        if (o.prop_r) ++b.n_prop_r;
    }
    // Sort by transmitted power, descending — the dominant channels first.
    std::sort(orders.begin(), orders.end(),
              [](const OrderEfficiency& a, const OrderEfficiency& c) {
                  return a.de_t > c.de_t;
              });
    b.orders = std::move(orders);
    return b;
}

} // namespace celeris
