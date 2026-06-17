#include "celeris/design/optimize.hpp"

#include "celeris/core.hpp"
#include "celeris/rcwa/rcwa2d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <future>

namespace celeris {
namespace {

double wrap_pi(double a) {
    while (a > pi) a -= 2.0 * pi;
    while (a <= -pi) a += 2.0 * pi;
    return a;
}

} // namespace

OptimizedPillar optimize_pillar(const Material& pillar, const Material& background,
                                const Material& incident, const Material& substrate,
                                double period_um, const PillarTarget& target,
                                int M, double fill0, double thickness0,
                                int max_iters) {
    // Design variables x = {fill, thickness_um}, with box constraints.
    const std::array<double, 2> lo{0.08, 0.20};
    const std::array<double, 2> hi{0.92, 1.00};
    std::array<double, 2> x{std::clamp(fill0, lo[0], hi[0]),
                            std::clamp(thickness0, lo[1], hi[1])};

    // Evaluate the unit cell: returns {phase, amplitude} of zeroth-order t.
    auto evaluate = [&](const std::array<double, 2>& v) -> std::array<double, 2> {
        Rcwa2DStack cell{period_um, period_um,
                         {RectCell2D{pillar, background, v[0], v[0], v[1]}}};
        auto r = solve_rcwa_2d(incident, cell, substrate, target.wavelength_um,
                               0.0, 0.0, /*Ex0=*/1.0, /*Ey0=*/0.0, M, M);
        return {std::arg(r.tx0), std::abs(r.tx0)};
    };
    auto loss = [&](const std::array<double, 2>& v) -> double {
        auto pa = evaluate(v);
        double dphi = wrap_pi(pa[0] - target.target_phase_rad);
        double damp = 1.0 - pa[1];
        return dphi * dphi + target.amplitude_weight * damp * damp;
    };

    // Adam optimizer with central-difference gradients.
    std::array<double, 2> m{0.0, 0.0}, v2{0.0, 0.0};
    const double b1 = 0.9, b2 = 0.999, eps = 1e-8;
    const std::array<double, 2> lr{0.04, 0.04};
    const std::array<double, 2> h{5e-3, 5e-3};  // FD step per variable

    std::array<double, 2> best = x;
    double best_loss = loss(x);

    for (int t = 1; t <= max_iters; ++t) {
        // Central-difference gradient: the 2*Ndim perturbed evaluations are
        // independent RCWA solves, so run them concurrently.
        std::array<double, 2> grad{0.0, 0.0};
        std::array<std::array<double, 2>, 2> xp{x, x}, xm{x, x};
        for (int k = 0; k < 2; ++k) {
            xp[k][k] = std::clamp(x[k] + h[k], lo[k], hi[k]);
            xm[k][k] = std::clamp(x[k] - h[k], lo[k], hi[k]);
        }
        auto j0p = std::async(std::launch::async, [&] { return loss(xp[0]); });
        auto j0m = std::async(std::launch::async, [&] { return loss(xm[0]); });
        auto j1p = std::async(std::launch::async, [&] { return loss(xp[1]); });
        auto j1m = std::async(std::launch::async, [&] { return loss(xm[1]); });
        grad[0] = (j0p.get() - j0m.get()) / (xp[0][0] - xm[0][0]);
        grad[1] = (j1p.get() - j1m.get()) / (xp[1][1] - xm[1][1]);
        for (int k = 0; k < 2; ++k) {
            m[k] = b1 * m[k] + (1 - b1) * grad[k];
            v2[k] = b2 * v2[k] + (1 - b2) * grad[k] * grad[k];
            double mhat = m[k] / (1 - std::pow(b1, t));
            double vhat = v2[k] / (1 - std::pow(b2, t));
            x[k] = std::clamp(x[k] - lr[k] * mhat / (std::sqrt(vhat) + eps),
                              lo[k], hi[k]);
        }
        double l = loss(x);
        if (l < best_loss) { best_loss = l; best = x; }
    }

    auto pa = evaluate(best);
    return {best[0], best[1], pa[0], pa[1], best_loss, max_iters};
}

} // namespace celeris
