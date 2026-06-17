#include "celeris/analysis/polarization.hpp"

#include "celeris/core.hpp"
#include "celeris/rcwa/rcwa2d.hpp"

#include <algorithm>
#include <cmath>
#include <future>
#include <thread>

namespace celeris {
namespace {

// Wrap an angle (radians) to (-pi, pi], return degrees.
double wrap_deg(double a) {
    while (a > pi) a -= 2.0 * pi;
    while (a <= -pi) a += 2.0 * pi;
    return a * 180.0 / pi;
}

} // namespace

std::vector<BirefringencePoint> analyze_birefringence(
    const Material& pillar, const Material& background, const Material& incident,
    const Material& substrate, double period_um, double wavelength_um,
    double thickness_um, double fill_y, double fill_min, double fill_max,
    int n_samples, int M) {
    std::vector<BirefringencePoint> out(std::max(0, n_samples));

    auto solve_one = [&](int i) {
        double fx = fill_min + (fill_max - fill_min) * i / (n_samples - 1);
        Rcwa2DStack cell{period_um, period_um,
                         {RectCell2D{pillar, background, fx, fill_y, thickness_um}}};
        // x-polarized illumination -> co-pol response is tx0.
        auto rx = solve_rcwa_2d(incident, cell, substrate, wavelength_um, 0.0, 0.0,
                                /*Ex0=*/1.0, /*Ey0=*/0.0, M, M);
        // y-polarized illumination -> co-pol response is ty0.
        auto ry = solve_rcwa_2d(incident, cell, substrate, wavelength_um, 0.0, 0.0,
                                /*Ex0=*/0.0, /*Ey0=*/1.0, M, M);
        BirefringencePoint p;
        p.fill_x = fx;
        p.fill_y = fill_y;
        p.phase_x_deg = wrap_deg(std::arg(rx.tx0));
        p.phase_y_deg = wrap_deg(std::arg(ry.ty0));
        p.retardance_deg = wrap_deg(std::arg(rx.tx0) - std::arg(ry.ty0));
        p.tx = std::abs(rx.tx0);
        p.ty = std::abs(ry.ty0);
        out[i] = p;
    };

    unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    int workers = std::min<int>(static_cast<int>(hw), std::max(1, n_samples));
    std::vector<std::future<void>> jobs;
    for (int w = 0; w < workers; ++w)
        jobs.push_back(std::async(std::launch::async, [&, w] {
            for (int i = w; i < n_samples; i += workers) solve_one(i);
        }));
    for (auto& j : jobs) j.get();
    return out;
}

} // namespace celeris
