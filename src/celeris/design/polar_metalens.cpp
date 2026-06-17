#include "celeris/design/polar_metalens.hpp"

#include "celeris/core.hpp"
#include "celeris/rcwa/rcwa2d.hpp"

#include <algorithm>
#include <cmath>
#include <future>
#include <thread>

namespace celeris {
namespace {

double angle_diff(double a, double b) {
    double d = a - b;
    while (d > pi) d -= 2.0 * pi;
    while (d <= -pi) d += 2.0 * pi;
    return d;
}

} // namespace

int PolarizationLibrary::lookup(double target_x, double target_y,
                                double amplitude_weight) const {
    int best = 0;
    double best_cost = 1e300;
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        const PolarCell& c = cells[i];
        double dx = angle_diff(c.phase_x, target_x);
        double dy = angle_diff(c.phase_y, target_y);
        double ax = 1.0 - c.amp_x, ay = 1.0 - c.amp_y;
        double cost = dx * dx + dy * dy + amplitude_weight * (ax * ax + ay * ay);
        if (cost < best_cost) { best_cost = cost; best = i; }
    }
    return best;
}

PolarizationLibrary build_polarization_library(
    const Material& pillar, const Material& background, const Material& incident,
    const Material& substrate, double period_um, double wavelength_um,
    double thickness_um, double fill_min, double fill_max, int n_samples, int M) {
    PolarizationLibrary lib;
    lib.period_um = period_um;
    lib.wavelength_um = wavelength_um;
    lib.thickness_um = thickness_um;
    const int total = n_samples * n_samples;
    lib.cells.resize(total);

    auto solve_one = [&](int idx) {
        int ix = idx % n_samples, iy = idx / n_samples;
        double fx = fill_min + (fill_max - fill_min) * ix / (n_samples - 1);
        double fy = fill_min + (fill_max - fill_min) * iy / (n_samples - 1);
        Rcwa2DStack cell{period_um, period_um,
                         {RectCell2D{pillar, background, fx, fy, thickness_um}}};
        auto rx = solve_rcwa_2d(incident, cell, substrate, wavelength_um, 0.0, 0.0,
                                1.0, 0.0, M, M);
        auto ry = solve_rcwa_2d(incident, cell, substrate, wavelength_um, 0.0, 0.0,
                                0.0, 1.0, M, M);
        lib.cells[idx] = {fx, fy, std::arg(rx.tx0), std::arg(ry.ty0),
                          std::abs(rx.tx0), std::abs(ry.ty0)};
    };

    unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    int workers = std::min<int>(static_cast<int>(hw), std::max(1, total));
    std::vector<std::future<void>> jobs;
    for (int w = 0; w < workers; ++w)
        jobs.push_back(std::async(std::launch::async, [&, w] {
            for (int i = w; i < total; i += workers) solve_one(i);
        }));
    for (auto& j : jobs) j.get();
    return lib;
}

PolarMetalensDesign design_polarization_metalens(const PolarizationLibrary& lib,
                                                 double focal_x_um,
                                                 double focal_y_um,
                                                 double diameter_um,
                                                 double amplitude_weight) {
    const double p = lib.period_um;
    const double lambda = lib.wavelength_um;
    const int n = std::max(1, static_cast<int>(std::round(diameter_um / p)));
    const double center = (n - 1) / 2.0;

    PolarMetalensDesign d;
    d.n_cells = n;
    d.period_um = p;
    d.fill_x.resize(static_cast<std::size_t>(n) * n);
    d.fill_y.resize(static_cast<std::size_t>(n) * n);

    auto target_phase = [&](double r, double f) {
        return -(2.0 * pi / lambda) * (std::sqrt(r * r + f * f) - f);
    };

    double sqx = 0, sqy = 0, ampx = 0, ampy = 0;
    int count = 0;
    for (int iy = 0; iy < n; ++iy)
        for (int ix = 0; ix < n; ++ix) {
            double x = (ix - center) * p, y = (iy - center) * p;
            double r = std::sqrt(x * x + y * y);
            double tx = target_phase(r, focal_x_um);
            double ty = target_phase(r, focal_y_um);
            int k = lib.lookup(tx, ty, amplitude_weight);
            const PolarCell& c = lib.cells[k];
            std::size_t off = static_cast<std::size_t>(iy) * n + ix;
            d.fill_x[off] = c.fill_x;
            d.fill_y[off] = c.fill_y;
            double ex = angle_diff(c.phase_x, tx), ey = angle_diff(c.phase_y, ty);
            sqx += ex * ex; sqy += ey * ey;
            ampx += c.amp_x; ampy += c.amp_y;
            ++count;
        }
    if (count) {
        d.rms_phase_error_x_deg = std::sqrt(sqx / count) * 180.0 / pi;
        d.rms_phase_error_y_deg = std::sqrt(sqy / count) * 180.0 / pi;
        d.mean_amp_x = ampx / count;
        d.mean_amp_y = ampy / count;
    }
    return d;
}

} // namespace celeris
