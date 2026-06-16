#include "celeris/design/metalens.hpp"

#include "celeris/core.hpp"
#include "celeris/rcwa/rcwa2d.hpp"

#include <algorithm>
#include <cmath>
#include <future>
#include <thread>

namespace celeris {
namespace {

// Smallest signed difference between two angles, in (-pi, pi].
double angle_diff(double a, double b) {
    double d = a - b;
    while (d > pi) d -= 2.0 * pi;
    while (d <= -pi) d += 2.0 * pi;
    return d;
}

} // namespace

double UnitCellLibrary::phase_span() const {
    if (phase.empty()) return 0.0;
    auto [lo, hi] = std::minmax_element(phase.begin(), phase.end());
    return *hi - *lo;
}

int UnitCellLibrary::lookup(double target_phase_rad) const {
    int best = 0;
    double best_err = std::abs(angle_diff(phase[0], target_phase_rad));
    for (int i = 1; i < static_cast<int>(phase.size()); ++i) {
        double e = std::abs(angle_diff(phase[i], target_phase_rad));
        if (e < best_err) {
            best_err = e;
            best = i;
        }
    }
    return best;
}

UnitCellLibrary build_unit_cell_library(const Material& pillar,
                                        const Material& background,
                                        const Material& incident,
                                        const Material& substrate,
                                        double period_um, double wavelength_um,
                                        double thickness_um, double fill_min,
                                        double fill_max, int n_samples, int M) {
    UnitCellLibrary lib;
    lib.period_um = period_um;
    lib.wavelength_um = wavelength_um;
    lib.thickness_um = thickness_um;
    lib.fill.resize(n_samples);
    lib.phase.resize(n_samples);
    lib.amplitude.resize(n_samples);

    // Each pillar size is an independent RCWA solve — embarrassingly parallel.
    // Spread the samples across the hardware threads. (This same batch-over-
    // geometries structure is exactly what the GPU port will parallelize next.)
    auto solve_one = [&](int i) {
        double f = fill_min + (fill_max - fill_min) * i / (n_samples - 1);
        Rcwa2DStack cell{period_um, period_um,
                         {RectCell2D{pillar, background, f, f, thickness_um}}};
        auto r = solve_rcwa_2d(incident, cell, substrate, wavelength_um, 0.0,
                               0.0, /*Ex0=*/1.0, /*Ey0=*/0.0, M, M);
        lib.fill[i] = f;
        lib.phase[i] = std::arg(r.tx0);
        lib.amplitude[i] = std::abs(r.tx0);
    };

    unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    int workers = std::min<int>(static_cast<int>(hw), n_samples);
    std::vector<std::future<void>> jobs;
    for (int w = 0; w < workers; ++w) {
        jobs.push_back(std::async(std::launch::async, [&, w] {
            for (int i = w; i < n_samples; i += workers) solve_one(i);
        }));
    }
    for (auto& j : jobs) j.get();
    return lib;
}

MetalensDesign design_metalens(const UnitCellLibrary& lib,
                               double focal_length_um, double diameter_um) {
    const double p = lib.period_um;
    const double lambda = lib.wavelength_um;
    const int n = std::max(1, static_cast<int>(std::round(diameter_um / p)));
    const double center = (n - 1) / 2.0;

    MetalensDesign d;
    d.n_cells = n;
    d.period_um = p;
    d.fill_map.resize(static_cast<std::size_t>(n) * n);

    double sq_err = 0.0;
    double amp_sum = 0.0;
    int count = 0;
    for (int iy = 0; iy < n; ++iy) {
        for (int ix = 0; ix < n; ++ix) {
            const double x = (ix - center) * p;
            const double y = (iy - center) * p;
            const double r = std::sqrt(x * x + y * y);
            // Ideal focusing (hyperbolic) phase profile.
            const double target =
                -(2.0 * pi / lambda) * (std::sqrt(r * r + focal_length_um *
                                                          focal_length_um) -
                                        focal_length_um);
            const int idx = lib.lookup(target);
            d.fill_map[static_cast<std::size_t>(iy) * n + ix] = lib.fill[idx];

            const double err = angle_diff(lib.phase[idx], target);
            sq_err += err * err;
            amp_sum += lib.amplitude[idx];
            ++count;
        }
    }
    d.rms_phase_error_deg = std::sqrt(sq_err / count) * 180.0 / pi;
    d.mean_amplitude = amp_sum / count;
    return d;
}

} // namespace celeris
