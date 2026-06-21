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

double UnitCellLibrary::coverage() const {
    const int n = static_cast<int>(phase.size());
    if (n < 2) return 0.0;
    std::vector<double> p(phase.begin(), phase.end());
    for (double& v : p) {  // normalize to [0, 2pi)
        while (v < 0.0) v += 2.0 * pi;
        while (v >= 2.0 * pi) v -= 2.0 * pi;
    }
    std::sort(p.begin(), p.end());
    double max_gap = 0.0;
    for (int i = 1; i < n; ++i) max_gap = std::max(max_gap, p[i] - p[i - 1]);
    max_gap = std::max(max_gap, (p[0] + 2.0 * pi) - p[n - 1]);  // wrap-around gap
    return 2.0 * pi - max_gap;
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

int UnitCellLibrary::lookup_weighted(double target_phase_rad,
                                     double amplitude_weight) const {
    int best = 0;
    double best_cost = 1e300;
    for (int i = 0; i < static_cast<int>(phase.size()); ++i) {
        double dphi = angle_diff(phase[i], target_phase_rad);
        double damp = 1.0 - amplitude[i];
        double cost = dphi * dphi + amplitude_weight * damp * damp;
        if (cost < best_cost) { best_cost = cost; best = i; }
    }
    return best;
}

UnitCellLibrary build_unit_cell_library_stack(Rcwa2DStack stack, int active_layer,
                                              const Material& incident,
                                              const Material& substrate,
                                              double wavelength_um, double fill_min,
                                              double fill_max, int n_samples, int M) {
    UnitCellLibrary lib;
    lib.period_um = stack.period_x_um;
    lib.wavelength_um = wavelength_um;
    lib.thickness_um =
        (active_layer >= 0 && active_layer < static_cast<int>(stack.layers.size()))
            ? stack.layers[active_layer].thickness_um
            : 0.0;
    lib.fill.resize(n_samples);
    lib.phase.resize(n_samples);
    lib.amplitude.resize(n_samples);

    auto solve_one = [&](int i) {
        double f = fill_min + (fill_max - fill_min) * i / (n_samples - 1);
        Rcwa2DStack s = stack;  // private copy, set the active layer's fill
        if (active_layer >= 0 && active_layer < static_cast<int>(s.layers.size())) {
            s.layers[active_layer].fill_x = f;
            s.layers[active_layer].fill_y = f;
        }
        auto r = solve_rcwa_2d(incident, s, substrate, wavelength_um, 0.0, 0.0,
                               1.0, 0.0, M, M);
        lib.fill[i] = f;
        lib.phase[i] = std::arg(r.tx0);
        lib.amplitude[i] = std::abs(r.tx0);
    };

    unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    int workers = std::min<int>(static_cast<int>(hw), n_samples);
    std::vector<std::future<void>> jobs;
    for (int w = 0; w < workers; ++w)
        jobs.push_back(std::async(std::launch::async, [&, w] {
            for (int i = w; i < n_samples; i += workers) solve_one(i);
        }));
    for (auto& j : jobs) j.get();
    return lib;
}

cdouble UnitCellLibrary::transmission_for_fill(double f) const {
    int best = 0;
    double best_d = std::abs(fill[0] - f);
    for (int i = 1; i < static_cast<int>(fill.size()); ++i) {
        double d = std::abs(fill[i] - f);
        if (d < best_d) { best_d = d; best = i; }
    }
    return std::polar(amplitude[best], phase[best]);
}

UnitCellLibrary build_unit_cell_library(const Material& pillar,
                                        const Material& background,
                                        const Material& incident,
                                        const Material& substrate,
                                        double period_um, double wavelength_um,
                                        double thickness_um, double fill_min,
                                        double fill_max, int n_samples, int M,
                                        MetaShape shape, double shape_param) {
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
                         {RectCell2D{pillar, background, f, f, thickness_um,
                                     shape, shape_param}}};
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

HeightOptResult optimize_height_for_2pi(const Material& pillar,
                                        const Material& background,
                                        const Material& incident,
                                        const Material& substrate,
                                        double period_um, double wavelength_um,
                                        double thick_lo, double thick_hi,
                                        int n_heights, double fill_min,
                                        double fill_max, int fill_samples, int M,
                                        double coverage_target_deg,
                                        MetaShape shape, double shape_param) {
    n_heights = std::max(2, n_heights);
    HeightOptResult res{};
    res.coverage_target_deg = coverage_target_deg;

    auto mean_T = [](const UnitCellLibrary& l) {
        double s = 0.0;
        for (double a : l.amplitude) s += a * a;
        return l.amplitude.empty() ? 0.0 : s / l.amplitude.size();
    };

    // Build a library at each height; keep them all so selection can be a clean
    // two-phase decision (the libraries are tiny -- a few dozen doubles each).
    std::vector<UnitCellLibrary> libs;
    libs.reserve(n_heights);
    double best_cov = 0.0;
    for (int i = 0; i < n_heights; ++i) {
        double h = thick_lo + (thick_hi - thick_lo) * i / (n_heights - 1);
        auto lib = build_unit_cell_library(pillar, background, incident, substrate,
                                           period_um, wavelength_um, h, fill_min,
                                           fill_max, fill_samples, M, shape, shape_param);
        double cov_deg = lib.coverage() * 180.0 / pi;
        res.sweep.push_back({h, cov_deg, mean_T(lib)});
        best_cov = std::max(best_cov, cov_deg);
        libs.push_back(std::move(lib));
    }

    // Selection. Candidate heights are those whose coverage is "good enough":
    // either they clear the absolute target, or -- if nothing does -- they sit
    // within a tolerance of the best coverage achieved. Among the candidates the
    // winner is the HIGHEST transmittance, the lever on transmission-weighted
    // Strehl. (The old "just take max coverage" fallback could pick a high-
    // coverage / low-transmittance height, which is exactly wrong.)
    const double tol = 12.0;  // degrees
    bool any_clears = best_cov >= coverage_target_deg - 1e-9;
    double floor_deg = any_clears ? coverage_target_deg : best_cov - tol;

    int best_idx = 0;
    double best_T = -1.0;
    for (int i = 0; i < n_heights; ++i) {
        if (res.sweep[i].coverage_deg < floor_deg - 1e-9) continue;
        if (res.sweep[i].mean_transmittance > best_T + 1e-12) {
            best_T = res.sweep[i].mean_transmittance;
            best_idx = i;
        }
    }

    const auto& b = res.sweep[best_idx];
    res.best_thickness_um = b.thickness_um;
    res.coverage_deg = b.coverage_deg;
    res.mean_transmittance = b.mean_transmittance;
    res.reached_target = b.coverage_deg >= coverage_target_deg - 1e-9;
    res.best_library = std::move(libs[best_idx]);
    return res;
}

MetalensDesign design_metalens(const UnitCellLibrary& lib,
                               const PhaseProfile& profile, double diameter_um,
                               double amplitude_weight) {
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
            const double target = phase_profile_value(profile, x, y, lambda);
            const int idx = lib.lookup_weighted(target, amplitude_weight);
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

MetalensDesign design_metalens(const UnitCellLibrary& lib,
                               double focal_length_um, double diameter_um,
                               double amplitude_weight) {
    PhaseProfile profile;
    profile.kind = PhaseProfileKind::Focusing;
    profile.focal_length_um = focal_length_um;
    return design_metalens(lib, profile, diameter_um, amplitude_weight);
}

} // namespace celeris
