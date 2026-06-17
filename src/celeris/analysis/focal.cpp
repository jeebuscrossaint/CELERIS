#include "celeris/analysis/focal.hpp"

#include "celeris/core.hpp"

#include <algorithm>
#include <cmath>
#include <future>
#include <thread>
#include <vector>

namespace celeris {
namespace {
// Run fn(j) for j in [0,n) across the hardware threads (focal-plane rows are
// independent, so this is a clean speedup for large apertures).
template <class Fn>
void parallel_rows(int n, Fn fn) {
    unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    int workers = std::min<int>(static_cast<int>(hw), std::max(1, n));
    std::vector<std::future<void>> jobs;
    for (int w = 0; w < workers; ++w)
        jobs.push_back(std::async(std::launch::async, [=, &fn] {
            for (int j = w; j < n; j += workers) fn(j);
        }));
    for (auto& j : jobs) j.get();
}
} // namespace

FocalAnalysis analyze_focus(const MetalensDesign& lens,
                            const UnitCellLibrary& lib, double focal_length_um,
                            double wavelength_um, double diameter_um) {
    const double p = lens.period_um;
    const double center = (lens.n_cells - 1) / 2.0;
    const double k = 2.0 * pi / wavelength_um;
    const double R_ap = diameter_um / 2.0;

    // Build the aperture: cell position, designed transmission, and the ideal
    // (target) phase, for every cell inside the circular aperture.
    struct Cell { double x, y; cdouble t; double phi_ideal; };
    std::vector<Cell> cells;
    cells.reserve(static_cast<std::size_t>(lens.n_cells) * lens.n_cells);
    for (int iy = 0; iy < lens.n_cells; ++iy) {
        for (int ix = 0; ix < lens.n_cells; ++ix) {
            const double x = (ix - center) * p;
            const double y = (iy - center) * p;
            const double r = std::sqrt(x * x + y * y);
            if (r > R_ap) continue;  // circular aperture mask
            const double fill = lens.fill_map[static_cast<std::size_t>(iy) * lens.n_cells + ix];
            const double phi_ideal =
                -k * (std::sqrt(r * r + focal_length_um * focal_length_um) -
                      focal_length_um);
            cells.push_back({x, y, lib.transmission_for_fill(fill), phi_ideal});
        }
    }

    // Focal-plane window and sampling. The diffraction-limited spot is
    // ~lambda*f/D, so a few-micron window finely sampled resolves it.
    const double dl = wavelength_um * focal_length_um / diameter_um;
    const double W = std::max(6.0 * dl, 4.0);  // half-window (µm)
    const int NG = 121;                         // grid points per axis
    const double step = 2.0 * W / (NG - 1);

    auto field_at = [&](double fx, double fy, bool ideal) -> cdouble {
        cdouble E{0.0, 0.0};
        for (const Cell& c : cells) {
            const double dx = fx - c.x, dy = fy - c.y;
            const double Rr = std::sqrt(dx * dx + dy * dy + focal_length_um * focal_length_um);
            const cdouble prop = std::polar(1.0 / Rr, k * Rr);
            E += (ideal ? std::polar(1.0, c.phi_ideal) : c.t) * prop;
        }
        return E;
    };

    // Intensity grid for the design (rows computed in parallel).
    std::vector<double> I(static_cast<std::size_t>(NG) * NG, 0.0);
    int center_row = NG / 2;
    parallel_rows(NG, [&](int j) {
        double fy = -W + j * step;
        for (int i = 0; i < NG; ++i) {
            double fx = -W + i * step;
            I[static_cast<std::size_t>(j) * NG + i] =
                std::norm(field_at(fx, fy, /*ideal=*/false));
        }
    });
    double peak = *std::max_element(I.begin(), I.end());
    // Ideal peak (at focal center): perfect-phase lens.
    double peak_ideal = std::norm(field_at(0.0, 0.0, /*ideal=*/true));

    // FWHM along the central row.
    double half = peak / 2.0;
    int lo = -1, hi = -1;
    for (int i = 0; i < NG; ++i) {
        double v = I[static_cast<std::size_t>(center_row) * NG + i];
        if (v >= half) { if (lo < 0) lo = i; hi = i; }
    }
    double fwhm = (hi > lo && lo >= 0) ? (hi - lo) * step : 0.0;

    // Encircled energy within the first Airy null (~1.22 lambda f / D).
    double r_null = 1.22 * dl;
    double in = 0.0, total = 0.0;
    for (int j = 0; j < NG; ++j) {
        double fy = -W + j * step;
        for (int i = 0; i < NG; ++i) {
            double fx = -W + i * step;
            double v = I[static_cast<std::size_t>(j) * NG + i];
            total += v;
            if (std::sqrt(fx * fx + fy * fy) <= r_null) in += v;
        }
    }

    FocalAnalysis out;
    out.strehl = peak_ideal > 0 ? peak / peak_ideal : 0.0;
    out.fwhm_um = fwhm;
    out.diffraction_limit_um = dl;
    out.encircled_energy = total > 0 ? in / total : 0.0;
    return out;
}

PsfMap compute_psf(const MetalensDesign& lens, const UnitCellLibrary& lib,
                   double focal_length_um, double wavelength_um,
                   double diameter_um, int n, double half_window_um) {
    const double p = lens.period_um;
    const double center = (lens.n_cells - 1) / 2.0;
    const double R_ap = diameter_um / 2.0;
    const double k = 2.0 * pi / wavelength_um;

    struct Cell { double x, y; cdouble t; };
    std::vector<Cell> cells;
    for (int iy = 0; iy < lens.n_cells; ++iy)
        for (int ix = 0; ix < lens.n_cells; ++ix) {
            const double x = (ix - center) * p, y = (iy - center) * p;
            if (std::sqrt(x * x + y * y) > R_ap) continue;
            const double fill = lens.fill_map[static_cast<std::size_t>(iy) * lens.n_cells + ix];
            cells.push_back({x, y, lib.transmission_for_fill(fill)});
        }

    PsfMap m;
    m.n = n;
    m.half_window_um = half_window_um;
    m.intensity.assign(static_cast<std::size_t>(n) * n, 0.0);
    const double step = 2.0 * half_window_um / (n - 1);
    parallel_rows(n, [&](int j) {
        double fy = -half_window_um + j * step;
        for (int i = 0; i < n; ++i) {
            double fx = -half_window_um + i * step;
            cdouble E{0.0, 0.0};
            for (const Cell& c : cells) {
                double R = std::sqrt((fx - c.x) * (fx - c.x) +
                                     (fy - c.y) * (fy - c.y) +
                                     focal_length_um * focal_length_um);
                E += c.t * std::polar(1.0 / R, k * R);
            }
            m.intensity[static_cast<std::size_t>(j) * n + i] = std::norm(E);
        }
    });
    return m;
}

} // namespace celeris
