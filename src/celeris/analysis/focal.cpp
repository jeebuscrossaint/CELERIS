#include "celeris/analysis/focal.hpp"

#include "celeris/core.hpp"

#include <algorithm>
#include <cmath>
#include <future>
#include <thread>
#include <vector>

#ifdef CELERIS_USE_CUDA_KERNELS
#include "celeris/cuda/propagate.hpp"
#endif

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

    // Intensity grid for the design (rows computed in parallel, or on the GPU).
    std::vector<double> I(static_cast<std::size_t>(NG) * NG, 0.0);
    int center_row = NG / 2;
    bool grid_done = false;
#ifdef CELERIS_USE_CUDA_KERNELS
    {
        std::vector<double> px(cells.size()), py(cells.size());
        std::vector<cdouble> pt(cells.size());
        for (std::size_t c = 0; c < cells.size(); ++c) {
            px[c] = cells[c].x; py[c] = cells[c].y; pt[c] = cells[c].t;
        }
        grid_done = cuda::propagate_psf(px.data(), py.data(), pt.data(),
                                        static_cast<int>(cells.size()), 0.0, 0.0,
                                        focal_length_um, k, NG, W, I.data());
    }
#endif
    if (!grid_done)
        parallel_rows(NG, [&](int j) {
            double fy = -W + j * step;
            for (int i = 0; i < NG; ++i) {
                double fx = -W + i * step;
                I[static_cast<std::size_t>(j) * NG + i] =
                    std::norm(field_at(fx, fy, /*ideal=*/false));
            }
        });
    double peak = *std::max_element(I.begin(), I.end());
    // Ideal peak (at focal center): perfect-phase lens, |t|=1 everywhere.
    double peak_ideal = std::norm(field_at(0.0, 0.0, /*ideal=*/true));
    // Same-amplitude perfectly-phased reference: each cell contributes its
    // actual |t| but with the ideal focusing phase, so all add in phase at the
    // center. peak = (Σ |t|/R)^2. Dividing the design peak by this isolates the
    // phase/wavefront quality from the transmission loss (the optical Strehl).
    double amp_inphase = 0.0;
    for (const Cell& c : cells) {
        const double Rr = std::sqrt(c.x * c.x + c.y * c.y +
                                    focal_length_um * focal_length_um);
        amp_inphase += std::abs(c.t) / Rr;
    }
    double peak_same_amp = amp_inphase * amp_inphase;

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
    out.phase_strehl = peak_same_amp > 0 ? peak / peak_same_amp : 0.0;
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

#ifdef CELERIS_USE_CUDA_KERNELS
    // GPU far-field propagation (hundreds of x faster for large apertures).
    {
        std::vector<double> px(cells.size()), py(cells.size());
        std::vector<cdouble> pt(cells.size());
        for (std::size_t c = 0; c < cells.size(); ++c) {
            px[c] = cells[c].x; py[c] = cells[c].y; pt[c] = cells[c].t;
        }
        if (cuda::propagate_psf(px.data(), py.data(), pt.data(),
                                static_cast<int>(cells.size()), 0.0, 0.0,
                                focal_length_um, k, n, half_window_um,
                                m.intensity.data()))
            return m;  // CPU fallback below if this returns false
    }
#endif

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

PsfMap propagate_pillars(const std::vector<double>& px,
                         const std::vector<double>& py,
                         const std::vector<cdouble>& t, double cx, double cy,
                         double z, double wavelength_um, int n,
                         double half_window_um) {
    PsfMap m;
    m.n = n;
    m.half_window_um = half_window_um;
    m.intensity.assign(static_cast<std::size_t>(n) * n, 0.0);
    const double k = 2.0 * pi / wavelength_um;
    const int npil = static_cast<int>(px.size());
    if (npil == 0 || n <= 0) return m;

#ifdef CELERIS_USE_CUDA_KERNELS
    if (cuda::propagate_psf(px.data(), py.data(), t.data(), npil, cx, cy, z, k, n,
                            half_window_um, m.intensity.data()))
        return m;
#endif
    const double step = 2.0 * half_window_um / (n - 1);
    parallel_rows(n, [&](int j) {
        double fy = cy - half_window_um + j * step;
        for (int i = 0; i < n; ++i) {
            double fx = cx - half_window_um + i * step;
            cdouble E{0.0, 0.0};
            for (int p = 0; p < npil; ++p) {
                double R = std::sqrt((fx - px[p]) * (fx - px[p]) +
                                     (fy - py[p]) * (fy - py[p]) + z * z);
                E += t[p] * std::polar(1.0 / R, k * R);
            }
            m.intensity[static_cast<std::size_t>(j) * n + i] = std::norm(E);
        }
    });
    return m;
}

FieldPsf compute_psf_field(const MetalensDesign& lens, const UnitCellLibrary& lib,
                           double focal_length_um, double wavelength_um,
                           double diameter_um, double angle_deg, int n,
                           double half_window_um, double on_axis_peak) {
    const double p = lens.period_um;
    const double center = (lens.n_cells - 1) / 2.0;
    const double R_ap = diameter_um / 2.0;
    const double k = 2.0 * pi / wavelength_um;
    const double theta = angle_deg * pi / 180.0;
    const double sin_t = std::sin(theta);

    // Aperture cells, each carrying the designed transmission times the tilted
    // incident plane-wave phase exp(i k sin(theta) x).
    struct Cell { double x, y; cdouble t; };
    std::vector<Cell> cells;
    for (int iy = 0; iy < lens.n_cells; ++iy)
        for (int ix = 0; ix < lens.n_cells; ++ix) {
            const double x = (ix - center) * p, y = (iy - center) * p;
            if (std::sqrt(x * x + y * y) > R_ap) continue;
            const double fill = lens.fill_map[static_cast<std::size_t>(iy) * lens.n_cells + ix];
            cdouble t = lib.transmission_for_fill(fill) * std::polar(1.0, k * sin_t * x);
            cells.push_back({x, y, t});
        }

    const double cx = focal_length_um * std::tan(theta);  // chief-ray landing
    FieldPsf out;
    out.angle_deg = angle_deg;
    out.cx_um = cx;
    out.psf.n = n;
    out.psf.half_window_um = half_window_um;
    out.psf.intensity.assign(static_cast<std::size_t>(n) * n, 0.0);
    const double step = 2.0 * half_window_um / (n - 1);

#ifdef CELERIS_USE_CUDA_KERNELS
    {
        std::vector<double> px(cells.size()), py(cells.size());
        std::vector<cdouble> pt(cells.size());
        for (std::size_t c = 0; c < cells.size(); ++c) {
            px[c] = cells[c].x; py[c] = cells[c].y; pt[c] = cells[c].t;
        }
        if (cuda::propagate_psf(px.data(), py.data(), pt.data(),
                                static_cast<int>(cells.size()), cx, 0.0,
                                focal_length_um, k, n, half_window_um,
                                out.psf.intensity.data())) {
            double peak = 0.0;
            for (double v : out.psf.intensity) peak = std::max(peak, v);
            out.rel_strehl = on_axis_peak > 0 ? peak / on_axis_peak : 1.0;
            return out;
        }
    }
#endif

    parallel_rows(n, [&](int j) {
        double fy = -half_window_um + j * step;
        for (int i = 0; i < n; ++i) {
            double fx = cx - half_window_um + i * step;  // window centered on cx
            cdouble E{0.0, 0.0};
            for (const Cell& c : cells) {
                double R = std::sqrt((fx - c.x) * (fx - c.x) +
                                     (fy - c.y) * (fy - c.y) +
                                     focal_length_um * focal_length_um);
                E += c.t * std::polar(1.0 / R, k * R);
            }
            out.psf.intensity[static_cast<std::size_t>(j) * n + i] = std::norm(E);
        }
    });

    double peak = 0.0;
    for (double v : out.psf.intensity) peak = std::max(peak, v);
    out.rel_strehl = on_axis_peak > 0 ? peak / on_axis_peak : 1.0;
    return out;
}

} // namespace celeris
