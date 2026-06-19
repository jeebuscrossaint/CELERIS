#include "app_state.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <vector>

#include "celeris/analysis/chromatic.hpp"
#include "celeris/materials/database.hpp"
#include "celeris/rcwa/rcwa2d.hpp"

using namespace celeris;

namespace celeris::gui {

void set_phase(const char* msg, float progress) {
    std::lock_guard<std::mutex> lk(g_mtx);
    g_status = msg;
    g_progress = progress;
}

void run_design(Params p) {
    g_running = true;
    set_phase("Building unit-cell library (RCWA sweep)...", 0.05f);
    // Assemble the unit-cell stack: extra layers (caps/AR coatings) above the
    // patterned pillar layer, which is the active (fill-swept) layer.
    Rcwa2DStack stack;
    stack.period_x_um = stack.period_y_um = p.period;
    for (const auto& L : p.extra_layers)
        stack.layers.push_back(RectCell2D{
            Material::constant(cdouble{L.n, 0.0}, "layer"), materials::air(),
            L.fill, L.fill, L.thickness});
    const int active = static_cast<int>(stack.layers.size());
    stack.layers.push_back(RectCell2D{make_pillar(p), materials::air(),
                                      0.5, 0.5, p.thickness});
    auto lib = build_unit_cell_library_stack(stack, active, materials::air(),
                                             make_substrate(p), p.wavelength, 0.08,
                                             0.92, p.fill_samples, p.harmonics);
    set_phase("Assembling lens (phase profile -> pillar map)...", 0.45f);
    auto lens = design_metalens(lib, p.focal, p.diameter);
    set_phase("Analyzing focus (Rayleigh-Sommerfeld)...", 0.55f);
    auto foc = analyze_focus(lens, lib, p.focal, p.wavelength, p.diameter);
    double dl = p.wavelength * p.focal / p.diameter;
    set_phase("Rendering point-spread function...", 0.72f);
    auto psf = compute_psf(lens, lib, p.focal, p.wavelength, p.diameter, 161,
                           std::max(5.0 * dl, 4.0));
    set_phase("Chromatic sweep (re-solving meta-atoms per wavelength)...", 0.9f);
    auto build_lib_at = [p](double lam) {
        Rcwa2DStack s;
        s.period_x_um = s.period_y_um = p.period;
        for (const auto& L : p.extra_layers)
            s.layers.push_back(RectCell2D{
                Material::constant(cdouble{L.n, 0.0}, "layer"), materials::air(),
                L.fill, L.fill, L.thickness});
        const int act = static_cast<int>(s.layers.size());
        s.layers.push_back(RectCell2D{make_pillar(p), materials::air(),
                                      0.5, 0.5, p.thickness});
        return build_unit_cell_library_stack(s, act, materials::air(),
                                             make_substrate(p), lam, 0.08, 0.92,
                                             p.fill_samples, p.harmonics);
    };
    auto chrom = analyze_chromatic_dispersive(
        lens, build_lib_at, p.focal, p.wavelength, p.diameter,
        p.wavelength * 0.85, p.wavelength * 1.25, 11);

    Results r;
    r.strehl = foc.strehl;
    r.fwhm = foc.fwhm_um;
    r.dl = foc.diffraction_limit_um;
    r.encircled = foc.encircled_energy;
    r.rms = lens.rms_phase_error_deg;
    r.meanT = lens.mean_amplitude;
    r.coverage_deg = lib.phase_span() * 180.0 / 3.14159265358979;
    r.na = std::sin(std::atan((p.diameter / 2.0) / p.focal));
    r.n_cells = lens.n_cells;
    r.pillars = lens.n_cells * lens.n_cells;
    r.psf = std::move(psf);
    for (auto& c : chrom) {
        r.chrom_wl.push_back(static_cast<float>(c.wavelength_um * 1000.0));
        r.chrom_focus.push_back(static_cast<float>(c.focal_length_um));
    }
    r.wf = analyze_wavefront(lens, lib, p.focal, p.wavelength, p.diameter);
    r.mtf = analyze_mtf(lens, lib, p.focal, p.wavelength, p.diameter);
    r.tf = analyze_through_focus(lens, lib, p.focal, p.wavelength, p.diameter);
    r.design = lens;
    r.lib = std::move(lib);
    r.used = p;
    {
        std::lock_guard<std::mutex> lk(g_mtx);
        g_res = std::move(r);
        g_status = "Done.";
        g_progress = 1.0f;
    }
    g_pending = true;
    g_running = false;
}

void run_tolerance() {
    g_running = true;
    set_phase("Monte-Carlo fabrication tolerance...", 0.2f);
    MetalensDesign d; UnitCellLibrary lib; Params p;
    { std::lock_guard<std::mutex> lk(g_mtx); d = g_res.design; lib = g_res.lib; p = g_res.used; }
    auto tol = analyze_tolerance(d, lib, p.focal, p.wavelength, p.diameter,
                                 {0.0, 5.0, 10.0, 20.0}, 12, 12345);
    { std::lock_guard<std::mutex> lk(g_mtx); g_tol = std::move(tol);
      g_status = "Tolerance analysis done."; g_progress = 1.0f; }
    g_tol_pending = true;
    g_running = false;
}

void run_fov() {
    g_running = true;
    set_phase("Field-of-view (off-axis) analysis...", 0.2f);
    MetalensDesign d; UnitCellLibrary lib; Params p;
    { std::lock_guard<std::mutex> lk(g_mtx); d = g_res.design; lib = g_res.lib; p = g_res.used; }
    auto fov = analyze_field_of_view(d, lib, p.focal, p.wavelength, p.diameter,
                                     {0.0, 1.0, 2.0, 5.0, 10.0});
    { std::lock_guard<std::mutex> lk(g_mtx); g_fov = std::move(fov);
      g_status = "Field-of-view analysis done."; g_progress = 1.0f; }
    g_fov_pending = true;
    g_running = false;
}

void run_spotgrid() {
    g_running = true;
    set_phase("Spot-vs-field diagram (off-axis PSFs)...", 0.1f);
    MetalensDesign d; UnitCellLibrary lib; Params p;
    { std::lock_guard<std::mutex> lk(g_mtx); d = g_res.design; lib = g_res.lib; p = g_res.used; }
    const double dl = p.wavelength * p.focal / p.diameter;
    const double win = std::max(8.0 * dl, 5.0);
    const int ng = 81;
    const std::vector<double> angles = {0.0, 2.0, 4.0, 6.0, 8.0};
    std::vector<FieldPsf> spots;
    auto ax = compute_psf_field(d, lib, p.focal, p.wavelength, p.diameter, 0.0, ng, win);
    double peak = 0.0;
    for (double v : ax.psf.intensity) peak = std::max(peak, v);
    ax.rel_strehl = 1.0;
    spots.push_back(std::move(ax));
    for (std::size_t i = 1; i < angles.size(); ++i) {
        set_phase("Spot-vs-field diagram (off-axis PSFs)...",
                  0.1f + 0.85f * static_cast<float>(i) / angles.size());
        spots.push_back(compute_psf_field(d, lib, p.focal, p.wavelength, p.diameter,
                                          angles[i], ng, win, peak));
    }
    { std::lock_guard<std::mutex> lk(g_mtx); g_spot = std::move(spots);
      g_status = "Spot-vs-field diagram done."; g_progress = 1.0f; }
    g_spot_pending = true;
    g_running = false;
}

void run_polardesign(Params p) {
    g_running = true;
    set_phase("Polarization library (fill_x x fill_y, 2 solves/cell)...", 0.1f);
    auto lib = build_polarization_library(
        make_pillar(p), materials::air(), materials::air(), make_substrate(p),
        p.period, p.wavelength, p.thickness, 0.10, 0.90,
        std::max(6, p.fill_samples), p.harmonics);
    set_phase("Assigning rectangular pillars (dual phase profile)...", 0.7f);
    auto d = design_polarization_metalens(lib, p.focal, p.focal_y, p.diameter);

    set_phase("Propagating X/Y-pol focal spots...", 0.9f);
    std::vector<double> px, py;
    std::vector<cdouble> tx, ty;
    const double cen = (d.n_cells - 1) / 2.0, R_ap = p.diameter / 2.0;
    for (int iy = 0; iy < d.n_cells; ++iy)
        for (int ix = 0; ix < d.n_cells; ++ix) {
            double x = (ix - cen) * d.period_um, y = (iy - cen) * d.period_um;
            if (std::sqrt(x * x + y * y) > R_ap) continue;
            std::size_t off = (std::size_t)iy * d.n_cells + ix;
            px.push_back(x); py.push_back(y);
            tx.push_back(d.t_x[off]); ty.push_back(d.t_y[off]);
        }
    double dlx = p.wavelength * p.focal / p.diameter;
    double dly = p.wavelength * p.focal_y / p.diameter;
    auto psfx = propagate_pillars(px, py, tx, 0, 0, p.focal, p.wavelength, 161,
                                  std::max(5.0 * dlx, 4.0));
    auto psfy = propagate_pillars(px, py, ty, 0, 0, p.focal_y, p.wavelength, 161,
                                  std::max(5.0 * dly, 4.0));

    const double kk = 2.0 * pi / p.wavelength;
    auto on_axis = [&](const std::vector<cdouble>& t, double z) {
        cdouble E{0, 0};
        for (std::size_t q = 0; q < px.size(); ++q) {
            double r = std::sqrt(px[q] * px[q] + py[q] * py[q] + z * z);
            E += t[q] * std::polar(1.0 / r, kk * r);
        }
        return std::norm(E);
    };
    auto peak_z = [&](const std::vector<cdouble>& t) {
        double zlo = 0.5 * std::min(p.focal, p.focal_y), zhi = 1.5 * std::max(p.focal, p.focal_y);
        double bz = zlo, bi = -1;
        for (int j = 0; j < 200; ++j) { double z = zlo + (zhi - zlo) * j / 199.0;
            double I = on_axis(t, z); if (I > bi) { bi = I; bz = z; } }
        return bz;
    };
    double zx = peak_z(tx), zy = peak_z(ty);
    double isox = (std::abs(p.focal - p.focal_y) > 1e-6)
        ? 10.0 * std::log10(on_axis(tx, zx) / std::max(on_axis(ty, zx), 1e-30)) : 0.0;
    double isoy = (std::abs(p.focal - p.focal_y) > 1e-6)
        ? 10.0 * std::log10(on_axis(ty, zy) / std::max(on_axis(tx, zy), 1e-30)) : 0.0;
    {
        std::lock_guard<std::mutex> lk(g_mtx);
        g_polar = std::move(d);
        g_polar_psf_x = std::move(psfx);
        g_polar_psf_y = std::move(psfy);
        g_polar_fx = p.focal; g_polar_fy = p.focal_y;
        g_polar_zx = zx; g_polar_zy = zy; g_polar_iso_x = isox; g_polar_iso_y = isoy;
        g_status = std::format("Polarization design: X@{:.0f}um RMS {:.1f}deg, "
                               "Y@{:.0f}um RMS {:.1f}deg",
                               p.focal, g_polar.rms_phase_error_x_deg, p.focal_y,
                               g_polar.rms_phase_error_y_deg);
        g_progress = 1.0f;
    }
    g_polar_pending = true;
    g_running = false;
}

void run_optimize(Params p) {
    g_running = true;
    set_phase("Optimizing design (period x height search)...", 0.0f);
    auto res = optimize_system(
        make_pillar(p), materials::air(), materials::air(), make_substrate(p), p.focal,
        p.diameter, p.wavelength, 0.20, 0.45, 0.30, 1.00, /*grid=*/5, /*M=*/5,
        /*fill_samples=*/12, /*efficiency_weight=*/0.3,
        [](float fr) { set_phase("Optimizing design (period x height search)...", fr); });
    {
        std::lock_guard<std::mutex> lk(g_mtx);
        g_opt = res;
        g_status = std::format("Optimized: period {:.3f} um, height {:.3f} um "
                               "(Strehl~{:.3f}). Applied -- press F5 to run.",
                               res.period_um, res.thickness_um, res.strehl);
        g_progress = 1.0f;
    }
    g_opt_pending = true;
    g_running = false;
}

} // namespace celeris::gui
