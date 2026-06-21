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
    r.phase_strehl = foc.phase_strehl;
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

void run_achromatic(Params p) {
    g_running = true;
    // Band samples about the design wavelength (ascending). Kept modest (5 samples)
    // for GUI responsiveness -- the CLI exposes finer control. All three library
    // kinds share the same band, the same two-design (gd_weight 0 vs p.gd_weight)
    // comparison, and the same chromatic-focus verification, so the summary metrics
    // and the focus-vs-lambda plot are uniform regardless of which library built it.
    const double l0 = p.wavelength;
    const double bw = std::clamp(p.band_frac, 0.02f, 0.6f);
    const int nb = 5;
    std::vector<double> band;
    for (int i = 0; i < nb; ++i) {
        double t = (nb > 1) ? static_cast<double>(i) / (nb - 1) : 0.5;
        band.push_back(l0 * (1.0 - bw / 2.0 + bw * t));
    }

    // Fill the shared summary + plot globals from the two designs and their
    // verified chromatic curves. Generic over the design type (AchromaticDesign
    // and PbAchromaticDesign share the field names used here); cs/ca are always
    // std::vector<AchromaticFocalPoint>. The design object itself (different type
    // per kind) is stored by the caller before this runs.
    auto finalize = [&](const auto& sd, const auto& ad, const auto& cs, const auto& ca,
                        bool single_h, int lib_kind) {
        auto drift = [](const std::vector<AchromaticFocalPoint>& c) {
            double mn = 1e300, mx = -1e300;
            for (const auto& q : c) { mn = std::min(mn, q.focal_length_um); mx = std::max(mx, q.focal_length_um); }
            return (mx > mn) ? mx - mn : 0.0;
        };
        AchroSummary s;
        s.drift_std = drift(cs); s.drift_ach = drift(ca);
        s.gd_rms_std = sd.rms_group_delay_error_fs; s.gd_rms_ach = ad.rms_group_delay_error_fs;
        s.base_rms_std = sd.rms_phase_error_deg; s.base_rms_ach = ad.rms_phase_error_deg;
        s.gd_coverage = ad.gd_coverage; s.meanT = ad.mean_amplitude;
        s.center_wl = l0; s.focal = p.focal;
        s.single_height = single_h; s.n_cells = ad.n_cells; s.lib_kind = lib_kind;
        std::lock_guard<std::mutex> lk(g_mtx);
        g_achro_sum = s;
        g_achro_wl.clear(); g_achro_focus_std.clear(); g_achro_focus_ach.clear();
        for (const auto& q : cs) {
            g_achro_wl.push_back(static_cast<float>(q.wavelength_um * 1000.0));
            g_achro_focus_std.push_back(static_cast<float>(q.focal_length_um));
        }
        for (const auto& q : ca)
            g_achro_focus_ach.push_back(static_cast<float>(q.focal_length_um));
        g_status = std::format("Achromatic: focal drift {:.2f}->{:.2f} um, "
                               "group-delay RMS {:.2f}->{:.2f} fs",
                               s.drift_std, s.drift_ach, s.gd_rms_std, s.gd_rms_ach);
        g_progress = 1.0f;
    };

    if (p.achro_lib == 2) {
        // Achromatic Pancharatnam-Berry: a dispersive BIREFRINGENT library at one
        // etch depth (two RCWA solves per atom per band). The rotation stamps the
        // base phase exactly, so the atom is chosen purely for group delay.
        set_phase("Building dispersive birefringent PB library (2 solves/atom/band)...", 0.05f);
        auto lib = build_dispersive_pb_library(
            make_pillar(p), materials::air(), materials::air(), make_substrate(p),
            p.period, band, l0, 0.10, 0.90, /*n_fills=*/8, p.etch_height, p.harmonics);
        set_phase("Group-delay atom selection (standard + achromatic PB)...", 0.78f);
        auto sd = design_pb_achromatic_metalens(lib, p.focal, p.diameter, /*handedness=*/+1, /*gd_weight=*/0.0);
        auto ad = design_pb_achromatic_metalens(lib, p.focal, p.diameter, /*handedness=*/+1, p.gd_weight);
        set_phase("Verifying chromatic focus (stored band response)...", 0.9f);
        auto cs = verify_pb_achromatic_focus(lib, sd);
        auto ca = verify_pb_achromatic_focus(lib, ad);
        { std::lock_guard<std::mutex> lk(g_mtx); g_achro_pb = ad; }
        finalize(sd, ad, cs, ca, /*single_h=*/true, /*lib_kind=*/2);
    } else {
        // Propagation-phase achromat from a dispersive library. Two flavours:
        //   0 = fill x height grid (multi-DOF, may need a grayscale/multi-level etch)
        //   1 = single-etch shape-diverse (one height, fabricable in one etch step)
        DispersiveLibrary dl;
        if (p.achro_lib == 1) {
            set_phase("Building single-etch dispersive library (shape-diverse, one height)...", 0.05f);
            dl = build_single_etch_library(
                make_pillar(p), materials::air(), materials::air(), make_substrate(p),
                p.period, band, l0, 0.08, 0.92, /*n_fills=*/10, p.etch_height, p.harmonics);
        } else {
            set_phase("Building dispersive library (fill x height, RCWA per band)...", 0.05f);
            dl = build_dispersive_library(
                make_pillar(p), materials::air(), materials::air(), make_substrate(p),
                p.period, band, l0, 0.08, 0.92, /*n_fills=*/10, /*thick_lo=*/0.40,
                /*thick_hi=*/1.40, /*n_heights=*/5, p.harmonics);
        }
        set_phase("Two-objective atom selection (standard + achromatic)...", 0.78f);
        auto sd = design_achromatic_metalens(dl, p.focal, p.diameter, /*gd_weight=*/0.0);
        auto ad = design_achromatic_metalens(dl, p.focal, p.diameter, p.gd_weight);
        set_phase("Verifying chromatic focus (stored band response)...", 0.9f);
        auto cs = verify_achromatic_focus(dl, sd);
        auto ca = verify_achromatic_focus(dl, ad);
        { std::lock_guard<std::mutex> lk(g_mtx); g_achro = ad; }
        finalize(sd, ad, cs, ca, ad.single_height, p.achro_lib);
    }
    g_achro_pending = true;
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
