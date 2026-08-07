#include "cli.hpp"

// celeris widefov: design a HYPERBOLIC lens and a QUADRATIC (parabolic) lens for
// the SAME f/D from the SAME library, sweep the field (incidence) angle, and show
// the classic wide-FOV trade. The hyperbolic lens is perfect on-axis but develops
// coma off-axis; the quadratic lens is the same parabola recentered under tilt, so
// its focus stays sharp across a wide angular range (the spot just shifts to
// x ~ f*sin(theta)) at the cost of some on-axis quality and a curved focal surface.
int cmd_widefov(int argc, char** argv) {
    const double focal = std::atof(arg_value(argc, argv, "--focal", "30"));
    const double diameter = std::atof(arg_value(argc, argv, "--diameter", "70"));
    const double lambda = std::atof(arg_value(argc, argv, "--wavelength", "0.532"));
    const double period = std::atof(arg_value(argc, argv, "--period", "0.35"));
    const double thickness = std::atof(arg_value(argc, argv, "--thickness", "0.6"));
    const double pillar_n = std::atof(arg_value(argc, argv, "--pillar-n", "2.4"));
    const int samples = std::atoi(arg_value(argc, argv, "--fill-samples", "24"));
    const int M = std::atoi(arg_value(argc, argv, "--harmonics", "6"));
    const double max_angle = std::atof(arg_value(argc, argv, "--max-angle", "30"));
    const double angle_step = std::atof(arg_value(argc, argv, "--angle-step", "3"));
    // Aperture stop: the limiting aperture (sets resolution), placed in FRONT of
    // the metasurface so each field angle samples a decentered patch -- the
    // configuration that makes a quadratic-phase lens wide-FOV. The stop defaults
    // to the FRONT FOCAL PLANE (distance = f): there the moving patch stays
    // centered on the tilted parabola's recentered vertex, so the quadratic lens
    // sees a near-on-axis (coma-free) patch at every field angle.
    const double stop_D = std::atof(arg_value(argc, argv, "--stop-diameter", "20"));
    const double stop_d = std::atof(arg_value(argc, argv, "--stop-distance",
                                              std::to_string(focal).c_str()));
    const std::string sub_name = arg_value(argc, argv, "--substrate", "bk7");
    const Material& substrate = sub_name == "air"  ? materials::air()
                                : sub_name == "sio2" ? materials::fused_silica()
                                                     : materials::bk7();
    const char* pillar_csv = arg_value(argc, argv, "--pillar-csv", nullptr);
    const Material pillar =
        pillar_csv ? load_material_csv(pillar_csv, "pillar-csv")
                   : Material::constant(cdouble{pillar_n, 0.0}, "pillar");

    std::println("CELERIS wide-FOV comparison (hyperbolic vs quadratic phase)");
    std::println("  f={}µm  lens D={}µm  λ={}µm  Λ={}µm  h={}µm  n_pillar={}  substrate={}",
                 focal, diameter, lambda, period, thickness, pillar_n, sub_name);
    std::println("  aperture stop: D={}µm at {}µm in FRONT of the metasurface "
                 "(stop is the limiting aperture)", stop_D, stop_d);
    // The lens must be big enough to catch the off-axis walk-off of the stopped
    // beam, else the patch clips the lens edge and BOTH lenses vignette.
    const double walk = stop_d * std::tan(max_angle * pi / 180.0);
    if (walk + stop_D / 2.0 > diameter / 2.0)
        std::println("  WARNING: at {:.0f}° the stop patch (center {:.1f}µm + radius {:.1f}µm) "
                     "exceeds the lens radius {:.1f}µm -> edge vignetting; raise --diameter or "
                     "lower --stop-distance/--max-angle.", max_angle, walk, stop_D / 2.0,
                     diameter / 2.0);

    // One library, two designs -> the only difference is the target phase profile.
    bool auto_height = false;
    for (int i = 2; i < argc; ++i)
        if (std::string(argv[i]) == "--auto-height") auto_height = true;
    UnitCellLibrary lib;
    if (auto_height) {
        std::println("  auto-height: sweeping for best coverage at high |t|...");
        auto opt = optimize_height_for_2pi(pillar, materials::air(), materials::air(),
                                           substrate, period, lambda, 0.30, 1.20, 12,
                                           0.08, 0.92, std::max(samples, 32), M, 330.0);
        lib = std::move(opt.best_library);
        std::println("  chosen height {:.3f} µm -> coverage {:.0f}°, mean |t|² {:.3f}",
                     opt.best_thickness_um, opt.coverage_deg, opt.mean_transmittance);
    } else {
        std::println("  building unit-cell library ({} pillars, M={})...", samples, M);
        lib = build_unit_cell_library(pillar, materials::air(), materials::air(),
                                      substrate, period, lambda, thickness,
                                      0.08, 0.92, samples, M);
    }
    std::println("  library coverage {:.0f}°", lib.coverage() * 180.0 / pi);
    std::println("  resolution: stop sets spot ~λf/D_stop = {:.2f}µm "
                 "(the full lens would give {:.2f}µm -- the wide-FOV resolution cost)",
                 lambda * focal / stop_D, lambda * focal / diameter);

    PhaseProfile hyp; hyp.kind = PhaseProfileKind::Focusing;  hyp.focal_length_um = focal;
    PhaseProfile quad; quad.kind = PhaseProfileKind::Quadratic; quad.focal_length_um = focal;
    auto lens_h = design_metalens(lib, hyp, diameter);
    auto lens_q = design_metalens(lib, quad, diameter);

    std::vector<double> angles;
    for (double a = 0.0; a <= max_angle + 1e-9; a += angle_step) angles.push_back(a);
    auto fov_h = analyze_wide_fov(lens_h, lib, focal, lambda, diameter, stop_D, stop_d, angles);
    auto fov_q = analyze_wide_fov(lens_q, lib, focal, lambda, diameter, stop_D, stop_d, angles);

    std::println("");
    std::println("  field-of-view sweep (relative Strehl vs each lens's own on-axis peak):");
    std::println("      {:>9}  {:>11}  {:>11}   {:>12}  {:>12}", "angle(°)",
                 "hyperbolic", "quadratic", "spot_hyp(µm)", "spot_quad(µm)");
    for (std::size_t i = 0; i < angles.size(); ++i)
        std::println("      {:>9.0f}  {:>11.3f}  {:>11.3f}   {:>12.2f}  {:>12.2f}",
                     angles[i], fov_h[i].rel_strehl, fov_q[i].rel_strehl,
                     fov_h[i].spot_shift_um, fov_q[i].spot_shift_um);

    // First angle where the relative Strehl falls below `thresh`, linearly
    // interpolated; returns the last swept angle if it never drops below.
    auto fov_to = [](const std::vector<FieldPoint>& fp, double thresh) {
        for (std::size_t i = 1; i < fp.size(); ++i)
            if (fp[i].rel_strehl < thresh) {
                double a0 = fp[i - 1].angle_deg, a1 = fp[i].angle_deg;
                double s0 = fp[i - 1].rel_strehl, s1 = fp[i].rel_strehl;
                return s0 == s1 ? a0 : a0 + (a1 - a0) * (s0 - thresh) / (s0 - s1);
            }
        return fp.back().angle_deg;
    };
    // "≥" when a lens never drops below the threshold within the swept range.
    auto cap = [](const std::vector<FieldPoint>& fp, double thresh) {
        return fp.back().rel_strehl >= thresh ? "≥±" : "±";
    };
    const double h50 = fov_to(fov_h, 0.5), q50 = fov_to(fov_q, 0.5);
    const double h80 = fov_to(fov_h, 0.8), q80 = fov_to(fov_q, 0.8);
    std::println("");
    std::println("  FOV half-angle to Strehl>=0.5:  hyperbolic {}{:.1f}°   quadratic {}{:.1f}°  ->  "
                 "{:.1f}x wider", cap(fov_h, 0.5), h50, cap(fov_q, 0.5), q50,
                 h50 > 0 ? q50 / h50 : 0.0);
    std::println("  FOV half-angle to Strehl>=0.8:  hyperbolic {}{:.1f}°   quadratic {}{:.1f}°",
                 cap(fov_h, 0.8), h80, cap(fov_q, 0.8), q80);
    std::println("  spot landing: quadratic spot tracks ~f*sin(θ); hyperbolic chief ray ~f*tan(θ).");
    std::println("  NOTE: the quadratic lens is wide-FOV ONLY with this offset stop (each angle");
    std::println("        sees a decentered low-NA patch -> a recentered parabola = sharp shifted");
    std::println("        focus). The price: resolution set by the stop (above) and a curved");
    std::println("        (Petzval) focal surface. Without the stop both lenses degrade alike.");
    return 0;
}

// celeris efficiency: the per-meta-atom energy budget -- where does the incident
// power go? Reflection, transmission, material absorption (1-R-T), and the split
// between the useful zeroth order and higher (stray) diffraction orders. A
// subwavelength pitch has ONLY order 0 propagating (no stray light by design);
// a lossy pillar (metal) shows real absorption. The per-element analog of a
// Zemax surface-efficiency / ghost budget.
int cmd_efficiency(int argc, char** argv) {
    const double lambda = std::atof(arg_value(argc, argv, "--wavelength", "0.532"));
    const double period = std::atof(arg_value(argc, argv, "--period", "0.35"));
    const double thickness = std::atof(arg_value(argc, argv, "--thickness", "0.6"));
    const double pillar_n = std::atof(arg_value(argc, argv, "--pillar-n", "2.4"));
    const double fill = std::atof(arg_value(argc, argv, "--fill", "0.5"));
    const int M = std::atoi(arg_value(argc, argv, "--harmonics", "8"));
    const std::string pol = arg_value(argc, argv, "--pol", "x");

    const Material substrate = resolve_substrate(argc, argv);
    const Material pillar = resolve_pillar(argc, argv, pillar_n);

    const std::string shape_name = arg_value(argc, argv, "--shape", "square");
    const double shape_param = std::atof(arg_value(argc, argv, "--shape-param", "0.5"));
    const MetaShape shape = shape_name == "circle" || shape_name == "ellipse"
                                ? MetaShape::Ellipse
                            : shape_name == "cross" ? MetaShape::Cross
                            : shape_name == "ring"  ? MetaShape::Ring
                                                    : MetaShape::Rectangle;
    const cdouble Ex0 = pol == "y" ? cdouble{0.0, 0.0} : cdouble{1.0, 0.0};
    const cdouble Ey0 = pol == "y" ? cdouble{1.0, 0.0} : cdouble{0.0, 0.0};

    std::println("CELERIS meta-atom efficiency budget");
    std::println("  λ={}µm  Λ={}µm  h={}µm  fill={}  shape={}  pol={}  M={}",
                 lambda, period, thickness, fill, shape_name, pol, M);
    std::println("  pillar={} (n={:.3f}+{:.3f}i)  substrate={} (n={:.3f})",
                 pillar.name(), pillar.index(lambda).real(),
                 pillar.index(lambda).imag(), substrate.name(),
                 substrate.index(lambda).real());
    const bool sub_wl = period < lambda;  // subwavelength pitch?
    std::println("  pitch is {} (Λ/λ={:.3f}) -> {} propagating diffraction orders expected",
                 sub_wl ? "SUBWAVELENGTH" : "NOT subwavelength", period / lambda,
                 sub_wl ? "only the 0th" : "higher (stray)");

    auto b = analyze_efficiency(materials::air(), pillar, materials::air(),
                                substrate, period, period, fill, fill, thickness,
                                lambda, Ex0, Ey0, M, shape, shape_param);

    std::println("");
    std::println("  energy budget (fractions of incident power):");
    std::println("      transmission : {:.4f}   (0th order {:.4f} + stray {:.4f})",
                 b.transmission, b.t_zero, b.t_stray);
    std::println("      reflection   : {:.4f}   (0th order {:.4f} + stray {:.4f})",
                 b.reflection, b.r_zero, b.r_stray);
    std::println("      absorption   : {:.4f}   (= 1 - R - T; material loss)",
                 b.absorption);
    std::println("      sum check    : {:.6f}   (R+T+A)",
                 b.transmission + b.reflection + b.absorption);
    std::println("  propagating channels: {} transmitted, {} reflected",
                 b.n_prop_t, b.n_prop_r);

    // Per-order table for the propagating transmitted orders + any stray.
    std::println("");
    std::println("  per-order split (propagating orders, by transmitted power):");
    std::println("      {:>4} {:>4}  {:>10}  {:>10}  {}", "p", "q", "T_order",
                 "R_order", "channel");
    int shown = 0;
    for (const auto& o : b.orders) {
        if (!o.prop_t && !o.prop_r) continue;       // skip evanescent bookkeeping
        if (o.de_t < 1e-6 && o.de_r < 1e-6) continue;
        std::println("      {:>4} {:>4}  {:>10.5f}  {:>10.5f}  {}", o.p, o.q,
                     o.de_t, o.de_r,
                     (o.p == 0 && o.q == 0) ? "0th (useful)" : "stray");
        if (++shown >= 12) { std::println("      ... (further orders omitted)"); break; }
    }
    if (shown == 0)
        std::println("      (none above 1e-6)");

    // Interpretation line: a good metalens atom is high-T, low absorption, all
    // the transmitted power in the 0th order (subwavelength pitch suppresses stray).
    const double useful = b.transmission > 0 ? b.t_zero / b.transmission : 0.0;
    std::println("");
    std::println("  -> {:.1f}% of transmitted power is in the useful 0th order; "
                 "absorption {:.1f}%, reflection {:.1f}%.",
                 useful * 100.0, b.absorption * 100.0, b.reflection * 100.0);
    return 0;
}

// celeris fieldmap: full-field PSF/Strehl grid. Designs a focusing lens, then
// sweeps a 2D grid of incidence (field) angles (theta_x, theta_y) and reports
// the relative Strehl + tangential/sagittal spot FWHM at each -- the full-field
// quality map (not just the center-row spot-vs-field of `design --fov`). A
// single hyperbolic lens is sharp on-axis and aberrates (coma) off-axis, so the
// map falls off toward the corners.
int cmd_fieldmap(int argc, char** argv) {
    const double focal = std::atof(arg_value(argc, argv, "--focal", "50"));
    const double diameter = std::atof(arg_value(argc, argv, "--diameter", "20"));
    const double lambda = std::atof(arg_value(argc, argv, "--wavelength", "0.532"));
    const double period = std::atof(arg_value(argc, argv, "--period", "0.35"));
    const double thickness = std::atof(arg_value(argc, argv, "--thickness", "0.6"));
    const double pillar_n = std::atof(arg_value(argc, argv, "--pillar-n", "2.4"));
    const int samples = std::atoi(arg_value(argc, argv, "--fill-samples", "18"));
    const int M = std::atoi(arg_value(argc, argv, "--harmonics", "6"));
    const double max_angle = std::atof(arg_value(argc, argv, "--max-angle", "12"));
    const int n_half = std::atoi(arg_value(argc, argv, "--angle-steps", "4"));
    const int psf_n = std::atoi(arg_value(argc, argv, "--psf-n", "101"));
    const char* out_pgm = arg_value(argc, argv, "--out", nullptr);

    const Material substrate = resolve_substrate(argc, argv);
    const Material pillar = resolve_pillar(argc, argv, pillar_n);

    std::println("CELERIS field-resolved analysis (full-field Strehl/FWHM grid)");
    std::println("  f={}µm  D={}µm  λ={}µm  Λ={}µm  h={}µm  pillar={}  substrate={}",
                 focal, diameter, lambda, period, thickness, pillar.name(),
                 substrate.name());
    std::println("  building unit-cell library ({} pillars, M={})...", samples, M);
    auto lib = build_unit_cell_library(pillar, materials::air(), materials::air(),
                                       substrate, period, lambda, thickness,
                                       0.08, 0.92, samples, M);
    auto lens = design_metalens(lib, focal, diameter);
    std::println("  designed {0}x{0} pillars, RMS phase error {1:.1f}°, mean |t| {2:.3f}",
                 lens.n_cells, lens.rms_phase_error_deg, lens.mean_amplitude);
    std::println("  sweeping a {0}x{0} field grid to ±{1:.1f}° (psf {2}x{2})...",
                 2 * n_half + 1, max_angle, psf_n);

    auto g = analyze_field_grid(lens, lib, focal, lambda, diameter, max_angle,
                                n_half, psf_n);

    // Relative-Strehl grid (rows = theta_y, cols = theta_x).
    std::println("");
    std::println("  relative Strehl across the field (rows θy, cols θx, degrees):");
    std::print("        θx:");
    for (int jx = -n_half; jx <= n_half; ++jx)
        std::print(" {:>6.1f}", n_half > 0 ? max_angle * jx / n_half : 0.0);
    std::println("");
    for (int jy = 0; jy < g.n; ++jy) {
        std::print("    {:>6.1f} :", g.points[(std::size_t)jy * g.n].theta_y_deg);
        for (int jx = 0; jx < g.n; ++jx)
            std::print(" {:>6.3f}", g.points[(std::size_t)jy * g.n + jx].rel_strehl);
        std::println("");
    }

    // Diagonal FWHM trend (on-axis to corner) — shows the off-axis spot blowing up.
    std::println("");
    std::println("  spot FWHM along the diagonal (on-axis -> corner):");
    std::println("      {:>8} {:>8}  {:>10}  {:>10}  {:>10}", "θx(°)", "θy(°)",
                 "rel.Strehl", "FWHMx(µm)", "FWHMy(µm)");
    for (int d = 0; d <= n_half; ++d) {
        int jx = n_half + d, jy = n_half + d;  // center is (n_half,n_half)
        const auto& pnt = g.points[(std::size_t)jy * g.n + jx];
        std::println("      {:>8.1f} {:>8.1f}  {:>10.3f}  {:>10.3f}  {:>10.3f}",
                     pnt.theta_x_deg, pnt.theta_y_deg, pnt.rel_strehl,
                     pnt.fwhm_x_um, pnt.fwhm_y_um);
    }

    // Field half-angle to Strehl >= 0.8 / 0.5 along +x (a single-number FOV).
    auto fov_to = [&](double thr) {
        double last = 0.0;
        for (int jx = n_half; jx <= 2 * n_half; ++jx) {  // center..+x edge
            const auto& pnt = g.points[(std::size_t)n_half * g.n + jx];
            if (pnt.rel_strehl >= thr) last = pnt.theta_x_deg;
            else break;
        }
        return last;
    };
    std::println("");
    std::println("  on-axis-row FOV half-angle: to Strehl>=0.8 = {:.1f}°, "
                 ">=0.5 = {:.1f}° (sweep cap ±{:.1f}°)",
                 fov_to(0.8), fov_to(0.5), max_angle);

    if (out_pgm) {
        // Write the rel-Strehl grid as a small PGM heatmap (bright = sharp).
        std::vector<double> img(g.points.size());
        for (std::size_t i = 0; i < g.points.size(); ++i)
            img[i] = g.points[i].rel_strehl;
        if (write_pgm(out_pgm, g.n, g.n, img, 1.0))
            std::println("  wrote field Strehl map ({0}x{0}) -> {1}", g.n, out_pgm);
        else
            std::println("  ERROR: could not write {}", out_pgm);
    }
    return 0;
}

// celeris design: build a unit-cell library, design a focusing metalens, export
// GDSII, and report design fidelity + focal performance.
int cmd_design(int argc, char** argv) {
    const double focal = std::atof(arg_value(argc, argv, "--focal", "50"));
    const double diameter = std::atof(arg_value(argc, argv, "--diameter", "20"));
    const double lambda = std::atof(arg_value(argc, argv, "--wavelength", "0.532"));
    const double period = std::atof(arg_value(argc, argv, "--period", "0.35"));
    const double thickness = std::atof(arg_value(argc, argv, "--thickness", "0.6"));
    const double pillar_n = std::atof(arg_value(argc, argv, "--pillar-n", "2.4"));
    const int samples = std::atoi(arg_value(argc, argv, "--fill-samples", "18"));
    const int M = std::atoi(arg_value(argc, argv, "--harmonics", "6"));
    const std::string out = arg_value(argc, argv, "--out", "metalens.gds");

    const Material substrate = resolve_substrate(argc, argv);
    const Material pillar = resolve_pillar(argc, argv, pillar_n);

    // --profile: the target wavefront. Default focusing (the hyperbolic lens);
    // vortex/deflector/axicon/freeform stamp other profiles on the SAME
    // propagation-phase path (vary pillar size to hit phi(x,y) via the library).
    const std::string profile_name = arg_value(argc, argv, "--profile", "focusing");
    auto profile_opt = parse_phase_profile(argc, argv, focal, diameter);
    if (!profile_opt) return 1;
    const PhaseProfile profile = *profile_opt;

    std::println("CELERIS metalens design");
    std::println("  profile={}  f={}µm  D={}µm  λ={}µm  Λ={}µm  h={}µm  pillar={} (n={:.3f}+{:.3f}i)  substrate={}",
                 profile_name, focal, diameter, lambda, period, thickness,
                 pillar.name(), pillar.index(lambda).real(),
                 pillar.index(lambda).imag(), substrate.name());
    // --auto-height: instead of using the fixed --thickness, sweep pillar height
    // and pick the (shortest, highest-transmittance) single etch depth that
    // reaches full 2pi phase coverage. This lifts the transmission-weighted
    // Strehl that a too-short pillar caps. The chosen height is reported and used.
    bool auto_height = false;
    for (int i = 2; i < argc; ++i)
        if (std::string(argv[i]) == "--auto-height") auto_height = true;

    // --shape: meta-atom cross-section. square (default, analytic & fastest),
    // circle/ellipse, cross, or ring -- the non-rectangular ones use the grid
    // (Laurent) factorization. --shape-param tunes cross arm width / ring inner-r.
    const std::string shape_name = arg_value(argc, argv, "--shape", "square");
    const double shape_param = std::atof(arg_value(argc, argv, "--shape-param", "0.5"));
    const MetaShape shape = shape_name == "circle" || shape_name == "ellipse"
                                ? MetaShape::Ellipse
                            : shape_name == "cross" ? MetaShape::Cross
                            : shape_name == "ring"  ? MetaShape::Ring
                                                    : MetaShape::Rectangle;
    if (shape != MetaShape::Rectangle)
        std::println("  meta-atom shape: {} (param {}) -- grid Laurent factorization",
                     shape_name, shape_param);

    double used_thickness = thickness;
    UnitCellLibrary lib;
    if (auto_height) {
        const double h_lo = std::atof(arg_value(argc, argv, "--height-lo", "0.30"));
        const double h_hi = std::atof(arg_value(argc, argv, "--height-hi", "1.20"));
        const int n_h = std::atoi(arg_value(argc, argv, "--height-steps", "12"));
        // Coverage is sampling-limited when the phase wraps faster than the fill
        // grid resolves, so sweep fills at least moderately densely here.
        const int sweep_fills = std::max(samples, 32);
        std::println("  auto-height: sweeping {} heights in [{}, {}] µm "
                     "(fill samples {}) for best coverage at high |t|...",
                     n_h, h_lo, h_hi, sweep_fills);
        auto opt = optimize_height_for_2pi(pillar, materials::air(), materials::air(),
                                           substrate, period, lambda, h_lo, h_hi,
                                           n_h, 0.08, 0.92, sweep_fills, M, 330.0,
                                           shape, shape_param);
        std::println("      {:>10}  {:>12}  {:>14}", "height(µm)", "coverage", "mean |t|²");
        for (auto& e : opt.sweep)
            std::println("      {:>10.3f}  {:>11.0f}°  {:>14.3f}",
                         e.thickness_um, e.coverage_deg, e.mean_transmittance);
        used_thickness = opt.best_thickness_um;
        lib = std::move(opt.best_library);
        std::println("  chosen height {:.3f} µm -> coverage {:.0f}° ({}), mean |t|² {:.3f}",
                     used_thickness, opt.coverage_deg,
                     opt.reached_target ? "clears target" : "best available", opt.mean_transmittance);
    } else {
        std::println("  building unit-cell library ({} pillars, M={})...", samples, M);
        lib = build_unit_cell_library(pillar, materials::air(), materials::air(),
                                      substrate, period, lambda, thickness,
                                      0.08, 0.92, samples, M, shape, shape_param);
    }
    std::println("  library phase coverage: {:.0f}°  (effective {:.0f}°)",
                 lib.phase_span() * 180.0 / pi, lib.coverage() * 180.0 / pi);

    auto lens = design_metalens(lib, profile, diameter);
    std::println("  designed {0}x{0} pillars, RMS phase error {1:.1f}°, mean |t| {2:.3f}",
                 lens.n_cells, lens.rms_phase_error_deg, lens.mean_amplitude);

    int np = write_metalens_gds(lens, out);
    if (np < 0) { std::println("  ERROR: could not write {}", out); return 1; }
    std::println("  wrote {} pillars -> {}", np, out);
    if (shape != MetaShape::Rectangle)
        std::println("  NOTE: the {} library was used for the PHYSICS, but the GDS writes "
                     "square footprints; shape-aware (polygon) GDS export is a TODO.",
                     shape_name);

    // Non-focusing profiles on the propagation-phase path: the focal-spot battery
    // below assumes a focusing lens, so run the profile's own propagation proof
    // (the same one `pbdesign` uses) and finish here. The propagation path carries
    // a finite library phase error + non-uniform amplitude -- the geometric-phase
    // `pbdesign` path stamps these same profiles exactly.
    if (profile.kind != PhaseProfileKind::Focusing) {
        const double pp = lens.period_um;
        const double cen = (lens.n_cells - 1) / 2.0;
        const double R_ap = diameter / 2.0;
        std::vector<double> px, py;
        std::vector<cdouble> tc;
        for (int iy = 0; iy < lens.n_cells; ++iy)
            for (int ix = 0; ix < lens.n_cells; ++ix) {
                double x = (ix - cen) * pp, y = (iy - cen) * pp;
                if (std::sqrt(x * x + y * y) > R_ap) continue;
                px.push_back(x);
                py.push_back(y);
                tc.push_back(lib.transmission_for_fill(
                    lens.fill_map[(std::size_t)iy * lens.n_cells + ix]));
            }
        const double recon_z = std::atof(arg_value(argc, argv, "--recon-z",
                                                   std::to_string(focal).c_str()));
        ProfileProof proof = profile_optical_proof(px, py, tc, profile, lambda, focal,
                                                   diameter, recon_z);
        std::println("  optical check: {}", proof.summary);
        std::println("  NOTE: propagation-phase profile -- finite library => RMS phase error "
                     "{:.1f}° + non-uniform |t| (`pbdesign` stamps these profiles exactly).",
                     lens.rms_phase_error_deg);

        const char* psf_path = arg_value(argc, argv, "--psf", nullptr);
        const char* rprefix = arg_value(argc, argv, "--report", nullptr);
        if (psf_path || rprefix) {
            auto psf = propagate_pillars(px, py, tc, proof.psf_cx, proof.psf_cy,
                                         proof.psf_z, lambda, 201, proof.psf_hw);
            if (psf_path && write_pgm(psf_path, psf.n, psf.n, psf.intensity, 2.2))
                std::println("  wrote reconstruction image ({0}x{0}) -> {1}", psf.n, psf_path);
            if (rprefix) {
                std::string base = rprefix;
                std::ofstream f(base + "_report.txt");
                bool ok = static_cast<bool>(f);
                if (f) {
                    f << "CELERIS metalens design report (propagation phase)\n";
                    f << "==================================================\n\n";
                    f << std::format("profile           : {}\n", profile_name);
                    f << std::format("wavelength        : {} um\n", lambda);
                    f << std::format("aperture diameter : {} um\n", diameter);
                    f << std::format("period            : {} um\n", period);
                    f << std::format("pillar height     : {} um{}\n", used_thickness,
                                     auto_height ? " (auto-selected)" : "");
                    f << std::format("array             : {0} x {0} pillars\n", lens.n_cells);
                    f << std::format("RMS phase error   : {:.1f} deg\n", lens.rms_phase_error_deg);
                    f << std::format("mean transmission : {:.3f}\n", lens.mean_amplitude);
                    f << std::format("optical check     : {}\n", proof.summary);
                }
                ok &= write_pgm(base + "_recon.pgm", psf.n, psf.n, psf.intensity, 2.2);
                ok &= (write_metalens_gds(lens, base + "_layout.gds") >= 0);
                std::println("  report bundle -> {0}_report.txt (+ _recon.pgm, _layout.gds)  ok={1}",
                             base, ok);
            }
        }
        return 0;
    }

    auto foc = analyze_focus(lens, lib, focal, lambda, diameter);
    std::println("  focal: Strehl {:.3f}, FWHM {:.2f}µm (diffraction limit {:.2f}µm), "
                 "encircled {:.0f}%",
                 foc.strehl, foc.fwhm_um, foc.diffraction_limit_um,
                 foc.encircled_energy * 100.0);

    bool tolerance = false;
    for (int i = 2; i < argc; ++i)
        if (std::string(argv[i]) == "--tolerance") tolerance = true;
    if (tolerance) {
        std::println("  fabrication tolerance (Monte-Carlo CD error, Strehl):");
        std::println("      {:>8}  {:>10}  {:>10}  {:>10}", "σ(nm)", "mean", "std",
                     "worst");
        auto tol = analyze_tolerance(lens, lib, focal, lambda, diameter,
                                     {0.0, 5.0, 10.0, 20.0}, /*trials=*/12,
                                     /*seed=*/12345);
        for (auto& t : tol)
            std::println("      {:>8.0f}  {:>10.3f}  {:>10.3f}  {:>10.3f}",
                         t.sigma_nm, t.mean_strehl, t.std_strehl, t.worst_strehl);
    }

    bool fov = false;
    for (int i = 2; i < argc; ++i)
        if (std::string(argv[i]) == "--fov") fov = true;
    if (fov) {
        std::println("  field of view (off-axis focus quality):");
        std::println("      {:>8}  {:>12}  {:>12}", "angle(°)", "rel.Strehl",
                     "shift(µm)");
        for (auto& f : analyze_field_of_view(lens, lib, focal, lambda, diameter,
                                             {0.0, 1.0, 2.0, 5.0, 10.0}))
            std::println("      {:>8.0f}  {:>12.3f}  {:>12.2f}", f.angle_deg,
                         f.rel_strehl, f.spot_shift_um);
    }

    const char* psf_path = arg_value(argc, argv, "--psf", nullptr);
    if (psf_path) {
        const double dl = lambda * focal / diameter;
        auto psf = compute_psf(lens, lib, focal, lambda, diameter, /*n=*/201,
                               /*half_window=*/std::max(5.0 * dl, 4.0));
        // gamma 2.2 brightens the faint Airy rings for visibility.
        if (write_pgm(psf_path, psf.n, psf.n, psf.intensity, /*gamma=*/2.2))
            std::println("  wrote focal PSF image ({0}x{0}) -> {1}", psf.n, psf_path);
        else
            std::println("  ERROR: could not write {}", psf_path);
    }

    // --report <prefix>: write a full deliverable bundle (metrics txt + PSF and
    // caustic images + GDS), the same artifacts the GUI's Save Report produces.
    const char* report_prefix = arg_value(argc, argv, "--report", nullptr);
    if (report_prefix) {
        std::string base = report_prefix;
        const double dl = lambda * focal / diameter;
        auto wf = analyze_wavefront(lens, lib, focal, lambda, diameter);
        auto tf = analyze_through_focus(lens, lib, focal, lambda, diameter);
        auto psf = compute_psf(lens, lib, focal, lambda, diameter, 201,
                               std::max(5.0 * dl, 4.0));
        std::ofstream f(base + "_report.txt");
        bool okall = static_cast<bool>(f);
        if (f) {
            f << "CELERIS metalens design report\n==============================\n\n";
            f << std::format("focal length      : {} um\n", focal);
            f << std::format("aperture diameter : {} um\n", diameter);
            f << std::format("wavelength        : {} um\n", lambda);
            f << std::format("period            : {} um\n", period);
            f << std::format("pillar height     : {} um{}\n\n", used_thickness,
                             auto_height ? " (auto-selected for full 2pi)" : "");
            f << std::format("array             : {0} x {0} ({1} pillars)\n",
                             lens.n_cells, lens.n_cells * lens.n_cells);
            f << std::format("phase coverage    : {:.0f} deg\n", lib.phase_span() * 180.0 / pi);
            f << std::format("RMS phase error   : {:.1f} deg\n", lens.rms_phase_error_deg);
            f << std::format("mean transmission : {:.3f}\n\n", lens.mean_amplitude);
            f << std::format("Strehl ratio      : {:.3f}\n", foc.strehl);
            f << std::format("spot FWHM         : {:.3f} um\n", foc.fwhm_um);
            f << std::format("diffraction limit : {:.3f} um\n", foc.diffraction_limit_um);
            f << std::format("encircled energy  : {:.1f} %\n", foc.encircled_energy * 100.0);
            f << std::format("wavefront RMS     : {:.4f} waves\n", wf.rms_waves);
            f << std::format("depth of focus    : {:.2f} um\n", tf.dof_um);
        }
        okall &= write_pgm(base + "_psf.pgm", psf.n, psf.n, psf.intensity, 2.2);
        if (!tf.caustic.empty()) {
            std::vector<double> cd(tf.caustic.begin(), tf.caustic.end());
            okall &= write_pgm(base + "_caustic.pgm", tf.caustic_nx, tf.caustic_nz, cd, 1.6);
        }
        okall &= (write_metalens_gds(lens, base + "_layout.gds") >= 0);
        std::println("  report bundle -> {0}_report.txt (+ _psf.pgm, _caustic.pgm, "
                     "_layout.gds)  ok={1}", base, okall);
    }
    return 0;
}
