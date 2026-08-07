#include "cli.hpp"

// celeris birefringence: sweep a rectangular pillar's x-width and report the
// form birefringence (x- vs y-polarized phase + retardance). Demonstrates the
// polarization-optics building block from the vectorial solver.
int cmd_birefringence(int argc, char** argv) {
    const double lambda = std::atof(arg_value(argc, argv, "--wavelength", "0.532"));
    const double period = std::atof(arg_value(argc, argv, "--period", "0.35"));
    const double thickness = std::atof(arg_value(argc, argv, "--thickness", "0.6"));
    const double pillar_n = std::atof(arg_value(argc, argv, "--pillar-n", "2.4"));
    const double fill_y = std::atof(arg_value(argc, argv, "--fill-y", "0.5"));
    const int samples = std::atoi(arg_value(argc, argv, "--samples", "15"));
    const int M = std::atoi(arg_value(argc, argv, "--harmonics", "6"));
    const Material pillar = Material::constant(cdouble{pillar_n, 0.0}, "pillar");

    std::println("CELERIS form-birefringence sweep");
    std::println("  lambda={}um  period={}um  height={}um  n_pillar={}  fill_y={}",
                 lambda, period, thickness, pillar_n, fill_y);
    std::println("  {:>8}  {:>10}  {:>10}  {:>12}  {:>6}  {:>6}", "fill_x", "phase_x",
                 "phase_y", "retardance", "|tx|", "|ty|");
    auto pts = analyze_birefringence(pillar, materials::air(), materials::air(),
                                     materials::bk7(), period, lambda, thickness,
                                     fill_y, 0.10, 0.90, samples, M);
    for (const auto& p : pts)
        std::println("  {:>8.3f}  {:>9.1f}d  {:>9.1f}d  {:>11.1f}d  {:>6.3f}  {:>6.3f}",
                     p.fill_x, p.phase_x_deg, p.phase_y_deg, p.retardance_deg, p.tx,
                     p.ty);

    // Sanity: at fill_x == fill_y the pillar is square -> retardance ~ 0.
    double min_asym = 1e9;
    for (const auto& p : pts)
        min_asym = std::min(min_asym, std::abs(p.fill_x - fill_y));
    std::println("  (square-pillar retardance should be ~0; max retardance shows the "
                 "achievable waveplate range)");
    return 0;
}

// celeris polardesign: a polarization-multiplexed metalens -- x-polarized light
// focuses at one distance, y-polarized at another -- built on the rectangular-
// pillar (fill_x, fill_y) library. Writes a rectangular-pillar GDS.
int cmd_polardesign(int argc, char** argv) {
    const double focal_x = std::atof(arg_value(argc, argv, "--focal-x", "50"));
    const double focal_y = std::atof(arg_value(argc, argv, "--focal-y", "80"));
    const double diameter = std::atof(arg_value(argc, argv, "--diameter", "20"));
    const double lambda = std::atof(arg_value(argc, argv, "--wavelength", "0.532"));
    const double period = std::atof(arg_value(argc, argv, "--period", "0.35"));
    const double thickness = std::atof(arg_value(argc, argv, "--thickness", "0.6"));
    const double pillar_n = std::atof(arg_value(argc, argv, "--pillar-n", "2.4"));
    const int samples = std::atoi(arg_value(argc, argv, "--samples", "12"));
    const int M = std::atoi(arg_value(argc, argv, "--harmonics", "6"));
    const std::string out = arg_value(argc, argv, "--out", "polar_metalens.gds");
    const Material pillar = Material::constant(cdouble{pillar_n, 0.0}, "pillar");

    std::println("CELERIS polarization-multiplexed metalens");
    std::println("  focal_x(X-pol)={}um  focal_y(Y-pol)={}um  D={}um  lambda={}um",
                 focal_x, focal_y, diameter, lambda);
    std::println("  building (fill_x, fill_y) library ({0}x{0} = {1} pillars, "
                 "2 solves each)...", samples, samples * samples);

    auto lib = build_polarization_library(pillar, materials::air(), materials::air(),
                                          materials::bk7(), period, lambda,
                                          thickness, 0.10, 0.90, samples, M);
    auto d = design_polarization_metalens(lib, focal_x, focal_y, diameter);
    std::println("  designed {0}x{0} rectangular pillars", d.n_cells);
    std::println("  X-pol: RMS phase error {:.1f} deg, mean |t| {:.3f}",
                 d.rms_phase_error_x_deg, d.mean_amp_x);
    std::println("  Y-pol: RMS phase error {:.1f} deg, mean |t| {:.3f}",
                 d.rms_phase_error_y_deg, d.mean_amp_y);

    // Aperture pillar list with per-polarization transmission.
    const double pp = d.period_um;
    const double cen = (d.n_cells - 1) / 2.0;
    const double R_ap = diameter / 2.0;
    std::vector<double> px, py;
    std::vector<cdouble> tx, ty;
    for (int iy = 0; iy < d.n_cells; ++iy)
        for (int ix = 0; ix < d.n_cells; ++ix) {
            double x = (ix - cen) * pp, y = (iy - cen) * pp;
            if (std::sqrt(x * x + y * y) > R_ap) continue;
            std::size_t off = (std::size_t)iy * d.n_cells + ix;
            px.push_back(x); py.push_back(y);
            tx.push_back(d.t_x[off]); ty.push_back(d.t_y[off]);
        }

    // Optical proof: scan on-axis intensity along z for each polarization and
    // confirm the peaks land at the two target focal planes (the bifocal claim).
    const double k = 2.0 * pi / lambda;
    auto peak_z = [&](const std::vector<cdouble>& t) {
        double zlo = 0.5 * std::min(focal_x, focal_y);
        double zhi = 1.5 * std::max(focal_x, focal_y);
        const int NZ = 240;
        double best_z = zlo, best_I = -1.0;
        for (int j = 0; j < NZ; ++j) {
            double z = zlo + (zhi - zlo) * j / (NZ - 1);
            cdouble E{0, 0};
            for (std::size_t q = 0; q < px.size(); ++q) {
                double r = std::sqrt(px[q] * px[q] + py[q] * py[q] + z * z);
                E += t[q] * std::polar(1.0 / r, k * r);
            }
            double I = std::norm(E);
            if (I > best_I) { best_I = I; best_z = z; }
        }
        return best_z;
    };
    double zx = peak_z(tx), zy = peak_z(ty);
    std::println("  optical check (on-axis focus):");
    std::println("    X-pol focuses at z = {:.1f} um  (target {:.1f})", zx, focal_x);
    std::println("    Y-pol focuses at z = {:.1f} um  (target {:.1f})", zy, focal_y);

    // Focal isolation: how much stronger each channel focuses at its OWN plane
    // than the other channel's light does there. This is the key spec for a
    // polarization-multiplexed lens (channel cross-talk).
    auto on_axis_I = [&](const std::vector<cdouble>& t, double z) {
        cdouble E{0, 0};
        for (std::size_t q = 0; q < px.size(); ++q) {
            double r = std::sqrt(px[q] * px[q] + py[q] * py[q] + z * z);
            E += t[q] * std::polar(1.0 / r, k * r);
        }
        return std::norm(E);
    };
    if (std::abs(focal_x - focal_y) > 1e-6) {
        double iso_x = 10.0 * std::log10(on_axis_I(tx, zx) / std::max(on_axis_I(ty, zx), 1e-30));
        double iso_y = 10.0 * std::log10(on_axis_I(ty, zy) / std::max(on_axis_I(tx, zy), 1e-30));
        std::println("    focal isolation: X-plane {:.1f} dB, Y-plane {:.1f} dB "
                     "(channel separation)", iso_x, iso_y);
    }

    int np = write_rect_gds(out, d.n_cells, d.period_um, d.fill_x, d.fill_y);
    if (np < 0) { std::println("  ERROR: could not write {}", out); return 1; }
    std::println("  wrote {} rectangular pillars -> {}", np, out);

    // --report <prefix>: metrics txt + per-polarization focal PSF images + GDS.
    const char* rp = arg_value(argc, argv, "--report", nullptr);
    if (rp) {
        std::string base = rp;
        double dlx = lambda * focal_x / diameter, dly = lambda * focal_y / diameter;
        auto psfx = propagate_pillars(px, py, tx, 0, 0, focal_x, lambda, 201,
                                      std::max(5.0 * dlx, 4.0));
        auto psfy = propagate_pillars(px, py, ty, 0, 0, focal_y, lambda, 201,
                                      std::max(5.0 * dly, 4.0));
        bool ok = write_pgm(base + "_xpol_psf.pgm", psfx.n, psfx.n, psfx.intensity, 2.2);
        ok &= write_pgm(base + "_ypol_psf.pgm", psfy.n, psfy.n, psfy.intensity, 2.2);
        std::ofstream f(base + "_polar_report.txt");
        if (f) {
            f << "CELERIS polarization-multiplexed metalens report\n";
            f << "=================================================\n\n";
            f << std::format("wavelength        : {} um\n", lambda);
            f << std::format("aperture diameter : {} um\n", diameter);
            f << std::format("array             : {0} x {0} rectangular pillars\n\n", d.n_cells);
            f << std::format("X-pol focal target: {} um   measured focus: {:.1f} um\n", focal_x, zx);
            f << std::format("X-pol RMS phase   : {:.1f} deg   mean |t|: {:.3f}\n\n",
                             d.rms_phase_error_x_deg, d.mean_amp_x);
            f << std::format("Y-pol focal target: {} um   measured focus: {:.1f} um\n", focal_y, zy);
            f << std::format("Y-pol RMS phase   : {:.1f} deg   mean |t|: {:.3f}\n",
                             d.rms_phase_error_y_deg, d.mean_amp_y);
        } else ok = false;
        std::println("  report bundle -> {0}_polar_report.txt (+ _xpol_psf.pgm, "
                     "_ypol_psf.pgm)  ok={1}", base, ok);
    }
    return 0;
}

// celeris pbdesign: a Pancharatnam-Berry (geometric-phase) focusing metalens.
// One fixed half-wave-plate meta-atom, ROTATED per site, imprints the focusing
// phase as 2*theta on the spin-flipped (cross-circular) output. Geometric phase
// is exact and amplitude is uniform -> phase-error-free, diffraction-limited,
// capped only by the atom's polarization conversion. RCWA-verifies the 2*theta
// relation and writes a rotated-pillar GDS.
int cmd_pbdesign(int argc, char** argv) {
    const double focal = std::atof(arg_value(argc, argv, "--focal", "50"));
    const double diameter = std::atof(arg_value(argc, argv, "--diameter", "20"));
    const double lambda = std::atof(arg_value(argc, argv, "--wavelength", "0.532"));
    const double period = std::atof(arg_value(argc, argv, "--period", "0.35"));
    const double thickness = std::atof(arg_value(argc, argv, "--thickness", "0.6"));
    const double pillar_n = std::atof(arg_value(argc, argv, "--pillar-n", "2.4"));
    const int samples = std::atoi(arg_value(argc, argv, "--samples", "12"));
    const int M = std::atoi(arg_value(argc, argv, "--harmonics", "6"));
    const int handedness = std::atoi(arg_value(argc, argv, "--handedness", "1")) >= 0 ? 1 : -1;
    const std::string out = arg_value(argc, argv, "--out", "pb_metalens.gds");
    const Material pillar = Material::constant(cdouble{pillar_n, 0.0}, "pillar");

    // The geometric-phase profile to imprint. The rotation map is identical for
    // all of them -- only the target phase phi(x,y) changes.
    const std::string profile_name = arg_value(argc, argv, "--profile", "focusing");
    auto profile_opt = parse_phase_profile(argc, argv, focal, diameter);
    if (!profile_opt) return 1;
    PhaseProfile profile = *profile_opt;

    std::println("CELERIS Pancharatnam-Berry (geometric-phase) metasurface");
    std::println("  profile={}  D={}um  lambda={}um  illumination={}",
                 profile_name, diameter, lambda, handedness > 0 ? "RCP" : "LCP");
    switch (profile.kind) {
        case PhaseProfileKind::Focusing:
            std::println("  focusing lens: f={}um", focal); break;
        case PhaseProfileKind::Vortex:
            std::println("  vortex/OAM: charge l={}{}", profile.topological_charge,
                         focal > 0 ? std::format("  (focused at f={}um)", focal)
                                   : "  (collimated)"); break;
        case PhaseProfileKind::Deflector:
            std::println("  beam deflector: {}deg toward azimuth {}deg",
                         profile.deflect_deg, profile.deflect_azimuth_deg); break;
        case PhaseProfileKind::Axicon:
            std::println("  axicon: cone half-angle {}deg (Bessel/line focus)",
                         profile.axicon_deg); break;
        case PhaseProfileKind::Freeform:
            std::println("  freeform hologram: {0}x{0} loaded phase map spanning {1}um",
                         profile.freeform_n, profile.freeform_extent_um); break;
    }

    // 1. Find the half-wave-plate meta-atom (best spin-flip conversion).
    std::println("  searching for the HWP atom ({0}x{0} fill grid, 2 solves each)...",
                 samples);
    HwpAtom atom = find_hwp_atom(pillar, materials::air(), materials::air(),
                                 materials::bk7(), period, lambda, thickness,
                                 0.10, 0.90, samples, M);
    std::println("  HWP atom: fill_x={:.3f} fill_y={:.3f}  retardance={:.1f} deg  "
                 "conversion eff={:.3f}",
                 atom.fill_x, atom.fill_y, atom.retardance_deg,
                 atom.conversion_efficiency);

    // 2. Design: rotate the fixed atom per site to stamp phi(x,y).
    PbMetalensDesign d = design_pb_metalens(atom, period, lambda, profile, diameter,
                                            handedness);
    std::println("  designed {0}x{0} rotated pillars  RMS phase error={1:.2e} deg "
                 "(geometric phase is exact)",
                 d.n_cells, d.rms_phase_error_deg);

    // 3. RCWA-verify the geometric-phase relation: solve the ROTATED atom at a
    //    few angles and confirm the spin-flip phase tracks -handedness*2*theta.
    std::vector<double> test_rot;
    for (int j = 0; j < 7; ++j) test_rot.push_back(j * (pi / 6.0));  // 0..180 deg
    auto vpts = verify_pb_phase(pillar, materials::air(), materials::air(),
                                materials::bk7(), period, lambda, atom, test_rot, M);
    // The atom's absolute conversion phase is a global piston -- recover it as the
    // circular mean of (measured + handedness*2*theta) so we test the SLOPE
    // (-handedness*2*theta), not the irrelevant constant.
    double sx = 0, sy = 0;
    for (const auto& v : vpts) {
        double resid = (v.cross_phase_deg + handedness * 2.0 * v.rotation_deg) * pi / 180.0;
        sx += std::cos(resid); sy += std::sin(resid);
    }
    const double piston_deg = std::atan2(sy, sx) * 180.0 / pi;
    std::println("  RCWA check (rotate the atom, measure spin-flip output; "
                 "global piston {:.1f} deg removed):", piston_deg);
    std::println("    theta(deg)   expected phase   measured phase   conv.eff   leakage");
    double phase_track_err = 0;
    for (const auto& v : vpts) {
        double expected = std::remainder(-handedness * 2.0 * v.rotation_deg + piston_deg, 360.0);
        double meas = std::remainder(v.cross_phase_deg, 360.0);
        double e = std::remainder(meas - expected, 360.0);
        phase_track_err += e * e;
        std::println("    {:7.1f}     {:11.1f}      {:11.1f}       {:.3f}     {:.3f}",
                     v.rotation_deg, expected, meas, v.conversion_eff, v.copol_leakage);
    }
    phase_track_err = std::sqrt(phase_track_err / vpts.size());
    std::println("    -> geometric-phase tracking RMS = {:.2f} deg (RCWA vs -2*theta)",
                 phase_track_err);

    // 4. Optical proof: propagate the spin-flipped field and confirm the profile
    //    does what it should -- a different test per element.
    const double pp = d.period_um;
    const double cen = (d.n_cells - 1) / 2.0;
    const double R_ap = diameter / 2.0;
    std::vector<double> px, py;
    std::vector<cdouble> tc;
    for (int iy = 0; iy < d.n_cells; ++iy)
        for (int ix = 0; ix < d.n_cells; ++ix) {
            double x = (ix - cen) * pp, y = (iy - cen) * pp;
            if (std::sqrt(x * x + y * y) > R_ap) continue;
            std::size_t off = (std::size_t)iy * d.n_cells + ix;
            px.push_back(x); py.push_back(y);
            tc.push_back(d.t_cross[off]);
        }
    const double recon_z = std::atof(arg_value(argc, argv, "--recon-z",
                                               std::to_string(focal).c_str()));
    ProfileProof proof = profile_optical_proof(px, py, tc, profile, lambda, focal,
                                               diameter, recon_z);
    const std::string& optical_line = proof.summary;
    const double psf_cx = proof.psf_cx, psf_cy = proof.psf_cy;
    const double psf_z = proof.psf_z, psf_hw = proof.psf_hw;
    std::println("  optical check: {}; conversion-efficiency cap = {:.1f}%",
                 optical_line, 100.0 * d.conversion_efficiency);

    // 5. GDS: rotated pillars (the rotation IS the design).
    int np = write_pb_gds(out, d.n_cells, d.period_um, atom.fill_x, atom.fill_y,
                          d.rotation_rad);
    if (np < 0) { std::println("  ERROR: could not write {}", out); return 1; }
    std::println("  wrote {} rotated pillars -> {}", np, out);

    // --report <prefix>: metrics txt + focal PSF image + GDS.
    const char* rp = arg_value(argc, argv, "--report", nullptr);
    if (rp) {
        std::string base = rp;
        auto psf = propagate_pillars(px, py, tc, psf_cx, psf_cy, psf_z, lambda, 201,
                                     psf_hw);
        bool ok = write_pgm(base + "_pb_psf.pgm", psf.n, psf.n, psf.intensity, 2.2);
        std::ofstream f(base + "_pb_report.txt");
        if (f) {
            f << "CELERIS Pancharatnam-Berry (geometric-phase) metasurface report\n";
            f << "===============================================================\n\n";
            f << std::format("profile            : {}\n", profile_name);
            f << std::format("wavelength         : {} um\n", lambda);
            f << std::format("aperture diameter  : {} um\n", diameter);
            f << std::format("illumination       : {} (circular)\n", handedness > 0 ? "RCP" : "LCP");
            f << std::format("optical check      : {}\n\n", optical_line);
            f << std::format("HWP atom fill_x    : {:.3f}\n", atom.fill_x);
            f << std::format("HWP atom fill_y    : {:.3f}\n", atom.fill_y);
            f << std::format("HWP retardance     : {:.1f} deg (ideal 180)\n",
                             atom.retardance_deg);
            f << std::format("conversion eff     : {:.3f} (efficiency cap)\n",
                             atom.conversion_efficiency);
            f << std::format("array              : {0} x {0} rotated pillars\n", d.n_cells);
            f << std::format("design phase error : {:.2e} deg (geometric phase is exact)\n",
                             d.rms_phase_error_deg);
            f << std::format("RCWA 2*theta track : {:.2f} deg RMS\n", phase_track_err);
        } else ok = false;
        std::println("  report bundle -> {0}_pb_report.txt (+ _pb_psf.pgm)  ok={1}",
                     base, ok);
    }
    return 0;
}
