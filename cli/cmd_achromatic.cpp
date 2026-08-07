#include "cli.hpp"

// celeris achromatic: design a broadband (achromatic) focusing metalens by
// DISPERSION ENGINEERING. Builds a meta-atom library characterized across the
// band (phase + group delay per atom), then picks at each site the atom matching
// BOTH the base focusing phase AND the radius-dependent group delay. Compares the
// chromatic focal shift of the achromatic design against a standard single-
// wavelength design to show the flattening -- and honestly reports the group-
// delay budget (a single-etch square-pillar library supplies a limited GD range,
// so the achromatic band/aperture is fundamentally bounded). ROADMAP #2.
int cmd_achromatic(int argc, char** argv) {
    const double focal = std::atof(arg_value(argc, argv, "--focal", "30"));
    const double diameter = std::atof(arg_value(argc, argv, "--diameter", "10"));
    const double lambda0 = std::atof(arg_value(argc, argv, "--wavelength", "0.532"));
    const double frac_bw = std::atof(arg_value(argc, argv, "--bandwidth", "0.20"));
    const int nb = std::atoi(arg_value(argc, argv, "--band-samples", "7"));
    const double period = std::atof(arg_value(argc, argv, "--period", "0.35"));
    const double pillar_n = std::atof(arg_value(argc, argv, "--pillar-n", "2.4"));
    const int samples = std::atoi(arg_value(argc, argv, "--fill-samples", "24"));
    const int M = std::atoi(arg_value(argc, argv, "--harmonics", "6"));
    const double gd_weight = std::atof(arg_value(argc, argv, "--gd-weight", "1.0"));
    const std::string out = arg_value(argc, argv, "--out", "achromatic.gds");
    const std::string sub_name = arg_value(argc, argv, "--substrate", "sio2");

    // The 2-DOF library: a single geometric DOF (fill only) traces a 1-D curve in
    // the (phase, group-delay) plane and cannot set both independently, so sweep a
    // fill x height grid. n_heights==1 falls back to a single-etch (1-DOF) library.
    const double h_lo = std::atof(arg_value(argc, argv, "--height-lo", "0.40"));
    const double h_hi = std::atof(arg_value(argc, argv, "--height-hi", "1.40"));
    const int n_h = std::atoi(arg_value(argc, argv, "--height-steps", "10"));

    // --single-etch: span the (phase, group-delay) plane with SHAPE variety at ONE
    // height instead of a fill x height grid -- fabricable in a single etch (no
    // grayscale lithography). The honest tradeoff is a smaller group-delay span.
    bool single_etch = false;
    for (int i = 2; i < argc; ++i)
        if (std::string(argv[i]) == "--single-etch") single_etch = true;
    // Single-etch wants a TALLER pillar: group delay accumulates with optical path,
    // so a higher-aspect atom widens the GD span each shape can supply (h~1.1µm at
    // Λ=0.35 ~ AR 4, fabricable). 3.2× drift flattening here vs 1.3× at 0.6µm.
    const double etch_height = std::atof(arg_value(argc, argv, "--height", "1.10"));

    const Material& substrate = sub_name == "air"  ? materials::air()
                                : sub_name == "bk7" ? materials::bk7()
                                                    : materials::fused_silica();
    const char* pillar_csv = arg_value(argc, argv, "--pillar-csv", nullptr);
    const Material pillar =
        pillar_csv ? load_material_csv(pillar_csv, "pillar-csv")
                   : Material::constant(cdouble{pillar_n, 0.0}, "pillar");

    // Band samples, ascending, centered on lambda0 (fractional bandwidth frac_bw).
    const double lam_lo = lambda0 * (1.0 - 0.5 * frac_bw);
    const double lam_hi = lambda0 * (1.0 + 0.5 * frac_bw);
    std::vector<double> band(nb);
    for (int j = 0; j < nb; ++j)
        band[j] = lam_lo + (lam_hi - lam_lo) * j / (nb - 1);

    std::println("CELERIS achromatic (broadband) metalens design");
    std::println("  f={}µm  D={}µm  λ0={}µm  band=[{:.3f},{:.3f}]µm ({:.0f}% BW, {} samples)",
                 focal, diameter, lambda0, lam_lo, lam_hi, frac_bw * 100.0, nb);

    DispersiveLibrary dlib;
    if (single_etch) {
        std::println("  Λ={}µm  n_pillar={}  substrate={}  SINGLE-ETCH (one height "
                     "{:.2f}µm, shape-diverse)", period, pillar_n, sub_name, etch_height);
        std::println("  building single-etch dispersive library ({} fills × "
                     "7 shapes (square+circle+3 cross+2 ring) × {} wavelengths, M={})...",
                     samples, nb, M);
        dlib = build_single_etch_library(pillar, materials::air(), materials::air(),
                                         substrate, period, band, lambda0, 0.08, 0.92,
                                         samples, etch_height, M);
    } else {
        std::println("  Λ={}µm  n_pillar={}  substrate={}  height grid [{:.2f},{:.2f}]µm × {}",
                     period, pillar_n, sub_name, h_lo, h_hi, n_h);
        std::println("  building dispersive library ({} fills × {} heights × {} wavelengths, M={})...",
                     samples, n_h, nb, M);
        dlib = build_dispersive_library(pillar, materials::air(), materials::air(),
                                        substrate, period, band, lambda0, 0.08, 0.92,
                                        samples, h_lo, h_hi, n_h, M);
    }
    std::println("  group-delay range supplied by the library: [{:.2f}, {:.2f}] fs "
                 "(span {:.2f} fs, {} atoms)",
                 dlib.gd_min_fs, dlib.gd_max_fs, dlib.gd_max_fs - dlib.gd_min_fs,
                 static_cast<int>(dlib.atoms.size()));

    // Two designs from the SAME library isolate the effect of dispersion
    // engineering: gd_weight=0 is dispersion-blind (the standard baseline);
    // gd_weight>0 adds the group-delay objective (achromatic).
    auto std_des = design_achromatic_metalens(dlib, focal, diameter, /*gd_weight=*/0.0);
    auto ach_des = design_achromatic_metalens(dlib, focal, diameter, gd_weight);

    std::println("  standard   (gd_weight 0): base-phase RMS {:.1f}°, GD RMS {:.2f} fs, "
                 "mean |t| {:.3f}, heights [{:.2f},{:.2f}]µm",
                 std_des.rms_phase_error_deg, std_des.rms_group_delay_error_fs,
                 std_des.mean_amplitude, std_des.min_height_um, std_des.max_height_um);
    std::println("  achromatic (gd_weight {:g}): base-phase RMS {:.1f}°, GD RMS {:.2f} fs, "
                 "mean |t| {:.3f}, heights [{:.2f},{:.2f}]µm",
                 gd_weight, ach_des.rms_phase_error_deg, ach_des.rms_group_delay_error_fs,
                 ach_des.mean_amplitude, ach_des.min_height_um, ach_des.max_height_um);
    std::println("  group-delay budget: required span {:.2f} fs, available {:.2f} fs "
                 "-> coverage {:.2f}",
                 ach_des.required_gd_span_fs, ach_des.available_gd_span_fs, ach_des.gd_coverage);
    if (ach_des.gd_coverage < 1.0)
        std::println("  HONEST: GD coverage < 1 -> the library cannot supply the full "
                     "group-delay span; achromatic only over a reduced aperture/bandwidth "
                     "(taller/higher-aspect or coupled atoms widen the span).");
    if (!ach_des.single_height)
        std::println("  NOTE: the achromatic design uses MULTIPLE pillar heights -> a multi-"
                     "level / grayscale etch (a single GDS layer encodes only the in-plane "
                     "footprints, not the per-site depth).");
    else if (single_etch) {
        // Confirm the single-etch property and show which shapes the design leaned
        // on -- shape diversity (not depth) is what supplied the group delay.
        int n_sq = 0, n_ci = 0, n_cr = 0, n_ri = 0;
        for (int q : ach_des.atom_index) {
            switch (dlib.atoms[q].shape) {
                case MetaShape::Rectangle: ++n_sq; break;
                case MetaShape::Ellipse:   ++n_ci; break;
                case MetaShape::Cross:     ++n_cr; break;
                case MetaShape::Ring:      ++n_ri; break;
            }
        }
        std::println("  SINGLE-ETCH: every pillar shares one height ({:.2f}µm) -> one "
                     "lithography step, no grayscale.", ach_des.min_height_um);
        std::println("  shape mix chosen (achromatic): {} square, {} circle, {} cross, "
                     "{} ring", n_sq, n_ci, n_cr, n_ri);
    }

    // Rigorous chromatic response of BOTH designs from the library's stored per-
    // atom band data (each atom's true dispersion; no extra RCWA solves).
    auto chrom_std = verify_achromatic_focus(dlib, std_des);
    auto chrom_ach = verify_achromatic_focus(dlib, ach_des);
    std::println("  chromatic focal length across the band (rigorous; flat = achromatic):");
    std::println("      {:>9}  {:>16}  {:>16}", "λ(nm)", "standard f(µm)", "achromatic f(µm)");
    double std_lo = 1e300, std_hi = -1e300, ach_lo = 1e300, ach_hi = -1e300;
    for (int j = 0; j < nb; ++j) {
        std::println("      {:>9.0f}  {:>16.2f}  {:>16.2f}",
                     chrom_std[j].wavelength_um * 1000.0,
                     chrom_std[j].focal_length_um, chrom_ach[j].focal_length_um);
        std_lo = std::min(std_lo, chrom_std[j].focal_length_um);
        std_hi = std::max(std_hi, chrom_std[j].focal_length_um);
        ach_lo = std::min(ach_lo, chrom_ach[j].focal_length_um);
        ach_hi = std::max(ach_hi, chrom_ach[j].focal_length_um);
    }
    const double std_drift = std_hi - std_lo, ach_drift = ach_hi - ach_lo;
    std::println("  focal drift over the band: standard {:.2f} µm  ->  achromatic {:.2f} µm "
                 "({:.1f}× tighter)",
                 std_drift, ach_drift,
                 ach_drift > 1e-9 ? std_drift / ach_drift : 0.0);

    MetalensDesign ach_lens = to_metalens_design(ach_des);
    int np = write_metalens_gds(ach_lens, out);
    if (np < 0) { std::println("  ERROR: could not write {}", out); return 1; }
    std::println("  wrote {} pillars -> {}", np, out);

    const char* rprefix = arg_value(argc, argv, "--report", nullptr);
    if (rprefix) {
        std::string base = rprefix;
        std::ofstream f(base + "_achromatic_report.txt");
        bool ok = static_cast<bool>(f);
        if (f) {
            f << "CELERIS achromatic metalens design report\n";
            f << "=========================================\n\n";
            f << std::format("center wavelength : {} um\n", lambda0);
            f << std::format("band              : [{:.3f}, {:.3f}] um ({:.0f}% BW)\n",
                             lam_lo, lam_hi, frac_bw * 100.0);
            f << std::format("focal length      : {} um\n", focal);
            f << std::format("aperture diameter : {} um\n", diameter);
            f << std::format("period            : {} um\n", period);
            if (single_etch)
                f << std::format("library           : single-etch, shape-diverse @ {:.2f} um\n",
                                 etch_height);
            else
                f << std::format("height grid       : [{:.2f}, {:.2f}] um x {}\n", h_lo, h_hi, n_h);
            f << std::format("library atoms     : {}\n", static_cast<int>(dlib.atoms.size()));
            f << std::format("array             : {0} x {0} pillars\n\n", ach_des.n_cells);
            f << std::format("base-phase RMS    : {:.1f} deg (standard {:.1f})\n",
                             ach_des.rms_phase_error_deg, std_des.rms_phase_error_deg);
            f << std::format("group-delay RMS   : {:.2f} fs (standard {:.2f})\n",
                             ach_des.rms_group_delay_error_fs, std_des.rms_group_delay_error_fs);
            f << std::format("mean transmission : {:.3f}\n", ach_des.mean_amplitude);
            f << std::format("GD required span  : {:.2f} fs\n", ach_des.required_gd_span_fs);
            f << std::format("GD available span : {:.2f} fs\n", ach_des.available_gd_span_fs);
            f << std::format("GD coverage       : {:.2f}\n", ach_des.gd_coverage);
            f << std::format("pillar heights    : [{:.2f}, {:.2f}] um ({})\n\n",
                             ach_des.min_height_um, ach_des.max_height_um,
                             ach_des.single_height ? "single etch" : "multi-level etch");
            f << std::format("focal drift (std) : {:.2f} um\n", std_drift);
            f << std::format("focal drift (ach) : {:.2f} um\n", ach_drift);
            f << "\nlambda(nm)   standard_f(um)   achromatic_f(um)\n";
            for (int j = 0; j < nb; ++j)
                f << std::format("{:>9.0f}   {:>14.2f}   {:>16.2f}\n",
                                 chrom_std[j].wavelength_um * 1000.0,
                                 chrom_std[j].focal_length_um, chrom_ach[j].focal_length_um);
        }
        ok &= (write_metalens_gds(ach_lens, base + "_achromatic_layout.gds") >= 0);
        std::println("  report bundle -> {0}_achromatic_report.txt (+ _achromatic_layout.gds)  ok={1}",
                     base, ok);
    }
    return 0;
}

// celeris pbachromatic: the MODERN achromatic recipe -- geometric (PB) phase +
// dispersion engineering, in a SINGLE etch. A dispersive birefringent library
// (rectangles over a fill_x x fill_y grid at one height) supplies a range of group
// delays; at each site the atom is chosen PURELY for its group delay and then
// ROTATED so the base focusing phase is hit EXACTLY (geometric phase is exact and
// wavelength-independent). So the base-phase residual is ~0 -- unlike the
// propagation-phase achromat -- and the achromatic limit is only the library's
// group-delay coverage. Compares against a standard (chromatic) PB lens from the
// same library and writes a rotated-rectangle GDS (one etch depth). ROADMAP #2.
int cmd_pb_achromatic(int argc, char** argv) {
    const double focal = std::atof(arg_value(argc, argv, "--focal", "30"));
    const double diameter = std::atof(arg_value(argc, argv, "--diameter", "10"));
    const double lambda0 = std::atof(arg_value(argc, argv, "--wavelength", "0.532"));
    const double frac_bw = std::atof(arg_value(argc, argv, "--bandwidth", "0.20"));
    const int nb = std::atoi(arg_value(argc, argv, "--band-samples", "7"));
    const double period = std::atof(arg_value(argc, argv, "--period", "0.35"));
    const double pillar_n = std::atof(arg_value(argc, argv, "--pillar-n", "2.4"));
    const int samples = std::atoi(arg_value(argc, argv, "--fill-samples", "12"));
    const int M = std::atoi(arg_value(argc, argv, "--harmonics", "6"));
    const double gd_weight = std::atof(arg_value(argc, argv, "--gd-weight", "1.0"));
    // Birefringent atoms accumulate group delay with optical path, so a taller
    // pillar widens the GD span each footprint can supply (h~1.1µm at Λ=0.35 ~ AR 4).
    const double height = std::atof(arg_value(argc, argv, "--height", "1.10"));
    const int handedness = std::atoi(arg_value(argc, argv, "--handedness", "1")) >= 0 ? 1 : -1;
    const std::string out = arg_value(argc, argv, "--out", "pb_achromatic.gds");
    const std::string sub_name = arg_value(argc, argv, "--substrate", "bk7");

    const Material& substrate = sub_name == "air" ? materials::air()
                                : sub_name == "sio2" ? materials::fused_silica()
                                                     : materials::bk7();
    const char* pillar_csv = arg_value(argc, argv, "--pillar-csv", nullptr);
    const Material pillar =
        pillar_csv ? load_material_csv(pillar_csv, "pillar-csv")
                   : Material::constant(cdouble{pillar_n, 0.0}, "pillar");

    const double lam_lo = lambda0 * (1.0 - 0.5 * frac_bw);
    const double lam_hi = lambda0 * (1.0 + 0.5 * frac_bw);
    std::vector<double> band(nb);
    for (int j = 0; j < nb; ++j)
        band[j] = lam_lo + (lam_hi - lam_lo) * j / (nb - 1);

    std::println("CELERIS achromatic Pancharatnam-Berry (geometric phase + dispersion)");
    std::println("  f={}µm  D={}µm  λ0={}µm  band=[{:.3f},{:.3f}]µm ({:.0f}% BW, {} samples)",
                 focal, diameter, lambda0, lam_lo, lam_hi, frac_bw * 100.0, nb);
    std::println("  Λ={}µm  n_pillar={}  substrate={}  illumination={}  SINGLE-ETCH height {:.2f}µm",
                 period, pillar_n, sub_name, handedness > 0 ? "RCP" : "LCP", height);
    std::println("  building dispersive birefringent library ({0}×{0} fill grid × {1} "
                 "wavelengths, 2 solves each, M={2})...", samples, nb, M);

    DispersivePbLibrary lib = build_dispersive_pb_library(
        pillar, materials::air(), materials::air(), substrate, period, band, lambda0,
        0.10, 0.90, samples, height, M);
    std::println("  group-delay range supplied by the library: [{:.2f}, {:.2f}] fs "
                 "(span {:.2f} fs, {} atoms)",
                 lib.gd_min_fs, lib.gd_max_fs, lib.gd_max_fs - lib.gd_min_fs,
                 static_cast<int>(lib.atoms.size()));

    // Two designs from the SAME library isolate the dispersion engineering:
    // gd_weight=0 = the best-conversion atom everywhere = a STANDARD (chromatic) PB
    // lens; gd_weight>0 varies the atom per radius to match the group delay.
    auto std_des = design_pb_achromatic_metalens(lib, focal, diameter, handedness, 0.0);
    auto ach_des = design_pb_achromatic_metalens(lib, focal, diameter, handedness, gd_weight);

    std::println("  standard PB (gd_weight 0): base-phase RMS {:.2e}°, GD RMS {:.2f} fs, "
                 "mean |a_cross| {:.3f}, conv {:.3f}",
                 std_des.rms_phase_error_deg, std_des.rms_group_delay_error_fs,
                 std_des.mean_amplitude, std_des.mean_conversion);
    std::println("  achromatic (gd_weight {:g}): base-phase RMS {:.2e}°, GD RMS {:.2f} fs, "
                 "mean |a_cross| {:.3f}, conv {:.3f}",
                 gd_weight, ach_des.rms_phase_error_deg, ach_des.rms_group_delay_error_fs,
                 ach_des.mean_amplitude, ach_des.mean_conversion);
    std::println("  base-phase RMS ~0 for BOTH: the rotation stamps the base phase EXACTLY "
                 "(geometric phase) -- the atom is free to chase group delay.");
    std::println("  group-delay budget: required span {:.2f} fs, available {:.2f} fs "
                 "-> coverage {:.2f}",
                 ach_des.required_gd_span_fs, ach_des.available_gd_span_fs, ach_des.gd_coverage);
    if (ach_des.gd_coverage < 1.0)
        std::println("  HONEST: GD coverage < 1 -> the library cannot supply the full "
                     "group-delay span; achromatic only over a reduced aperture/bandwidth "
                     "(taller/higher-aspect or coupled atoms widen the span).");

    // Rigorous chromatic response of BOTH designs from the library's stored band data.
    auto chrom_std = verify_pb_achromatic_focus(lib, std_des);
    auto chrom_ach = verify_pb_achromatic_focus(lib, ach_des);
    std::println("  chromatic focal length across the band (rigorous; flat = achromatic):");
    std::println("      {:>9}  {:>16}  {:>16}", "λ(nm)", "standard f(µm)", "achromatic f(µm)");
    double std_lo = 1e300, std_hi = -1e300, ach_lo = 1e300, ach_hi = -1e300;
    for (int j = 0; j < nb; ++j) {
        std::println("      {:>9.0f}  {:>16.2f}  {:>16.2f}",
                     chrom_std[j].wavelength_um * 1000.0,
                     chrom_std[j].focal_length_um, chrom_ach[j].focal_length_um);
        std_lo = std::min(std_lo, chrom_std[j].focal_length_um);
        std_hi = std::max(std_hi, chrom_std[j].focal_length_um);
        ach_lo = std::min(ach_lo, chrom_ach[j].focal_length_um);
        ach_hi = std::max(ach_hi, chrom_ach[j].focal_length_um);
    }
    const double std_drift = std_hi - std_lo, ach_drift = ach_hi - ach_lo;
    std::println("  focal drift over the band: standard {:.2f} µm  ->  achromatic {:.2f} µm "
                 "({:.1f}× tighter)",
                 std_drift, ach_drift, ach_drift > 1e-9 ? std_drift / ach_drift : 0.0);

    // GDS: per-site rotated rectangles, all one etch depth (a single mask layer).
    int np = write_pb_rect_gds(out, ach_des.n_cells, ach_des.period_um,
                               ach_des.fill_x_map, ach_des.fill_y_map,
                               ach_des.rotation_rad);
    if (np < 0) { std::println("  ERROR: could not write {}", out); return 1; }
    std::println("  wrote {} rotated rectangles (one etch depth) -> {}", np, out);

    const char* rprefix = arg_value(argc, argv, "--report", nullptr);
    if (rprefix) {
        std::string base = rprefix;
        std::ofstream f(base + "_pb_achromatic_report.txt");
        bool ok = static_cast<bool>(f);
        if (f) {
            f << "CELERIS achromatic Pancharatnam-Berry metalens design report\n";
            f << "===========================================================\n\n";
            f << std::format("recipe            : geometric (PB) phase + dispersion (single etch)\n");
            f << std::format("center wavelength : {} um\n", lambda0);
            f << std::format("band              : [{:.3f}, {:.3f}] um ({:.0f}% BW)\n",
                             lam_lo, lam_hi, frac_bw * 100.0);
            f << std::format("focal length      : {} um\n", focal);
            f << std::format("aperture diameter : {} um\n", diameter);
            f << std::format("period            : {} um\n", period);
            f << std::format("etch depth        : {:.2f} um (single etch)\n", height);
            f << std::format("illumination      : {} (circular)\n", handedness > 0 ? "RCP" : "LCP");
            f << std::format("library atoms     : {}\n", static_cast<int>(lib.atoms.size()));
            f << std::format("array             : {0} x {0} rotated rectangles\n\n", ach_des.n_cells);
            f << std::format("base-phase RMS    : {:.2e} deg (geometric phase is exact)\n",
                             ach_des.rms_phase_error_deg);
            f << std::format("group-delay RMS   : {:.2f} fs (standard {:.2f})\n",
                             ach_des.rms_group_delay_error_fs, std_des.rms_group_delay_error_fs);
            f << std::format("mean |a_cross|    : {:.3f}\n", ach_des.mean_amplitude);
            f << std::format("conversion (cap)  : {:.3f}\n", ach_des.mean_conversion);
            f << std::format("GD required span  : {:.2f} fs\n", ach_des.required_gd_span_fs);
            f << std::format("GD available span : {:.2f} fs\n", ach_des.available_gd_span_fs);
            f << std::format("GD coverage       : {:.2f}\n", ach_des.gd_coverage);
            f << std::format("focal drift (std) : {:.2f} um\n", std_drift);
            f << std::format("focal drift (ach) : {:.2f} um\n", ach_drift);
            f << "\nlambda(nm)   standard_f(um)   achromatic_f(um)\n";
            for (int j = 0; j < nb; ++j)
                f << std::format("{:>9.0f}   {:>14.2f}   {:>16.2f}\n",
                                 chrom_std[j].wavelength_um * 1000.0,
                                 chrom_std[j].focal_length_um, chrom_ach[j].focal_length_um);
        }
        ok &= (write_pb_rect_gds(base + "_pb_achromatic_layout.gds", ach_des.n_cells,
                                 ach_des.period_um, ach_des.fill_x_map, ach_des.fill_y_map,
                                 ach_des.rotation_rad) >= 0);
        std::println("  report bundle -> {0}_pb_achromatic_report.txt (+ _pb_achromatic_layout.gds)  ok={1}",
                     base, ok);
    }
    return 0;
}
