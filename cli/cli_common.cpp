#include "cli.hpp"

const char* arg_value(int argc, char** argv, const std::string& key,
                      const char* fallback) {
    for (int i = 2; i + 1 < argc; ++i)
        if (key == argv[i]) return argv[i + 1];
    return fallback;
}

// Resolve the substrate material from --substrate <name>. Accepts any registry
// name (air/sio2/bk7/sapphire/gan/...); the old air|sio2|bk7 spellings still
// work since they are registry names. Falls back to the default name on error.
Material resolve_substrate(int argc, char** argv,
                                  const std::string& def) {
    const std::string name = arg_value(argc, argv, "--substrate", def.c_str());
    try {
        return materials::by_name(name);
    } catch (const std::exception& e) {
        std::println("  WARNING: substrate {}; using {}", e.what(), def);
        return materials::by_name(def);
    }
}

// Resolve the pillar material. Precedence: --pillar-csv <file> (explicit data)
// > --pillar <name> (registry name, e.g. tio2/gan/a-si — real n,k incl. loss)
// > --pillar-n <n> (a lossless constant index, the legacy default). Reports the
// resolved choice. `def_n` is the constant-index fallback.
Material resolve_pillar(int argc, char** argv, double def_n) {
    if (const char* csv = arg_value(argc, argv, "--pillar-csv", nullptr))
        return load_material_csv(csv, "pillar-csv");
    if (const char* name = arg_value(argc, argv, "--pillar", nullptr)) {
        try {
            Material m = materials::by_name(name);
            std::println("  pillar material: {} (registry)", m.name());
            return m;
        } catch (const std::exception& e) {
            std::println("  WARNING: pillar {}; using constant n={}", e.what(), def_n);
        }
    }
    const double pillar_n = std::atof(arg_value(argc, argv, "--pillar-n",
                                                std::to_string(def_n).c_str()));
    return Material::constant(cdouble{pillar_n, 0.0}, "pillar");
}

// Parse --profile (+ its parameters) into a PhaseProfile, shared by the
// propagation-phase (`design`) and geometric-phase (`pbdesign`) paths. Returns
// nullopt and prints an error on a bad spec. `focal`/`diameter` supply the
// defaults for the focusing term and the freeform map extent.
std::optional<PhaseProfile> parse_phase_profile(int argc, char** argv,
                                                       double focal, double diameter) {
    const std::string name = arg_value(argc, argv, "--profile", "focusing");
    PhaseProfile p;
    if (name == "focusing") {
        p.kind = PhaseProfileKind::Focusing;
        p.focal_length_um = focal;
    } else if (name == "quadratic") {
        // Parabolic wide-FOV lens (focuses on-axis at z=f like focusing, but its
        // off-axis focus stays sharp -- see `celeris widefov`).
        p.kind = PhaseProfileKind::Quadratic;
        p.focal_length_um = focal;
    } else if (name == "vortex") {
        p.kind = PhaseProfileKind::Vortex;
        p.focal_length_um = focal;  // focused vortex (donut focal spot); <=0 => collimated
        p.topological_charge = std::atoi(arg_value(argc, argv, "--charge", "1"));
    } else if (name == "deflector") {
        p.kind = PhaseProfileKind::Deflector;
        p.deflect_deg = std::atof(arg_value(argc, argv, "--deflect-deg", "10"));
        p.deflect_azimuth_deg = std::atof(arg_value(argc, argv, "--deflect-azimuth", "0"));
    } else if (name == "axicon") {
        p.kind = PhaseProfileKind::Axicon;
        p.axicon_deg = std::atof(arg_value(argc, argv, "--axicon-deg", "5"));
    } else if (name == "freeform") {
        // Arbitrary loaded phi(x,y) map (a computer-generated hologram). The map
        // spans --freeform-extent um (default: the full aperture diameter).
        const char* ff = arg_value(argc, argv, "--freeform-file", nullptr);
        if (!ff) {
            std::println("ERROR: --profile freeform needs --freeform-file <grid.txt>");
            return std::nullopt;
        }
        const double extent = std::atof(arg_value(argc, argv, "--freeform-extent",
                                                  std::to_string(diameter).c_str()));
        try {
            p = load_freeform_phase(ff, extent);
        } catch (const std::exception& e) {
            std::println("ERROR: {}", e.what());
            return std::nullopt;
        }
    } else {
        std::println("ERROR: unknown --profile '{}' "
                     "(use focusing|quadratic|vortex|deflector|axicon|freeform)", name);
        return std::nullopt;
    }
    return p;
}

// Propagate a realized aperture field {(px[q], py[q]) -> tc[q]} and run the
// profile's optical proof -- a different physical test per element (focus z-scan,
// donut null, beam deflection, line focus, freeform reconstruction). Shared by
// the geometric-phase (`pbdesign`) and propagation-phase (`design`) paths so both
// verify the SAME way. `recon_z` is where a freeform hologram is reconstructed.
ProfileProof profile_optical_proof(const std::vector<double>& px,
                                          const std::vector<double>& py,
                                          const std::vector<cdouble>& tc,
                                          const PhaseProfile& profile, double lambda,
                                          double focal, double diameter,
                                          double recon_z) {
    const double k = 2.0 * pi / lambda;
    const double dl = lambda * std::max(focal, 1.0) / diameter;  // ~spot scale
    const double R_ap = diameter / 2.0;

    // On-axis Rayleigh-Sommerfeld intensity at distance z (focusing/axicon).
    auto on_axis_I = [&](double z) {
        cdouble E{0, 0};
        for (std::size_t q = 0; q < px.size(); ++q) {
            double r = std::sqrt(px[q] * px[q] + py[q] * py[q] + z * z);
            E += tc[q] * std::polar(1.0 / r, k * r);
        }
        return std::norm(E);
    };

    ProfileProof out;
    out.psf_z = focal;
    out.psf_hw = std::max(5.0 * dl, 4.0);

    if (profile.kind == PhaseProfileKind::Focusing ||
        profile.kind == PhaseProfileKind::Quadratic) {
        const int NZ = 240;
        double best_z = 0.5 * focal, best_I = -1.0;
        for (int j = 0; j < NZ; ++j) {
            double z = 0.5 * focal + focal * j / (NZ - 1);
            double I = on_axis_I(z);
            if (I > best_I) { best_I = I; best_z = z; }
        }
        out.summary = std::format("on-axis focus at z = {:.1f} um (target {:.1f})",
                                  best_z, focal);
    } else if (profile.kind == PhaseProfileKind::Vortex) {
        // Focused vortex -> a DONUT focal spot: deep on-axis null, bright ring.
        out.psf_hw = std::max(8.0 * dl, 5.0);
        auto psf = propagate_pillars(px, py, tc, 0, 0, out.psf_z, lambda, 201, out.psf_hw);
        int nn = psf.n, pk = 0; double peak = 0;
        for (int i = 0; i < nn * nn; ++i)
            if (psf.intensity[i] > peak) { peak = psf.intensity[i]; pk = i; }
        double center_I = psf.intensity[(nn / 2) * nn + nn / 2];
        double pix = 2.0 * out.psf_hw / (nn - 1);
        double ring_r = std::hypot((double)(pk % nn - nn / 2), (double)(pk / nn - nn / 2)) * pix;
        out.summary = std::format(
            "focal-plane donut (l={}): on-axis/peak = {:.3f} (null), ring radius {:.2f} um",
            profile.topological_charge, peak > 0 ? center_I / peak : 0.0, ring_r);
    } else if (profile.kind == PhaseProfileKind::Deflector) {
        // The beam lands displaced by z*tan(a) along the azimuth on a screen z=focal.
        const double a = profile.deflect_deg * pi / 180.0;
        const double az = profile.deflect_azimuth_deg * pi / 180.0;
        out.psf_cx = out.psf_z * std::tan(a) * std::cos(az);
        out.psf_cy = out.psf_z * std::tan(a) * std::sin(az);
        out.psf_hw = std::max(8.0 * dl, 5.0);
        auto psf = propagate_pillars(px, py, tc, out.psf_cx, out.psf_cy, out.psf_z,
                                     lambda, 201, out.psf_hw);
        int nn = psf.n, pk = 0; double peak = 0;
        for (int i = 0; i < nn * nn; ++i)
            if (psf.intensity[i] > peak) { peak = psf.intensity[i]; pk = i; }
        double pix = 2.0 * out.psf_hw / (nn - 1);
        double mx = out.psf_cx + (pk % nn - nn / 2) * pix;
        double my = out.psf_cy + (pk / nn - nn / 2) * pix;
        double disp = std::hypot(mx, my);
        double meas_deg = std::atan2(disp, out.psf_z) * 180.0 / pi;
        out.summary = std::format(
            "beam at screen z={:.0f}um lands {:.2f}um off-axis -> {:.2f}deg (target {:.2f}deg)",
            out.psf_z, disp, meas_deg, profile.deflect_deg);
    } else if (profile.kind == PhaseProfileKind::Axicon) {
        // Bessel beam with an extended on-axis line focus ~ R/tan(b).
        const double b = profile.axicon_deg * pi / 180.0;
        const double dof_pred = R_ap / std::tan(b);
        const double zlo = 0.1 * dof_pred, zhi = 1.4 * dof_pred;
        const int NZ = 400;
        std::vector<double> Iz(NZ);
        double maxI = 0;
        for (int j = 0; j < NZ; ++j) {
            Iz[j] = on_axis_I(zlo + (zhi - zlo) * j / (NZ - 1));
            if (Iz[j] > maxI) maxI = Iz[j];
        }
        int jlo = -1, jhi = -1;
        for (int j = 0; j < NZ; ++j)
            if (Iz[j] >= 0.5 * maxI) { if (jlo < 0) jlo = j; jhi = j; }
        double zspan = jhi > jlo ? (jhi - jlo) * (zhi - zlo) / (NZ - 1) : 0.0;
        out.psf_z = 0.5 * dof_pred;
        out.psf_hw = std::max(6.0 * dl, 4.0);
        // The FWHM span is a fraction of the geometric max range R/tan(b) -- the
        // on-axis intensity ramps up then falls, it isn't flat over the full range.
        out.summary = std::format(
            "extended on-axis line focus: FWHM {:.1f}um (geometric max range R/tan b ~{:.1f}um)",
            zspan, dof_pred);
    } else {  // Freeform: an arbitrary loaded phi(x,y); no single crisp metric.
        // The realized field reconstructs whatever wavefront the loaded map encodes;
        // render it at z=recon_z. Per-site phase fidelity is reported by the caller.
        out.psf_z = recon_z;
        out.psf_hw = std::max(8.0 * dl, 5.0);
        out.summary = std::format(
            "freeform {0}x{0} hologram: reconstruction rendered at z={1:.0f}um "
            "(per-site phase fidelity reported above)",
            profile.freeform_n, recon_z);
    }
    return out;
}
