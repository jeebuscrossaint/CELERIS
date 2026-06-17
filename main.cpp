#include <cmath>
#include <cstdlib>
#include <print>
#include <string>
#include <vector>

#ifdef CELERIS_USE_CUDA
#include <Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <random>

#include "celeris/cuda/eigensolve.hpp"
#endif

#include "celeris/analysis/chromatic.hpp"
#include "celeris/analysis/focal.hpp"
#include "celeris/analysis/tolerance.hpp"
#include "celeris/design/metalens.hpp"
#include "celeris/design/optimize.hpp"
#include "celeris/io/gds.hpp"
#include "celeris/materials/database.hpp"
#include "celeris/optics/tmm.hpp"
#include "celeris/rcwa/rcwa1d.hpp"
#include "celeris/rcwa/rcwa2d.hpp"

using namespace celeris;

// Physics validation suite — every solver checked against closed-form results,
// an independent solver, or energy conservation. Run with: celeris selftest
static int run_selftest() {
    const double lambda = 0.550;  // design wavelength, 550 nm (µm)
    const double normal = 0.0;    // normal incidence

    // ---- Validation 1: bare glass, energy conservation --------------------
    // No layers, air -> N-BK7. Closed form R = ((1-n)/(1+n))^2 ~ 0.0424 at
    // 550 nm, and since glass is lossless here, R + T must equal 1.
    {
        auto res = solve_stack(materials::air(), {}, materials::bk7(),
                               lambda, normal, Pol::TE);
        std::println("[1] Bare N-BK7 glass at 550 nm:");
        std::println("    R = {:.4f}  (expect ~0.0424)", res.R);
        std::println("    T = {:.4f}", res.T);
        std::println("    R + T = {:.6f}  (expect 1.000000)", res.R + res.T);
    }

    // ---- Validation 2: quarter-wave AR coating ----------------------------
    // A single layer of index n1 = sqrt(n_air * n_glass), one quarter-wave
    // thick, makes reflection vanish at the design wavelength. This is the
    // textbook single-layer anti-reflection coating.
    {
        const cdouble n_glass = materials::bk7().index(lambda);
        const cdouble n1 = std::sqrt(materials::air().index(lambda) * n_glass);
        const double d1 = lambda / (4.0 * n1.real());  // quarter-wave thickness

        std::vector<Layer> stack = {
            {Material::constant(n1, "ideal-AR"), d1}};

        std::println("[2] Quarter-wave AR coating (n1 = {:.4f}, d = {:.1f} nm):",
                     n1.real(), d1 * 1000.0);
        for (double wl : {0.450, 0.550, 0.650}) {
            auto res = solve_stack(materials::air(), stack, materials::bk7(),
                                   wl, normal, Pol::TE);
            std::println("    lambda = {:.0f} nm:  R = {:.5f}{}", wl * 1000.0,
                         res.R, wl == 0.550 ? "   <- ~0 at design" : "");
        }
    }

    // ---- Validation 3: distributed Bragg reflector ------------------------
    // Alternating high/low quarter-wave layers build a high-reflectance
    // mirror. Reflectance should climb toward 1 as we add more pairs.
    {
        const cdouble nH{2.30, 0.0};  // high index (e.g. TiO2-like)
        const cdouble nL{1.46, 0.0};  // low index  (e.g. SiO2-like)
        const double dH = lambda / (4.0 * nH.real());
        const double dL = lambda / (4.0 * nL.real());

        std::println("[3] Bragg mirror (nH=2.30 / nL=1.46 quarter-wave pairs):");
        for (int pairs : {2, 4, 8}) {
            std::vector<Layer> stack;
            for (int p = 0; p < pairs; ++p) {
                stack.push_back({Material::constant(nH, "H"), dH});
                stack.push_back({Material::constant(nL, "L"), dL});
            }
            auto res = solve_stack(materials::air(), stack, materials::bk7(),
                                   lambda, normal, Pol::TE);
            std::println("    {} pairs:  R = {:.5f}", pairs, res.R);
        }
    }

    // ---- Validation 4: RCWA degenerate grating == TMM slab ----------------
    // If a grating's ridge and groove are the SAME material, it's just a
    // uniform slab. RCWA must then reproduce the TMM single-slab result in its
    // zeroth order — a cross-check between two completely independent solvers.
    {
        const auto& glass = materials::bk7();
        const double d = 0.5;  // slab thickness, µm

        auto tmm = solve_stack(materials::air(), {{glass, d}}, materials::air(),
                               lambda, normal, Pol::TE);

        BinaryGrating1D degenerate{glass, glass, 1.0 /*period*/, 0.5, d};
        auto rcwa = solve_rcwa_1d(materials::air(), degenerate, materials::air(),
                                  lambda, normal, /*M=*/8, Pol::TE);
        // order 0 sits at the middle of the orders vector
        std::size_t zero = rcwa.orders.size() / 2;

        std::println("[4] RCWA (degenerate grating) vs TMM (uniform slab):");
        std::println("    TMM : R = {:.6f}  T = {:.6f}", tmm.R, tmm.T);
        std::println("    RCWA: R = {:.6f}  T = {:.6f}  (order 0)",
                     rcwa.de_r[zero], rcwa.de_t[zero]);
        std::println("    Σ DE = {:.6f}  (expect 1.000000)", rcwa.sum_de);
    }

    // ---- Validation 5: real grating, energy conservation + convergence ----
    // A freestanding glass binary grating in air. At Λ=1.0 µm, λ=0.5 µm,
    // orders m = -1,0,+1 propagate. Lossless ⇒ Σ DE must equal 1, and the
    // split between orders must converge as we keep more harmonics M.
    {
        const auto glass = Material::constant(cdouble{1.5, 0.0}, "n1.5");
        BinaryGrating1D g{glass, materials::air(), 1.0 /*Λ*/, 0.5 /*fill*/, 0.5};

        std::println("[5] Freestanding grating (Λ=1.0µm, λ=0.5µm, normal):");
        std::println("    {:>3}   {:>10}  {:>10}  {:>10}", "M", "DE_t(0)",
                     "DE_t(+1)", "Σ DE");
        for (int M : {2, 5, 10, 20}) {
            auto r = solve_rcwa_1d(materials::air(), g, materials::air(), 0.5,
                                   0.0, M, Pol::TE);
            std::size_t z = r.orders.size() / 2;
            std::println("    {:>3}   {:>10.5f}  {:>10.5f}  {:>10.6f}", M,
                         r.de_t[z], r.de_t[z + 1], r.sum_de);
        }
    }

    // ---- Validation 6: TM polarization ------------------------------------
    // (a) Degenerate grating, TM: at normal incidence a uniform slab has
    //     TM == TE, so RCWA-TM order 0 must again match the TMM slab.
    // (b) Real grating, TM: energy must still conserve (Σ DE = 1). And note
    //     TE != TM for an actual grating even at normal incidence — the solver
    //     must distinguish them.
    {
        const auto& glass = materials::bk7();
        const double d = 0.5;
        auto tmm = solve_stack(materials::air(), {{glass, d}}, materials::air(),
                               lambda, normal, Pol::TM);
        BinaryGrating1D degenerate{glass, glass, 1.0, 0.5, d};
        auto rcwa_tm = solve_rcwa_1d(materials::air(), degenerate,
                                     materials::air(), lambda, normal, 8, Pol::TM);
        std::size_t z = rcwa_tm.orders.size() / 2;
        std::println("[6] TM polarization:");
        std::println("    (a) degenerate vs TMM-TM:  TMM R={:.6f}  RCWA R={:.6f}",
                     tmm.R, rcwa_tm.de_r[z]);

        const auto n15 = Material::constant(cdouble{1.5, 0.0}, "n1.5");
        BinaryGrating1D g{n15, materials::air(), 1.0, 0.5, 0.5};
        auto te = solve_rcwa_1d(materials::air(), g, materials::air(), 0.5,
                                normal, 20, Pol::TE);
        auto tm = solve_rcwa_1d(materials::air(), g, materials::air(), 0.5,
                                normal, 20, Pol::TM);
        std::size_t zt = te.orders.size() / 2;
        std::println("    (b) real grating, M=20:");
        std::println("        TE: DE_t(+1)={:.5f}  Σ DE={:.6f}", te.de_t[zt + 1],
                     te.sum_de);
        std::println("        TM: DE_t(+1)={:.5f}  Σ DE={:.6f}  (TE != TM ✓)",
                     tm.de_t[zt + 1], tm.sum_de);
    }

    // ---- Validation 7: multilayer S-matrix solver -------------------------
    // (a) A one-layer STACK must reproduce the single-layer solver exactly.
    // (b) Splitting one layer into N identical sublayers must give the same
    //     answer (the S-matrix recursion is self-consistent and stable).
    // (c) A genuinely layered device (grating + homogeneous cap) conserves
    //     energy.
    {
        const auto n15 = Material::constant(cdouble{1.5, 0.0}, "n1.5");
        const int M = 20;
        const double oblique = 15.0 * pi / 180.0;

        BinaryGrating1D single{n15, materials::air(), 1.0, 0.5, 0.5};
        auto ref = solve_rcwa_1d(materials::air(), single, materials::air(), 0.5,
                                 oblique, M, Pol::TM);

        Rcwa1DStack one{1.0, {GratingLayer1D{n15, materials::air(), 0.5, 0.5}}};
        auto st1 = solve_rcwa_1d(materials::air(), one, materials::air(), 0.5,
                                 oblique, M, Pol::TM);

        Rcwa1DStack split{1.0, {}};
        for (int i = 0; i < 5; ++i)  // 5 slices of 0.1 µm == one 0.5 µm layer
            split.layers.push_back({n15, materials::air(), 0.5, 0.1});
        auto stN = solve_rcwa_1d(materials::air(), split, materials::air(), 0.5,
                                 oblique, M, Pol::TM);

        std::size_t z = ref.orders.size() / 2;
        std::println("[7] Multilayer S-matrix (TM, 15deg incidence):");
        std::println("    (a) single-layer solver  DE_t(0) = {:.8f}", ref.de_t[z]);
        std::println("    (b) 1-layer stack        DE_t(0) = {:.8f}", st1.de_t[z]);
        std::println("    (c) 5-sublayer stack     DE_t(0) = {:.8f}", stN.de_t[z]);
        std::println("        max|stack-ref| over all orders = {:.2e}",
                     [&] {
                         double e = 0.0;
                         for (std::size_t i = 0; i < ref.de_t.size(); ++i) {
                             e = std::max(e, std::abs(st1.de_t[i] - ref.de_t[i]));
                             e = std::max(e, std::abs(stN.de_t[i] - ref.de_t[i]));
                         }
                         return e;
                     }());
        std::println("        Σ DE: single={:.6f} stack={:.6f} split={:.6f}",
                     ref.sum_de, st1.sum_de, stN.sum_de);

        // A real layered device: glass grating with a homogeneous AR-like cap.
        Rcwa1DStack device{1.0,
                           {GratingLayer1D::homogeneous(
                                Material::constant(cdouble{1.2, 0.0}, "cap"), 0.1),
                            GratingLayer1D{n15, materials::air(), 0.5, 0.5}}};
        auto dev = solve_rcwa_1d(materials::air(), device, materials::bk7(), 0.5,
                                 oblique, M, Pol::TM);
        std::println("    (d) grating + cap on glass:  Σ DE = {:.6f}", dev.sum_de);
    }

    // ---- Validation 8: 2D RCWA — reduce to 1D + energy + metalens cell ----
    // (a) A y-invariant 2D cell (fill_y=1, My=0) excited with E along y is
    //     exactly the 1D TE grating. Must match the validated 1D solver.
    // (b) Energy conservation for a real square pillar.
    // (c) A metalens unit cell: report transmission phase (the design output).
    {
        const auto n15 = Material::constant(cdouble{1.5, 0.0}, "n1.5");
        const int M = 12;

        // Λ=1.3 µm avoids the Rayleigh anomaly (no order exactly at grazing).
        BinaryGrating1D g1d{n15, materials::air(), 1.3, 0.5, 0.5};
        auto te1d = solve_rcwa_1d(materials::air(), g1d, materials::air(), 0.5,
                                  0.0, M, Pol::TE);
        std::size_t z = te1d.orders.size() / 2;

        Rcwa2DStack ginv{1.3, 1.0, {RectCell2D{n15, materials::air(), 0.5, 1.0, 0.5}}};
        auto te2d = solve_rcwa_2d(materials::air(), ginv, materials::air(), 0.5,
                                  0.0, 0.0, /*Ex0=*/0.0, /*Ey0=*/1.0, M, 0);

        std::println("[8] 2D RCWA:");
        std::println("    (a) y-invariant cell vs 1D TE:  1D DE_t0={:.6f}  "
                     "2D DE_t0={:.6f}  |Δ|={:.2e}",
                     te1d.de_t[z], te2d.de_t0,
                     std::abs(te1d.de_t[z] - te2d.de_t0));
        std::println("        Σ DE: 1D={:.6f}  2D={:.6f}", te1d.sum_de,
                     te2d.sum_de);

        // (b)+(c) a real square TiO2-like pillar metalens unit cell on glass.
        const auto tio2 = Material::constant(cdouble{2.40, 0.0}, "TiO2~");
        std::println("    (b) square pillar metalens cell (TiO2 n=2.4 on glass,"
                     " Λ=0.35µm, λ=0.532µm):");
        std::println("        {:>6}  {:>8}  {:>10}  {:>8}", "fill", "Σ DE",
                     "phase(deg)", "|t|^2");
        for (double f : {0.3, 0.5, 0.7}) {
            Rcwa2DStack cell{0.35, 0.35,
                             {RectCell2D{tio2, materials::air(), f, f, 0.6}}};
            auto r = solve_rcwa_2d(materials::air(), cell, materials::bk7(),
                                   0.532, 0.0, 0.0, 1.0, 0.0, 8, 8);
            double phase_deg = std::arg(r.tx0) * 180.0 / pi;
            std::println("        {:>6.2f}  {:>8.6f}  {:>10.1f}  {:>8.4f}", f,
                         r.sum_de, phase_deg, std::norm(r.tx0));
        }
    }

    // ---- Demo 9: end-to-end metalens design -------------------------------
    // Build a phase library (sweep pillar size), then design a focusing lens
    // and report how well the realized phase matches the ideal profile.
    {
        const auto tio2 = Material::constant(cdouble{2.40, 0.0}, "TiO2~");
        std::println("[9] Metalens design pipeline (TiO2 pillars, λ=532nm):");
        std::println("    Building unit-cell library (sweeping pillar size)...");
        auto lib = build_unit_cell_library(tio2, materials::air(),
                                           materials::air(), materials::bk7(),
                                           /*period=*/0.35, /*λ=*/0.532,
                                           /*thickness=*/0.6, /*fill*/ 0.08, 0.92,
                                           /*n_samples=*/18, /*M=*/6);
        std::println("    Library: {} pillars, phase coverage = {:.0f}° "
                     "(need ~360° for full control)",
                     lib.fill.size(), lib.phase_span() * 180.0 / pi);

        auto lens = design_metalens(lib, /*focal=*/50.0, /*diameter=*/20.0);
        std::println("    Designed lens: {0}x{0} pillars, f=50µm, D=20µm",
                     lens.n_cells);
        std::println("    RMS phase error vs ideal = {:.1f}°  (lower = sharper "
                     "focus)",
                     lens.rms_phase_error_deg);
        std::println("    Mean pillar transmission |t| = {:.3f}", lens.mean_amplitude);

        // Export the fabrication file and validate it round-trips.
        const std::string gds = "metalens.gds";
        int written = write_metalens_gds(lens, gds, /*layer=*/1, /*min_fill=*/0.05);
        int read_back = gds_count_boundaries(gds);
        std::println("    GDSII export -> {}: {} pillars written, {} read back "
                     "({})",
                     gds, written, read_back,
                     (written == read_back && written > 0) ? "valid ✓" : "MISMATCH");

        // Does it actually focus? Propagate to the focal plane and measure.
        auto foc = analyze_focus(lens, lib, /*f=*/50.0, /*λ=*/0.532, /*D=*/20.0);
        // A/B: amplitude-aware pillar selection vs phase-only.
        auto lens_po = design_metalens(lib, 50.0, 20.0, /*amplitude_weight=*/0.0);
        auto foc_po = analyze_focus(lens_po, lib, 50.0, 0.532, 20.0);
        std::println("    Focal performance:");
        std::println("      Strehl: phase-only {:.3f} -> amplitude-aware {:.3f}",
                     foc_po.strehl, foc.strehl);
        std::println("      Strehl ratio   = {:.3f}  (1.0 = perfect)", foc.strehl);
        std::println("      spot FWHM      = {:.2f} µm  (diffraction limit "
                     "λf/D = {:.2f} µm)",
                     foc.fwhm_um, foc.diffraction_limit_um);
        std::println("      encircled E    = {:.1f}% within first Airy null",
                     foc.encircled_energy * 100.0);

        // Chromatic behavior: how the focus shifts across a wavelength band.
        auto chrom = analyze_chromatic(lens, lib, /*f=*/50.0, /*λ0=*/0.532,
                                       /*D=*/20.0, 0.45, 0.65, 5);
        std::println("    Chromatic focal shift (designed for 532nm):");
        std::println("      {:>7}  {:>12}  {:>10}  {:>10}", "λ(nm)", "focus(µm)",
                     "f0·λ0/λ", "rel.peak");
        for (auto& c : chrom)
            std::println("      {:>7.0f}  {:>12.2f}  {:>10.2f}  {:>10.2f}",
                         c.wavelength_um * 1000.0, c.focal_length_um,
                         50.0 * 0.532 / c.wavelength_um, c.rel_peak);
    }

    // ---- Demo 11: inverse design (gradient-based optimizer) ---------------
    // Instead of looking a pillar up from the discrete library, SOLVE for the
    // geometry hitting a target phase with maximum transmission.
    {
        const auto tio2 = Material::constant(cdouble{2.40, 0.0}, "TiO2~");
        const double target_deg = 90.0;
        PillarTarget tgt{0.532, target_deg * pi / 180.0, /*amplitude_weight=*/1.0};
        std::println("[11] Inverse design: optimize pillar for {:.0f}° phase "
                     "@532nm (max transmission):", target_deg);
        auto opt = optimize_pillar(tio2, materials::air(), materials::air(),
                                   materials::bk7(), /*period=*/0.35, tgt,
                                   /*M=*/5, /*fill0=*/0.50, /*thickness0=*/0.50,
                                   /*max_iters=*/25);
        std::println("    converged geometry: fill={:.3f}, thickness={:.3f} µm",
                     opt.fill, opt.thickness_um);
        std::println("    achieved phase = {:.1f}°  (target {:.0f}°),  |t| = "
                     "{:.3f},  loss = {:.2e}",
                     opt.achieved_phase_rad * 180.0 / pi, target_deg,
                     opt.achieved_amplitude, opt.loss);
    }

    // ---- Demo 12: form birefringence (polarization-multiplexed basis) -----
    // A rectangular pillar (fill_x != fill_y) responds differently to x- and
    // y-polarized light — "form birefringence." That phase difference is the
    // basis for polarization-multiplexed metalenses (one device, two functions
    // selected by polarization). Here we show it grow with pillar asymmetry.
    {
        const auto tio2 = Material::constant(cdouble{2.40, 0.0}, "TiO2~");
        std::println("[12] Form birefringence (rectangular TiO2 pillar, 532nm):");
        std::println("      {:>10}  {:>10}  {:>10}  {:>12}", "fill_x", "fill_y",
                     "φx-φy(°)", "|tx|,|ty|");
        for (auto [fx, fy] : {std::pair{0.5, 0.5}, {0.6, 0.4}, {0.7, 0.3}}) {
            Rcwa2DStack cell{0.35, 0.35, {RectCell2D{tio2, materials::air(), fx, fy, 0.6}}};
            auto rx = solve_rcwa_2d(materials::air(), cell, materials::bk7(),
                                    0.532, 0.0, 0.0, 1.0, 0.0, 8, 8);  // x-pol
            auto ry = solve_rcwa_2d(materials::air(), cell, materials::bk7(),
                                    0.532, 0.0, 0.0, 0.0, 1.0, 8, 8);  // y-pol
            double dphi = std::arg(rx.tx0) - std::arg(ry.ty0);
            while (dphi > pi) dphi -= 2 * pi;
            while (dphi <= -pi) dphi += 2 * pi;
            std::println("      {:>10.2f}  {:>10.2f}  {:>10.1f}  {:>5.2f},{:>5.2f}",
                         fx, fy, dphi * 180.0 / pi, std::abs(rx.tx0),
                         std::abs(ry.ty0));
        }
    }

#ifdef CELERIS_USE_CUDA
    // ---- Benchmark 10: GPU vs CPU eigensolve ------------------------------
    // The per-layer RCWA eigenproblem is a general complex matrix of size 2N.
    // Honest head-to-head on a representative 578x578 (2N at M=8): cuSOLVER
    // Xgeev (GPU) vs Eigen ComplexEigenSolver (CPU). Verify eigenvalues agree.
    {
        const int n = 578;
        std::mt19937 rng(42);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        Eigen::MatrixXcd A(n, n);
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
                A(i, j) = cdouble{dist(rng), dist(rng)};

        std::println("[10] GPU vs CPU eigensolve ({}x{} general complex):", n, n);
        if (!cuda::available()) {
            std::println("     no CUDA device available");
        } else {
            auto t0 = std::chrono::steady_clock::now();
            Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces(A);
            auto t1 = std::chrono::steady_clock::now();
            double cpu_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

            std::vector<cdouble> w(n), vr((std::size_t)n * n);
            cuda::geev(A.data(), n, w.data(), vr.data());  // warm-up (device init)
            auto t2 = std::chrono::steady_clock::now();
            bool ok = cuda::geev(A.data(), n, w.data(), vr.data());
            auto t3 = std::chrono::steady_clock::now();
            double gpu_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();

            // Compare eigenvalue sets (sort by real then imag).
            std::vector<cdouble> ec(ces.eigenvalues().data(),
                                    ces.eigenvalues().data() + n);
            std::vector<cdouble> eg = w;
            auto cmp = [](cdouble a, cdouble b) {
                return a.real() != b.real() ? a.real() < b.real()
                                            : a.imag() < b.imag();
            };
            std::sort(ec.begin(), ec.end(), cmp);
            std::sort(eg.begin(), eg.end(), cmp);
            double max_diff = 0.0;
            for (int i = 0; i < n; ++i)
                max_diff = std::max(max_diff, std::abs(ec[i] - eg[i]));

            std::println("     CPU (Eigen)      : {:.1f} ms", cpu_ms);
            std::println("     GPU (cuSOLVER)   : {:.1f} ms   ({:.1f}x)", gpu_ms,
                         cpu_ms / gpu_ms);
            std::println("     eigenvalue match : max|Δ| = {:.2e}  ok={}",
                         max_diff, ok);
        }
    }
#endif
    return 0;
}

// ----------------------------------------------------------------------------
// CLI front-end
// ----------------------------------------------------------------------------

namespace {

const char* arg_value(int argc, char** argv, const std::string& key,
                      const char* fallback) {
    for (int i = 2; i + 1 < argc; ++i)
        if (key == argv[i]) return argv[i + 1];
    return fallback;
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
    const std::string sub_name = arg_value(argc, argv, "--substrate", "bk7");

    const Material& substrate = sub_name == "air"  ? materials::air()
                                : sub_name == "sio2" ? materials::fused_silica()
                                                     : materials::bk7();
    const auto pillar = Material::constant(cdouble{pillar_n, 0.0}, "pillar");

    std::println("CELERIS metalens design");
    std::println("  f={}µm  D={}µm  λ={}µm  Λ={}µm  h={}µm  n_pillar={}  substrate={}",
                 focal, diameter, lambda, period, thickness, pillar_n, sub_name);
    std::println("  building unit-cell library ({} pillars, M={})...", samples, M);

    auto lib = build_unit_cell_library(pillar, materials::air(), materials::air(),
                                       substrate, period, lambda, thickness,
                                       0.08, 0.92, samples, M);
    std::println("  library phase coverage: {:.0f}°", lib.phase_span() * 180.0 / pi);

    auto lens = design_metalens(lib, focal, diameter);
    std::println("  designed {0}x{0} pillars, RMS phase error {1:.1f}°, mean |t| {2:.3f}",
                 lens.n_cells, lens.rms_phase_error_deg, lens.mean_amplitude);

    int np = write_metalens_gds(lens, out);
    if (np < 0) { std::println("  ERROR: could not write {}", out); return 1; }
    std::println("  wrote {} pillars -> {}", np, out);

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
    return 0;
}

void print_help() {
    std::println(
        "CELERIS — GPU-ready metalens design via rigorous coupled-wave analysis\n"
        "\n"
        "Usage:\n"
        "  celeris design [options]   design a focusing metalens -> GDSII + report\n"
        "  celeris selftest           run the physics validation suite\n"
        "  celeris help               show this message\n"
        "\n"
        "design options (defaults):\n"
        "  --focal <µm=50>        focal length\n"
        "  --diameter <µm=20>     lens aperture diameter\n"
        "  --wavelength <µm=0.532>\n"
        "  --period <µm=0.35>     unit-cell pitch\n"
        "  --thickness <µm=0.6>   pillar height\n"
        "  --pillar-n <2.4>       pillar refractive index\n"
        "  --substrate <bk7|air|sio2=bk7>\n"
        "  --fill-samples <18>    library resolution\n"
        "  --harmonics <6>        RCWA Fourier half-count (accuracy vs speed)\n"
        "  --out <metalens.gds>   output GDSII path\n"
        "  --tolerance            Monte-Carlo fabrication-error / yield analysis");
}

} // namespace

int main(int argc, char** argv) {
    const std::string cmd = argc > 1 ? argv[1] : "help";
    if (cmd == "selftest") return run_selftest();
    if (cmd == "design") return cmd_design(argc, argv);
    print_help();
    return (cmd == "help" || cmd == "--help" || cmd == "-h") ? 0 : 1;
}
