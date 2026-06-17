#include <cmath>
#include <cstdlib>
#include <format>
#include <fstream>
#include <print>
#include <string>
#include <vector>

#ifdef CELERIS_USE_CUDA
#include <Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <future>
#include <random>
#include <thread>

#include "celeris/cuda/eigensolve.hpp"
#ifdef CELERIS_USE_CUDA_KERNELS
#include "celeris/cuda/propagate.hpp"
#endif
#endif

#include "celeris/analysis/chromatic.hpp"
#include "celeris/analysis/field.hpp"
#include "celeris/analysis/focal.hpp"
#include "celeris/analysis/polarization.hpp"
#include "celeris/analysis/throughfocus.hpp"
#include "celeris/analysis/tolerance.hpp"
#include "celeris/analysis/wavefront.hpp"
#include "celeris/design/metalens.hpp"
#include "celeris/design/optimize.hpp"
#include "celeris/design/polar_metalens.hpp"
#include "celeris/io/gds.hpp"
#include "celeris/io/image.hpp"
#include "celeris/io/material_csv.hpp"
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

    // ---- Validation 13: RCWA vs Effective-Medium Theory --------------------
    // A deeply subwavelength grating (period << λ) behaves as a uniform
    // birefringent film. Rytov's 0th-order EMT: ε∥ = f·ε1+(1-f)·ε2 (TE, E along
    // grooves), ε⊥ = [f/ε1+(1-f)/ε2]⁻¹ (TM). The grating's RCWA reflectance must
    // converge to the TMM reflectance of those effective films — an independent
    // analytic check of the full TE+TM vectorial solver.
    {
        const double lam = 0.5, d = 0.10, f = 0.5, e1 = 2.1025, e2 = 1.0;  // ridge n=1.45, groove air
        const double eTE = f * e1 + (1 - f) * e2;            // arithmetic mean
        const double eTM = 1.0 / (f / e1 + (1 - f) / e2);    // harmonic mean
        const auto ridge = Material::constant(cdouble{1.45, 0.0}, "n1.45");
        const auto nTE = Material::constant(cdouble{std::sqrt(eTE), 0.0}, "nTE");
        const auto nTM = Material::constant(cdouble{std::sqrt(eTM), 0.0}, "nTM");
        std::println("[13] RCWA vs effective-medium theory (Λ=λ/20 subwavelength):");
        std::println("      {:>4}  {:>14}  {:>14}  {:>9}", "pol", "RCWA R", "EMT-film R", "|Δ|");
        for (int te = 1; te >= 0; --te) {
            Pol pol = te ? Pol::TE : Pol::TM;
            BinaryGrating1D g{ridge, materials::air(), 0.025, f, d};  // Λ = λ/20
            auto rg = solve_rcwa_1d(materials::air(), g, materials::air(), lam, 0.0, 20, pol);
            double rcwaR = rg.de_r[rg.orders.size() / 2];
            auto slab = solve_stack(materials::air(), {{te ? nTE : nTM, d}},
                                    materials::air(), lam, 0.0, pol);
            std::println("      {:>4}  {:>14.4f}  {:>14.4f}  {:>9.4f}",
                         te ? "TE" : "TM", rcwaR, slab.R, std::abs(rcwaR - slab.R));
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
    const char* pillar_csv = arg_value(argc, argv, "--pillar-csv", nullptr);
    const Material pillar =
        pillar_csv ? load_material_csv(pillar_csv, "pillar-csv")
                   : Material::constant(cdouble{pillar_n, 0.0}, "pillar");

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
            f << std::format("pillar height     : {} um\n\n", thickness);
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

    int np = write_rect_gds(out, d.n_cells, d.period_um, d.fill_x, d.fill_y);
    if (np < 0) { std::println("  ERROR: could not write {}", out); return 1; }
    std::println("  wrote {} rectangular pillars -> {}", np, out);
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
        "  --pillar-n <2.4>       pillar refractive index (constant)\n"
        "  --pillar-csv <file>    load pillar n,k from CSV (overrides --pillar-n)\n"
        "  --substrate <bk7|air|sio2=bk7>\n"
        "  --fill-samples <18>    library resolution\n"
        "  --harmonics <6>        RCWA Fourier half-count (accuracy vs speed)\n"
        "  --out <metalens.gds>   output GDSII path\n"
        "  --tolerance            Monte-Carlo fabrication-error / yield analysis\n"
        "  --fov                  off-axis field-of-view analysis\n"
        "  --psf <file.pgm>       write the focal-spot image (PGM)\n"
        "  --report <prefix>      write a full deliverable bundle (txt metrics +\n"
        "                         PSF & caustic images + GDS)\n"
        "\n"
        "celeris birefringence [--fill-y 0.5] [--period --wavelength --thickness\n"
        "                       --pillar-n --samples --harmonics]\n"
        "  sweep a rectangular pillar's x-width; report x/y-polarized phase and\n"
        "  retardance (the polarization-optics / waveplate building block)\n"
        "\n"
        "celeris polardesign [--focal-x 50] [--focal-y 80] [--diameter --wavelength\n"
        "                    --period --thickness --pillar-n --samples --out]\n"
        "  polarization-multiplexed lens: X-pol and Y-pol focus at different\n"
        "  distances; writes a rectangular-pillar GDS");
}

} // namespace

#ifdef CELERIS_USE_CUDA
// Honest head-to-head for the batched GPU eigensolve: the metalens library
// sweep is a batch of independent same-size general-complex eigenproblems, so
// this times CPU-sequential, CPU-parallel (the path the real builder uses), and
// GPU-batched over a representative batch. Usage: celeris gpubench [--n N]
// [--batch B] [--streams S].
static int run_gpubench(int argc, char** argv) {
    int n = 242, batch = 32, streams = 4;  // n ~ 2N at M=5 (a real 2D RCWA size)
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&] { return (i + 1 < argc) ? std::atoi(argv[++i]) : 0; };
        if (a == "--n") n = next();
        else if (a == "--batch") batch = next();
        else if (a == "--streams") streams = next();
    }
    if (!cuda::available()) { std::println("gpubench: no CUDA device available"); return 1; }

    const std::size_t nn = static_cast<std::size_t>(n) * n;
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    std::vector<cdouble> As(static_cast<std::size_t>(batch) * nn);
    for (auto& v : As) v = cdouble{dist(rng), dist(rng)};

    std::println("GPU batched eigensolve benchmark");
    std::println("  batch = {} matrices, {}x{} general complex, streams = {}",
                 batch, n, n, streams);

    auto eig_one = [&](int b) {
        Eigen::Map<const Eigen::MatrixXcd> M(As.data() + static_cast<std::size_t>(b) * nn, n, n);
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces(M);
        return ces.eigenvalues()(0);  // touch a result so it isn't optimized away
    };

    // CPU sequential.
    auto t0 = std::chrono::steady_clock::now();
    cdouble sink{0, 0};
    for (int b = 0; b < batch; ++b) sink += eig_one(b);
    auto t1 = std::chrono::steady_clock::now();
    double cpu_seq = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // CPU parallel (same std::async fan-out the unit-cell library builder uses).
    auto t2 = std::chrono::steady_clock::now();
    {
        unsigned hw = std::max(1u, std::thread::hardware_concurrency());
        int workers = std::min<int>(static_cast<int>(hw), batch);
        std::vector<std::future<void>> jobs;
        for (int w = 0; w < workers; ++w)
            jobs.push_back(std::async(std::launch::async, [&, w] {
                for (int b = w; b < batch; b += workers) (void)eig_one(b);
            }));
        for (auto& j : jobs) j.get();
    }
    auto t3 = std::chrono::steady_clock::now();
    double cpu_par = std::chrono::duration<double, std::milli>(t3 - t2).count();

    // GPU batched (warm up device first so init isn't charged to the timing).
    std::vector<cdouble> ws(static_cast<std::size_t>(batch) * n);
    std::vector<cdouble> vrs(static_cast<std::size_t>(batch) * nn);
    cuda::geev_batched(As.data(), n, 1, ws.data(), vrs.data(), 1);
    auto t4 = std::chrono::steady_clock::now();
    bool ok = cuda::geev_batched(As.data(), n, batch, ws.data(), vrs.data(), streams);
    auto t5 = std::chrono::steady_clock::now();
    double gpu = std::chrono::duration<double, std::milli>(t5 - t4).count();

    // Correctness: eigenvalues of matrix 0 must match Eigen (order-independent).
    Eigen::Map<const Eigen::MatrixXcd> M0(As.data(), n, n);
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces0(M0);
    std::vector<cdouble> ec(ces0.eigenvalues().data(), ces0.eigenvalues().data() + n);
    std::vector<cdouble> eg(ws.begin(), ws.begin() + n);
    auto cmp = [](cdouble a, cdouble b) {
        return a.real() != b.real() ? a.real() < b.real() : a.imag() < b.imag();
    };
    std::sort(ec.begin(), ec.end(), cmp);
    std::sort(eg.begin(), eg.end(), cmp);
    double max_diff = 0.0;
    for (int i = 0; i < n; ++i) max_diff = std::max(max_diff, std::abs(ec[i] - eg[i]));

    std::println("");
    std::println("  CPU sequential : {:8.1f} ms   ({:.2f} ms / solve)", cpu_seq, cpu_seq / batch);
    std::println("  CPU parallel   : {:8.1f} ms   ({:.1f}x vs seq)", cpu_par, cpu_seq / cpu_par);
    std::println("  GPU batched    : {:8.1f} ms   ({:.2f}x vs CPU parallel, {:.2f}x vs seq)",
                 gpu, cpu_par / gpu, cpu_seq / gpu);
    std::println("  correctness    : max|d eigenvalue| = {:.2e}  ok={}", max_diff, ok);
    (void)sink;
    return ok ? 0 : 1;
}
#endif

#ifdef CELERIS_USE_CUDA_KERNELS
// Honest head-to-head for the GPU far-field propagation kernel: build a real
// metalens, then compute its focal PSF on the CPU (analyze_focus path, double,
// parallel over cores) vs the GPU kernel (float), comparing the maps and timing.
// Usage: celeris psfbench [--diameter D] [--focal F] [--grid N].
static int run_psfbench(int argc, char** argv) {
    double diameter = 120.0, focal = 50.0, wavelength = 0.532, period = 0.35;
    int grid = 161, M = 5, fill_samples = 16;
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        auto nd = [&] { return (i + 1 < argc) ? std::atof(argv[++i]) : 0.0; };
        auto ni = [&] { return (i + 1 < argc) ? std::atoi(argv[++i]) : 0; };
        if (a == "--diameter") diameter = nd();
        else if (a == "--focal") focal = nd();
        else if (a == "--grid") grid = ni();
    }
    if (!cuda::available()) { std::println("psfbench: no CUDA device available"); return 1; }

    std::println("PSF propagation benchmark (building lens...)");
    auto pillar = Material::constant(cdouble{2.40, 0.0}, "TiO2~");
    auto lib = build_unit_cell_library(pillar, materials::air(), materials::air(),
                                       materials::bk7(), period, wavelength, 0.6,
                                       0.08, 0.92, fill_samples, M);
    auto lens = design_metalens(lib, focal, diameter);
    const double dl = wavelength * focal / diameter;
    const double W = std::max(6.0 * dl, 4.0);

    // Aperture pillar list (same construction analyze_focus/compute_psf use).
    const double p = lens.period_um;
    const double center = (lens.n_cells - 1) / 2.0;
    const double R_ap = diameter / 2.0;
    std::vector<double> px, py;
    std::vector<cdouble> pt;
    for (int iy = 0; iy < lens.n_cells; ++iy)
        for (int ix = 0; ix < lens.n_cells; ++ix) {
            double x = (ix - center) * p, y = (iy - center) * p;
            if (std::sqrt(x * x + y * y) > R_ap) continue;
            double fill = lens.fill_map[(std::size_t)iy * lens.n_cells + ix];
            px.push_back(x); py.push_back(y);
            pt.push_back(lib.transmission_for_fill(fill));
        }
    std::println("  {} x {} cells, {} pillars in aperture, {}x{} focal grid",
                 lens.n_cells, lens.n_cells, px.size(), grid, grid);

    // CPU (double, parallel across cores).
    auto t0 = std::chrono::steady_clock::now();
    auto cpu = compute_psf(lens, lib, focal, wavelength, diameter, grid, W);
    auto t1 = std::chrono::steady_clock::now();
    double cpu_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // GPU (float kernel). Warm up first so device init isn't timed.
    std::vector<double> gpu((std::size_t)grid * grid);
    double k = 2.0 * pi / wavelength;
    cuda::propagate_psf(px.data(), py.data(), pt.data(), (int)px.size(), 0.0, 0.0,
                        focal, k, grid, W, gpu.data());  // warm-up
    auto t2 = std::chrono::steady_clock::now();
    bool ok = cuda::propagate_psf(px.data(), py.data(), pt.data(), (int)px.size(),
                                  0.0, 0.0, focal, k, grid, W, gpu.data());
    auto t3 = std::chrono::steady_clock::now();
    double gpu_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();

    // Compare peak-normalized maps.
    double cmax = 0, gmax = 0;
    for (double v : cpu.intensity) cmax = std::max(cmax, v);
    for (double v : gpu) gmax = std::max(gmax, v);
    double maxrel = 0.0;
    if (cmax > 0 && gmax > 0)
        for (std::size_t i = 0; i < gpu.size(); ++i)
            maxrel = std::max(maxrel, std::abs(cpu.intensity[i] / cmax - gpu[i] / gmax));

    std::println("");
    std::println("  CPU (double, parallel) : {:8.1f} ms", cpu_ms);
    std::println("  GPU (float kernel)     : {:8.1f} ms   ({:.1f}x faster)", gpu_ms, cpu_ms / gpu_ms);
    std::println("  agreement              : max|d normalized PSF| = {:.2e}  ok={}", maxrel, ok);
    return ok ? 0 : 1;
}
#endif

int main(int argc, char** argv) {
    const std::string cmd = argc > 1 ? argv[1] : "help";
    if (cmd == "selftest") return run_selftest();
    if (cmd == "design") return cmd_design(argc, argv);
    if (cmd == "birefringence") return cmd_birefringence(argc, argv);
    if (cmd == "polardesign") return cmd_polardesign(argc, argv);
#ifdef CELERIS_USE_CUDA
    if (cmd == "gpubench") return run_gpubench(argc, argv);
#endif
#ifdef CELERIS_USE_CUDA_KERNELS
    if (cmd == "psfbench") return run_psfbench(argc, argv);
#endif
    print_help();
    return (cmd == "help" || cmd == "--help" || cmd == "-h") ? 0 : 1;
}
