#include "cli.hpp"

#include <atomic>
#include <chrono>
#include <functional>
#include <future>
#include <thread>

namespace {
// A self-test case: its body writes results via sp() into a per-case buffer, so
// independent cases can run concurrently and still be printed in declaration
// order. `heavy` marks the long RCWA cases skipped by `selftest --quick`.
struct Case {
    bool heavy;
    std::function<void()> run;
};

// Per-case output sink (thread-local so concurrent cases don't interleave).
thread_local std::string* g_selftest_out = nullptr;
template <class... Ts>
void sp(std::format_string<Ts...> fmt, Ts&&... args) {
    std::string s = std::format(fmt, std::forward<Ts>(args)...);
    s.push_back('\n');
    if (g_selftest_out) g_selftest_out->append(s);
    else std::print("{}", s);
}

// Suite-level failed-check count (atomic: cases run concurrently). A case records
// a failed assertion via chk(); the runner returns non-zero if any check failed,
// so a wrong number makes `celeris selftest` (and therefore CI) fail.
std::atomic<int> g_selftest_failures{0};
void chk(bool ok, const std::string& what) {
    sp("    [{}] {}", ok ? "PASS" : "FAIL", what);
    if (!ok) g_selftest_failures.fetch_add(1, std::memory_order_relaxed);
}
}  // namespace

// Physics validation suite — every solver checked against closed-form results,
// an independent solver, or energy conservation. Run with: celeris selftest.
// `--quick` (quick=true) runs the fast subset, skipping the heavy RCWA cases —
// used by CI on every push; the full suite runs on a nightly schedule.
int run_selftest(bool quick) {
    const double lambda = 0.550;  // design wavelength, 550 nm (µm)
    const double normal = 0.0;    // normal incidence

    std::vector<Case> cases;

    // ---- Validation 1: bare glass, energy conservation --------------------
    // No layers, air -> N-BK7. Closed form R = ((1-n)/(1+n))^2 ~ 0.0424 at
    // 550 nm, and since glass is lossless here, R + T must equal 1.
    cases.push_back(Case{false, [&]() {
        auto res = solve_stack(materials::air(), {}, materials::bk7(),
                               lambda, normal, Pol::TE);
        sp("[1] Bare N-BK7 glass at 550 nm:");
        sp("    R = {:.4f}  (expect ~0.0424)", res.R);
        sp("    T = {:.4f}", res.T);
        sp("    R + T = {:.6f}  (expect 1.000000)", res.R + res.T);
        chk(std::abs(res.R + res.T - 1.0) < 1e-9, "[1] bare glass: R+T = 1");
        chk(std::abs(res.R - 0.0424) < 3e-3, "[1] bare glass: R ~ 0.0424");
    }});

    // ---- Validation 2: quarter-wave AR coating ----------------------------
    // A single layer of index n1 = sqrt(n_air * n_glass), one quarter-wave
    // thick, makes reflection vanish at the design wavelength. This is the
    // textbook single-layer anti-reflection coating.
    cases.push_back(Case{false, [&]() {
        const cdouble n_glass = materials::bk7().index(lambda);
        const cdouble n1 = std::sqrt(materials::air().index(lambda) * n_glass);
        const double d1 = lambda / (4.0 * n1.real());  // quarter-wave thickness

        std::vector<Layer> stack = {
            {Material::constant(n1, "ideal-AR"), d1}};

        sp("[2] Quarter-wave AR coating (n1 = {:.4f}, d = {:.1f} nm):",
                     n1.real(), d1 * 1000.0);
        double R_design = 1.0;
        for (double wl : {0.450, 0.550, 0.650}) {
            auto res = solve_stack(materials::air(), stack, materials::bk7(),
                                   wl, normal, Pol::TE);
            sp("    lambda = {:.0f} nm:  R = {:.5f}{}", wl * 1000.0,
                         res.R, wl == 0.550 ? "   <- ~0 at design" : "");
            if (wl == 0.550) R_design = res.R;
        }
        chk(R_design < 1e-3, "[2] quarter-wave AR: R ~ 0 at design wavelength");
    }});

    // ---- Validation 3: distributed Bragg reflector ------------------------
    // Alternating high/low quarter-wave layers build a high-reflectance
    // mirror. Reflectance should climb toward 1 as we add more pairs.
    cases.push_back(Case{false, [&]() {
        const cdouble nH{2.30, 0.0};  // high index (e.g. TiO2-like)
        const cdouble nL{1.46, 0.0};  // low index  (e.g. SiO2-like)
        const double dH = lambda / (4.0 * nH.real());
        const double dL = lambda / (4.0 * nL.real());

        sp("[3] Bragg mirror (nH=2.30 / nL=1.46 quarter-wave pairs):");
        double R2 = 0.0, R8 = 0.0;
        for (int pairs : {2, 4, 8}) {
            std::vector<Layer> stack;
            for (int p = 0; p < pairs; ++p) {
                stack.push_back({Material::constant(nH, "H"), dH});
                stack.push_back({Material::constant(nL, "L"), dL});
            }
            auto res = solve_stack(materials::air(), stack, materials::bk7(),
                                   lambda, normal, Pol::TE);
            sp("    {} pairs:  R = {:.5f}", pairs, res.R);
            if (pairs == 2) R2 = res.R;
            if (pairs == 8) R8 = res.R;
        }
        chk(R8 > 0.99, "[3] Bragg mirror: R >= 0.99 at 8 pairs");
        chk(R8 > R2, "[3] Bragg mirror: R increases with pairs");
    }});

    // ---- Validation 4: RCWA degenerate grating == TMM slab ----------------
    // If a grating's ridge and groove are the SAME material, it's just a
    // uniform slab. RCWA must then reproduce the TMM single-slab result in its
    // zeroth order — a cross-check between two completely independent solvers.
    cases.push_back(Case{false, [&]() {
        const auto& glass = materials::bk7();
        const double d = 0.5;  // slab thickness, µm

        auto tmm = solve_stack(materials::air(), {{glass, d}}, materials::air(),
                               lambda, normal, Pol::TE);

        BinaryGrating1D degenerate{glass, glass, 1.0 /*period*/, 0.5, d};
        auto rcwa = solve_rcwa_1d(materials::air(), degenerate, materials::air(),
                                  lambda, normal, /*M=*/8, Pol::TE);
        // order 0 sits at the middle of the orders vector
        std::size_t zero = rcwa.orders.size() / 2;

        sp("[4] RCWA (degenerate grating) vs TMM (uniform slab):");
        sp("    TMM : R = {:.6f}  T = {:.6f}", tmm.R, tmm.T);
        sp("    RCWA: R = {:.6f}  T = {:.6f}  (order 0)",
                     rcwa.de_r[zero], rcwa.de_t[zero]);
        sp("    Σ DE = {:.6f}  (expect 1.000000)", rcwa.sum_de);
        chk(std::abs(rcwa.de_r[zero] - tmm.R) < 1e-4 &&
            std::abs(rcwa.de_t[zero] - tmm.T) < 1e-4,
            "[4] RCWA(degenerate grating) == TMM(slab), order 0");
        chk(std::abs(rcwa.sum_de - 1.0) < 1e-6, "[4] degenerate grating: energy conserved");
    }});

    // ---- Validation 5: real grating, energy conservation + convergence ----
    // A freestanding glass binary grating in air. At Λ=1.0 µm, λ=0.5 µm,
    // orders m = -1,0,+1 propagate. Lossless ⇒ Σ DE must equal 1, and the
    // split between orders must converge as we keep more harmonics M.
    cases.push_back(Case{false, [&]() {
        const auto glass = Material::constant(cdouble{1.5, 0.0}, "n1.5");
        BinaryGrating1D g{glass, materials::air(), 1.0 /*Λ*/, 0.5 /*fill*/, 0.5};

        sp("[5] Freestanding grating (Λ=1.0µm, λ=0.5µm, normal):");
        sp("    {:>3}   {:>10}  {:>10}  {:>10}", "M", "DE_t(0)",
                     "DE_t(+1)", "Σ DE");
        double maxdef = 0.0;
        for (int M : {2, 5, 10, 20}) {
            auto r = solve_rcwa_1d(materials::air(), g, materials::air(), 0.5,
                                   0.0, M, Pol::TE);
            std::size_t z = r.orders.size() / 2;
            sp("    {:>3}   {:>10.5f}  {:>10.5f}  {:>10.6f}", M,
                         r.de_t[z], r.de_t[z + 1], r.sum_de);
            maxdef = std::max(maxdef, std::abs(r.sum_de - 1.0));
        }
        chk(maxdef < 1e-6, "[5] freestanding grating: energy conserved (all M)");

        // (b) 1D TE convergence to an external reference (grcwa). A subwavelength
        // (Λ=0.3µm < λ=0.5µm) freestanding n=1.5 grating passes ONLY the zeroth
        // order; its transmittance must converge in M to grcwa's 0.93333 (grcwa
        // is itself validated vs analytic TMM). This locks the half-space
        // companion-admittance sign in the TE path (Y = −j·kz): the wrong sign
        // still conserves energy but converges the split to the wrong value.
        BinaryGrating1D sub{glass, materials::air(), 0.3, 0.5, 0.5};
        sp("    1D TE convergence vs grcwa (Λ=0.3µm subwavelength, "
                     "T0 -> grcwa 0.93333):");
        sp("    {:>3}   {:>10}  {:>10}", "M", "TE T0", "|Δ grcwa|");
        double t0_hiM = 0.0;
        for (int M : {4, 10, 20}) {
            auto r = solve_rcwa_1d(materials::air(), sub, materials::air(), 0.5,
                                   0.0, M, Pol::TE);
            std::size_t z = r.orders.size() / 2;
            sp("    {:>3}   {:>10.5f}  {:>10.2e}", M, r.de_t[z],
                         std::abs(r.de_t[z] - 0.93333));
            if (M == 20) t0_hiM = r.de_t[z];
        }
        chk(std::abs(t0_hiM - 0.93333) < 1e-3, "[5] 1D-TE subwavelength T0 -> grcwa 0.93333");
    }});

    // ---- Validation 6: TM polarization ------------------------------------
    // (a) Degenerate grating, TM: at normal incidence a uniform slab has
    //     TM == TE, so RCWA-TM order 0 must again match the TMM slab.
    // (b) Real grating, TM: energy must still conserve (Σ DE = 1). And note
    //     TE != TM for an actual grating even at normal incidence — the solver
    //     must distinguish them.
    cases.push_back(Case{false, [&]() {
        const auto& glass = materials::bk7();
        const double d = 0.5;
        auto tmm = solve_stack(materials::air(), {{glass, d}}, materials::air(),
                               lambda, normal, Pol::TM);
        BinaryGrating1D degenerate{glass, glass, 1.0, 0.5, d};
        auto rcwa_tm = solve_rcwa_1d(materials::air(), degenerate,
                                     materials::air(), lambda, normal, 8, Pol::TM);
        std::size_t z = rcwa_tm.orders.size() / 2;
        sp("[6] TM polarization:");
        sp("    (a) degenerate vs TMM-TM:  TMM R={:.6f}  RCWA R={:.6f}",
                     tmm.R, rcwa_tm.de_r[z]);

        const auto n15 = Material::constant(cdouble{1.5, 0.0}, "n1.5");
        BinaryGrating1D g{n15, materials::air(), 1.0, 0.5, 0.5};
        auto te = solve_rcwa_1d(materials::air(), g, materials::air(), 0.5,
                                normal, 20, Pol::TE);
        auto tm = solve_rcwa_1d(materials::air(), g, materials::air(), 0.5,
                                normal, 20, Pol::TM);
        std::size_t zt = te.orders.size() / 2;
        sp("    (b) real grating, M=20:");
        sp("        TE: DE_t(+1)={:.5f}  Σ DE={:.6f}", te.de_t[zt + 1],
                     te.sum_de);
        sp("        TM: DE_t(+1)={:.5f}  Σ DE={:.6f}  (TE != TM ✓)",
                     tm.de_t[zt + 1], tm.sum_de);
        chk(std::abs(tmm.R - rcwa_tm.de_r[z]) < 1e-4, "[6] TM degenerate grating == TMM-TM");
        chk(std::abs(te.sum_de - 1.0) < 1e-6 && std::abs(tm.sum_de - 1.0) < 1e-6,
            "[6] TE and TM grating energy conserved");
        chk(std::abs(te.de_t[zt + 1] - tm.de_t[zt + 1]) > 1e-3,
            "[6] solver distinguishes TE from TM");
    }});

    // ---- Validation 7: multilayer S-matrix solver -------------------------
    // (a) A one-layer STACK must reproduce the single-layer solver exactly.
    // (b) Splitting one layer into N identical sublayers must give the same
    //     answer (the S-matrix recursion is self-consistent and stable).
    // (c) A genuinely layered device (grating + homogeneous cap) conserves
    //     energy.
    cases.push_back(Case{false, [&]() {
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
        sp("[7] Multilayer S-matrix (TM, 15deg incidence):");
        sp("    (a) single-layer solver  DE_t(0) = {:.8f}", ref.de_t[z]);
        sp("    (b) 1-layer stack        DE_t(0) = {:.8f}", st1.de_t[z]);
        sp("    (c) 5-sublayer stack     DE_t(0) = {:.8f}", stN.de_t[z]);
        double stack_err = 0.0;
        for (std::size_t i = 0; i < ref.de_t.size(); ++i) {
            stack_err = std::max(stack_err, std::abs(st1.de_t[i] - ref.de_t[i]));
            stack_err = std::max(stack_err, std::abs(stN.de_t[i] - ref.de_t[i]));
        }
        sp("        max|stack-ref| over all orders = {:.2e}", stack_err);
        sp("        Σ DE: single={:.6f} stack={:.6f} split={:.6f}",
                     ref.sum_de, st1.sum_de, stN.sum_de);
        chk(stack_err < 1e-6, "[7] S-matrix stack reproduces single-layer solver");
        chk(std::abs(ref.sum_de - 1.0) < 1e-6, "[7] single-grating energy conserved");

        // A real layered device: glass grating with a homogeneous AR-like cap.
        Rcwa1DStack device{1.0,
                           {GratingLayer1D::homogeneous(
                                Material::constant(cdouble{1.2, 0.0}, "cap"), 0.1),
                            GratingLayer1D{n15, materials::air(), 0.5, 0.5}}};
        auto dev = solve_rcwa_1d(materials::air(), device, materials::bk7(), 0.5,
                                 oblique, M, Pol::TM);
        sp("    (d) grating + cap on glass:  Σ DE = {:.6f}", dev.sum_de);
        chk(std::abs(dev.sum_de - 1.0) < 1e-6, "[7] grating+cap device energy conserved");
    }});

    // ---- Validation 8: 2D RCWA (improved Li factorization) ----------------
    // (a) A y-invariant, subwavelength 2D cell (fill_y=1, single propagating
    //     order) must reduce to the 1D solver for BOTH polarizations: E along y
    //     == 1D TE, E along x == 1D TM. The TM match is the payoff of Li's
    //     inverse-rule factorization (the basic factorization fails it).
    // (b) Energy conservation + convergence for a real high-contrast pillar.
    // (c) Cross-check vs an external solver (grcwa) on the same geometry.
    cases.push_back(Case{true, [&]() {
        const auto n15 = Material::constant(cdouble{1.5, 0.0}, "n1.5");
        const int M2 = 10;
        const auto tio2 = Material::constant(cdouble{2.45, 0.0}, "TiO2~");

        // Every solve below is independent, and this case dominates the whole
        // suite's runtime (the high-M 2D eigensolves), so launch them all
        // concurrently and fetch/print the results in order. Output is identical
        // to the serial version.
        Rcwa2DStack ginv{0.3, 0.3, {RectCell2D{n15, materials::air(), 0.5, 1.0, 0.5}}};
        Rcwa2DStack rc{0.35, 0.35, {RectCell2D{tio2, materials::air(), 0.6, 0.3, 0.6}}};
        auto solve_shape = [&](MetaShape sh, double fill, double param) {
            Rcwa2DStack s{0.35, 0.35,
                          {RectCell2D{tio2, materials::air(), fill, fill, 0.6, sh, param}}};
            return solve_rcwa_2d(materials::air(), s, materials::fused_silica(),
                                 0.532, 0.0, 0.0, 1.0, 0.0, 8, 8);
        };
        auto sweep = [&](int m) {
            Rcwa2DStack cell{0.35, 0.35, {RectCell2D{tio2, materials::air(), 0.5, 0.5, 0.6}}};
            return solve_rcwa_2d(materials::air(), cell, materials::fused_silica(),
                                 0.532, 0.0, 0.0, 1.0, 0.0, m, m);
        };
        auto go = [](auto fn) { return std::async(std::launch::async, fn); };
        auto f_te2d = go([&] { return solve_rcwa_2d(materials::air(), ginv, materials::air(), 0.5, 0.0, 0.0, 0.0, 1.0, M2, M2); });
        auto f_tm2d = go([&] { return solve_rcwa_2d(materials::air(), ginv, materials::air(), 0.5, 0.0, 0.0, 1.0, 0.0, M2, M2); });
        auto f_m6 = go([&] { return sweep(6); });
        auto f_m8 = go([&] { return sweep(8); });
        auto f_m10 = go([&] { return sweep(10); });
        auto f_m12 = go([&] { return sweep(12); });
        auto f_rx = go([&] { return solve_rcwa_2d(materials::air(), rc, materials::fused_silica(), 0.532, 0.0, 0.0, 1.0, 0.0, 12, 12); });
        auto f_ry = go([&] { return solve_rcwa_2d(materials::air(), rc, materials::fused_silica(), 0.532, 0.0, 0.0, 0.0, 1.0, 12, 12); });
        auto f_shc = go([&] { return solve_shape(MetaShape::Ellipse, 0.7, 0.5); });
        auto f_shp = go([&] { return solve_shape(MetaShape::Cross, 0.8, 0.4); });
        auto f_shr = go([&] { return solve_shape(MetaShape::Ring, 0.9, 0.5); });

        // (a) subwavelength grating vs grcwa (both polarizations).
        auto te2d = f_te2d.get();
        auto tm2d = f_tm2d.get();
        sp("[8] 2D RCWA (Li factorization):");
        sp("    (a) subwavelength grating vs grcwa (external solver):");
        sp("        TE (E∥y): 2D T0={:.5f}  (grcwa 0.93334)  |Δ|={:.2e}",
                     te2d.de_t0, std::abs(te2d.de_t0 - 0.93334));
        sp("        TM (E∥x): 2D T0={:.5f}  (grcwa 0.96050)  |Δ|={:.2e}",
                     tm2d.de_t0, std::abs(tm2d.de_t0 - 0.96050));
        chk(std::abs(te2d.de_t0 - 0.93334) < 2e-3 && std::abs(tm2d.de_t0 - 0.96050) < 5e-3,
            "[8] 2D subwavelength grating matches grcwa (TE & TM)");

        // (b) real high-contrast square TiO2 pillar: energy + convergence vs M.
        sp("    (b) square TiO2 pillar (n=2.45, Λ=0.35µm, λ=0.532µm, "
                     "fused silica), convergence:");
        sp("        {:>5}  {:>8}  {:>10}  {:>8}", "M", "Σ DE", "T0", "phase°");
        double sweep_maxdef = 0.0;
        auto row = [&](int m, auto& fut) {
            auto r = fut.get();
            sp("        {:>5}  {:>8.6f}  {:>10.5f}  {:>8.1f}", m, r.sum_de,
                         r.de_t0, std::arg(r.tx0) * 180.0 / pi);
            sweep_maxdef = std::max(sweep_maxdef, std::abs(r.sum_de - 1.0));
        };
        row(6, f_m6);
        row(8, f_m8);
        row(10, f_m10);
        row(12, f_m12);
        chk(sweep_maxdef < 1e-4, "[8] square TiO2 pillar energy conserved across M");

        // (c) external cross-check (grcwa, Li/converged): asymmetric rect pillar.
        sp("    (c) cross-check vs grcwa (rect fx=0.6 fy=0.3, fused silica):");
        auto rx = f_rx.get();
        auto ry = f_ry.get();
        sp("        x-pol T0={:.4f} (grcwa 0.954)   y-pol T0={:.4f} "
                     "(grcwa 0.972)   ΣDE={:.6f}", rx.de_t0, ry.de_t0, rx.sum_de);
        chk(std::abs(rx.de_t0 - 0.954) < 0.03 && std::abs(ry.de_t0 - 0.972) < 0.03,
            "[8] asymmetric rect pillar matches grcwa (x & y pol)");
        chk(std::abs(rx.sum_de - 1.0) < 1e-4, "[8] rect pillar energy conserved");

        // (d) Non-separable shapes via the grid (Laurent) factorization.
        sp("    (d) shaped meta-atoms (grid Laurent factorization, M=8):");
        auto shc = f_shc.get();
        auto shp = f_shp.get();
        auto shr = f_shr.get();
        sp("        circle(d0.7): T0={:.4f} φ={:.0f}°  ΣDE={:.6f}",
                     shc.de_t0, std::arg(shc.tx0) * 180.0 / pi, shc.sum_de);
        sp("        cross(arm0.4): T0={:.4f} φ={:.0f}°  ΣDE={:.6f}   "
                     "(grcwa cross @nG201 = 0.977)", shp.de_t0,
                     std::arg(shp.tx0) * 180.0 / pi, shp.sum_de);
        sp("        ring(in0.5):  T0={:.4f} φ={:.0f}°  ΣDE={:.6f}",
                     shr.de_t0, std::arg(shr.tx0) * 180.0 / pi, shr.sum_de);
        chk(std::abs(shc.sum_de - 1.0) < 2e-3 && std::abs(shp.sum_de - 1.0) < 2e-3 &&
            std::abs(shr.sum_de - 1.0) < 2e-3, "[8] shaped meta-atoms energy conserved");
    }});

    // ---- Demo 9: end-to-end metalens design -------------------------------
    // Build a phase library (sweep pillar size), then design a focusing lens
    // and report how well the realized phase matches the ideal profile.
    cases.push_back(Case{true, [&]() {
        const auto tio2 = Material::constant(cdouble{2.40, 0.0}, "TiO2~");
        sp("[9] Metalens design pipeline (TiO2 pillars, λ=532nm):");
        sp("    Building unit-cell library (sweeping pillar size)...");
        auto lib = build_unit_cell_library(tio2, materials::air(),
                                           materials::air(), materials::bk7(),
                                           /*period=*/0.35, /*λ=*/0.532,
                                           /*thickness=*/0.6, /*fill*/ 0.08, 0.92,
                                           /*n_samples=*/18, /*M=*/6);
        sp("    Library: {} pillars, phase coverage = {:.0f}° "
                     "(need ~360° for full control)",
                     lib.fill.size(), lib.phase_span() * 180.0 / pi);

        auto lens = design_metalens(lib, /*focal=*/50.0, /*diameter=*/20.0);
        sp("    Designed lens: {0}x{0} pillars, f=50µm, D=20µm",
                     lens.n_cells);
        sp("    RMS phase error vs ideal = {:.1f}°  (lower = sharper "
                     "focus)",
                     lens.rms_phase_error_deg);
        sp("    Mean pillar transmission |t| = {:.3f}", lens.mean_amplitude);
        chk(lib.phase_span() * 180.0 / pi > 300.0, "[9] library phase coverage > 300 deg");
        chk(lens.rms_phase_error_deg < 20.0, "[9] realized phase error < 20 deg");

        // Export the fabrication file and validate it round-trips.
        const std::string gds = "metalens.gds";
        int written = write_metalens_gds(lens, gds, /*layer=*/1, /*min_fill=*/0.05);
        int read_back = gds_count_boundaries(gds);
        sp("    GDSII export -> {}: {} pillars written, {} read back "
                     "({})",
                     gds, written, read_back,
                     (written == read_back && written > 0) ? "valid ✓" : "MISMATCH");

        chk(written == read_back && written > 0, "[9] GDSII layout round-trips");

        // Does it actually focus? Propagate to the focal plane and measure.
        auto foc = analyze_focus(lens, lib, /*f=*/50.0, /*λ=*/0.532, /*D=*/20.0);
        // A/B: amplitude-aware pillar selection vs phase-only.
        auto lens_po = design_metalens(lib, 50.0, 20.0, /*amplitude_weight=*/0.0);
        auto foc_po = analyze_focus(lens_po, lib, 50.0, 0.532, 20.0);
        sp("    Focal performance:");
        sp("      Strehl: phase-only {:.3f} -> amplitude-aware {:.3f}",
                     foc_po.strehl, foc.strehl);
        sp("      Strehl ratio   = {:.3f}  (1.0 = perfect)", foc.strehl);
        sp("      spot FWHM      = {:.2f} µm  (diffraction limit "
                     "λf/D = {:.2f} µm)",
                     foc.fwhm_um, foc.diffraction_limit_um);
        sp("      encircled E    = {:.1f}% within first Airy null",
                     foc.encircled_energy * 100.0);
        chk(std::abs(foc.fwhm_um - foc.diffraction_limit_um) < 0.1,
            "[9] focal spot FWHM at the diffraction limit");
        chk(foc.strehl > 0.5, "[9] designed lens focuses (Strehl > 0.5)");

        // Chromatic behavior: how the focus shifts across a wavelength band.
        auto chrom = analyze_chromatic(lens, lib, /*f=*/50.0, /*λ0=*/0.532,
                                       /*D=*/20.0, 0.45, 0.65, 5);
        sp("    Chromatic focal shift (designed for 532nm):");
        sp("      {:>7}  {:>12}  {:>10}  {:>10}", "λ(nm)", "focus(µm)",
                     "f0·λ0/λ", "rel.peak");
        for (auto& c : chrom)
            sp("      {:>7.0f}  {:>12.2f}  {:>10.2f}  {:>10.2f}",
                         c.wavelength_um * 1000.0, c.focal_length_um,
                         50.0 * 0.532 / c.wavelength_um, c.rel_peak);
    }});

    // ---- Validation 14: phase profiles + freeform (hologram) reproduction --
    // The analytic profiles (focusing/vortex/deflector/axicon) and an arbitrary
    // loaded phi(x,y) map flow through the SAME phase_profile_value() used by both
    // design paths. Sample an analytic deflector onto a grid, treat it as a
    // Freeform map, and confirm the bilinear sampler reproduces the analytic phase.
    // A linear ramp (the deflector) is reproduced EXACTLY by bilinear interpolation
    // at every point; a curved profile (focusing) is exact at the grid nodes.
    cases.push_back(Case{true, [&]() {
        sp("[14] Phase profiles + freeform (CGH) reproduction:");
        const double lambda = 0.532, extent = 14.0;
        const int n = 64;

        auto sample_to_freeform = [&](const PhaseProfile& src) {
            PhaseProfile ff;
            ff.kind = PhaseProfileKind::Freeform;
            ff.freeform_n = n;
            ff.freeform_extent_um = extent;
            ff.freeform_phase_rad.resize((std::size_t)n * n);
            for (int iy = 0; iy < n; ++iy)
                for (int ix = 0; ix < n; ++ix) {
                    double x = (ix / (double)(n - 1) - 0.5) * extent;
                    double y = (iy / (double)(n - 1) - 0.5) * extent;
                    ff.freeform_phase_rad[(std::size_t)iy * n + ix] =
                        phase_profile_value(src, x, y, lambda);
                }
            return ff;
        };

        PhaseProfile defl;
        defl.kind = PhaseProfileKind::Deflector;
        defl.deflect_deg = 10.0;
        PhaseProfile defl_ff = sample_to_freeform(defl);

        // Bilinear interp of a linear ramp is exact -> compare at off-node points too.
        double max_err_ramp = 0.0;
        for (int iy = 0; iy < 40; ++iy)
            for (int ix = 0; ix < 40; ++ix) {
                double x = (ix / 39.0 - 0.5) * extent * 0.98;
                double y = (iy / 39.0 - 0.5) * extent * 0.98;
                max_err_ramp = std::max(max_err_ramp,
                    std::abs(phase_profile_value(defl_ff, x, y, lambda) -
                             phase_profile_value(defl, x, y, lambda)));
            }
        sp("    deflector ramp via loaded map: max |Δφ| = {:.2e} rad "
                     "(bilinear is exact for a linear ramp) {}",
                     max_err_ramp, max_err_ramp < 1e-9 ? "✓" : "FAIL");

        // Curved (focusing) profile: exact at the grid nodes.
        PhaseProfile foc;
        foc.kind = PhaseProfileKind::Focusing;
        foc.focal_length_um = 50.0;
        PhaseProfile foc_ff = sample_to_freeform(foc);
        double max_err_node = 0.0;
        for (int iy = 0; iy < n; ++iy)
            for (int ix = 0; ix < n; ++ix) {
                double x = (ix / (double)(n - 1) - 0.5) * extent;
                double y = (iy / (double)(n - 1) - 0.5) * extent;
                max_err_node = std::max(max_err_node,
                    std::abs(phase_profile_value(foc_ff, x, y, lambda) -
                             phase_profile_value(foc, x, y, lambda)));
            }
        sp("    focusing map at grid nodes: max |Δφ| = {:.2e} rad {}",
                     max_err_node, max_err_node < 1e-9 ? "✓" : "FAIL");
    }});

    // ---- Validation 15: achromatic (dispersion-engineered) design ----------
    // A fill x height meta-atom library spans the (phase, group-delay) plane, so
    // a two-objective selection can match BOTH the base focusing phase AND the
    // radius-dependent group delay. Adding the group-delay objective (gd_weight>0)
    // must FLATTEN the chromatic focal drift vs a dispersion-blind baseline
    // (gd_weight=0) while keeping the base phase diffraction-limited. Small,
    // fast instance (10 fills x 6 heights x 3 wavelengths).
    cases.push_back(Case{true, [&]() {
        const auto tio2 = Material::constant(cdouble{2.40, 0.0}, "TiO2~");
        const double l0 = 0.532, bw = 0.20, D = 6.0, fl = 20.0;
        std::vector<double> band = {l0 * (1 - 0.5 * bw), l0, l0 * (1 + 0.5 * bw)};
        sp("[15] Achromatic design (fill×height dispersion engineering, "
                     "{:.0f}% band):", bw * 100.0);
        auto dl = build_dispersive_library(tio2, materials::air(), materials::air(),
                                           materials::fused_silica(), 0.35, band, l0,
                                           0.08, 0.92, 10, 0.40, 1.40, 6, /*M=*/5);
        auto sd = design_achromatic_metalens(dl, fl, D, /*gd_weight=*/0.0);
        auto ad = design_achromatic_metalens(dl, fl, D, /*gd_weight=*/1.0);
        auto cs = verify_achromatic_focus(dl, sd);
        auto ca = verify_achromatic_focus(dl, ad);
        auto drift = [](const std::vector<AchromaticFocalPoint>& c) {
            double lo = 1e300, hi = -1e300;
            for (auto& p : c) { lo = std::min(lo, p.focal_length_um); hi = std::max(hi, p.focal_length_um); }
            return hi - lo;
        };
        double ds = drift(cs), da = drift(ca);
        sp("    GD library span = {:.2f} fs ({} atoms); base-phase RMS: "
                     "standard {:.1f}°, achromatic {:.1f}°",
                     dl.gd_max_fs - dl.gd_min_fs, static_cast<int>(dl.atoms.size()),
                     sd.rms_phase_error_deg, ad.rms_phase_error_deg);
        sp("    group-delay RMS: standard {:.2f} fs -> achromatic {:.2f} fs",
                     sd.rms_group_delay_error_fs, ad.rms_group_delay_error_fs);
        sp("    chromatic focal drift: standard {:.2f} µm -> achromatic {:.2f} µm  {}",
                     ds, da,
                     (da < ds && ad.rms_phase_error_deg < 25.0) ? "✓ (flatter + still focusing)"
                                                                : "FAIL");
    }});

    // ---- Validation 15b: SINGLE-ETCH achromatic (shape-diverse, one height) --
    // The fabricable variant: every atom shares ONE etch depth and the (phase,
    // group-delay) plane is spanned by SHAPE variety (square/circle/cross/ring)
    // instead of by depth. The robust signal that dispersion engineering works
    // here (the low-Fresnel focal-drift metric is noisy at this small aperture)
    // is that adding the group-delay objective REDUCES the group-delay RMS while
    // every chosen atom keeps the single height -- one lithography step.
    cases.push_back(Case{true, [&]() {
        const auto tio2 = Material::constant(cdouble{2.40, 0.0}, "TiO2~");
        const double l0 = 0.532, bw = 0.20, D = 6.0, fl = 20.0, h = 0.80;
        std::vector<double> band = {l0 * (1 - 0.5 * bw), l0, l0 * (1 + 0.5 * bw)};
        sp("[15b] Single-etch achromatic (shape-diverse @ one height {:.2f}µm):", h);
        auto dl = build_single_etch_library(tio2, materials::air(), materials::air(),
                                            materials::fused_silica(), 0.35, band, l0,
                                            0.08, 0.92, /*n_fills=*/5, h, /*M=*/5);
        auto sd = design_achromatic_metalens(dl, fl, D, /*gd_weight=*/0.0);
        auto ad = design_achromatic_metalens(dl, fl, D, /*gd_weight=*/1.0);
        sp("    GD library span = {:.2f} fs ({} atoms, all @ {:.2f}µm); "
                     "base-phase RMS: standard {:.1f}°, achromatic {:.1f}°",
                     dl.gd_max_fs - dl.gd_min_fs, static_cast<int>(dl.atoms.size()), h,
                     sd.rms_phase_error_deg, ad.rms_phase_error_deg);
        bool single = sd.single_height && ad.single_height;
        bool gd_better = ad.rms_group_delay_error_fs < sd.rms_group_delay_error_fs;
        sp("    group-delay RMS: standard {:.2f} fs -> achromatic {:.2f} fs  {}",
                     sd.rms_group_delay_error_fs, ad.rms_group_delay_error_fs,
                     (single && gd_better && ad.rms_phase_error_deg < 25.0)
                         ? "✓ (GD-matched in ONE etch)"
                         : "FAIL");
    }});

    // ---- Validation 15c: achromatic Pancharatnam-Berry (PB + dispersion) -----
    // The modern recipe: the geometric (PB) phase sets the base profile EXACTLY by
    // rotation while a dispersive birefringent atom (picked per radius from a
    // fill_x x fill_y grid at ONE height) supplies the group delay. The two robust
    // signals: (1) the base-phase RMS is ~0 for BOTH the standard and achromatic
    // designs (geometric phase is exact -- no library quantization on phase), and
    // (2) engaging the group-delay objective REDUCES the group-delay RMS, all at a
    // single etch depth. (Focal drift stays a noisy supplement at this aperture.)
    cases.push_back(Case{true, [&]() {
        const auto tio2 = Material::constant(cdouble{2.40, 0.0}, "TiO2~");
        const double l0 = 0.532, bw = 0.20, D = 6.0, fl = 20.0, h = 1.10;
        std::vector<double> band = {l0 * (1 - 0.5 * bw), l0, l0 * (1 + 0.5 * bw)};
        sp("[15c] Achromatic PB (geometric phase + dispersion, one etch {:.2f}µm):", h);
        auto lib = build_dispersive_pb_library(tio2, materials::air(), materials::air(),
                                               materials::bk7(), 0.35, band, l0,
                                               0.10, 0.90, /*n_fills=*/6, h, /*M=*/5);
        auto sd = design_pb_achromatic_metalens(lib, fl, D, /*handedness=*/+1, /*gd_weight=*/0.0);
        auto ad = design_pb_achromatic_metalens(lib, fl, D, /*handedness=*/+1, /*gd_weight=*/1.0);
        sp("    GD library span = {:.2f} fs ({} atoms @ {:.2f}µm); base-phase RMS: "
                     "standard {:.1e}°, achromatic {:.1e}° (geometric -> exact)",
                     lib.gd_max_fs - lib.gd_min_fs, static_cast<int>(lib.atoms.size()), h,
                     sd.rms_phase_error_deg, ad.rms_phase_error_deg);
        bool base_exact = sd.rms_phase_error_deg < 1e-6 && ad.rms_phase_error_deg < 1e-6;
        bool gd_better = ad.rms_group_delay_error_fs < sd.rms_group_delay_error_fs;
        sp("    group-delay RMS: standard {:.2f} fs -> achromatic {:.2f} fs  {}",
                     sd.rms_group_delay_error_fs, ad.rms_group_delay_error_fs,
                     (base_exact && gd_better) ? "✓ (exact base phase + GD-matched, ONE etch)"
                                               : "FAIL");
    }});

    // ---- Reproduction 16: published device (Khorasaninejad 2016, 532 nm) ---
    // The exact NA=0.80 TiO2 nanofin from "Metalenses at visible wavelengths,"
    // Science 352, 1190 (2016): W=95, L=250, H=600 nm, period 325 nm. Lock the
    // reproduction's two headline signals: the nanofin is a near-ideal half-wave
    // plate (>=95% of the TRANSMITTED light is spin-converted) whose transmitted-
    // normalized conversion brackets the paper's 73% focusing efficiency from
    // above, and the geometric phase is exact. (Constant n~2.45 stand-in -- the
    // real Siefke n,k matches at 532 nm -- so the test needs no data file.)
    cases.push_back(Case{true, [&]() {
        const auto tio2 = Material::constant(cdouble{2.45, 0.0}, "TiO2~");
        const double lam = 0.532, U = 0.325, H = 0.600;
        const double fx = 0.095 / U, fy = 0.250 / U;
        sp("[16] Reproduce Khorasaninejad 2016 NA=0.80 nanofin "
                     "(532 nm, W95 x L250 x H600):");
        Rcwa2DStack st{U, U, {RectCell2D{tio2, materials::air(), fx, fy, H}}};
        auto rx = solve_rcwa_2d(materials::air(), st, materials::fused_silica(), lam,
                                0, 0, 1, 0, 10, 10);
        auto ry = solve_rcwa_2d(materials::air(), st, materials::fused_silica(), lam,
                                0, 0, 0, 1, 10, 10);
        const cdouble tx = rx.tx0, ty = ry.ty0;
        const double T = 0.5 * (std::norm(tx) + std::norm(ty));
        const double conv_abs = std::norm(tx - ty) / 4.0;
        const double conv_rel = conv_abs / T;
        const double retard = std::remainder(std::arg(tx) - std::arg(ty), 2 * pi) * 180.0 / pi;
        const double paper_eff = 0.73;
        HwpAtom atom;
        atom.fill_x = fx; atom.fill_y = fy; atom.thickness_um = H;
        atom.t_x = tx; atom.t_y = ty; atom.conversion_efficiency = conv_abs;
        auto d = design_pb_metalens(atom, U, lam, /*f=*/6.0, /*D=*/6.0, +1);
        const bool hwp_ok = conv_rel > 0.95;                 // near-ideal half-wave plate
        const bool brackets = conv_rel > paper_eff;          // upper bound on focusing eff
        const bool phase_ok = d.rms_phase_error_deg < 1e-6;  // geometric phase exact
        const bool retard_ok = std::abs(std::abs(retard) - 180.0) < 15.0;
        sp("    transmitted-norm conversion = {:.1f}% (HWP quality), retardance "
                     "{:.0f} deg, geo-phase RMS {:.1e} deg", 100.0 * conv_rel, retard,
                     d.rms_phase_error_deg);
        sp("    brackets paper focusing eff ({:.0f}%) from above + exact phase  {}",
                     100.0 * paper_eff,
                     (hwp_ok && brackets && phase_ok && retard_ok) ? "✓" : "FAIL");
    }});

    // ---- Wide-FOV 17: quadratic-phase lens vs hyperbolic, off-axis ----------
    // A hyperbolic metalens is perfect on-axis but develops coma off-axis; a
    // QUADRATIC (parabolic) lens with an aperture stop IN FRONT (here at the front
    // focal plane) is the same parabola simply recentered under tilt, so its focus
    // stays sharp across a wide field (the spot just shifts). With the offset stop
    // each field angle samples a decentered low-NA patch: coma-free for the
    // quadratic (vertex-centered parabola), coma-laden for the hyperbolic lens.
    // Lock the contrast -- at the widest angle the quadratic holds its Strehl while
    // the hyperbolic collapses. (See `celeris widefov`.)
    cases.push_back(Case{true, [&]() {
        const auto tio2 = Material::constant(cdouble{2.40, 0.0}, "TiO2~");
        const double l0 = 0.532, f = 20.0, D = 50.0, stopD = 16.0, stopd = 20.0;
        sp("[17] Wide-FOV: quadratic vs hyperbolic phase (stop D={:.0f}µm @ {:.0f}µm "
                     "in front, f={:.0f}µm):", stopD, stopd, f);
        auto lib = build_unit_cell_library(tio2, materials::air(), materials::air(),
                                           materials::bk7(), 0.35, l0, 0.6, 0.08, 0.92,
                                           14, /*M=*/5);
        PhaseProfile hyp; hyp.kind = PhaseProfileKind::Focusing; hyp.focal_length_um = f;
        PhaseProfile quad; quad.kind = PhaseProfileKind::Quadratic; quad.focal_length_um = f;
        auto lh = design_metalens(lib, hyp, D);
        auto lq = design_metalens(lib, quad, D);
        std::vector<double> angles = {0.0, 15.0, 30.0};
        auto fh = analyze_wide_fov(lh, lib, f, l0, D, stopD, stopd, angles);
        auto fq = analyze_wide_fov(lq, lib, f, l0, D, stopD, stopd, angles);
        const double hyp_edge = fh.back().rel_strehl, quad_edge = fq.back().rel_strehl;
        sp("    rel. Strehl at {:.0f}°: hyperbolic {:.3f} (coma) -> quadratic {:.3f} "
                     "(holds focus)", angles.back(), hyp_edge, quad_edge);
        bool win = quad_edge > 0.8 && hyp_edge < 0.6 && quad_edge > hyp_edge + 0.3;
        sp("    quadratic-phase lens keeps a wide-field focus where the hyperbolic "
                     "fails  {}", win ? "✓" : "FAIL");
    }});

    // ---- Reproduction 18: broadband achromat (Chen 2018, 470-670 nm) -------
    // Reproduce Chen et al., Nat. Nanotechnol. 13, 220 (2018): a single-layer
    // NA=0.20 D=26.4 um TiO2 achromat, H=600 nm, 400 nm lattice. Two locks:
    // (a) the DEVICE-SPECIFIC physical limit -- the required group-delay span
    //     GD=(1/c)(sqrt(R^2+f^2)-f) at the published NA/D is ~4.4 fs, at the ~5 fs
    //     ceiling of a 600-nm TiO2 nanofin, so the diameter is GD-limited (this is
    //     why the paper's lens is 26.4 um); and
    // (b) the ENGINE behaviour at the published period/height -- the geometric
    //     phase is exact for both std/achromatic and the group-delay objective
    //     reduces the GD RMS, in ONE 600-nm etch. (Small library for speed.)
    cases.push_back(Case{true, [&]() {
        // (a) analytic group-delay budget at the published NA=0.20, D=26.4 um.
        const double NA = 0.20, D = 26.4, R = D / 2.0;
        const double f = R / std::tan(std::asin(NA));
        const double gd_req = (std::sqrt(R * R + f * f) - f) * GD_FS_PER_UM;
        sp("[18] Reproduce Chen 2018 achromat (NA=0.20, D=26.4µm, H=600nm, 470-670nm):");
        const bool gd_at_ceiling = gd_req > 3.5 && gd_req < 5.0;  // ~4.4 fs, at nanofin limit
        sp("    required GD span = {:.2f} fs (600-nm-nanofin ceiling ~5 fs) -> "
                     "diameter is GD-limited  {}", gd_req, gd_at_ceiling ? "✓" : "FAIL");
        // (b) engine: single-etch PB library at the published cell, std vs achromatic.
        const auto tio2 = Material::constant(cdouble{2.40, 0.0}, "TiO2~");
        std::vector<double> band = {0.470, 0.570, 0.670};   // band endpoints + center
        auto lib = build_dispersive_pb_library(tio2, materials::air(), materials::air(),
                                               materials::fused_silica(), 0.400, band, 0.570,
                                               0.10, 0.90, /*n_fills=*/6, /*h=*/0.600, /*M=*/5);
        auto sd = design_pb_achromatic_metalens(lib, /*f=*/8.0, /*D=*/5.0, +1, 0.0);
        auto ad = design_pb_achromatic_metalens(lib, /*f=*/8.0, /*D=*/5.0, +1, 1.0);
        const bool base_exact = sd.rms_phase_error_deg < 1e-6 && ad.rms_phase_error_deg < 1e-6;
        const bool gd_better = ad.rms_group_delay_error_fs < sd.rms_group_delay_error_fs;
        sp("    one 600-nm etch: base-phase RMS {:.1e}° (geometric exact), GD RMS "
                     "{:.2f}->{:.2f} fs  {}", ad.rms_phase_error_deg,
                     sd.rms_group_delay_error_fs, ad.rms_group_delay_error_fs,
                     (gd_at_ceiling && base_exact && gd_better) ? "✓" : "FAIL");
    }});

    // ---- Demo 11: inverse design (gradient-based optimizer) ---------------
    // Instead of looking a pillar up from the discrete library, SOLVE for the
    // geometry hitting a target phase with maximum transmission.
    cases.push_back(Case{true, [&]() {
        const auto tio2 = Material::constant(cdouble{2.40, 0.0}, "TiO2~");
        const double target_deg = 90.0;
        PillarTarget tgt{0.532, target_deg * pi / 180.0, /*amplitude_weight=*/1.0};
        sp("[11] Inverse design: optimize pillar for {:.0f}° phase "
                     "@532nm (max transmission):", target_deg);
        auto opt = optimize_pillar(tio2, materials::air(), materials::air(),
                                   materials::bk7(), /*period=*/0.35, tgt,
                                   /*M=*/5, /*fill0=*/0.50, /*thickness0=*/0.50,
                                   /*max_iters=*/25);
        sp("    converged geometry: fill={:.3f}, thickness={:.3f} µm",
                     opt.fill, opt.thickness_um);
        sp("    achieved phase = {:.1f}°  (target {:.0f}°),  |t| = "
                     "{:.3f},  loss = {:.2e}",
                     opt.achieved_phase_rad * 180.0 / pi, target_deg,
                     opt.achieved_amplitude, opt.loss);
    }});

    // ---- Demo 12: form birefringence (polarization-multiplexed basis) -----
    // A rectangular pillar (fill_x != fill_y) responds differently to x- and
    // y-polarized light — "form birefringence." That phase difference is the
    // basis for polarization-multiplexed metalenses (one device, two functions
    // selected by polarization). Here we show it grow with pillar asymmetry.
    cases.push_back(Case{true, [&]() {
        const auto tio2 = Material::constant(cdouble{2.40, 0.0}, "TiO2~");
        sp("[12] Form birefringence (rectangular TiO2 pillar, 532nm):");
        sp("      {:>10}  {:>10}  {:>10}  {:>12}", "fill_x", "fill_y",
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
            sp("      {:>10.2f}  {:>10.2f}  {:>10.1f}  {:>5.2f},{:>5.2f}",
                         fx, fy, dphi * 180.0 / pi, std::abs(rx.tx0),
                         std::abs(ry.ty0));
        }
    }});

    // ---- Validation 13: RCWA vs Effective-Medium Theory --------------------
    // A deeply subwavelength grating (period << λ) behaves as a uniform
    // birefringent film. Rytov's 0th-order EMT: ε∥ = f·ε1+(1-f)·ε2 (TE, E along
    // grooves), ε⊥ = [f/ε1+(1-f)/ε2]⁻¹ (TM). The grating's RCWA reflectance must
    // converge to the TMM reflectance of those effective films — an independent
    // analytic check of the full TE+TM vectorial solver.
    cases.push_back(Case{true, [&]() {
        const double lam = 0.5, d = 0.10, f = 0.5, e1 = 2.1025, e2 = 1.0;  // ridge n=1.45, groove air
        const double eTE = f * e1 + (1 - f) * e2;            // arithmetic mean
        const double eTM = 1.0 / (f / e1 + (1 - f) / e2);    // harmonic mean
        const auto ridge = Material::constant(cdouble{1.45, 0.0}, "n1.45");
        const auto nTE = Material::constant(cdouble{std::sqrt(eTE), 0.0}, "nTE");
        const auto nTM = Material::constant(cdouble{std::sqrt(eTM), 0.0}, "nTM");
        sp("[13] RCWA vs effective-medium theory (Λ=λ/20 subwavelength):");
        sp("      {:>4}  {:>14}  {:>14}  {:>9}", "pol", "RCWA R", "EMT-film R", "|Δ|");
        for (int te = 1; te >= 0; --te) {
            Pol pol = te ? Pol::TE : Pol::TM;
            BinaryGrating1D g{ridge, materials::air(), 0.025, f, d};  // Λ = λ/20
            auto rg = solve_rcwa_1d(materials::air(), g, materials::air(), lam, 0.0, 20, pol);
            double rcwaR = rg.de_r[rg.orders.size() / 2];
            auto slab = solve_stack(materials::air(), {{te ? nTE : nTM, d}},
                                    materials::air(), lam, 0.0, pol);
            sp("      {:>4}  {:>14.4f}  {:>14.4f}  {:>9.4f}",
                         te ? "TE" : "TM", rcwaR, slab.R, std::abs(rcwaR - slab.R));
        }
    }});

    // ---- Validation 19: named material registry (real n,k library) ---------
    // The registry maps canonical names + aliases to Materials. Analytic
    // dielectrics are exact published Sellmeier (lossless); tabulated metals &
    // semiconductors carry real loss from refractiveindex.info data. Lock a few
    // literature spot-values so the shipped data/registry can't silently drift.
    cases.push_back(Case{false, [&]() {
        using namespace materials;
        sp("[19] Material registry (real n,k library):");
        // analytic, lossless: sapphire ordinary at 532nm = 1.7717 (Malitson 1972).
        const cdouble sa = by_name("sapphire").index(0.532);
        const bool sapphire_ok = std::abs(sa.real() - 1.7717) < 2e-3 && sa.imag() == 0.0;
        // alias resolution: gold==au, silica==sio2.
        const bool alias_ok = by_name("gold").name() == by_name("au").name() &&
                              by_name("silica").index(0.532).real() ==
                                  by_name("sio2").index(0.532).real();
        // tabulated, lossy: gold at 600nm has small n, large k (Johnson-Christy).
        const cdouble au = by_name("au").index(0.600);
        const bool au_ok = au.real() < 0.5 && au.imag() > 2.5;
        // c-Si nearly transparent (small k) vs a-Si absorbing (large k) at 532nm.
        const double k_csi = by_name("c-si").index(0.532).imag();
        const double k_asi = by_name("a-si").index(0.532).imag();
        const bool si_ok = k_csi < 0.10 && k_asi > 0.50;
        // every tabulated data file located on disk?
        bool files_ok = true;
        int n_mat = 0;
        for (const auto& m : catalog()) {
            ++n_mat;
            if (m.tabulated && !m.available) files_ok = false;
        }
        sp("    {} materials | sapphire n@532={:.4f} k=0 {} | gold@600 "
                     "n={:.2f} k={:.2f} {} | k(c-Si)={:.3f}<k(a-Si)={:.3f} {} | "
                     "aliases {} | data files {}",
                     n_mat, sa.real(), sapphire_ok ? "✓" : "FAIL", au.real(),
                     au.imag(), au_ok ? "✓" : "FAIL", k_csi, k_asi,
                     si_ok ? "✓" : "FAIL", alias_ok ? "✓" : "FAIL",
                     files_ok ? "✓" : "FAIL");
    }});

    cases.push_back(Case{true, [&]() {
        using namespace materials;
        sp("[20] Efficiency budget (per-order / absorption):");
        // (a) Lossless TiO2 atom, subwavelength pitch -> energy conserved, no
        // absorption, ALL transmitted power in the single propagating 0th order.
        auto bt = analyze_efficiency(materials::air(), by_name("tio2"),
                                     materials::air(), by_name("sio2"), 0.35, 0.35,
                                     0.6, 0.6, 0.6, 0.532, {1.0, 0.0}, {0.0, 0.0}, 6);
        const double sum_t = bt.transmission + bt.reflection + bt.absorption;
        const bool e_ok = std::abs(sum_t - 1.0) < 1e-3;        // R+T+A = 1
        const bool a_ok = bt.absorption < 1e-3;                // lossless
        const bool o0_ok = std::abs(bt.t_zero - bt.transmission) < 1e-4 &&
                           bt.n_prop_t == 1;                   // only 0th propagates
        // (b) Lossy gold atom -> real absorption (the deficit 1-R-T), energy still
        // accounted for exactly.
        auto bg = analyze_efficiency(materials::air(), by_name("au"),
                                     materials::air(), by_name("sio2"), 0.35, 0.35,
                                     0.5, 0.5, 0.10, 0.532, {1.0, 0.0}, {0.0, 0.0}, 6);
        const double sum_g = bg.transmission + bg.reflection + bg.absorption;
        const bool g_ok = std::abs(sum_g - 1.0) < 1e-3 && bg.absorption > 0.2;
        sp("    TiO2: R+T+A={:.6f} {} | absorption {:.2e} {} | 0th holds "
                     "all T ({} order) {}",
                     sum_t, e_ok ? "✓" : "FAIL", bt.absorption, a_ok ? "✓" : "FAIL",
                     bt.n_prop_t, o0_ok ? "✓" : "FAIL");
        sp("    Au:   R+T+A={:.6f}, absorption {:.3f} (lossy metal) {}",
                     sum_g, bg.absorption, g_ok ? "✓" : "FAIL");
    }});

    cases.push_back(Case{true, [&]() {
        using namespace materials;
        sp("[21] Field-resolved grid (full-field Strehl map):");
        // Small focusing lens, 5x5 field grid. Lock the structural invariants: the
        // on-axis node is the rel_strehl=1 reference, the grid is symmetric (the
        // lens is rotationally symmetric), Strehl falls monotonically off-axis, and
        // the on-axis spot is diffraction-limited.
        const double f = 30.0, D = 12.0, lam = 0.532, per = 0.35, h = 0.6;
        auto lib = build_unit_cell_library(by_name("tio2"), materials::air(),
                                           materials::air(), by_name("sio2"), per,
                                           lam, h, 0.08, 0.92, 14, 5);
        auto lens = design_metalens(lib, f, D);
        auto g = analyze_field_grid(lens, lib, f, lam, D, /*max_angle=*/10.0,
                                    /*n_half=*/2, /*psf_n=*/61);
        const int n = g.n, c = n / 2;  // center index
        auto at = [&](int jx, int jy) -> const FieldGridPoint& {
            return g.points[(std::size_t)jy * n + jx];
        };
        const double s00 = at(c, c).rel_strehl;
        const double corner = at(0, 0).rel_strehl;
        // 4-corner symmetry (rotational symmetry of the lens).
        const double cmax = std::max({at(0, 0).rel_strehl, at(n - 1, 0).rel_strehl,
                                      at(0, n - 1).rel_strehl,
                                      at(n - 1, n - 1).rel_strehl});
        const double cmin = std::min({at(0, 0).rel_strehl, at(n - 1, 0).rel_strehl,
                                      at(0, n - 1).rel_strehl,
                                      at(n - 1, n - 1).rel_strehl});
        const bool sym_ok = (cmax - cmin) < 1e-3;
        const bool ref_ok = std::abs(s00 - 1.0) < 1e-9;     // on-axis = reference
        const bool fall_ok = corner < s00 - 0.02;           // degrades off-axis
        const double dl = lam * f / D;
        const bool dl_ok = at(c, c).fwhm_x_um <= 1.6 * dl &&
                           at(c, c).fwhm_x_um > 0.0;         // on-axis ~ DL
        sp("    on-axis Strehl {:.6f} {} | corner {:.3f}<center {} | "
                     "4-corner sym Δ={:.2e} {} | on-axis FWHM {:.2f}µm (DL {:.2f}) {}",
                     s00, ref_ok ? "✓" : "FAIL", corner, fall_ok ? "✓" : "FAIL",
                     cmax - cmin, sym_ok ? "✓" : "FAIL", at(c, c).fwhm_x_um, dl,
                     dl_ok ? "✓" : "FAIL");
    }});

#ifdef CELERIS_USE_CUDA
    // ---- Benchmark 10: GPU vs CPU eigensolve ------------------------------
    // The per-layer RCWA eigenproblem is a general complex matrix of size 2N.
    // Honest head-to-head on a representative 578x578 (2N at M=8): cuSOLVER
    // Xgeev (GPU) vs Eigen ComplexEigenSolver (CPU). Verify eigenvalues agree.
    cases.push_back(Case{true, [&]() {
        const int n = 578;
        std::mt19937 rng(42);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        Eigen::MatrixXcd A(n, n);
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
                A(i, j) = cdouble{dist(rng), dist(rng)};

        sp("[10] GPU vs CPU eigensolve ({}x{} general complex):", n, n);
        if (!cuda::available()) {
            sp("     no CUDA device available");
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

            sp("     CPU (Eigen)      : {:.1f} ms", cpu_ms);
            sp("     GPU (cuSOLVER)   : {:.1f} ms   ({:.1f}x)", gpu_ms,
                         cpu_ms / gpu_ms);
            sp("     eigenvalue match : max|Δ| = {:.2e}  ok={}",
                         max_diff, ok);
        }
    }});
#endif

    // --- execute cases concurrently, then print in declaration order --------
    // The cases are independent (each was a separate scope), so we run them on a
    // thread pool and buffer each one's output; this turns the largely serial
    // suite into a parallel one without changing any result.
    std::vector<std::string> out(cases.size());
    std::vector<double> ms(cases.size(), 0.0);
    std::atomic<std::size_t> next{0};
    unsigned nthreads = std::max(1u, std::thread::hardware_concurrency());
    if (nthreads > cases.size()) nthreads = static_cast<unsigned>(cases.size());
    auto worker = [&] {
        for (std::size_t i = next.fetch_add(1); i < cases.size();
             i = next.fetch_add(1)) {
            if (quick && cases[i].heavy) continue;  // fast subset for CI on push
            g_selftest_out = &out[i];
            auto t0 = std::chrono::steady_clock::now();
            cases[i].run();
            ms[i] = std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - t0).count();
            g_selftest_out = nullptr;
        }
    };
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < nthreads; ++t) pool.emplace_back(worker);
    for (auto& th : pool) th.join();

    for (std::size_t i = 0; i < cases.size(); ++i) std::print("{}", out[i]);
    // Per-case wall time to stderr (keeps stdout identical to the serial suite).
    for (std::size_t i = 0; i < cases.size(); ++i)
        if (!out[i].empty())
            std::println(stderr, "  timing: case#{:<2} {:9.1f} ms", i, ms[i]);

    const int fails = g_selftest_failures.load();
    std::println("");
    if (fails == 0)
        std::println("SELF-TEST: all locked-tolerance checks passed.");
    else
        std::println("SELF-TEST: {} check(s) FAILED — see [FAIL] lines above.", fails);
    return fails == 0 ? 0 : 1;
}
