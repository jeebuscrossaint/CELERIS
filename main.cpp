#include <cmath>
#include <print>
#include <vector>

#include "celeris/materials/database.hpp"
#include "celeris/optics/tmm.hpp"
#include "celeris/rcwa/rcwa1d.hpp"

using namespace celeris;

int main() {
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
}
