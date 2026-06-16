#include "celeris/rcwa/rcwa1d.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

namespace celeris {

Rcwa1DResult solve_rcwa_1d(const Material& incident,
                           const BinaryGrating1D& grating,
                           const Material& substrate,
                           double wavelength_um,
                           double theta0_rad,
                           int M,
                           Pol pol) {
    if (M < 0) throw std::invalid_argument("solve_rcwa_1d: M must be >= 0");

    using Eigen::MatrixXcd;
    using Eigen::VectorXcd;
    using Eigen::ComplexEigenSolver;

    const int n = 2 * M + 1;          // number of retained orders
    const double k0 = 2.0 * pi / wavelength_um;

    const cdouble n_inc = incident.index(wavelength_um);
    const cdouble n_sub = substrate.index(wavelength_um);

    // --- Order wavevectors (normalized by k0) ------------------------------
    // Floquet condition: kx_m/k0 = n_inc·sinθ − m·(λ/Λ).
    VectorXcd kx(n);
    for (int m = -M; m <= M; ++m) {
        kx(m + M) = n_inc * std::sin(theta0_rad) -
                    static_cast<double>(m) * (wavelength_um / grating.period_um);
    }

    // Longitudinal (z) wavevectors in each half-space, normalized by k0.
    // std::sqrt's principal branch gives Re>=0 (propagating) or Im>0
    // (evanescent, decaying away from the grating) — exactly what we want.
    auto kz_of = [&](cdouble nreg, int idx) -> cdouble {
        return std::sqrt(nreg * nreg - kx(idx) * kx(idx));
    };
    VectorXcd kzI(n), kzII(n);  // longitudinal wavevectors kz/k0 in each region
    for (int idx = 0; idx < n; ++idx) {
        kzI(idx) = kz_of(n_inc, idx);
        kzII(idx) = kz_of(n_sub, idx);
    }

    // --- Permittivity Toeplitz matrix E_{a,b} = ε_{(a−b)} ------------------
    MatrixXcd E(n, n);
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b)
            E(a, b) = grating.eps_fourier((a - M) - (b - M), wavelength_um);

    const MatrixXcd Kx = kx.asDiagonal();
    const MatrixXcd Id = MatrixXcd::Identity(n, n);

    // --- Layer eigenproblem  d²ψ/dz'² = A ψ  and region "admittances" ------
    // Both polarizations share the same boundary-matching structure; only the
    // eigen-operator A, the companion modes V, and the region admittance Y
    // (relating the continuous companion field to the order amplitudes) differ.
    MatrixXcd A;
    MatrixXcd Vmodes;
    VectorXcd Yreg_I(n), Yreg_II(n);

    ComplexEigenSolver<MatrixXcd> ces;
    if (pol == Pol::TE) {
        // E_y is primary; companion H_x. A = Kx² − E,  V = W·Q,  Y = j·kz/k0.
        A = Kx * Kx - E;
        ces.compute(A);
    } else {
        // TM: H_y primary; companion E_x. Li's rule -> middle uses Toeplitz of
        // 1/ε; companion uses the matrix inverse of [ε].  A = E·(Kx·Einv·Kx−I).
        MatrixXcd Einv(n, n);  // Toeplitz of reciprocal permittivity 1/ε
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b)
                Einv(a, b) =
                    grating.eps_inv_fourier((a - M) - (b - M), wavelength_um);
        A = E * (Kx * Einv * Kx - Id);
        ces.compute(A);
    }

    VectorXcd q = ces.eigenvalues().cwiseSqrt();  // q = kz_layer/k0 per mode
    for (int i = 0; i < n; ++i)
        if (q(i).real() < 0) q(i) = -q(i);  // decaying / outgoing branch
    MatrixXcd W = ces.eigenvectors();

    const cdouble j_unit{0.0, 1.0};
    if (pol == Pol::TE) {
        Vmodes = W * q.asDiagonal();
        Yreg_I = j_unit * kzI;
        Yreg_II = j_unit * kzII;
    } else {
        Vmodes = j_unit * (E.inverse() * W) * q.asDiagonal();
        Yreg_I = kzI.array() / (n_inc * n_inc);
        Yreg_II = kzII.array() / (n_sub * n_sub);
    }
    const MatrixXcd& V = Vmodes;

    // Stable exponential of the layer: X = exp(−k0·q·d), |entries| <= 1.
    VectorXcd xdiag = (-k0 * grating.thickness_um * q.array()).exp();
    MatrixXcd X = xdiag.asDiagonal();

    // --- Boundary-condition linear system ----------------------------------
    // Unknowns u = [R; T; a; b], each length n. Field/derivative continuity at
    // z=0 and z=d (enhanced-transmittance referencing keeps it numerically
    // stable). See header reference.
    const MatrixXcd YI = Yreg_I.asDiagonal();
    const MatrixXcd YII = Yreg_II.asDiagonal();
    const MatrixXcd WX = W * X;
    const MatrixXcd VX = V * X;

    VectorXcd delta = VectorXcd::Zero(n);
    delta(M) = 1.0;  // unit-amplitude incident wave in order 0

    MatrixXcd S = MatrixXcd::Zero(4 * n, 4 * n);
    VectorXcd rhs = VectorXcd::Zero(4 * n);

    // z=0, field:        δ + R = W a + WX b
    S.block(0, 0, n, n) = Id;
    S.block(0, 2 * n, n, n) = -W;
    S.block(0, 3 * n, n, n) = -WX;
    rhs.segment(0, n) = -delta;
    // z=0, companion:    Y_I(δ − R) = V a − VX b
    S.block(n, 0, n, n) = -YI;
    S.block(n, 2 * n, n, n) = -V;
    S.block(n, 3 * n, n, n) = VX;
    rhs.segment(n, n) = -YI * delta;
    // z=d, field:        WX a + W b = T
    S.block(2 * n, n, n, n) = -Id;
    S.block(2 * n, 2 * n, n, n) = WX;
    S.block(2 * n, 3 * n, n, n) = W;
    // z=d, companion:    VX a − V b = Y_II T
    S.block(3 * n, n, n, n) = -YII;
    S.block(3 * n, 2 * n, n, n) = VX;
    S.block(3 * n, 3 * n, n, n) = -V;

    VectorXcd sol = S.partialPivLu().solve(rhs);
    VectorXcd R = sol.segment(0, n);
    VectorXcd T = sol.segment(n, n);

    // --- Diffraction efficiencies ------------------------------------------
    // Each order carries z-directed power ∝ Re(kz/ε_flux), where ε_flux = 1 for
    // TE and n² for TM. Normalize by the incident order's flux. In region I the
    // ε_flux cancels (same medium), so reflected DE is identical in form for
    // both polarizations; transmitted DE picks up the ε_inc/ε_sub ratio for TM.
    const double kz_inc = (n_inc * std::cos(cdouble{theta0_rad, 0.0})).real();
    const cdouble eps_inc_flux = (pol == Pol::TE) ? cdouble{1.0, 0.0} : n_inc * n_inc;
    const cdouble eps_sub_flux = (pol == Pol::TE) ? cdouble{1.0, 0.0} : n_sub * n_sub;

    Rcwa1DResult out;
    out.orders.reserve(n);
    out.de_r.reserve(n);
    out.de_t.reserve(n);
    double total = 0.0;
    for (int idx = 0; idx < n; ++idx) {
        double der = std::norm(R(idx)) * kzI(idx).real() / kz_inc;
        double det = std::norm(T(idx)) *
                     (kzII(idx) / eps_sub_flux).real() *
                     eps_inc_flux.real() / kz_inc;
        out.orders.push_back(idx - M);
        out.de_r.push_back(der);
        out.de_t.push_back(det);
        total += der + det;
    }
    out.sum_de = total;
    return out;
}

} // namespace celeris
