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
    if (pol != Pol::TE) {
        throw std::invalid_argument(
            "solve_rcwa_1d: only TE is implemented so far (TM needs Li's "
            "inverse-rule factorization — next milestone)");
    }
    if (M < 0) throw std::invalid_argument("solve_rcwa_1d: M must be >= 0");

    using Eigen::MatrixXcd;
    using Eigen::VectorXcd;
    using Eigen::ComplexEigenSolver;

    const int n = 2 * M + 1;          // number of retained orders
    const cdouble j{0.0, 1.0};
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
    VectorXcd Yi(n), Yii(n);  // for TE the matching "admittance" is just kz/k0
    for (int idx = 0; idx < n; ++idx) {
        Yi(idx) = kz_of(n_inc, idx);
        Yii(idx) = kz_of(n_sub, idx);
    }

    // --- Permittivity Toeplitz matrix E_{a,b} = ε_{(a−b)} ------------------
    MatrixXcd E(n, n);
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b)
            E(a, b) = grating.eps_fourier((a - M) - (b - M), wavelength_um);

    // --- Eigenproblem inside the grating: d²S/dz'² = A S, A = Kx² − E ------
    MatrixXcd Kx = kx.asDiagonal();
    MatrixXcd A = Kx * Kx - E;

    ComplexEigenSolver<MatrixXcd> ces(A);
    VectorXcd q = ces.eigenvalues().cwiseSqrt();  // q = kz_layer/k0 per mode
    // Force decaying branch (Re(q) >= 0).
    for (int i = 0; i < n; ++i)
        if (q(i).real() < 0) q(i) = -q(i);
    MatrixXcd W = ces.eigenvectors();
    MatrixXcd V = W * q.asDiagonal();             // companion (field-derivative) modes

    // Stable exponential of the layer: X = exp(−k0·q·d), |entries| <= 1.
    VectorXcd xdiag = (-k0 * grating.thickness_um * q.array()).exp();
    MatrixXcd X = xdiag.asDiagonal();

    // --- Boundary-condition linear system ----------------------------------
    // Unknowns u = [R; T; a; b], each length n. Field/derivative continuity at
    // z=0 and z=d (enhanced-transmittance referencing keeps it numerically
    // stable). See header reference.
    const MatrixXcd I = MatrixXcd::Identity(n, n);
    const MatrixXcd jYi = (j * Yi).asDiagonal();
    const MatrixXcd jYii = (j * Yii).asDiagonal();
    const MatrixXcd WX = W * X;
    const MatrixXcd VX = V * X;

    VectorXcd delta = VectorXcd::Zero(n);
    delta(M) = 1.0;  // unit-amplitude incident wave in order 0

    MatrixXcd S = MatrixXcd::Zero(4 * n, 4 * n);
    VectorXcd rhs = VectorXcd::Zero(4 * n);

    // z=0, field:        δ + R = W a + WX b
    S.block(0, 0, n, n) = I;
    S.block(0, 2 * n, n, n) = -W;
    S.block(0, 3 * n, n, n) = -WX;
    rhs.segment(0, n) = -delta;
    // z=0, derivative:   jYi(δ − R) = V a − VX b
    S.block(n, 0, n, n) = -jYi;
    S.block(n, 2 * n, n, n) = -V;
    S.block(n, 3 * n, n, n) = VX;
    rhs.segment(n, n) = -jYi * delta;
    // z=d, field:        WX a + W b = T
    S.block(2 * n, n, n, n) = -I;
    S.block(2 * n, 2 * n, n, n) = WX;
    S.block(2 * n, 3 * n, n, n) = W;
    // z=d, derivative:   VX a − V b = jYii T
    S.block(3 * n, n, n, n) = -jYii;
    S.block(3 * n, 2 * n, n, n) = VX;
    S.block(3 * n, 3 * n, n, n) = -V;

    VectorXcd sol = S.partialPivLu().solve(rhs);
    VectorXcd R = sol.segment(0, n);
    VectorXcd T = sol.segment(n, n);

    // --- Diffraction efficiencies ------------------------------------------
    // Normalize by the longitudinal wavevector of the incident order.
    const double kz_inc = (n_inc * std::cos(cdouble{theta0_rad, 0.0})).real();

    Rcwa1DResult out;
    out.orders.reserve(n);
    out.de_r.reserve(n);
    out.de_t.reserve(n);
    double total = 0.0;
    for (int idx = 0; idx < n; ++idx) {
        double der = std::norm(R(idx)) * Yi(idx).real() / kz_inc;
        double det = std::norm(T(idx)) * Yii(idx).real() / kz_inc;
        out.orders.push_back(idx - M);
        out.de_r.push_back(der);
        out.de_t.push_back(det);
        total += der + det;
    }
    out.sum_de = total;
    return out;
}

} // namespace celeris
