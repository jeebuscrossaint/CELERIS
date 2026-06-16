// Multilayer 1D RCWA via the scattering-matrix (Redheffer) method.
//
// Each layer is reduced to its Bloch eigenmodes (same as the single-layer
// solver). Interfaces between adjacent mode bases and propagation within each
// layer are expressed as 2x2-block scattering matrices, then combined with the
// Redheffer star product. Because propagation only ever uses the DECAYING
// exponential exp(-k0·q·d) (|entries| <= 1), the recursion is numerically
// stable for arbitrarily thick or evanescent layers — the whole reason to use
// S-matrices instead of a transfer-matrix product.
//
// Reference: Li, "Formulation and comparison of two recursive matrix
// algorithms for modeling layered diffraction gratings," JOSA A 13 (1996).

#include "celeris/rcwa/rcwa1d.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

namespace celeris {
namespace {

using Eigen::MatrixXcd;
using Eigen::VectorXcd;

// Eigenmodes of one layer: field modes W, companion modes V, and the stable
// propagation factor X = exp(-k0·q·d).
struct LayerModes {
    MatrixXcd W;
    MatrixXcd V;
    MatrixXcd X;
};

LayerModes compute_layer_modes(const GratingLayer1D& layer,
                               const MatrixXcd& Kx, int M, double wavelength_um,
                               double k0, Pol pol) {
    const int n = 2 * M + 1;
    const MatrixXcd Id = MatrixXcd::Identity(n, n);

    MatrixXcd E(n, n);
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b)
            E(a, b) = layer.eps_fourier((a - M) - (b - M), wavelength_um);

    MatrixXcd A;
    if (pol == Pol::TE) {
        A = Kx * Kx - E;
    } else {
        MatrixXcd Einv(n, n);
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b)
                Einv(a, b) =
                    layer.eps_inv_fourier((a - M) - (b - M), wavelength_um);
        A = E * (Kx * Einv * Kx - Id);
    }

    Eigen::ComplexEigenSolver<MatrixXcd> ces(A);
    VectorXcd q = ces.eigenvalues().cwiseSqrt();
    for (int i = 0; i < n; ++i)
        if (q(i).real() < 0) q(i) = -q(i);  // decaying branch

    LayerModes m;
    m.W = ces.eigenvectors();
    const cdouble j_unit{0.0, 1.0};
    m.V = (pol == Pol::TE) ? MatrixXcd(m.W * q.asDiagonal())
                           : MatrixXcd(j_unit * (E.inverse() * m.W) *
                                       q.asDiagonal());
    m.X = (-k0 * layer.thickness_um * q.array()).exp().matrix().asDiagonal();
    return m;
}

// A 2x2-block scattering matrix mapping incoming -> outgoing amplitudes.
struct SMatrix {
    MatrixXcd S11, S12, S21, S22;
};

SMatrix identity_smatrix(int n) {
    SMatrix s;
    s.S11 = MatrixXcd::Zero(n, n);
    s.S12 = MatrixXcd::Identity(n, n);
    s.S21 = MatrixXcd::Identity(n, n);
    s.S22 = MatrixXcd::Zero(n, n);
    return s;
}

// Redheffer star product  A ⋆ B  (A above, B below).
SMatrix star(const SMatrix& A, const SMatrix& B) {
    const int n = A.S11.rows();
    const MatrixXcd Id = MatrixXcd::Identity(n, n);
    const MatrixXcd D = (Id - B.S11 * A.S22).inverse();
    const MatrixXcd F = (Id - A.S22 * B.S11).inverse();
    SMatrix S;
    S.S11 = A.S11 + A.S12 * D * B.S11 * A.S21;
    S.S12 = A.S12 * D * B.S12;
    S.S21 = B.S21 * F * A.S21;
    S.S22 = B.S22 + B.S21 * F * A.S22 * B.S12;
    return S;
}

// Scattering matrix of the interface between region "a" (above, modes Wa,Va)
// and region "b" (below, modes Wb,Vb), with all amplitudes referenced at the
// interface. Derived from continuity of the field (W) and companion (V):
//   Wa(u_a+ + u_a-) = Wb(u_b+ + u_b-)
//   Va(u_a+ - u_a-) = Vb(u_b+ - u_b-)
SMatrix interface_smatrix(const MatrixXcd& Wa, const MatrixXcd& Va,
                          const MatrixXcd& Wb, const MatrixXcd& Vb) {
    const int n = Wa.rows();
    // [ Wa  -Wb ] [u_a-]   [ -Wa   Wb ] [u_a+]
    // [ Va   Vb ] [u_b+] = [  Va   Vb ] [u_b-]
    MatrixXcd Mo(2 * n, 2 * n), Mi(2 * n, 2 * n);
    Mo.block(0, 0, n, n) = Wa;
    Mo.block(0, n, n, n) = -Wb;
    Mo.block(n, 0, n, n) = Va;
    Mo.block(n, n, n, n) = Vb;
    Mi.block(0, 0, n, n) = -Wa;
    Mi.block(0, n, n, n) = Wb;
    Mi.block(n, 0, n, n) = Va;
    Mi.block(n, n, n, n) = Vb;

    MatrixXcd Sfull = Mo.partialPivLu().solve(Mi);  // outgoing = S * incoming
    SMatrix s;
    s.S11 = Sfull.block(0, 0, n, n);
    s.S12 = Sfull.block(0, n, n, n);
    s.S21 = Sfull.block(n, 0, n, n);
    s.S22 = Sfull.block(n, n, n, n);
    return s;
}

// Propagation across a layer with modal factor X (decaying both ways).
SMatrix propagation_smatrix(const MatrixXcd& X) {
    const int n = X.rows();
    SMatrix s;
    s.S11 = MatrixXcd::Zero(n, n);
    s.S12 = X;
    s.S21 = X;
    s.S22 = MatrixXcd::Zero(n, n);
    return s;
}

} // namespace

Rcwa1DResult solve_rcwa_1d(const Material& incident,
                           const Rcwa1DStack& stack,
                           const Material& substrate,
                           double wavelength_um,
                           double theta0_rad,
                           int M,
                           Pol pol) {
    if (M < 0) throw std::invalid_argument("solve_rcwa_1d: M must be >= 0");
    if (stack.layers.empty())
        throw std::invalid_argument("solve_rcwa_1d: stack has no layers");

    const int n = 2 * M + 1;
    const double k0 = 2.0 * pi / wavelength_um;
    const cdouble n_inc = incident.index(wavelength_um);
    const cdouble n_sub = substrate.index(wavelength_um);

    // Floquet order wavevectors (normalized by k0).
    VectorXcd kx(n);
    for (int m = -M; m <= M; ++m)
        kx(m + M) = n_inc * std::sin(theta0_rad) -
                    static_cast<double>(m) * (wavelength_um / stack.period_um);
    const MatrixXcd Kx = kx.asDiagonal();

    auto kz_of = [&](cdouble nreg, int idx) {
        return std::sqrt(nreg * nreg - kx(idx) * kx(idx));
    };
    VectorXcd kzI(n), kzII(n);
    for (int idx = 0; idx < n; ++idx) {
        kzI(idx) = kz_of(n_inc, idx);
        kzII(idx) = kz_of(n_sub, idx);
    }

    // Half-space "modes": plane waves, so field modes are the identity; the
    // companion modes are the region admittances (matching the single-layer
    // solver's convention exactly so the two agree for a 1-layer stack).
    const cdouble j_unit{0.0, 1.0};
    const MatrixXcd Wh = MatrixXcd::Identity(n, n);
    MatrixXcd Vinc, Vsub;
    if (pol == Pol::TE) {
        Vinc = (j_unit * kzI).matrix().asDiagonal();
        Vsub = (j_unit * kzII).matrix().asDiagonal();
    } else {
        Vinc = (kzI.array() / (n_inc * n_inc)).matrix().asDiagonal();
        Vsub = (kzII.array() / (n_sub * n_sub)).matrix().asDiagonal();
    }

    // Build the global scattering matrix top -> bottom.
    SMatrix global = identity_smatrix(n);
    MatrixXcd W_above = Wh, V_above = Vinc;  // start in the incident half-space
    for (const GratingLayer1D& layer : stack.layers) {
        LayerModes lm = compute_layer_modes(layer, Kx, M, wavelength_um, k0, pol);
        global = star(global, interface_smatrix(W_above, V_above, lm.W, lm.V));
        global = star(global, propagation_smatrix(lm.X));
        W_above = lm.W;
        V_above = lm.V;
    }
    // Final interface into the substrate half-space.
    global = star(global, interface_smatrix(W_above, V_above, Wh, Vsub));

    // Incident: unit amplitude in order 0, nothing entering from below.
    VectorXcd delta = VectorXcd::Zero(n);
    delta(M) = 1.0;
    VectorXcd R = global.S11 * delta;  // reflected order amplitudes
    VectorXcd T = global.S21 * delta;  // transmitted order amplitudes

    // Diffraction efficiencies (same flux bookkeeping as the single-layer path).
    const double kz_inc = (n_inc * std::cos(cdouble{theta0_rad, 0.0})).real();
    const cdouble eps_inc = (pol == Pol::TE) ? cdouble{1.0, 0.0} : n_inc * n_inc;
    const cdouble eps_sub = (pol == Pol::TE) ? cdouble{1.0, 0.0} : n_sub * n_sub;

    Rcwa1DResult out;
    double total = 0.0;
    for (int idx = 0; idx < n; ++idx) {
        double der = std::norm(R(idx)) * kzI(idx).real() / kz_inc;
        double det = std::norm(T(idx)) * (kzII(idx) / eps_sub).real() *
                     eps_inc.real() / kz_inc;
        out.orders.push_back(idx - M);
        out.de_r.push_back(der);
        out.de_t.push_back(det);
        total += der + det;
    }
    out.sum_de = total;
    return out;
}

} // namespace celeris
