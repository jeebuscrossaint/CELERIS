#include "celeris/rcwa/rcwa2d.hpp"
#include "celeris/rcwa/smatrix.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

namespace celeris {

cdouble RectCell2D::eps_fourier(int p, int q, double wavelength_um) const {
    const cdouble ep = pillar.permittivity(wavelength_um);
    const cdouble eb = background.permittivity(wavelength_um);
    // Separable rectangle: ε(x,y) = eb + (ep-eb)·Rect_x·Rect_y, centered.
    // Fourier: ε_{pq} = eb·δ_{p0}δ_{q0} + (ep-eb)·Sx(p)·Sy(q).
    const double Sx = (p == 0) ? fill_x : std::sin(pi * p * fill_x) / (pi * p);
    const double Sy = (q == 0) ? fill_y : std::sin(pi * q * fill_y) / (pi * q);
    const cdouble base = (p == 0 && q == 0) ? eb : cdouble{0.0, 0.0};
    return base + (ep - eb) * (Sx * Sy);
}

namespace {

using Eigen::MatrixXcd;
using Eigen::VectorXcd;
using detail::SMatrix;

// Flattened 2D order index for (p,q), p in [-Mx,Mx], q in [-My,My].
inline int order_index(int p, int q, int Mx, int My) {
    return (p + Mx) * (2 * My + 1) + (q + My);
}

struct Modes {
    MatrixXcd W;  // tangential-E field modes  [Ex;Ey], 2N×2N
    MatrixXcd V;  // tangential-H companion modes, 2N×2N
};

} // namespace

Rcwa2DResult solve_rcwa_2d(const Material& incident,
                           const Rcwa2DStack& stack,
                           const Material& substrate,
                           double wavelength_um,
                           double theta_rad,
                           double phi_rad,
                           cdouble Ex0,
                           cdouble Ey0,
                           int Mx,
                           int My) {
    if (Mx < 0 || My < 0) throw std::invalid_argument("M must be >= 0");
    if (stack.layers.empty())
        throw std::invalid_argument("solve_rcwa_2d: stack has no layers");

    const int N = (2 * Mx + 1) * (2 * My + 1);  // number of 2D orders
    const double k0 = 2.0 * pi / wavelength_um;
    const cdouble n_inc = incident.index(wavelength_um);
    const cdouble n_sub = substrate.index(wavelength_um);
    const cdouble j_unit{0.0, 1.0};

    // --- Floquet order wavevectors (normalized by k0) ----------------------
    const cdouble kx0 = n_inc * std::sin(theta_rad) * std::cos(phi_rad);
    const cdouble ky0 = n_inc * std::sin(theta_rad) * std::sin(phi_rad);
    VectorXcd kxv(N), kyv(N);
    for (int p = -Mx; p <= Mx; ++p)
        for (int q = -My; q <= My; ++q) {
            const int i = order_index(p, q, Mx, My);
            kxv(i) = kx0 - static_cast<double>(p) * (wavelength_um / stack.period_x_um);
            kyv(i) = ky0 - static_cast<double>(q) * (wavelength_um / stack.period_y_um);
        }
    const MatrixXcd Kx = kxv.asDiagonal();
    const MatrixXcd Ky = kyv.asDiagonal();
    const MatrixXcd Id = MatrixXcd::Identity(N, N);

    auto kz_region = [&](cdouble nreg, int i) {
        cdouble kz = std::sqrt(nreg * nreg - kxv(i) * kxv(i) - kyv(i) * kyv(i));
        // Regularize orders sitting exactly at grazing (kz == 0, the Rayleigh
        // anomaly): nudge to a tiny evanescent value so reciprocals stay finite.
        // Such orders carry negligible power, so this does not affect results.
        if (std::abs(kz) < 1e-8) kz = cdouble{0.0, 1e-8};
        return kz;
    };

    // --- Homogeneous half-space modes: W = I, V = Q_h·Ω_h⁻¹ ----------------
    // For a homogeneous region the natural basis is plane waves (W = I); the
    // companion follows from the homogeneous Q with Ω = i·Kz.
    auto half_space_modes = [&](cdouble nreg) -> Modes {
        VectorXcd kz(N);
        for (int i = 0; i < N; ++i) kz(i) = kz_region(nreg, i);
        const cdouble er = nreg * nreg;
        MatrixXcd Qh(2 * N, 2 * N);
        Qh.block(0, 0, N, N) = Kx * Ky;
        Qh.block(0, N, N, N) = er * Id - Kx * Kx;
        Qh.block(N, 0, N, N) = Ky * Ky - er * Id;
        Qh.block(N, N, N, N) = -Ky * Kx;
        // Ω_h = i·Kz on each block  ⇒  Ω_h⁻¹ = -i·Kz⁻¹.
        VectorXcd omega_inv(2 * N);
        for (int i = 0; i < N; ++i) {
            omega_inv(i) = -j_unit / kz(i);
            omega_inv(N + i) = -j_unit / kz(i);
        }
        Modes m;
        m.W = MatrixXcd::Identity(2 * N, 2 * N);
        m.V = Qh * omega_inv.asDiagonal();
        return m;
    };

    // --- One patterned layer: eigenmodes of P·Q ----------------------------
    auto layer_modes = [&](const RectCell2D& layer) {
        MatrixXcd Erc(N, N);
        for (int pi = -Mx; pi <= Mx; ++pi)
            for (int qi = -My; qi <= My; ++qi) {
                const int row = order_index(pi, qi, Mx, My);
                for (int pj = -Mx; pj <= Mx; ++pj)
                    for (int qj = -My; qj <= My; ++qj) {
                        const int col = order_index(pj, qj, Mx, My);
                        Erc(row, col) =
                            layer.eps_fourier(pi - pj, qi - qj, wavelength_um);
                    }
            }
        const MatrixXcd Erci = Erc.inverse();

        // P and Q block operators (μ = 1).
        MatrixXcd P(2 * N, 2 * N), Q(2 * N, 2 * N);
        P.block(0, 0, N, N) = Kx * Erci * Ky;
        P.block(0, N, N, N) = Id - Kx * Erci * Kx;
        P.block(N, 0, N, N) = Ky * Erci * Ky - Id;
        P.block(N, N, N, N) = -Ky * Erci * Kx;
        Q.block(0, 0, N, N) = Kx * Ky;
        Q.block(0, N, N, N) = Erc - Kx * Kx;
        Q.block(N, 0, N, N) = Ky * Ky - Erc;
        Q.block(N, N, N, N) = -Ky * Kx;

        Eigen::ComplexEigenSolver<MatrixXcd> ces(P * Q);
        VectorXcd omega = ces.eigenvalues().cwiseSqrt();
        for (int i = 0; i < 2 * N; ++i)
            if (omega(i).real() < 0) omega(i) = -omega(i);  // decaying branch
        MatrixXcd W = ces.eigenvectors();
        VectorXcd omega_inv = omega.cwiseInverse();
        MatrixXcd V = Q * W * omega_inv.asDiagonal();
        MatrixXcd X = (-k0 * layer.thickness_um * omega.array()).exp().matrix().asDiagonal();
        return std::make_tuple(W, V, X);
    };

    // --- Chain the scattering matrices top -> bottom -----------------------
    using namespace detail;
    SMatrix global = identity_smatrix(2 * N);
    Modes inc_modes = half_space_modes(n_inc);
    MatrixXcd W_above = inc_modes.W, V_above = inc_modes.V;
    for (const RectCell2D& layer : stack.layers) {
        auto [W, V, X] = layer_modes(layer);
        global = star(global, interface_smatrix(W_above, V_above, W, V));
        global = star(global, propagation_smatrix(X));
        W_above = W;
        V_above = V;
    }
    Modes sub_modes = half_space_modes(n_sub);
    global = star(global, interface_smatrix(W_above, V_above, sub_modes.W, sub_modes.V));

    // --- Excite order 0 with the given polarization ------------------------
    const int i0 = order_index(0, 0, Mx, My);
    VectorXcd delta = VectorXcd::Zero(2 * N);
    delta(i0) = Ex0;       // Ex block
    delta(N + i0) = Ey0;   // Ey block
    VectorXcd Rvec = global.S11 * delta;
    VectorXcd Tvec = global.S21 * delta;

    // --- Diffraction efficiencies ------------------------------------------
    // Per order, recover Ez from transversality k·E = 0, then the z-directed
    // power is ∝ Re(kz)·|E|². Normalize by the incident order-0 flux.
    const double kz_inc = kz_region(n_inc, i0).real();
    Rcwa2DResult out{};
    double sumR = 0.0, sumT = 0.0;
    for (int i = 0; i < N; ++i) {
        const cdouble rx = Rvec(i), ry = Rvec(N + i);
        const cdouble tx = Tvec(i), ty = Tvec(N + i);
        const cdouble kzr = kz_region(n_inc, i);
        const cdouble kzt = kz_region(n_sub, i);
        const cdouble rz = (kxv(i) * rx + kyv(i) * ry) / kzr;
        const cdouble tz = (kxv(i) * tx + kyv(i) * ty) / kzt;
        const double der = (kzr.real() / kz_inc) *
                           (std::norm(rx) + std::norm(ry) + std::norm(rz));
        const double det = (kzt.real() / kz_inc) *
                           (std::norm(tx) + std::norm(ty) + std::norm(tz));
        sumR += der;
        sumT += det;
        if (i == i0) {
            out.de_r0 = der;
            out.de_t0 = det;
            out.tx0 = tx;
            out.ty0 = ty;
        }
    }
    out.R = sumR;
    out.T = sumT;
    out.sum_de = sumR + sumT;
    return out;
}

} // namespace celeris
