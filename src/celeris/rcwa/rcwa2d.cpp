#include "celeris/rcwa/rcwa2d.hpp"
#include "celeris/rcwa/eig.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <vector>

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

cdouble RectCell2D::inv_eps_fourier(int p, int q, double wavelength_um) const {
    // 1/ε is the same separable rectangle with ε -> 1/ε (the rectangle takes two
    // values, ep inside and eb outside, so its reciprocal is also two-valued):
    //   (1/ε)_{pq} = (1/eb)·δ_{p0}δ_{q0} + (1/ep - 1/eb)·Sx(p)·Sy(q).
    const cdouble ipe = 1.0 / pillar.permittivity(wavelength_um);
    const cdouble ibe = 1.0 / background.permittivity(wavelength_um);
    const double Sx = (p == 0) ? fill_x : std::sin(pi * p * fill_x) / (pi * p);
    const double Sy = (q == 0) ? fill_y : std::sin(pi * q * fill_y) / (pi * q);
    const cdouble base = (p == 0 && q == 0) ? ibe : cdouble{0.0, 0.0};
    return base + (ipe - ibe) * (Sx * Sy);
}

namespace {

using Eigen::MatrixXcd;
using Eigen::VectorXcd;

// Flattened 2D order index for (p,q), p in [-Mx,Mx], q in [-My,My].
inline int order_index(int p, int q, int Mx, int My) {
    return (p + Mx) * (2 * My + 1) + (q + My);
}

// 1D Fourier coefficient of a centered rect of fractional width `fill`.
inline double rect1d(int m, double fill) {
    return (m == 0) ? fill : std::sin(pi * m * fill) / (pi * m);
}

// Modes of one region (Liu-Fan / grcwa convention, physical units):
//   phi  : eigenvectors = tangential-H Fourier modes [Hx; Hy]  (2N x 2N)
//   q    : longitudinal kz per mode (length 2N)
//   kp   : ω²I − Jk·⟦1/ε⟧·Jkᵀ, the operator that maps H-modes to E-modes
//   thickness : layer thickness (0 for the bounding half-spaces)
struct Region {
    VectorXcd q;
    MatrixXcd phi;
    MatrixXcd kp;
    double thickness;
};

// Branch choice (grcwa): keep Im(kz) >= 0 so exp(i kz z) decays for evanescent
// orders and the scattering-matrix recursion stays stable.
inline cdouble branch(cdouble q) { return q.imag() < 0.0 ? -q : q; }

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
    const int N2 = 2 * N;
    const double k0 = 2.0 * pi / wavelength_um;
    const cdouble omega = k0;            // c = 1 units, ω = k0
    const cdouble omega2 = omega * omega;
    const cdouble n_inc = incident.index(wavelength_um);
    const cdouble n_sub = substrate.index(wavelength_um);
    const cdouble jj{0.0, 1.0};

    // --- Floquet order wavevectors (PHYSICAL, grcwa units) -----------------
    const cdouble kx0 = n_inc * std::sin(theta_rad) * std::cos(phi_rad) * omega;
    const cdouble ky0 = n_inc * std::sin(theta_rad) * std::sin(phi_rad) * omega;
    VectorXcd kxv(N), kyv(N);
    for (int p = -Mx; p <= Mx; ++p)
        for (int q = -My; q <= My; ++q) {
            const int i = order_index(p, q, Mx, My);
            kxv(i) = kx0 - 2.0 * pi * static_cast<double>(p) / stack.period_x_um;
            kyv(i) = ky0 - 2.0 * pi * static_cast<double>(q) / stack.period_y_um;
        }
    const MatrixXcd Kx = kxv.asDiagonal();
    const MatrixXcd Ky = kyv.asDiagonal();
    const MatrixXcd Id = MatrixXcd::Identity(N, N);

    // kkT = [[Kx Kx, Kx Ky],[Ky Kx, Ky Ky]]  (shared by every region).
    MatrixXcd kkT(N2, N2);
    kkT.block(0, 0, N, N) = Kx * Kx;
    kkT.block(0, N, N, N) = Kx * Ky;
    kkT.block(N, 0, N, N) = Ky * Kx;
    kkT.block(N, N, N, N) = Ky * Ky;

    // --- Homogeneous region: analytic plane-wave modes ---------------------
    auto uniform_region = [&](cdouble nreg, double thick) -> Region {
        const cdouble er = nreg * nreg, ier = 1.0 / er;
        VectorXcd q(N2);
        for (int i = 0; i < N; ++i) {
            cdouble kz = std::sqrt(er * omega2 - kxv(i) * kxv(i) - kyv(i) * kyv(i));
            if (std::abs(kz) < 1e-8) kz = cdouble{0.0, 1e-8};  // grazing guard
            kz = branch(kz);
            q(i) = kz;
            q(N + i) = kz;
        }
        MatrixXcd kp(N2, N2);
        kp.block(0, 0, N, N) = omega2 * Id - ier * (Ky * Ky);
        kp.block(0, N, N, N) = ier * (Ky * Kx);
        kp.block(N, 0, N, N) = ier * (Kx * Ky);
        kp.block(N, N, N, N) = omega2 * Id - ier * (Kx * Kx);
        return {q, MatrixXcd::Identity(N2, N2), kp, thick};
    };

    // --- Patterned region: improved (Li) factorization + Liu-Fan eigsystem -
    auto patterned_region = [&](const RectCell2D& layer) -> Region {
        const cdouble ep = layer.pillar.permittivity(wavelength_um);
        const cdouble eb = layer.background.permittivity(wavelength_um);
        const cdouble ipe = 1.0 / ep, ibe = 1.0 / eb;
        const double fx = layer.fill_x, fy = layer.fill_y;
        const int nx = 2 * Mx + 1, ny = 2 * My + 1;

        // ⟦1/ε⟧ : Toeplitz of the reciprocal permittivity (normal component).
        MatrixXcd Einv(N, N);
        for (int pi_ = -Mx; pi_ <= Mx; ++pi_)
            for (int qi = -My; qi <= My; ++qi) {
                const int row = order_index(pi_, qi, Mx, My);
                for (int pj = -Mx; pj <= Mx; ++pj)
                    for (int qj = -My; qj <= My; ++qj)
                        Einv(row, order_index(pj, qj, Mx, My)) =
                            layer.inv_eps_fourier(pi_ - pj, qi - qj, wavelength_um);
            }

        // Exx: inverse rule in x (invert the in-band 1D reciprocal Toeplitz),
        // Laurent in y. Eyy: symmetric. These multiply Ex / Ey respectively.
        MatrixXcd Hx(nx, nx);
        for (int a = -Mx; a <= Mx; ++a)
            for (int b = -Mx; b <= Mx; ++b)
                Hx(a + Mx, b + Mx) =
                    ((a == b) ? ibe : cdouble{0, 0}) + (ipe - ibe) * rect1d(a - b, fx);
        const MatrixXcd Axb = Hx.inverse();
        MatrixXcd Exx(N, N);
        for (int p = -Mx; p <= Mx; ++p)
            for (int q = -My; q <= My; ++q) {
                const int row = order_index(p, q, Mx, My);
                for (int pp = -Mx; pp <= Mx; ++pp)
                    for (int qq = -My; qq <= My; ++qq) {
                        const double sy = rect1d(q - qq, fy);
                        cdouble v = Axb(p + Mx, pp + Mx) * sy;
                        if (p == pp) v += eb * (((q == qq) ? 1.0 : 0.0) - sy);
                        Exx(row, order_index(pp, qq, Mx, My)) = v;
                    }
            }

        MatrixXcd Hy(ny, ny);
        for (int a = -My; a <= My; ++a)
            for (int b = -My; b <= My; ++b)
                Hy(a + My, b + My) =
                    ((a == b) ? ibe : cdouble{0, 0}) + (ipe - ibe) * rect1d(a - b, fy);
        const MatrixXcd Ayb = Hy.inverse();
        MatrixXcd Eyy(N, N);
        for (int p = -Mx; p <= Mx; ++p)
            for (int q = -My; q <= My; ++q) {
                const int row = order_index(p, q, Mx, My);
                for (int pp = -Mx; pp <= Mx; ++pp)
                    for (int qq = -My; qq <= My; ++qq) {
                        const double sx = rect1d(p - pp, fx);
                        cdouble v = Ayb(q + My, qq + My) * sx;
                        if (q == qq) v += eb * (((p == pp) ? 1.0 : 0.0) - sx);
                        Eyy(row, order_index(pp, qq, Mx, My)) = v;
                    }
            }

        // kp = ω²I − Jk·⟦1/ε⟧·Jkᵀ   (Jk = [−Ky; Kx])
        MatrixXcd kp(N2, N2);
        kp.block(0, 0, N, N) = omega2 * Id - Ky * Einv * Ky;
        kp.block(0, N, N, N) = Ky * Einv * Kx;
        kp.block(N, 0, N, N) = Kx * Einv * Ky;
        kp.block(N, N, N, N) = omega2 * Id - Kx * Einv * Kx;

        // ep2 = diag(Eyy, Exx); eigenproblem M = ep2·kp − kkT. Block order fixed
        // by the My=0 limit: the first (Hx/s-pol) block must reduce to the TE
        // operator ⟦ε⟧ω²−Kx² (Laurent → Eyy here) and the second to TM (Exx).
        MatrixXcd ep2 = MatrixXcd::Zero(N2, N2);
        ep2.block(0, 0, N, N) = Eyy;
        ep2.block(N, N, N, N) = Exx;

        VectorXcd evals;
        MatrixXcd phi;
        eig_general(ep2 * kp - kkT, evals, phi);
        VectorXcd q(N2);
        for (int i = 0; i < N2; ++i) q(i) = branch(std::sqrt(evals(i)));
        return {q, phi, kp, layer.thickness_um};
    };

    // --- Assemble regions: incident half-space, pillars, substrate ---------
    std::vector<Region> R;
    R.reserve(stack.layers.size() + 2);
    R.push_back(uniform_region(n_inc, 0.0));
    for (const RectCell2D& layer : stack.layers) R.push_back(patterned_region(layer));
    R.push_back(uniform_region(n_sub, 0.0));
    const int NL = static_cast<int>(R.size());

    // --- Global scattering matrix (grcwa GetSMatrix recursion) -------------
    MatrixXcd S11 = MatrixXcd::Identity(N2, N2);
    MatrixXcd S12 = MatrixXcd::Zero(N2, N2);
    MatrixXcd S21 = MatrixXcd::Zero(N2, N2);
    MatrixXcd S22 = MatrixXcd::Identity(N2, N2);
    for (int li = 0; li < NL - 1; ++li) {
        const Region& a = R[li];
        const Region& b = R[li + 1];
        const MatrixXcd Q = a.phi.inverse() * b.phi;
        const MatrixXcd P1 =
            a.q.asDiagonal().toDenseMatrix() * (a.kp * a.phi).inverse();
        const MatrixXcd P2 =
            (b.kp * b.phi) * b.q.cwiseInverse().asDiagonal();
        const MatrixXcd P = P1 * P2;
        const MatrixXcd T11 = 0.5 * (Q + P);
        const MatrixXcd T12 = 0.5 * (Q - P);
        const MatrixXcd d1 = (jj * a.q.array() * a.thickness).exp().matrix().asDiagonal();
        const MatrixXcd d2 = (jj * b.q.array() * b.thickness).exp().matrix().asDiagonal();

        const MatrixXcd Pinv = (T11 - d1 * S12 * T12).inverse();
        const MatrixXcd nS11 = Pinv * d1 * S11;
        const MatrixXcd nS12 = Pinv * (d1 * S12 * T11 - T12) * d2;
        const MatrixXcd nS21 = S21 + S22 * T12 * nS11;
        const MatrixXcd nS22 = S22 * T12 * nS12 + S22 * T11 * d2;
        S11 = nS11; S12 = nS12; S21 = nS21; S22 = nS22;
    }

    // --- Excite order 0; solve exterior amplitudes -------------------------
    // grcwa normal-incidence mapping: a0[i0] = −s, a0[N+i0] = p. The H-mode
    // amplitudes here are proportional to the incident E-field; we normalize the
    // reported t by the incident order-0 Ex below, so the absolute mapping
    // constant cancels.
    const int i0 = order_index(0, 0, Mx, My);
    VectorXcd a0 = VectorXcd::Zero(N2);
    a0(i0) = -Ey0;
    a0(N + i0) = Ex0;
    const VectorXcd bN = VectorXcd::Zero(N2);
    const VectorXcd aN = S11 * a0 + S12 * bN;  // forward amps in substrate
    const VectorXcd b0 = S21 * a0 + S22 * bN;  // backward amps in incident region

    // --- z-directed Poynting flux per order (grcwa GetZPoyntingFlux) -------
    // Returns the per-order forward and backward fluxes (length N each).
    auto flux = [&](const VectorXcd& ai, const VectorXcd& bi, const Region& reg,
                    VectorXcd& fwd, VectorXcd& bwd) {
        const MatrixXcd A =
            reg.kp * reg.phi * (omega * reg.q).cwiseInverse().asDiagonal();
        const VectorXcd pa = reg.phi * ai, pb = reg.phi * bi;
        const VectorXcd Aa = A * ai, Ab = A * bi;
        VectorXcd diff =
            0.5 * (pb.conjugate().cwiseProduct(Aa) - Ab.conjugate().cwiseProduct(pa));
        VectorXcd fxy = Aa.conjugate().cwiseProduct(pa).real().cast<cdouble>() + diff;
        VectorXcd bxy =
            -(Ab.conjugate().cwiseProduct(pb).real().cast<cdouble>()) + diff.conjugate();
        fwd = fxy.head(N) + fxy.tail(N);
        bwd = bxy.head(N) + bxy.tail(N);
    };

    VectorXcd fwd0, bwd0, fwdN, bwdN;
    flux(a0, b0, R.front(), fwd0, bwd0);          // incident region
    flux(aN, bN, R.back(), fwdN, bwdN);           // substrate

    const double norm = (n_inc.real()) / std::cos(theta_rad);  // grcwa normalization

    // --- Complex transmission (for the metalens phase library) -------------
    // E-field = kp·phi·(ai−bi)/(ω q), reordered: Ex = E[N+i], Ey = −E[i].
    auto efield = [&](const VectorXcd& ai, const VectorXcd& bi, const Region& reg) {
        const MatrixXcd A =
            reg.kp * reg.phi * (omega * reg.q).cwiseInverse().asDiagonal();
        return VectorXcd(A * (ai - bi));
    };
    const VectorXcd Esub = efield(aN, bN, R.back());
    const VectorXcd Einc = efield(a0, VectorXcd::Zero(N2), R.front());  // incident only
    const cdouble ex_inc = Einc(N + i0);
    const cdouble denom = (std::abs(ex_inc) > 1e-30) ? ex_inc : cdouble{1.0, 0.0};

    Rcwa2DResult out{};
    double sumR = 0.0, sumT = 0.0;
    for (int i = 0; i < N; ++i) {
        sumR += (-bwd0(i)).real() * norm;
        sumT += (fwdN(i)).real() * norm;
    }
    out.R = sumR;
    out.T = sumT;
    out.sum_de = sumR + sumT;
    out.de_r0 = (-bwd0(i0)).real() * norm;
    out.de_t0 = (fwdN(i0)).real() * norm;
    out.tx0 = Esub(N + i0) / denom;   // transmitted Ex / incident Ex
    out.ty0 = -Esub(i0) / denom;      // transmitted Ey / incident Ex
    return out;
}

} // namespace celeris
