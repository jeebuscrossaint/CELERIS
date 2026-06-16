#include "celeris/optics/tmm.hpp"

#include <Eigen/Dense>
#include <cmath>

namespace celeris {
namespace {

// Optical admittance η of a medium for a given polarization, given the complex
// cosine of the propagation angle inside it. Free-space admittance is
// normalized to 1, which cancels in the R/T ratios below.
//   TE (s): η = n·cosθ
//   TM (p): η = n/cosθ
cdouble admittance(cdouble n, cdouble cos_theta, Pol pol) {
    return pol == Pol::TE ? n * cos_theta : n / cos_theta;
}

// cosθ inside a medium of index n, given the conserved tangential wavevector
// from Snell's law: n0·sinθ0 = n·sinθ  ⇒  cosθ = sqrt(1 − (n0 sinθ0 / n)²).
cdouble cos_theta_in(cdouble n, cdouble n0_sin0) {
    cdouble s = n0_sin0 / n;
    return std::sqrt(cdouble{1.0, 0.0} - s * s);
}

} // namespace

TmmResult solve_stack(const Material& incident,
                      const std::vector<Layer>& layers,
                      const Material& substrate,
                      double wavelength_um,
                      double theta0_rad,
                      Pol pol) {
    using Mat2 = Eigen::Matrix2cd;

    const cdouble n0 = incident.index(wavelength_um);
    const cdouble ns = substrate.index(wavelength_um);

    // Conserved across every interface (Snell's law invariant).
    const cdouble n0_sin0 = n0 * std::sin(theta0_rad);
    const double k0 = 2.0 * pi / wavelength_um;  // vacuum wavenumber

    const cdouble cos0 = std::cos(cdouble{theta0_rad, 0.0});
    const cdouble eta0 = admittance(n0, cos0, pol);

    // Accumulate the total characteristic matrix M = M_1 · M_2 · … · M_N.
    Mat2 M = Mat2::Identity();
    for (const Layer& layer : layers) {
        const cdouble n = layer.material.index(wavelength_um);
        const cdouble cos_t = cos_theta_in(n, n0_sin0);
        const cdouble eta = admittance(n, cos_t, pol);

        // Phase thickness of the layer: δ = k0 · n · d · cosθ.
        const cdouble delta = k0 * n * layer.thickness_um * cos_t;
        const cdouble cos_d = std::cos(delta);
        const cdouble sin_d = std::sin(delta);
        const cdouble i{0.0, 1.0};

        Mat2 Ml;
        Ml << cos_d,            i * sin_d / eta,
              i * eta * sin_d,  cos_d;
        M = M * Ml;
    }

    // Substrate admittance and the [B; C] vector: [B;C] = M · [1; η_sub].
    const cdouble cos_s = cos_theta_in(ns, n0_sin0);
    const cdouble etas = admittance(ns, cos_s, pol);

    const cdouble B = M(0, 0) * 1.0 + M(0, 1) * etas;
    const cdouble C = M(1, 0) * 1.0 + M(1, 1) * etas;

    // Reflection/transmission from the standard Macleod relations
    // (valid for a lossless incident medium, so η0 is real).
    const cdouble denom = eta0 * B + C;
    const cdouble r = (eta0 * B - C) / denom;
    const cdouble t = 2.0 * eta0 / denom;

    TmmResult out;
    out.r = r;
    out.t = t;
    out.R = std::norm(r);
    out.T = 4.0 * eta0.real() * etas.real() / std::norm(denom);
    return out;
}

} // namespace celeris
