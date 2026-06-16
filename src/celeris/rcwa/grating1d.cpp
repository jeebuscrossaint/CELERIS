#include "celeris/rcwa/grating1d.hpp"

#include <cmath>

namespace celeris {

cdouble BinaryGrating1D::eps_fourier(int h, double wavelength_um) const {
    const cdouble er = ridge.permittivity(wavelength_um);
    const cdouble eg = groove.permittivity(wavelength_um);
    if (h == 0) {
        return fill * er + (1.0 - fill) * eg;
    }
    const double x = pi * h * fill;
    return (er - eg) * (std::sin(x) / (pi * h));
}

} // namespace celeris
