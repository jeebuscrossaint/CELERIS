#include "celeris/rcwa/grating1d.hpp"

#include <cmath>

namespace celeris {

namespace {
// Fourier coefficient of a binary function that takes value `a` over the ridge
// (fraction `fill`, centered at x=0) and `b` over the groove.
cdouble binary_fourier(cdouble a, cdouble b, double fill, int h) {
    if (h == 0) return fill * a + (1.0 - fill) * b;
    const double x = pi * h * fill;
    return (a - b) * (std::sin(x) / (pi * h));
}
} // namespace

cdouble BinaryGrating1D::eps_fourier(int h, double wavelength_um) const {
    return binary_fourier(ridge.permittivity(wavelength_um),
                          groove.permittivity(wavelength_um), fill, h);
}

cdouble BinaryGrating1D::eps_inv_fourier(int h, double wavelength_um) const {
    return binary_fourier(1.0 / ridge.permittivity(wavelength_um),
                          1.0 / groove.permittivity(wavelength_um), fill, h);
}

} // namespace celeris
