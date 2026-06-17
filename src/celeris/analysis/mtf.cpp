#include "celeris/analysis/mtf.hpp"

#include "celeris/core.hpp"

#include <cmath>
#include <vector>

namespace celeris {

MtfCurve analyze_mtf(const MetalensDesign& lens, const UnitCellLibrary& lib,
                     double focal_length_um, double wavelength_um,
                     double diameter_um) {
    const int n = lens.n_cells;
    const double p = lens.period_um;
    const double center = (n - 1) / 2.0;
    const double R = diameter_um / 2.0;
    const double k = 2.0 * pi / wavelength_um;

    // Build the complex exit pupil P = t·exp(-i·phi_ideal) (i.e. |t|·e^{i2π·OPD})
    // and the perfect reference pupil P0 = 1 inside the aperture.
    std::vector<cdouble> P(static_cast<std::size_t>(n) * n, cdouble{0, 0});
    std::vector<double> P0(static_cast<std::size_t>(n) * n, 0.0);
    for (int iy = 0; iy < n; ++iy)
        for (int ix = 0; ix < n; ++ix) {
            const double x = (ix - center) * p, y = (iy - center) * p;
            const double r = std::sqrt(x * x + y * y);
            if (r > R) continue;
            const double fill = lens.fill_map[static_cast<std::size_t>(iy) * n + ix];
            const double phi_ideal =
                -k * (std::sqrt(r * r + focal_length_um * focal_length_um) - focal_length_um);
            const std::size_t idx = static_cast<std::size_t>(iy) * n + ix;
            P[idx] = lib.transmission_for_fill(fill) * std::polar(1.0, -phi_ideal);
            P0[idx] = 1.0;
        }

    // OTF(shift) = autocorrelation of the pupil along x, normalized by shift 0.
    // A pupil shift of s cells maps to image-plane spatial frequency
    //   nu = s*p / (lambda*f)  [cycles/um]  ->  *1000 = cycles/mm.
    auto autocorr = [&](auto& pupil, int s, auto normRef) -> double {
        cdouble acc{0, 0};
        for (int iy = 0; iy < n; ++iy)
            for (int ix = 0; ix + s < n; ++ix) {
                acc += cdouble(pupil[static_cast<std::size_t>(iy) * n + ix]) *
                       std::conj(cdouble(pupil[static_cast<std::size_t>(iy) * n + ix + s]));
            }
        return std::abs(acc) / normRef;
    };

    double normP = autocorr(P, 0, 1.0);
    double normP0 = autocorr(P0, 0, 1.0);
    if (normP <= 0) normP = 1;
    if (normP0 <= 0) normP0 = 1;

    MtfCurve c;
    const int max_s = static_cast<int>(diameter_um / p) + 1;  // to the cutoff
    for (int s = 0; s <= max_s && s < n; ++s) {
        double nu_cyc_mm = (s * p) / (wavelength_um * focal_length_um) * 1000.0;
        c.freq_cyc_mm.push_back(static_cast<float>(nu_cyc_mm));
        c.mtf.push_back(static_cast<float>(autocorr(P, s, normP)));
        c.mtf_ideal.push_back(static_cast<float>(autocorr(P0, s, normP0)));
    }
    c.cutoff_cyc_mm = diameter_um / (wavelength_um * focal_length_um) * 1000.0;
    return c;
}

} // namespace celeris
