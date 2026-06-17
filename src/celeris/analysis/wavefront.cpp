#include "celeris/analysis/wavefront.hpp"

#include "celeris/core.hpp"

#include <array>
#include <cmath>

namespace celeris {
namespace {

double wrap_waves(double phase_rad) {
    while (phase_rad > pi) phase_rad -= 2 * pi;
    while (phase_rad <= -pi) phase_rad += 2 * pi;
    return phase_rad / (2 * pi);  // -> (-0.5, 0.5] waves
}

// Low-order Zernike polynomials (Noll order, OSA normalization) on the unit
// disk, evaluated at (rho, theta). Returns Z_noll for noll in 1..11.
double zernike(int noll, double r, double t) {
    const double r2 = r * r, r3 = r2 * r, r4 = r2 * r2;
    switch (noll) {
        case 1:  return 1.0;                                   // piston
        case 2:  return 2.0 * r * std::cos(t);                 // tip
        case 3:  return 2.0 * r * std::sin(t);                 // tilt
        case 4:  return std::sqrt(3.0) * (2 * r2 - 1);         // defocus
        case 5:  return std::sqrt(6.0) * r2 * std::sin(2 * t); // astig 45
        case 6:  return std::sqrt(6.0) * r2 * std::cos(2 * t); // astig 0/90
        case 7:  return std::sqrt(8.0) * (3 * r3 - 2 * r) * std::sin(t); // coma y
        case 8:  return std::sqrt(8.0) * (3 * r3 - 2 * r) * std::cos(t); // coma x
        case 9:  return std::sqrt(8.0) * r3 * std::sin(3 * t); // trefoil
        case 10: return std::sqrt(8.0) * r3 * std::cos(3 * t); // trefoil
        case 11: return std::sqrt(5.0) * (6 * r4 - 6 * r2 + 1); // spherical
        default: return 0.0;
    }
}

const char* zernike_name(int noll) {
    switch (noll) {
        case 2:  return "tip";
        case 3:  return "tilt";
        case 4:  return "defocus";
        case 5:  return "astigmatism 45";
        case 6:  return "astigmatism 0/90";
        case 7:  return "coma (y)";
        case 8:  return "coma (x)";
        case 9:  return "trefoil (y)";
        case 10: return "trefoil (x)";
        case 11: return "spherical";
        default: return "piston";
    }
}

} // namespace

WavefrontAnalysis analyze_wavefront(const MetalensDesign& lens,
                                    const UnitCellLibrary& lib,
                                    double focal_length_um, double wavelength_um,
                                    double diameter_um) {
    const int n = lens.n_cells;
    const double p = lens.period_um;
    const double center = (n - 1) / 2.0;
    const double R = diameter_um / 2.0;
    const double k = 2.0 * pi / wavelength_um;

    WavefrontAnalysis w;
    w.n = n;
    w.opd.assign(static_cast<std::size_t>(n) * n, 0.0);
    w.mask.assign(static_cast<std::size_t>(n) * n, 0);

    // Build the OPD map (realized minus ideal focusing phase), in waves.
    std::vector<double> rho, theta, val;  // per in-aperture cell
    for (int iy = 0; iy < n; ++iy) {
        for (int ix = 0; ix < n; ++ix) {
            const double x = (ix - center) * p, y = (iy - center) * p;
            const double r = std::sqrt(x * x + y * y);
            if (r > R) continue;
            const double fill = lens.fill_map[static_cast<std::size_t>(iy) * n + ix];
            const double phi_real = std::arg(lib.transmission_for_fill(fill));
            const double phi_ideal =
                -k * (std::sqrt(r * r + focal_length_um * focal_length_um) - focal_length_um);
            const double opd = wrap_waves(phi_real - phi_ideal);
            const std::size_t idx = static_cast<std::size_t>(iy) * n + ix;
            w.opd[idx] = opd;
            w.mask[idx] = 1;
            rho.push_back(r / R);
            theta.push_back(std::atan2(y, x));
            val.push_back(opd);
        }
    }
    if (val.empty()) { w.rms_waves = w.pv_waves = 0; w.strehl_marechal = 1; return w; }

    // Remove piston (mean), then RMS / peak-to-valley.
    double mean = 0.0;
    for (double v : val) mean += v;
    mean /= val.size();
    double sq = 0.0, lo = 1e300, hi = -1e300;
    for (double v : val) {
        double d = v - mean;
        sq += d * d;
        lo = std::min(lo, v);
        hi = std::max(hi, v);
    }
    w.rms_waves = std::sqrt(sq / val.size());
    w.pv_waves = hi - lo;
    w.strehl_marechal = std::exp(-std::pow(2 * pi * w.rms_waves, 2));

    // Zernike decomposition by discrete least-squares projection.
    for (int j = 2; j <= 11; ++j) {
        double num = 0.0, den = 0.0;
        for (std::size_t i = 0; i < val.size(); ++i) {
            double z = zernike(j, rho[i], theta[i]);
            num += (val[i] - mean) * z;  // piston-removed
            den += z * z;
        }
        if (den > 0) w.zernike.push_back({j, zernike_name(j), num / den});
    }

    // Re-bias the OPD map to piston-removed for display.
    for (std::size_t i = 0; i < w.opd.size(); ++i)
        if (w.mask[i]) w.opd[i] -= mean;
    return w;
}

} // namespace celeris
