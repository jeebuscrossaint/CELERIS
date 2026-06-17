#include "celeris/analysis/throughfocus.hpp"

#include "celeris/core.hpp"

#include <algorithm>
#include <cmath>

namespace celeris {

ThroughFocus analyze_through_focus(const MetalensDesign& lens,
                                   const UnitCellLibrary& lib,
                                   double focal_length_um, double wavelength_um,
                                   double diameter_um) {
    const int n = lens.n_cells;
    const double p = lens.period_um;
    const double center = (n - 1) / 2.0;
    const double R = diameter_um / 2.0;
    const double k = 2.0 * pi / wavelength_um;

    struct Cell { double x, y; cdouble t; };
    std::vector<Cell> cells;
    for (int iy = 0; iy < n; ++iy)
        for (int ix = 0; ix < n; ++ix) {
            const double x = (ix - center) * p, y = (iy - center) * p;
            if (std::sqrt(x * x + y * y) > R) continue;
            const double fill = lens.fill_map[static_cast<std::size_t>(iy) * n + ix];
            cells.push_back({x, y, lib.transmission_for_fill(fill)});
        }

    auto on_axis = [&](double z) {
        cdouble E{0, 0};
        for (const Cell& c : cells) {
            double r = std::sqrt(c.x * c.x + c.y * c.y + z * z);
            E += c.t * std::polar(1.0 / r, k * r);
        }
        return std::norm(E);
    };

    // Scan a window around the nominal focus; widen it for low NA (long DOF).
    const double na = std::sin(std::atan(R / focal_length_um));
    const double half = std::max(2.0 * wavelength_um / (na * na + 1e-6), 0.25 * focal_length_um);
    const double zlo = std::max(0.1 * focal_length_um, focal_length_um - half);
    const double zhi = focal_length_um + half;
    const int N = 161;

    ThroughFocus tf;
    double peak = 0.0;
    std::vector<double> I(N);
    for (int i = 0; i < N; ++i) {
        double z = zlo + (zhi - zlo) * i / (N - 1);
        I[i] = on_axis(z);
        peak = std::max(peak, I[i]);
        tf.z_um.push_back(static_cast<float>(z));
    }
    if (peak <= 0) peak = 1;

    int ipk = 0;
    for (int i = 0; i < N; ++i) {
        tf.intensity.push_back(static_cast<float>(I[i] / peak));
        if (I[i] > I[ipk]) ipk = i;
    }
    tf.z_peak_um = tf.z_um[ipk];

    // Depth of focus = axial full width at half maximum of the peak.
    double half_lvl = peak / 2.0;
    int lo = -1, hi = -1;
    for (int i = 0; i < N; ++i)
        if (I[i] >= half_lvl) { if (lo < 0) lo = i; hi = i; }
    tf.dof_um = (hi > lo && lo >= 0) ? (tf.z_um[hi] - tf.z_um[lo]) : 0.0;
    return tf;
}

} // namespace celeris
