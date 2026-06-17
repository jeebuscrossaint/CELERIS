#include "celeris/analysis/field.hpp"

#include "celeris/core.hpp"

#include <cmath>
#include <vector>

namespace celeris {

std::vector<FieldPoint> analyze_field_of_view(const MetalensDesign& lens,
                                              const UnitCellLibrary& lib,
                                              double focal_length_um,
                                              double wavelength_um,
                                              double diameter_um,
                                              const std::vector<double>& angles_deg) {
    const double p = lens.period_um;
    const double center = (lens.n_cells - 1) / 2.0;
    const double R_ap = diameter_um / 2.0;
    const double k = 2.0 * pi / wavelength_um;

    struct Cell { double x, y; cdouble t; };
    std::vector<Cell> cells;
    for (int iy = 0; iy < lens.n_cells; ++iy)
        for (int ix = 0; ix < lens.n_cells; ++ix) {
            const double x = (ix - center) * p, y = (iy - center) * p;
            if (std::sqrt(x * x + y * y) > R_ap) continue;
            const double fill = lens.fill_map[static_cast<std::size_t>(iy) * lens.n_cells + ix];
            cells.push_back({x, y, lib.transmission_for_fill(fill)});
        }

    // Focal-plane field for an incident tilt theta (in x-z): the illumination
    // adds a linear phase exp(i k sin(theta) x) across the aperture.
    auto peak_near = [&](double theta, double xc) {
        const double s = std::sin(theta);
        const double dl = wavelength_um * focal_length_um / diameter_um;
        const double W = std::max(3.0 * dl, 2.0);  // half-window around xc
        const int NG = 61;
        double best = 0.0, best_x = xc;
        for (int j = 0; j < NG; ++j) {
            double fy = -W + 2.0 * W * j / (NG - 1);
            for (int i = 0; i < NG; ++i) {
                double fx = xc - W + 2.0 * W * i / (NG - 1);
                cdouble E{0.0, 0.0};
                for (const Cell& c : cells) {
                    double R = std::sqrt((fx - c.x) * (fx - c.x) +
                                         (fy - c.y) * (fy - c.y) +
                                         focal_length_um * focal_length_um);
                    E += c.t * std::polar(1.0, k * s * c.x) * std::polar(1.0 / R, k * R);
                }
                double I = std::norm(E);
                if (I > best) { best = I; best_x = fx; }
            }
        }
        return std::pair<double, double>{best, best_x};
    };

    double on_axis_peak = peak_near(0.0, 0.0).first;

    std::vector<FieldPoint> out;
    for (double ang : angles_deg) {
        double theta = ang * pi / 180.0;
        double xc = focal_length_um * std::tan(theta);  // paraxial spot position
        auto [pk, px] = peak_near(theta, xc);
        out.push_back({ang, on_axis_peak > 0 ? pk / on_axis_peak : 0.0, px});
    }
    return out;
}

} // namespace celeris
