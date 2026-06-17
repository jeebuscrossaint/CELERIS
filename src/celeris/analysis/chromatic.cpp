#include "celeris/analysis/chromatic.hpp"

#include "celeris/core.hpp"

#include <cmath>

namespace celeris {

std::vector<ChromaticPoint> analyze_chromatic(const MetalensDesign& lens,
                                              const UnitCellLibrary& lib,
                                              double design_focal,
                                              double design_wavelength,
                                              double diameter, double lambda_min,
                                              double lambda_max, int n) {
    const double p = lens.period_um;
    const double center = (lens.n_cells - 1) / 2.0;
    const double R_ap = diameter / 2.0;

    // Aperture cells (fixed designed transmission; geometric-phase model).
    struct Cell { double x, y; cdouble t; };
    std::vector<Cell> cells;
    for (int iy = 0; iy < lens.n_cells; ++iy)
        for (int ix = 0; ix < lens.n_cells; ++ix) {
            const double x = (ix - center) * p, y = (iy - center) * p;
            if (std::sqrt(x * x + y * y) > R_ap) continue;
            const double fill = lens.fill_map[static_cast<std::size_t>(iy) * lens.n_cells + ix];
            cells.push_back({x, y, lib.transmission_for_fill(fill)});
        }

    auto on_axis = [&](double z, double k) {
        cdouble E{0.0, 0.0};
        for (const Cell& c : cells) {
            double R = std::sqrt(c.x * c.x + c.y * c.y + z * z);
            E += c.t * std::polar(1.0 / R, k * R);
        }
        return std::norm(E);
    };

    // Reference peak at the design wavelength/focus for relative efficiency.
    double design_peak = on_axis(design_focal, 2.0 * pi / design_wavelength);

    std::vector<ChromaticPoint> out;
    for (int i = 0; i < n; ++i) {
        double lam = lambda_min + (lambda_max - lambda_min) * i / (n - 1);
        double k = 2.0 * pi / lam;
        double expected = design_focal * design_wavelength / lam;  // ~ analytic shift
        double zlo = 0.55 * expected, zhi = 1.6 * expected;
        const int NZ = 120;
        double best_z = expected, best_I = -1.0;
        for (int j = 0; j < NZ; ++j) {
            double z = zlo + (zhi - zlo) * j / (NZ - 1);
            double I = on_axis(z, k);
            if (I > best_I) { best_I = I; best_z = z; }
        }
        out.push_back({lam, best_z, design_peak > 0 ? best_I / design_peak : 0.0});
    }
    return out;
}

std::vector<ChromaticPoint> analyze_chromatic_dispersive(
    const MetalensDesign& lens,
    const std::function<UnitCellLibrary(double)>& build_lib_at,
    double design_focal, double design_wavelength, double diameter,
    double lambda_min, double lambda_max, int n) {
    const double p = lens.period_um;
    const double center = (lens.n_cells - 1) / 2.0;
    const double R_ap = diameter / 2.0;

    // Aperture cell geometry is fixed by the design; only the transmission each
    // cell imparts changes with wavelength (recomputed from the per-λ library).
    struct Cell { double x, y, fill; };
    std::vector<Cell> cells;
    for (int iy = 0; iy < lens.n_cells; ++iy)
        for (int ix = 0; ix < lens.n_cells; ++ix) {
            const double x = (ix - center) * p, y = (iy - center) * p;
            if (std::sqrt(x * x + y * y) > R_ap) continue;
            const double fill =
                lens.fill_map[static_cast<std::size_t>(iy) * lens.n_cells + ix];
            cells.push_back({x, y, fill});
        }

    // On-axis intensity at distance z, given each cell's transmission at this λ.
    auto on_axis = [&](double z, double k, const std::vector<cdouble>& t) {
        cdouble E{0.0, 0.0};
        for (std::size_t i = 0; i < cells.size(); ++i) {
            const Cell& c = cells[i];
            double R = std::sqrt(c.x * c.x + c.y * c.y + z * z);
            E += t[i] * std::polar(1.0 / R, k * R);
        }
        return std::norm(E);
    };

    auto transmissions_at = [&](const UnitCellLibrary& lib) {
        std::vector<cdouble> t(cells.size());
        for (std::size_t i = 0; i < cells.size(); ++i)
            t[i] = lib.transmission_for_fill(cells[i].fill);
        return t;
    };

    // Reference peak at the design wavelength using its own (rigorous) library.
    UnitCellLibrary design_lib = build_lib_at(design_wavelength);
    std::vector<cdouble> design_t = transmissions_at(design_lib);
    double design_peak =
        on_axis(design_focal, 2.0 * pi / design_wavelength, design_t);

    std::vector<ChromaticPoint> out;
    for (int i = 0; i < n; ++i) {
        double lam = lambda_min + (lambda_max - lambda_min) * i / (n - 1);
        double k = 2.0 * pi / lam;
        // Rebuild the meta-atom library at this wavelength: pillar index, fill
        // response and waveguide dispersion are all re-solved here.
        UnitCellLibrary lib = build_lib_at(lam);
        std::vector<cdouble> t = transmissions_at(lib);

        double expected = design_focal * design_wavelength / lam;  // analytic guess
        double zlo = 0.55 * expected, zhi = 1.6 * expected;
        const int NZ = 120;
        double best_z = expected, best_I = -1.0;
        for (int j = 0; j < NZ; ++j) {
            double z = zlo + (zhi - zlo) * j / (NZ - 1);
            double I = on_axis(z, k, t);
            if (I > best_I) { best_I = I; best_z = z; }
        }
        out.push_back({lam, best_z, design_peak > 0 ? best_I / design_peak : 0.0});
    }
    return out;
}

} // namespace celeris
