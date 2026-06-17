#include "celeris/design/system_opt.hpp"

#include "celeris/analysis/focal.hpp"
#include "celeris/design/metalens.hpp"

#include <algorithm>

namespace celeris {

SystemOptResult optimize_system(const Material& pillar, const Material& background,
                                const Material& incident, const Material& substrate,
                                double focal_um, double diameter_um,
                                double wavelength_um, double period_lo,
                                double period_hi, double thick_lo, double thick_hi,
                                int grid, int M, int fill_samples,
                                double efficiency_weight,
                                const std::function<void(float)>& progress) {
    grid = std::max(2, grid);
    SystemOptResult best{period_lo, thick_lo, 0, 0, -1};
    const int total = grid * grid;
    int done = 0;

    for (int ip = 0; ip < grid; ++ip) {
        double period = period_lo + (period_hi - period_lo) * ip / (grid - 1);
        for (int it = 0; it < grid; ++it) {
            double thick = thick_lo + (thick_hi - thick_lo) * it / (grid - 1);

            auto lib = build_unit_cell_library(pillar, background, incident, substrate,
                                               period, wavelength_um, thick, 0.08,
                                               0.92, fill_samples, M);
            auto lens = design_metalens(lib, focal_um, diameter_um);
            auto foc = analyze_focus(lens, lib, focal_um, wavelength_um, diameter_um);
            double merit = foc.strehl + efficiency_weight * lens.mean_amplitude;
            if (merit > best.merit)
                best = {period, thick, foc.strehl, lens.mean_amplitude, merit};

            if (progress) progress(static_cast<float>(++done) / total);
        }
    }
    return best;
}

} // namespace celeris
