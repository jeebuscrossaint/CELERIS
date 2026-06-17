#include "celeris/analysis/tolerance.hpp"

#include "celeris/analysis/focal.hpp"

#include <algorithm>
#include <cmath>
#include <future>
#include <random>

namespace celeris {

std::vector<ToleranceResult> analyze_tolerance(
    const MetalensDesign& lens, const UnitCellLibrary& lib, double focal_length_um,
    double wavelength_um, double diameter_um, const std::vector<double>& sigmas_nm,
    int trials, unsigned seed) {
    std::vector<ToleranceResult> out;

    for (double sigma_nm : sigmas_nm) {
        // CD error in nm -> error in fill fraction (fill = pillar_size / period).
        const double sigma_fill = (sigma_nm / 1000.0) / lens.period_um;
        const int n = std::max(1, trials);

        // Independent trials run concurrently.
        std::vector<std::future<double>> futs;
        for (int t = 0; t < n; ++t) {
            futs.push_back(std::async(std::launch::async, [=, &lib] {
                std::mt19937 rng(seed + static_cast<unsigned>(t) * 7919u +
                                 static_cast<unsigned>(sigma_nm * 13.0));
                std::normal_distribution<double> g(0.0, sigma_fill);
                MetalensDesign perturbed = lens;  // copy, then jitter each pillar
                for (double& f : perturbed.fill_map)
                    f = std::clamp(f + g(rng), 0.05, 0.95);
                return analyze_focus(perturbed, lib, focal_length_um,
                                     wavelength_um, diameter_um)
                    .strehl;
            }));
        }

        double sum = 0.0, sum2 = 0.0, worst = 1e300;
        for (auto& f : futs) {
            double s = f.get();
            sum += s;
            sum2 += s * s;
            worst = std::min(worst, s);
        }
        double mean = sum / n;
        double var = std::max(0.0, sum2 / n - mean * mean);
        out.push_back({sigma_nm, mean, std::sqrt(var), worst});
    }
    return out;
}

} // namespace celeris
