#pragma once
// System-level design optimizer — the "merit function" idea from real optical
// design tools, applied to a metalens. Instead of the user guessing the period
// and pillar height, search that design space and return the combination that
// maximizes a merit (focal Strehl, optionally weighted by transmission
// efficiency). Each trial builds a full unit-cell library, assembles the lens,
// and measures the actual focal performance — so it optimizes the real metric,
// not a proxy.

#include "celeris/materials/material.hpp"

#include <functional>

namespace celeris {

struct SystemOptResult {
    double period_um;
    double thickness_um;
    double strehl;
    double mean_amplitude;
    double merit;
};

// Grid-search period × pillar-height; return the best by merit =
// strehl + efficiency_weight * mean_transmission. `progress` is called with a
// 0..1 fraction. Uses a coarse library (small M / few samples) for speed; the
// winning parameters should then be re-run at full quality.
SystemOptResult optimize_system(const Material& pillar, const Material& background,
                                const Material& incident, const Material& substrate,
                                double focal_um, double diameter_um,
                                double wavelength_um, double period_lo,
                                double period_hi, double thick_lo, double thick_hi,
                                int grid, int M, int fill_samples,
                                double efficiency_weight,
                                const std::function<void(float)>& progress);

} // namespace celeris
