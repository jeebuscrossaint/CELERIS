#pragma once
// Polarization / form-birefringence analysis. A rectangular nanopillar is
// anisotropic: an x-polarized wave and a y-polarized wave see different
// effective structure, so the pillar imparts different phase to each. That
// retardance (phi_x - phi_y) is the basis of polarization optics built from
// metasurfaces -- waveplates, polarization-multiplexed metalenses (two
// independent phase profiles on the two polarizations), the core of
// polarization imaging (e.g. Metalenz). A square pillar is symmetric, so its
// retardance is exactly zero; stretching it one way opens a controllable
// birefringence.
//
// This sweeps the pillar's x-width at a fixed y-width and reports both
// co-polarized phases (from two RCWA solves, x- and y-polarized illumination)
// and their retardance -- a quantitative map of the achievable birefringence.

#include "celeris/materials/material.hpp"

#include <vector>

namespace celeris {

struct BirefringencePoint {
    double fill_x;          // swept pillar width (fraction of period) along x
    double fill_y;          // fixed pillar width along y
    double phase_x_deg;     // phase imparted to x-polarized light (deg)
    double phase_y_deg;     // phase imparted to y-polarized light (deg)
    double retardance_deg;  // phi_x - phi_y, wrapped to (-180, 180]
    double tx;              // |t| for x polarization (transmission amplitude)
    double ty;              // |t| for y polarization
};

// Sweep fill_x over [fill_min, fill_max] (n_samples points) at fixed fill_y,
// solving the 2D RCWA twice per point (x- and y-polarized normal incidence).
// pillar/background are the meta-atom materials; incident/substrate the
// half-spaces. Solves run in parallel across CPU cores.
std::vector<BirefringencePoint> analyze_birefringence(
    const Material& pillar, const Material& background, const Material& incident,
    const Material& substrate, double period_um, double wavelength_um,
    double thickness_um, double fill_y, double fill_min, double fill_max,
    int n_samples, int M);

} // namespace celeris
