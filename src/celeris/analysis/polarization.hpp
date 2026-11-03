#pragma once
/// @file
/// @brief Polarization / form-birefringence analysis of rectangular pillars.
///
/// A rectangular nanopillar is anisotropic: an x-polarized wave and a y-polarized
/// wave see different effective structure, so the pillar imparts different phase
/// to each. That retardance (phi_x - phi_y) is the basis of polarization optics
/// built from metasurfaces -- waveplates, polarization-multiplexed metalenses (two
/// independent phase profiles on the two polarizations), the core of
/// polarization imaging (e.g. Metalenz). A square pillar is symmetric, so its
/// retardance is exactly zero; stretching it one way opens a controllable
/// birefringence.
///
/// This sweeps the pillar's x-width at a fixed y-width and reports both
/// co-polarized phases (from two RCWA solves, x- and y-polarized illumination)
/// and their retardance -- a quantitative map of the achievable birefringence.

#include "celeris/materials/material.hpp"

#include <vector>

namespace celeris {

/// @brief Co-polarized response of one rectangular pillar.
struct BirefringencePoint {
    double fill_x;          ///< Swept pillar width along x (fraction of the period).
    double fill_y;          ///< Fixed pillar width along y (fraction of the period).
    double phase_x_deg;     ///< Phase imparted to x-polarized light (degrees).
    double phase_y_deg;     ///< Phase imparted to y-polarized light (degrees).
    double retardance_deg;  ///< phi_x - phi_y, wrapped to (-180, 180] (degrees).
    double tx;              ///< |t| for x polarization (transmission amplitude).
    double ty;              ///< |t| for y polarization.
};

/// @brief Sweep pillar x-width at fixed y-width and report the birefringence.
///
/// Solves the 2D RCWA twice per point (x- and y-polarized normal incidence).
/// Solves run in parallel across CPU cores.
/// @param pillar,background  Meta-atom materials.
/// @param incident,substrate Half-space materials above and below.
/// @param period_um          Square unit-cell period (µm).
/// @param wavelength_um      Vacuum wavelength (µm).
/// @param thickness_um       Pillar height (µm).
/// @param fill_y             Fixed y-width (fraction of the period).
/// @param fill_min,fill_max  Range of the swept x-width (fractions of the period).
/// @param n_samples          Number of sweep points.
/// @param M                  Fourier half-count in each direction.
/// @return One BirefringencePoint per sample, in sweep order.
std::vector<BirefringencePoint> analyze_birefringence(
    const Material& pillar, const Material& background, const Material& incident,
    const Material& substrate, double period_um, double wavelength_um,
    double thickness_um, double fill_y, double fill_min, double fill_max,
    int n_samples, int M);

} // namespace celeris
