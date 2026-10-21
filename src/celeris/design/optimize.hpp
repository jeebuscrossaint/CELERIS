#pragma once
/// @file
/// @brief Inverse design of a single pillar: solve for the geometry that meets a target.
///
/// Instead of looking up a pillar from a discrete library, this searches for
/// the geometry that best meets a target ("tell me the goal, I'll find the
/// structure").
///
/// The method is gradient descent (Adam) with central-difference numerical
/// gradients, which is exact and well-suited to few-parameter unit cells (here:
/// pillar fill and height). A true analytic adjoint — differentiating through
/// the RCWA eigensolve/S-matrix — would cut the cost from O(params) solves per
/// step to O(1); it is planned, and not needed for correctness.

#include "celeris/materials/material.hpp"

namespace celeris {

/// @brief What the optimized pillar should achieve.
struct PillarTarget {
    double wavelength_um;      ///< Design wavelength (µm).
    double target_phase_rad;   ///< Desired zeroth-order transmission phase (rad).
    double amplitude_weight;   ///< Weight on the (1 - |t|)² transmission term.
};

/// @brief Result of optimize_pillar().
struct OptimizedPillar {
    double fill;                ///< Optimized square-pillar fill fraction.
    double thickness_um;        ///< Optimized pillar height (µm).
    double achieved_phase_rad;  ///< Transmission phase of the result (rad).
    double achieved_amplitude;  ///< |t| of the result.
    double loss;                ///< Final loss value.
    int iterations;             ///< Iterations run.
};

/// @brief Optimize a square pillar's fill and height to meet a target.
///
/// Each step runs a few 2D RCWA solves for the numerical gradient.
/// @param pillar,background  Pillar and surrounding materials.
/// @param incident,substrate Half-spaces above and below.
/// @param period_um          Lattice period (µm).
/// @param target             Target phase and transmission weighting.
/// @param M                  Harmonic half-count (Mx = My = M).
/// @param fill0,thickness0   Starting fill and height (µm).
/// @param max_iters          Iteration cap.
OptimizedPillar optimize_pillar(const Material& pillar, const Material& background,
                                const Material& incident, const Material& substrate,
                                double period_um, const PillarTarget& target,
                                int M, double fill0, double thickness0,
                                int max_iters);

} // namespace celeris
