#pragma once
// Inverse design: instead of looking up a pillar from a discrete library, SOLVE
// for the geometry that best meets a target. This is the Tier-2 differentiator
// — "tell me the goal, I'll find the structure."
//
// This first version is gradient descent (Adam) with numerical gradients, which
// is exact and well-suited to few-parameter unit cells (here: pillar fill and
// height). A true analytic adjoint — differentiating through the RCWA
// eigensolve/S-matrix — would cut the cost from O(params) solves per step to
// O(1) and is the natural next optimization; it is not needed for correctness.

#include "celeris/materials/material.hpp"

namespace celeris {

struct PillarTarget {
    double wavelength_um;
    double target_phase_rad;   // desired zeroth-order transmission phase
    double amplitude_weight;   // weight on the (1 - |t|)^2 transmission term
};

struct OptimizedPillar {
    double fill;                // optimized square-pillar fill fraction
    double thickness_um;        // optimized pillar height
    double achieved_phase_rad;
    double achieved_amplitude;  // |t| of the result
    double loss;
    int iterations;
};

// Optimize a square pillar (design variables: fill and thickness) to meet the
// target, starting from (fill0, thickness0). Each step runs a few 2D RCWA
// solves for the numerical gradient. Mx = My = M.
OptimizedPillar optimize_pillar(const Material& pillar, const Material& background,
                                const Material& incident, const Material& substrate,
                                double period_um, const PillarTarget& target,
                                int M, double fill0, double thickness0,
                                int max_iters);

} // namespace celeris
