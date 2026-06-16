#pragma once
// Rigorous Coupled-Wave Analysis (RCWA) for a single 1D grating layer.
//
// This is the core engine, in its simplest form. The method:
//   1. Expand the periodic permittivity ε(x) in a Fourier series (the orders).
//   2. Inside the grating, the coupled-wave equations become an eigenvalue
//      problem; its eigenmodes are the "Bloch" modes of the structure.
//   3. Match tangential fields at the top/bottom interfaces to the incident,
//      reflected, and transmitted plane-wave orders.
//   4. Read off the diffraction efficiency carried by each order.
//
// This first cut implements TE (s) polarization for one lamellar layer between
// two homogeneous half-spaces. TM and multi-layer S-matrix stacking come next.
//
// References: Moharam & Gaylord, "Formulation for stable and efficient
// implementation of the rigorous coupled-wave analysis," JOSA A 12 (1995).

#include "celeris/core.hpp"
#include "celeris/materials/material.hpp"
#include "celeris/rcwa/grating1d.hpp"

#include <vector>

namespace celeris {

struct Rcwa1DResult {
    std::vector<int> orders;     // diffraction order index m for each entry
    std::vector<double> de_r;    // reflected diffraction efficiency per order
    std::vector<double> de_t;    // transmitted diffraction efficiency per order
    double sum_de;               // Σ(de_r + de_t); == 1 for a lossless stack
};

// Solve one binary grating layer.
//   incident/substrate : the half-spaces above/below the grating
//   wavelength_um      : vacuum wavelength
//   theta0_rad         : incidence angle from normal (planar/classical mount)
//   n_harmonics (M)    : orders retained run m = -M..+M (so 2M+1 total).
//                        Larger M = more accurate, slower (eigenproblem is
//                        (2M+1)-dimensional). Convergence in M is expected.
//   pol                : TE only for now (TM throws).
Rcwa1DResult solve_rcwa_1d(const Material& incident,
                           const BinaryGrating1D& grating,
                           const Material& substrate,
                           double wavelength_um,
                           double theta0_rad,
                           int n_harmonics,
                           Pol pol);

// A multilayer stack: a shared period and an ordered list of layers, top
// (incident side) to bottom (substrate side).
struct Rcwa1DStack {
    double period_um;
    std::vector<GratingLayer1D> layers;
};

// Multilayer RCWA via stable scattering-matrix (Redheffer) recursion. Unlike a
// naive transfer-matrix product, the S-matrix stays bounded for thick or
// strongly evanescent layers. Reduces exactly to the single-layer solver above
// when the stack has one layer.
Rcwa1DResult solve_rcwa_1d(const Material& incident,
                           const Rcwa1DStack& stack,
                           const Material& substrate,
                           double wavelength_um,
                           double theta0_rad,
                           int n_harmonics,
                           Pol pol);

} // namespace celeris
