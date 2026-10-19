#pragma once
/// @file
/// @brief Rigorous coupled-wave analysis (RCWA) for 1D gratings and stacks.
///
/// The method:
///   1. Expand the periodic permittivity ε(x) in a Fourier series (the orders).
///   2. Inside the grating, the coupled-wave equations become an eigenvalue
///      problem; its eigenmodes are the "Bloch" modes of the structure.
///   3. Match tangential fields at the top/bottom interfaces to the incident,
///      reflected, and transmitted plane-wave orders.
///   4. Read off the diffraction efficiency carried by each order.
///
/// Both polarizations are supported; TM uses Li's inverse-rule factorization.
///
/// References: Moharam & Gaylord, "Formulation for stable and efficient
/// implementation of the rigorous coupled-wave analysis," JOSA A 12 (1995);
/// L. Li, JOSA A 13 (1996).

#include "celeris/core.hpp"
#include "celeris/materials/material.hpp"
#include "celeris/rcwa/grating1d.hpp"

#include <vector>

namespace celeris {

/// @brief Diffraction efficiencies of a 1D RCWA solve, one entry per order.
struct Rcwa1DResult {
    std::vector<int> orders;     ///< Diffraction order index m for each entry.
    std::vector<double> de_r;    ///< Reflected diffraction efficiency per order.
    std::vector<double> de_t;    ///< Transmitted diffraction efficiency per order.
    double sum_de;               ///< Σ(de_r + de_t); equals 1 for a lossless stack.
};

/// @brief Solve one binary grating layer between two half-spaces.
/// @param incident       Half-space above the grating (light comes from here).
/// @param grating        The grating layer (carries the period).
/// @param substrate      Half-space below the grating.
/// @param wavelength_um  Vacuum wavelength (µm).
/// @param theta0_rad     Incidence angle from the normal (planar/classical mount).
/// @param n_harmonics    M: orders m = -M..+M are retained (2M+1 total). Larger M
///                       is more accurate and slower; results converge in M.
/// @param pol            TE or TM.
/// @return Per-order reflected/transmitted efficiencies.
/// @throws std::invalid_argument if `n_harmonics < 0`.
Rcwa1DResult solve_rcwa_1d(const Material& incident,
                           const BinaryGrating1D& grating,
                           const Material& substrate,
                           double wavelength_um,
                           double theta0_rad,
                           int n_harmonics,
                           Pol pol);

/// @brief A 1D multilayer stack: one shared period and an ordered list of layers.
struct Rcwa1DStack {
    double period_um;                    ///< Period shared by every layer (µm).
    std::vector<GratingLayer1D> layers;  ///< Layers, top (incident side) to bottom.
};

/// @brief Multilayer 1D RCWA via a stable scattering-matrix (Redheffer) recursion.
///
/// Unlike a naive transfer-matrix product, the S-matrix stays bounded for thick
/// or strongly evanescent layers. Reduces exactly to the single-layer overload
/// when the stack has one layer. Parameters are as for the single-layer overload.
/// @param incident       Half-space above the stack.
/// @param stack          The layers and their shared period.
/// @param substrate      Half-space below the stack.
/// @param wavelength_um  Vacuum wavelength (µm).
/// @param theta0_rad     Incidence angle from the normal.
/// @param n_harmonics    M: orders m = -M..+M are retained.
/// @param pol            TE or TM.
/// @return Per-order reflected/transmitted efficiencies.
Rcwa1DResult solve_rcwa_1d(const Material& incident,
                           const Rcwa1DStack& stack,
                           const Material& substrate,
                           double wavelength_um,
                           double theta0_rad,
                           int n_harmonics,
                           Pol pol);

} // namespace celeris
