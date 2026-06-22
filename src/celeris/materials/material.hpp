#pragma once
// A material's optical response: the complex refractive index n(λ) = n + i*k
// as a function of vacuum wavelength λ (in micrometers, the optics convention).
//
// Design note: Material is *value-semantic*. It holds its dispersion model as a
// std::function, so you can copy it freely and store it by value in a Layer.
// No inheritance, no pointers, no lifetime headaches — important because the
// engine passes materials around constantly.

#include "celeris/core.hpp"

#include <array>
#include <functional>
#include <string>
#include <vector>

namespace celeris {

class Material {
public:
    Material(std::function<cdouble(double)> dispersion, std::string name);

    // Complex refractive index n + i*k at the given vacuum wavelength (µm).
    cdouble index(double wavelength_um) const { return dispersion_(wavelength_um); }

    // Relative permittivity ε = n². (RCWA works in ε; TMM works in n.)
    cdouble permittivity(double wavelength_um) const {
        cdouble n = index(wavelength_um);
        return n * n;
    }

    const std::string& name() const { return name_; }

    // ---- Factories --------------------------------------------------------

    // Non-dispersive material with a fixed index (use cdouble{n, k} for loss).
    static Material constant(cdouble n, std::string name = "constant");

    // Sellmeier dispersion (lossless, k = 0):
    //   n² = 1 + Σ_i  B_i · λ² / (λ² − C_i),   λ in µm.
    // Each term is {B_i, C_i}. This is how glass catalogs specify dispersion.
    static Material sellmeier(std::vector<std::array<double, 2>> bc_terms,
                              std::string name = "sellmeier");

    // refractiveindex.info "formula 1" (Sellmeier, lossless):
    //   n² = 1 + c₀ + Σ_i  c_{2i-1} · λ² / (λ² − c_{2i}²),   λ in µm.
    // `coeffs` = the raw coefficient list straight off a refractiveindex.info
    // YAML (c₀ followed by (strength, resonance-wavelength) pairs — note the
    // resonance is a WAVELENGTH that gets SQUARED here, unlike sellmeier()'s
    // pre-squared C). Lets canonical published dispersions be pasted verbatim.
    static Material sellmeier_ri_formula1(std::vector<double> coeffs,
                                          std::string name = "formula1");

    // Tabulated n and k vs wavelength (µm), linearly interpolated.
    // `wavelength_um` must be strictly ascending; out-of-range queries clamp
    // to the nearest endpoint. This is what you'll load from
    // refractiveindex.info CSVs for real materials (TiO2, Si, …).
    static Material tabulated(std::vector<double> wavelength_um,
                              std::vector<double> n,
                              std::vector<double> k,
                              std::string name = "tabulated");

private:
    std::function<cdouble(double)> dispersion_;
    std::string name_;
};

} // namespace celeris
