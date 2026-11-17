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

/// @brief A material's complex refractive index n(λ) + i·k(λ) vs vacuum wavelength (µm).
///
/// Value-semantic: the dispersion model is held as a `std::function`, so a
/// Material can be copied and stored by value. Build one with the static
/// factories or get a built-in from `celeris::materials`.
class Material {
public:
    /// @brief Wrap an arbitrary dispersion model.
    /// @param dispersion Callable returning n + i·k at a vacuum wavelength (µm).
    /// @param name       Label used in reports and `__repr__`.
    Material(std::function<cdouble(double)> dispersion, std::string name);

    /// @brief Complex refractive index n + i·k at the given vacuum wavelength (µm).
    cdouble index(double wavelength_um) const { return dispersion_(wavelength_um); }

    /// @brief Relative permittivity ε = n² at the given vacuum wavelength (µm).
    ///
    /// RCWA works in ε; TMM works in n.
    cdouble permittivity(double wavelength_um) const {
        cdouble n = index(wavelength_um);
        return n * n;
    }

    /// @brief The material's label.
    const std::string& name() const { return name_; }

    // ---- Factories --------------------------------------------------------

    /// @brief Non-dispersive material with a fixed index.
    /// @param n    Complex index; use cdouble{n, k} with k > 0 for loss.
    /// @param name Label.
    static Material constant(cdouble n, std::string name = "constant");

    /// @brief Lossless Sellmeier dispersion, n² = 1 + Σ_i B_i·λ² / (λ² − C_i), λ in µm.
    ///
    /// This is how glass catalogs specify dispersion.
    /// @param bc_terms One {B_i, C_i} pair per term (C_i in µm², already squared).
    /// @param name     Label.
    static Material sellmeier(std::vector<std::array<double, 2>> bc_terms,
                              std::string name = "sellmeier");

    /// @brief refractiveindex.info "formula 1" (lossless Sellmeier).
    ///
    /// n² = 1 + c₀ + Σ_i c_{2i−1}·λ² / (λ² − c_{2i}²), λ in µm. Unlike
    /// `sellmeier`, the resonance is a WAVELENGTH that is squared here, so the
    /// coefficient list can be pasted verbatim from a refractiveindex.info YAML.
    /// @param coeffs c₀ followed by (strength, resonance wavelength) pairs.
    /// @param name   Label.
    static Material sellmeier_ri_formula1(std::vector<double> coeffs,
                                          std::string name = "formula1");

    /// @brief Tabulated n and k vs wavelength, linearly interpolated.
    ///
    /// Queries outside the table CLAMP to the nearest endpoint (no
    /// extrapolation, no error), so check the table covers your band.
    /// @param wavelength_um Sample wavelengths (µm), strictly ascending.
    /// @param n             Real index at each sample.
    /// @param k             Extinction coefficient at each sample.
    /// @param name          Label.
    static Material tabulated(std::vector<double> wavelength_um,
                              std::vector<double> n,
                              std::vector<double> k,
                              std::string name = "tabulated");

private:
    std::function<cdouble(double)> dispersion_;
    std::string name_;
};

} // namespace celeris
