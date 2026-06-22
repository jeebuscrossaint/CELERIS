#include "celeris/materials/material.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace celeris {

Material::Material(std::function<cdouble(double)> dispersion, std::string name)
    : dispersion_(std::move(dispersion)), name_(std::move(name)) {}

Material Material::constant(cdouble n, std::string name) {
    return Material([n](double) { return n; }, std::move(name));
}

Material Material::sellmeier(std::vector<std::array<double, 2>> bc_terms,
                             std::string name) {
    return Material(
        [terms = std::move(bc_terms)](double wl) -> cdouble {
            const double l2 = wl * wl;
            double n2 = 1.0;
            for (const auto& [B, C] : terms) {
                n2 += B * l2 / (l2 - C);
            }
            return cdouble{std::sqrt(n2), 0.0};
        },
        std::move(name));
}

Material Material::sellmeier_ri_formula1(std::vector<double> coeffs,
                                         std::string name) {
    if (coeffs.empty() || coeffs.size() % 2 == 0) {
        throw std::invalid_argument(
            "Material::sellmeier_ri_formula1: need c0 + (strength, "
            "resonance) pairs (an odd-length coefficient list)");
    }
    return Material(
        [c = std::move(coeffs)](double wl) -> cdouble {
            const double l2 = wl * wl;
            double n2 = 1.0 + c[0];
            for (std::size_t i = 1; i + 1 < c.size(); i += 2) {
                const double res2 = c[i + 1] * c[i + 1];  // resonance λ², squared here
                n2 += c[i] * l2 / (l2 - res2);
            }
            return cdouble{std::sqrt(n2), 0.0};
        },
        std::move(name));
}

Material Material::tabulated(std::vector<double> wavelength_um,
                             std::vector<double> n,
                             std::vector<double> k,
                             std::string name) {
    if (wavelength_um.size() < 2 || wavelength_um.size() != n.size() ||
        n.size() != k.size()) {
        throw std::invalid_argument(
            "Material::tabulated: need matching arrays of length >= 2");
    }
    if (!std::is_sorted(wavelength_um.begin(), wavelength_um.end())) {
        throw std::invalid_argument(
            "Material::tabulated: wavelengths must be ascending");
    }

    return Material(
        [wl = std::move(wavelength_um), n = std::move(n),
         k = std::move(k)](double query) -> cdouble {
            // Clamp out-of-range queries to the table endpoints.
            if (query <= wl.front()) return cdouble{n.front(), k.front()};
            if (query >= wl.back()) return cdouble{n.back(), k.back()};

            // Find the bracketing interval [hi-1, hi] and lerp within it.
            auto hi = std::lower_bound(wl.begin(), wl.end(), query);
            std::size_t i = static_cast<std::size_t>(hi - wl.begin());
            double t = (query - wl[i - 1]) / (wl[i] - wl[i - 1]);
            double nn = n[i - 1] + t * (n[i] - n[i - 1]);
            double kk = k[i - 1] + t * (k[i] - k[i - 1]);
            return cdouble{nn, kk};
        },
        std::move(name));
}

} // namespace celeris
