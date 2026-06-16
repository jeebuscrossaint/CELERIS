#include "celeris/materials/database.hpp"

namespace celeris::materials {

const Material& air() {
    static const Material m = Material::constant(cdouble{1.0, 0.0}, "air");
    return m;
}

const Material& bk7() {
    // Schott N-BK7, Sellmeier coefficients {B_i, C_i} (C in µm²).
    static const Material m = Material::sellmeier(
        {{1.03961212, 0.00600069867},
         {0.231792344, 0.0200179144},
         {1.01046945, 103.560653}},
        "N-BK7");
    return m;
}

const Material& fused_silica() {
    // Malitson 1965 for fused silica. C_i are the squares of the published
    // resonance wavelengths (µm²): 0.0684043², 0.1162414², 9.896161².
    static const Material m = Material::sellmeier(
        {{0.6961663, 0.0684043 * 0.0684043},
         {0.4079426, 0.1162414 * 0.1162414},
         {0.8974794, 9.896161 * 9.896161}},
        "fused-silica");
    return m;
}

} // namespace celeris::materials
