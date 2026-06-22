#include "celeris/materials/database.hpp"

#include "celeris/io/material_csv.hpp"

#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>

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

const Material& silicon_nitride() {
    // Si3N4 (stoichiometric LPCVD), Luke et al., Opt. Lett. 40, 4823 (2015).
    // n^2 = 1 + 3.0249 l^2/(l^2 - 0.1353406^2) + 40314 l^2/(l^2 - 1239.842^2).
    static const Material m = Material::sellmeier(
        {{3.0249, 0.1353406 * 0.1353406},
         {40314.0, 1239.842 * 1239.842}},
        "Si3N4");
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

const Material& sapphire() {
    // Synthetic sapphire (Al2O3), ordinary ray, Malitson & Dodge, J. Opt. Soc.
    // Am. 62, 1405 (1972). refractiveindex.info formula-1 coefficients pasted
    // verbatim: c0=0, then (strength, resonance-λ) pairs. Valid 0.20–5.0 µm.
    static const Material m = Material::sellmeier_ri_formula1(
        {0.0, 1.4313493, 0.0726631, 0.65054713, 0.1193242, 5.3414021, 18.028251},
        "Al2O3-sapphire");
    return m;
}

const Material& gallium_nitride() {
    // α-GaN (wurtzite), ordinary ray, Barker & Ilegems, Phys. Rev. B 7, 743
    // (1973). refractiveindex.info formula-1: c0=2.60, pairs (1.75, 0.256),
    // (4.1, 17.86). Valid 0.35–10 µm.
    static const Material m = Material::sellmeier_ri_formula1(
        {2.60, 1.75, 0.256, 4.1, 17.86}, "GaN");
    return m;
}

namespace {

// Canonicalize a name for lookup: lowercase, treat '_' as '-'.
std::string norm(std::string s) {
    for (char& c : s) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (c == '_') c = '-';
    }
    return s;
}

// One catalog entry. Either analytic (a function pointer) or tabulated (a CSV
// basename shipped in data/). `names[0]` is canonical; the rest are aliases.
struct Entry {
    std::vector<std::string> names;
    std::string description;
    bool tabulated;
    double lo, hi;            // validity range (µm); 0,0 for broadband analytic
    std::string data_file;    // tabulated only
    const Material& (*analytic)();  // analytic only, else nullptr
};

const std::vector<Entry>& table() {
    static const std::vector<Entry> t = {
        // --- analytic, lossless (transparent dielectrics) ---
        {{"air"}, "Air (n=1, lossless)", false, 0, 0, "", &air},
        {{"sio2", "fused-silica", "silica"},
         "Fused silica SiO2 (Malitson 1965)", false, 0, 0, "", &fused_silica},
        {{"bk7", "n-bk7"}, "Schott N-BK7 crown glass (Sellmeier)", false, 0, 0,
         "", &bk7},
        {{"si3n4", "sin", "silicon-nitride"},
         "Si3N4 LPCVD (Luke 2015)", false, 0, 0, "", &silicon_nitride},
        {{"al2o3", "sapphire", "alumina"},
         "Al2O3 sapphire, ordinary ray (Malitson-Dodge 1972)", false, 0.20, 5.0,
         "", &sapphire},
        {{"gan", "gallium-nitride"},
         "alpha-GaN, ordinary ray (Barker-Ilegems 1973)", false, 0.35, 10.0, "",
         &gallium_nitride},
        // --- tabulated real n,k (may be lossy), shipped in data/ ---
        {{"tio2", "titania"},
         "TiO2 ALD amorphous (Siefke 2016) -- Capasso-group visible metalens",
         true, 0.38, 0.80, "TiO2_Siefke.csv", nullptr},
        {{"c-si", "csi", "silicon", "crystalline-silicon"},
         "c-Si crystalline silicon (Green 2008)", true, 0.25, 1.45,
         "cSi_Green2008.csv", nullptr},
        {{"a-si", "asi", "amorphous-silicon"},
         "a-Si amorphous silicon (Pierce-Spicer 1972)", true, 0.10, 2.07,
         "aSi_Pierce.csv", nullptr},
        {{"au", "gold"}, "Au gold (Johnson-Christy 1972)", true, 0.19, 1.94,
         "Au_JohnsonChristy.csv", nullptr},
        {{"ag", "silver"}, "Ag silver (Johnson-Christy 1972)", true, 0.19, 1.94,
         "Ag_JohnsonChristy.csv", nullptr},
        {{"al", "aluminum", "aluminium"}, "Al aluminum (McPeak 2015)", true,
         0.15, 1.70, "Al_McPeak.csv", nullptr},
    };
    return t;
}

const Entry* find_entry(const std::string& name) {
    const std::string q = norm(name);
    for (const auto& e : table())
        for (const auto& n : e.names)
            if (n == q) return &e;
    return nullptr;
}

} // namespace

std::string find_data_file(const std::string& basename) {
    namespace fs = std::filesystem;
    std::vector<std::string> roots;
    if (const char* env = std::getenv("CELERIS_DATA_DIR"); env && *env)
        roots.emplace_back(env);
    roots.insert(roots.end(),
                 {"data", "../data", "../../data", "../../../data"});
    std::error_code ec;
    for (const auto& r : roots) {
        fs::path p = fs::path(r) / basename;
        if (fs::exists(p, ec)) return p.string();
    }
    if (fs::exists(basename, ec)) return basename;  // already a usable path
    return "";
}

std::vector<MaterialInfo> catalog() {
    std::vector<MaterialInfo> out;
    for (const auto& e : table()) {
        bool avail = true;
        if (e.tabulated) avail = !find_data_file(e.data_file).empty();
        out.push_back({e.names.front(), e.description, e.tabulated, e.lo, e.hi,
                       e.data_file, avail});
    }
    return out;
}

Material by_name(const std::string& name) {
    const Entry* e = find_entry(name);
    if (!e)
        throw std::runtime_error("unknown material '" + name +
                                 "' (try `celeris materials` to list them)");
    if (!e->tabulated) return e->analytic();  // value copy of the singleton
    const std::string path = find_data_file(e->data_file);
    if (path.empty())
        throw std::runtime_error("material '" + e->names.front() +
                                 "': data file " + e->data_file +
                                 " not found (set CELERIS_DATA_DIR or run from "
                                 "the repo root)");
    return load_material_csv(path, e->names.front());
}

} // namespace celeris::materials
