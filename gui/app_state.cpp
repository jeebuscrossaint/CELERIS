#include "app_state.hpp"

#include <fstream>
#include <sstream>

#include "celeris/materials/database.hpp"

using namespace celeris;

namespace celeris::gui {

// --- global definitions -----------------------------------------------------
std::optional<Material> g_loaded_material;
std::string g_loaded_name = "(none)";

std::atomic<bool> g_running{false};
std::atomic<bool> g_pending{false};
std::atomic<bool> g_tol_pending{false};
std::atomic<bool> g_fov_pending{false};
std::atomic<bool> g_spot_pending{false};
std::atomic<bool> g_polar_pending{false};
std::atomic<bool> g_opt_pending{false};
std::atomic<bool> g_achro_pending{false};
std::atomic<float> g_progress{0.0f};
std::mutex g_mtx;
std::string g_status = "Ready — set parameters and click Design.";
Results g_res;

std::vector<ToleranceResult> g_tol;
std::vector<FieldPoint> g_fov;
std::vector<FieldPsf> g_spot;

PolarMetalensDesign g_polar;
double g_polar_fx = 0, g_polar_fy = 0, g_polar_zx = 0, g_polar_zy = 0,
       g_polar_iso_x = 0, g_polar_iso_y = 0;
PsfMap g_polar_psf_x, g_polar_psf_y;

SystemOptResult g_opt;

AchromaticDesign g_achro;
PbAchromaticDesign g_achro_pb;
AchroSummary g_achro_sum;
std::vector<float> g_achro_wl, g_achro_focus_std, g_achro_focus_ach;

// --- model helpers ----------------------------------------------------------
Material make_pillar(const Params& p) {
    switch (p.pillar_mat) {
        case 1: return materials::silicon_nitride();
        case 2: return materials::fused_silica();
        case 3: return Material::constant(cdouble{2.40, 0.0}, "TiO2~");
        case 4: return Material::constant(cdouble{3.50, 0.0}, "a-Si~");
        case 5: return Material::constant(cdouble{2.35, 0.0}, "GaN~");
        case 6:
            if (g_loaded_material) return *g_loaded_material;
            return Material::constant(cdouble{p.pillar_n, 0.0}, "custom");
        default: return Material::constant(cdouble{p.pillar_n, 0.0}, "custom");
    }
}

const Material& make_substrate(const Params& p) {
    switch (p.substrate_mat) {
        case 1: return materials::air();
        case 2: return materials::fused_silica();
        default: return materials::bk7();
    }
}

// --- project persistence (plain key/value .celeris file) --------------------
bool save_project(const std::string& path, const Params& p) {
    std::ofstream f(path);
    if (!f) return false;
    f << "celeris_project 1\n";
    f << "focal " << p.focal << "\n";
    f << "diameter " << p.diameter << "\n";
    f << "wavelength " << p.wavelength << "\n";
    f << "period " << p.period << "\n";
    f << "thickness " << p.thickness << "\n";
    f << "pillar_n " << p.pillar_n << "\n";
    f << "harmonics " << p.harmonics << "\n";
    f << "fill_samples " << p.fill_samples << "\n";
    f << "pillar_mat " << p.pillar_mat << "\n";
    f << "substrate_mat " << p.substrate_mat << "\n";
    f << "focal_y " << p.focal_y << "\n";
    f << "band_frac " << p.band_frac << "\n";
    f << "gd_weight " << p.gd_weight << "\n";
    f << "achro_lib " << p.achro_lib << "\n";
    f << "etch_height " << p.etch_height << "\n";
    for (const auto& L : p.extra_layers)
        f << "layer " << L.n << " " << L.fill << " " << L.thickness << "\n";
    return static_cast<bool>(f);
}

bool load_project(const std::string& path, Params& p) {
    std::ifstream f(path);
    if (!f) return false;
    Params np;
    np.extra_layers.clear();
    std::string line;
    bool header = false;
    while (std::getline(f, line)) {
        std::istringstream ss(line);
        std::string key;
        if (!(ss >> key)) continue;
        if (key == "celeris_project") { header = true; }
        else if (key == "focal") ss >> np.focal;
        else if (key == "diameter") ss >> np.diameter;
        else if (key == "wavelength") ss >> np.wavelength;
        else if (key == "period") ss >> np.period;
        else if (key == "thickness") ss >> np.thickness;
        else if (key == "pillar_n") ss >> np.pillar_n;
        else if (key == "harmonics") ss >> np.harmonics;
        else if (key == "fill_samples") ss >> np.fill_samples;
        else if (key == "pillar_mat") ss >> np.pillar_mat;
        else if (key == "substrate_mat") ss >> np.substrate_mat;
        else if (key == "focal_y") ss >> np.focal_y;
        else if (key == "band_frac") ss >> np.band_frac;
        else if (key == "gd_weight") ss >> np.gd_weight;
        else if (key == "achro_lib") ss >> np.achro_lib;
        else if (key == "etch_height") ss >> np.etch_height;
        else if (key == "layer") {
            LayerRow L; ss >> L.n >> L.fill >> L.thickness;
            np.extra_layers.push_back(L);
        }
    }
    if (!header) return false;
    p = std::move(np);
    return true;
}

// --- session auto-restore ---------------------------------------------------
// The session file is the project format plus a UI-preferences line (dark mode).
// Written on exit, loaded silently on startup, so the app reopens where you left
// it. load_project ignores the extra dark_mode key, so we re-scan for it here.
bool save_session(const std::string& path, const Params& p, bool dark) {
    if (!save_project(path, p)) return false;
    std::ofstream f(path, std::ios::app);
    if (!f) return false;
    f << "dark_mode " << (dark ? 1 : 0) << "\n";
    return static_cast<bool>(f);
}

bool load_session(const std::string& path, Params& p, bool& dark) {
    if (!load_project(path, p)) return false;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream ss(line);
        std::string key;
        if (ss >> key && key == "dark_mode") {
            int d = 0;
            if (ss >> d) dark = (d != 0);
            break;
        }
    }
    return true;
}

} // namespace celeris::gui
